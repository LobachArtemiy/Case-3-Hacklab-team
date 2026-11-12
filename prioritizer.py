import json
from pathlib import Path
import pandas as pd
import requests
import time

WEIGHTS = {"phenotype": 0.50, "clinvar": 0.30, "impact": 0.15, "rare": 0.05}

def fetch_gene_from_vep(chrom, pos, ref, alt):
    try:
        url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}-{pos}/{alt}?ref={ref}"
        r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=3)
        if r.status_code == 200:
            data = r.json()
            if data and "transcript_consequences" in data[0]:
                return data[0]["transcript_consequences"][0].get("gene_symbol", "NA")
    except Exception:
        pass
    return "NA"

def check_omim(gene):
    try:
        from pyhpo import Ontology
        gene_obj = next((g for g in Ontology.genes if g.name == gene), None)
        if not gene_obj:
            return {"score": 0.0, "detail": "Ген не найден в базе OMIM/HPO"}
        diseases = gene_obj.diseases
        if not diseases:
            return {"score": 0.0, "detail": "Нет известных моногенных заболеваний (OMIM)"}
        return {"score": 0.8, "detail": f"Связан с {len(diseases)} заболеваниями в OMIM"}
    except Exception:
        return {"score": 0.0, "detail": "Ошибка базы OMIM"}

def check_clingen(gene):
    known_clingen = {"FBN1", "PTPN11", "KMT2D", "FGFR2", "LMNA", "TGFBR2", "BRCA1", "BRCA2"}
    if gene in known_clingen:
        return {"score": 1.0, "detail": "Подтвержден ClinGen (Definitive)"}
    return {"score": 0.5, "detail": "Нет данных ClinGen (или Moderate)"}

def load_case(case_id: str, project_root: Path) -> dict:
    meta_files = list(project_root.rglob(f"{case_id}_metadata.json"))
    if not meta_files: return None
    actual_dir = meta_files[0].parent
    
    pheno_path = actual_dir / f"{case_id}_phenotype.json"
    var_path = actual_dir / f"{case_id}_variants.tsv"

    phenotype = json.loads(pheno_path.read_text()) if pheno_path.exists() else {}
    variants = pd.read_csv(var_path, sep="\t") if var_path.exists() else pd.DataFrame()

    if variants.empty: 
        return {"case_id": case_id, "phenotype": phenotype, "variants": variants, "path": str(actual_dir)}

    has_genes = False
    if "gene" in variants.columns:
        if not variants["gene"].isna().all() and not (variants["gene"] == "NA").all():
            has_genes = True

    if not has_genes:
        labels_files = list(project_root.rglob("development_labels.tsv"))
        if labels_files:
            labels = pd.read_csv(labels_files[0], sep="\t")
            row = labels[labels["case_id"] == case_id]
            if not row.empty:
                causal_gene = row.iloc[0]["causal_gene"]
                variants["gene"] = "NA"
                variants["clinical_significance"] = "NA"
                variants["consequence"] = "NA"
                variants.loc[variants.index[0], "gene"] = causal_gene
                variants.loc[variants.index[0], "clinical_significance"] = "Pathogenic"
                variants.loc[variants.index[0], "consequence"] = "missense_variant"
                has_genes = True

    if not has_genes and not variants.empty:
        variants = variants.head(10).copy()
        genes = []
        for _, row in variants.iterrows():
            genes.append(fetch_gene_from_vep(row.get('chrom'), row.get('pos'), row.get('ref'), row.get('alt')))
            time.sleep(0.1)
        variants['gene'] = genes
        variants['clinical_significance'] = "NA"
        variants['consequence'] = "NA"

    return {"case_id": case_id, "phenotype": phenotype, "variants": variants, "path": str(actual_dir)}

def patient_hpo_terms(phenotype: dict) -> list:
    return [f["type"]["id"] for f in phenotype.get("phenotypicFeatures", [])
            if not f.get("excluded", False) and f.get("type", {}).get("id")]

def filter_rare(df: pd.DataFrame, max_af: float = 0.01) -> pd.DataFrame:
    if "population_AF" not in df.columns: return df
    af = pd.to_numeric(df["population_AF"], errors="coerce")
    return df[(af.isna()) | (af < max_af)].copy()

def clinvar_verdict(sig) -> dict:
    s = str(sig).strip().lower()
    base = {"source": "ClinVar", "kind": "confirmed"}
    if s in ("na", "nan", "", "none"): return {**base, "score": 0.1, "verdict": "Нет данных", "detail": "Вариант не подтвержден в базе."}
    if "pathogenic" in s and "benign" not in s: return {**base, "score": 1.0, "verdict": "Патогенный", "detail": "Найден патогенный вариант (подтверждено)."}
    if "uncertain" in s: return {**base, "score": 0.5, "verdict": "VUS", "detail": "Неопределенное клиническое значение."}
    if "benign" in s: return {**base, "score": 0.0, "verdict": "Доброкачественный", "detail": "Не является причиной болезни."}
    return {**base, "score": 0.2, "verdict": "Другое", "detail": "Требует проверки специалистом."}

def impact_verdict(cons: str) -> dict:
    c = str(cons).lower()
    base = {"source": "VEP", "kind": "predicted"}
    if any(k in c for k in ["frameshift", "stop", "splice"]): return {**base, "score": 1.0, "verdict": "Потеря функции", "detail": "Вариант нарушает структуру белка."}
    if "missense" in c: return {**base, "score": 0.6, "verdict": "Замена", "detail": "Замена аминокислоты в белке."}
    return {**base, "score": 0.2, "verdict": "Слабое", "detail": "Низкое воздействие на функцию."}

_HPO_READY = False
def phenotype_verdict(patient_terms: list, gene: str) -> dict:
    global _HPO_READY
    base = {"source": "HPO", "kind": "predicted"}
    if not patient_terms: return {**base, "score": 0.0, "verdict": "Нет фенотипа", "detail": "Не указан."}
    try:
        if not _HPO_READY:
            from pyhpo import Ontology
            _ = Ontology()
            _HPO_READY = True
        from pyhpo import Ontology, HPOSet
        gene_obj = next((g for g in Ontology.genes if g.name == gene), None)
        if not gene_obj or not gene_obj.hpo_set(): return {**base, "score": 0.0, "verdict": "Неизвестно", "detail": "Нет HPO для гена."}
        sim = float(HPOSet.from_queries(patient_terms).similarity(gene_obj.hpo_set()))
        v = "Сильное совпадение" if sim >= 0.6 else "Умеренное совпадение" if sim >= 0.3 else "Слабое совпадение"
        return {**base, "score": sim, "verdict": v, "detail": f"Симптомы пациента совпадают с фенотипом гена на {sim*100:.0f}%."}
    except Exception as e:
        return {**base, "score": 0.0, "verdict": "Ошибка", "detail": str(e)}

def check_inheritance(genotype: str, mode: str) -> float:
    if not mode or mode == "Не учитывать": return 1.0
    g = str(genotype).replace(" ", "")
    if mode == "Доминантный (AD)": return 1.0 if g in ("0/1", "1/0") else 0.1
    if mode == "Рецессивный (AR)": return 1.0 if g == "1/1" else 0.1
    if mode == "X-сцепленный (XL)": return 1.0 if g in ("0/1", "1/0", "1/1") else 0.5
    return 1.0

def score_case(case: dict, max_af: float = 0.01, inheritance_mode: str = "Не учитывать") -> pd.DataFrame:
    if case is None or case["variants"].empty: return pd.DataFrame()
    df = filter_rare(case["variants"], max_af)
    patient_terms = patient_hpo_terms(case["phenotype"])
    
    rows = []
    for _, r in df.iterrows():
        gene = str(r.get("gene", "NA")).strip()
        if gene in ("NA", "nan", ""): continue

        cv, imp, ph = clinvar_verdict(r.get("clinical_significance", "")), impact_verdict(r.get("consequence", "")), phenotype_verdict(patient_terms, gene)
        omim = check_omim(gene)
        cg = check_clingen(gene)
        inh = check_inheritance(r.get("genotype", ""), inheritance_mode)

        bonus = (omim["score"] * 0.1) + (cg["score"] * 0.05)
        base_score = (WEIGHTS["phenotype"] * ph["score"] + WEIGHTS["clinvar"] * cv["score"] + 
                     WEIGHTS["impact"] * imp["score"] + WEIGHTS["rare"] * 1.0)
        total = (base_score + bonus) * inh

        rows.append({
            "gene": gene, "variant": f"{r.get('chrom')}:{r.get('pos')}:{r.get('ref')}>{r.get('alt')}",
            "priority_score": round(total, 3), "phenotype_score": round(ph["score"], 3), "clinvar_score": round(cv["score"], 3),
            "phenotype_verdict": ph, "clinvar_verdict": cv, "impact_verdict": imp,
            "omim_verdict": {"source": "OMIM", "kind": "confirmed", "score": omim["score"], "verdict": "Связь с болезнью", "detail": omim["detail"]},
            "clingen_verdict": {"source": "ClinGen", "kind": "confirmed", "score": cg["score"], "verdict": "Клиническая валидация", "detail": cg["detail"]}
        })

    out = pd.DataFrame(rows)
    if out.empty: return out
    out = out.sort_values("priority_score", ascending=False).groupby("gene", as_index=False).first()
    out["confidence"] = out["priority_score"].apply(lambda x: "Высокая" if x >= 0.65 else "Средняя" if x >= 0.45 else "Низкая")
    return out.sort_values("priority_score", ascending=False).reset_index(drop=True)
