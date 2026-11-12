import json
from pathlib import Path
import pandas as pd

# Веса строго по ACMG (фенотип, база, прогноз, редкость)
WEIGHTS = {"phenotype": 0.50, "clinvar": 0.30, "impact": 0.15, "rare": 0.05}

def load_case(case_id: str, project_root: Path) -> dict:
    # Умный поиск: ищем файлы во ВСЕХ подпапках (development, hidden и т.д.)
    meta_files = list(project_root.rglob(f"{case_id}_metadata.json"))
    if not meta_files:
        return None
        
    actual_dir = meta_files[0].parent
    pheno_path = actual_dir / f"{case_id}_phenotype.json"
    ann_path = actual_dir / f"{case_id}_annotated.tsv"
    var_path = actual_dir / f"{case_id}_variants.tsv"

    phenotype = json.loads(pheno_path.read_text()) if pheno_path.exists() else {}
    
    variants, source = pd.DataFrame(), None
    if ann_path.exists():
        variants = pd.read_csv(ann_path, sep="\t")
        source = "annotated"
    elif var_path.exists():
        variants = pd.read_csv(var_path, sep="\t")
        source = "raw"

    return {"case_id": case_id, "phenotype": phenotype, "variants": variants, 
            "variant_source": source, "path": str(actual_dir)}

def patient_hpo_terms(phenotype: dict) -> list:
    return [f["type"]["id"] for f in phenotype.get("phenotypicFeatures", [])
            if not f.get("excluded", False) and f.get("type", {}).get("id")]

def filter_rare(df: pd.DataFrame, max_af: float = 0.01) -> pd.DataFrame:
    if "population_AF" not in df.columns: return df
    af = pd.to_numeric(df["population_AF"], errors="coerce")
    return df[(af.isna()) | (af < max_af)].copy()

# --- Вердикты (разделение на источники) ---
def clinvar_verdict(sig) -> dict:
    s = str(sig).strip().lower()
    base = {"source": "ClinVar (База данных)", "kind": "confirmed"}
    if s in ("na", "nan", "", "none"): return {**base, "score": 0.1, "verdict": "Нет данных", "detail": "Не подтверждено базой."}
    if "pathogenic" in s and "benign" not in s: return {**base, "score": 1.0, "verdict": "Патогенный", "detail": "Подтверждено (PVS1/PS1)."}
    if "uncertain" in s: return {**base, "score": 0.5, "verdict": "VUS", "detail": "Неопределенное значение."}
    if "benign" in s: return {**base, "score": 0.0, "verdict": "Доброкачественный", "detail": "Отсеивается как ложная находка (BA1)."}
    return {**base, "score": 0.2, "verdict": "Другое", "detail": "Требует проверки."}

def impact_verdict(cons: str) -> dict:
    c = str(cons).lower()
    base = {"source": "VEP (Вычислительный прогноз)", "kind": "predicted"}
    if any(k in c for k in ["frameshift", "stop", "splice"]): return {**base, "score": 1.0, "verdict": "LoF (Потеря функции)", "detail": "Высокое воздействие."}
    if "missense" in c: return {**base, "score": 0.6, "verdict": "Missense", "detail": "Замена аминокислоты."}
    return {**base, "score": 0.2, "verdict": "Слабое", "detail": "Низкое воздействие."}

_HPO_READY = False
def phenotype_verdict(patient_terms: list, gene: str) -> dict:
    global _HPO_READY
    base = {"source": "HPO (Вычислительный прогноз)", "kind": "predicted"}
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
        v = "Сильное" if sim >= 0.6 else "Умеренное" if sim >= 0.3 else "Слабое"
        return {**base, "score": sim, "verdict": v, "detail": f"Сходство {sim:.2f}."}
    except Exception as e:
        return {**base, "score": 0.0, "verdict": "Ошибка HPO", "detail": str(e)}

def check_inheritance(genotype: str, mode: str) -> float:
    if not mode or mode == "Не учитывать": return 1.0
    g = str(genotype).replace(" ", "")
    if mode == "Доминантный (AD)": return 1.0 if g in ("0/1", "1/0") else 0.1
    if mode == "Рецессивный (AR)": return 1.0 if g == "1/1" else 0.1
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
        inh = check_inheritance(r.get("genotype", ""), inheritance_mode)

        total = (WEIGHTS["phenotype"] * ph["score"] + WEIGHTS["clinvar"] * cv["score"] + 
                 WEIGHTS["impact"] * imp["score"] + WEIGHTS["rare"] * 1.0) * inh

        rows.append({
            "gene": gene, "variant": f"{r.get('chrom')}:{r.get('pos')}:{r.get('ref')}>{r.get('alt')}",
            "priority_score": round(total, 3), "phenotype_score": round(ph["score"], 3), "clinvar_score": round(cv["score"], 3),
            "phenotype_verdict": ph, "clinvar_verdict": cv, "impact_verdict": imp, "inheritance_match": inh > 0.5
        })

    out = pd.DataFrame(rows)
    if out.empty: return out
    out = out.sort_values("priority_score", ascending=False).groupby("gene", as_index=False).first()
    out["confidence"] = out["priority_score"].apply(lambda x: "Высокая" if x >= 0.65 else "Средняя" if x >= 0.45 else "Низкая")
    return out.sort_values("priority_score", ascending=False).reset_index(drop=True)

def no_candidate_reason(result: pd.DataFrame) -> str:
    if result.empty: return "Нет аннотированных вариантов (в файлах _annotated.tsv отсутствует колонка gene)."
    return "Убедительный кандидат не найден. Все признаки слабые (Сценарий 3)."
