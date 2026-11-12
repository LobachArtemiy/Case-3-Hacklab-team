import json
from pathlib import Path
import pandas as pd

# Фиксированные веса (научно обоснованные, не меняются пользователем)
WEIGHTS = {"phenotype": 0.50, "clinvar": 0.30, "impact": 0.15, "rare": 0.05}

def load_case(case_id: str, data_dir: Path) -> dict:
    # Ищем файлы, даже если они лежат в подпапках (development / hidden)
    meta_files = list(data_dir.rglob(f"{case_id}_metadata.json"))
    actual_dir = meta_files[0].parent if meta_files else data_dir

    meta_path = actual_dir / f"{case_id}_metadata.json"
    pheno_path = actual_dir / f"{case_id}_phenotype.json"
    ann_path = actual_dir / f"{case_id}_annotated.tsv"
    var_path = actual_dir / f"{case_id}_variants.tsv"

    metadata, phenotype = {}, {}
    if meta_path.exists(): metadata = json.loads(meta_path.read_text())
    if pheno_path.exists(): phenotype = json.loads(pheno_path.read_text())

    variants, source = pd.DataFrame(), None
    if ann_path.exists():
        variants = pd.read_csv(ann_path, sep="\t")
        source = "annotated"
    elif var_path.exists():
        variants = pd.read_csv(var_path, sep="\t")
        source = "raw"

    return {"case_id": case_id, "metadata": metadata, "phenotype": phenotype,
            "variants": variants, "variant_source": source}

def patient_hpo_terms(phenotype: dict) -> list:
    return [f["type"]["id"] for f in phenotype.get("phenotypicFeatures", [])
            if not f.get("excluded", False) and f.get("type", {}).get("id")]

def filter_rare(df: pd.DataFrame, max_af: float = 0.01) -> pd.DataFrame:
    if "population_AF" not in df.columns: return df
    af = pd.to_numeric(df["population_AF"], errors="coerce")
    return df[(af.isna()) | (af < max_af)].copy()

# --- Вердикты ---
def clinvar_verdict(sig) -> dict:
    s = str(sig).strip().lower()
    base = {"source": "ClinVar", "link": "https://www.ncbi.nlm.nih.gov/clinvar/"}
    if s in ("na", "nan", "", "none"): return {**base, "score": 0.1, "verdict": "Нет данных", "kind": "missing", "detail": "Не подтверждено."}
    if "pathogenic" in s and "benign" not in s: return {**base, "score": 1.0, "verdict": "Патогенный", "kind": "confirmed", "detail": "Подтверждено (PVS1/PS1)."}
    if "uncertain" in s: return {**base, "score": 0.5, "verdict": "VUS", "kind": "uncertain", "detail": "Неопределенное значение."}
    if "benign" in s: return {**base, "score": 0.0, "verdict": "Доброкачественный", "kind": "benign", "detail": "Отсеивается (BA1)."}
    return {**base, "score": 0.2, "verdict": "Другое", "kind": "other", "detail": "Требует проверки."}

def impact_verdict(cons: str) -> dict:
    c = str(cons).lower()
    base = {"source": "VEP", "link": "https://www.ensembl.org/vep"}
    if any(k in c for k in ["frameshift", "stop", "splice"]): return {**base, "score": 1.0, "verdict": "LoF (потеря функции)", "kind": "predicted", "detail": "Высокое воздействие."}
    if "missense" in c: return {**base, "score": 0.6, "verdict": "Missense", "kind": "predicted", "detail": "Замена аминокислоты."}
    return {**base, "score": 0.2, "verdict": "Слабое", "kind": "predicted", "detail": "Низкое воздействие."}

# --- pyhpo ---
_HPO_READY = False
def _init_hpo():
    global _HPO_READY
    if not _HPO_READY:
        from pyhpo import Ontology
        _ = Ontology()
        _HPO_READY = True

def phenotype_verdict(patient_terms: list, gene: str) -> dict:
    base = {"source": "HPO (pyhpo)", "link": "https://hpo.jax.org/"}
    if not patient_terms: return {**base, "score": 0.0, "verdict": "Нет фенотипа", "kind": "missing", "detail": "Не указан."}
    
    _init_hpo()
    from pyhpo import Ontology, HPOSet
    gene_obj = next((g for g in Ontology.genes if g.name == gene), None)
    if not gene_obj or not gene_obj.hpo_set(): return {**base, "score": 0.0, "verdict": "Связь неизвестна", "kind": "missing", "detail": "Нет HPO для гена."}
    
    ps, gs = HPOSet.from_queries(patient_terms), gene_obj.hpo_set()
    sim = float(ps.similarity(gs))
    if sim >= 0.6: v = "Сильное соответствие"
    elif sim >= 0.3: v = "Умеренное соответствие"
    else: v = "Слабое соответствие"
    return {**base, "score": sim, "verdict": v, "kind": "computed", "detail": f"Сходство {sim:.2f}."}

def check_inheritance(genotype: str, inheritance_mode: str) -> float:
    """Проверяет, подходит ли генотип под выбранный тип наследования."""
    if not inheritance_mode or inheritance_mode == "Не учитывать": return 1.0
    g = str(genotype).replace(" ", "")
    
    if inheritance_mode == "Доминантный (AD)":
        return 1.0 if g in ("0/1", "1/0") else 0.1
    if inheritance_mode == "Рецессивный (AR)":
        return 1.0 if g == "1/1" else 0.1
    if inheritance_mode == "X-сцепленный (XL)":
        return 1.0 if g in ("0/1", "1/1", "1/0") else 0.5
    return 1.0

# --- Главная функция скоринга ---
def score_case(case: dict, max_af: float = 0.01, inheritance_mode: str = "Не учитывать") -> pd.DataFrame:
    df = case["variants"].copy()
    if df.empty: return pd.DataFrame()

    patient_terms = patient_hpo_terms(case["phenotype"])
    df = filter_rare(df, max_af)

    rows = []
    for _, r in df.iterrows():
        gene = str(r.get("gene", "NA")).strip()
        if gene in ("NA", "nan", ""): continue

        cv, imp, ph = clinvar_verdict(r.get("clinical_significance", "")), impact_verdict(r.get("consequence", "")), phenotype_verdict(patient_terms, gene)
        inh_score = check_inheritance(r.get("genotype", ""), inheritance_mode)

        # Если генотип не подходит под наследование, резко снижаем общий вес
        total = (WEIGHTS["phenotype"] * ph["score"] + WEIGHTS["clinvar"] * cv["score"] + 
                 WEIGHTS["impact"] * imp["score"] + WEIGHTS["rare"] * 1.0) * inh_score

        rows.append({
            "gene": gene, "variant": f"{r.get('chrom')}:{r.get('pos')}:{r.get('ref')}>{r.get('alt')}",
            "genotype": r.get("genotype"), "priority_score": round(total, 3),
            "phenotype_score": round(ph["score"], 3), "clinvar_score": round(cv["score"], 3),
            "phenotype_verdict": ph, "clinvar_verdict": cv, "impact_verdict": imp,
            "inheritance_match": inh_score > 0.5
        })

    out = pd.DataFrame(rows)
    if out.empty: return out

    # Оставляем лучший вариант на ген
    out = out.sort_values("priority_score", ascending=False).groupby("gene", as_index=False).first()

    def get_conf(row):
        if row["priority_score"] >= 0.65 and row["clinvar_verdict"]["kind"] == "confirmed": return "Высокая"
        if row["priority_score"] >= 0.45: return "Средняя"
        return "Низкая"

    out["confidence"] = out.apply(get_conf, axis=1)
    return out.sort_values("priority_score", ascending=False).reset_index(drop=True)

def no_candidate_reason(result: pd.DataFrame) -> str:
    if result.empty: return "Не осталось аннотированных вариантов."
    top = result.iloc[0]
    problems = []
    if top["phenotype_score"] < 0.3: problems.append("низкое соответствие фенотипу")
    if top["clinvar_verdict"]["kind"] == "missing": problems.append("нет данных ClinVar")
    if not top.get("inheritance_match", True): problems.append("генотип не подходит под наследование")
    if not problems: problems.append("все признаки слабые")
    return "Убедительный кандидат не найден. Причины: " + ", ".join(problems) + "."
