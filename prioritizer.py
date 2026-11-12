
import json
from pathlib import Path
import pandas as pd


def load_case(case_id: str, data_dir: Path) -> dict:
    meta_path = data_dir / f"{case_id}_metadata.json"
    pheno_path = data_dir / f"{case_id}_phenotype.json"
    ann_path = data_dir / f"{case_id}_annotated.tsv"
    var_path = data_dir / f"{case_id}_variants.tsv"

    metadata, phenotype = {}, {}
    if meta_path.exists():
        metadata = json.loads(meta_path.read_text())
    if pheno_path.exists():
        phenotype = json.loads(pheno_path.read_text())

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


def excluded_hpo_terms(phenotype: dict) -> list:
    return [f["type"]["id"] for f in phenotype.get("phenotypicFeatures", [])
            if f.get("excluded", False) and f.get("type", {}).get("id")]


def filter_rare(df: pd.DataFrame, max_af: float = 0.01) -> pd.DataFrame:
    if "population_AF" not in df.columns:
        return df
    af = pd.to_numeric(df["population_AF"], errors="coerce")
    return df[(af.isna()) | (af < max_af)].copy()


def clinvar_verdict(sig) -> dict:
    s = str(sig).strip()
    low = s.lower()
    base = {"source": "ClinVar", "link": "https://www.ncbi.nlm.nih.gov/clinvar/"}
    
    if low in ("na", "nan", "", "none", "not provided"):
        return {**base, "score": 0.1, "verdict": "Данных нет", "kind": "missing",
                "detail": "Клиническая значимость не подтверждена базой."}
    if "pathogenic" in low and "conflicting" not in low and "benign" not in low:
        return {**base, "score": 1.0, "verdict": "Патогенный / вероятно патогенный", "kind": "confirmed",
                "detail": f"«{s}». Подтверждено базой (ACMG: PVS1/PS1)."}
    if "conflicting" in low:
        return {**base, "score": 0.4, "verdict": "Противоречивые данные", "kind": "conflicting",
                "detail": f"«{s}». Лаборатории расходятся."}
    if "uncertain" in low:
        return {**base, "score": 0.5, "verdict": "VUS (неопределённое значение)", "kind": "uncertain",
                "detail": f"«{s}». Данных недостаточно."}
    if "benign" in low:
        return {**base, "score": 0.0, "verdict": "Доброкачественный", "kind": "benign",
                "detail": f"«{s}». Отсеивается как ложная находка (BA1/BS1)."}
    return {**base, "score": 0.2, "verdict": f"Другая значимость: {s}", "kind": "other",
            "detail": "Требует ручной интерпретации."}


def impact_verdict(consequence: str) -> dict:
    c = str(consequence).lower()
    base = {"source": "VEP (прогноз)", "link": "https://www.ensembl.org/info/docs/tools/vep/"}
    
    if any(k in c for k in ["frameshift", "stop_gained", "splice", "nonsense"]):
        return {**base, "score": 1.0, "verdict": "Высокое воздействие (loss-of-function)", "kind": "predicted",
                "detail": f"Тип: {consequence}. ACMG: PVS1."}
    if "missense" in c or "inframe" in c:
        return {**base, "score": 0.6, "verdict": "Умеренное воздействие (missense)", "kind": "predicted",
                "detail": f"Тип: {consequence}. ACMG: PM5/PP3."}
    if "synonymous" in c or "intron" in c:
        return {**base, "score": 0.2, "verdict": "Низкое воздействие", "kind": "predicted",
                "detail": f"Тип: {consequence}."}
    return {**base, "score": 0.3, "verdict": "Не классифицировано", "kind": "unknown",
            "detail": f"Тип: {consequence or 'NA'}."}


_HPO_READY = False
_gene_cache = {}


def _init_hpo():
    global _HPO_READY
    if not _HPO_READY:
        from pyhpo import Ontology
        _ = Ontology()
        _HPO_READY = True


def gene_hpo_set(gene: str):
    if gene in _gene_cache:
        return _gene_cache[gene]
    _init_hpo()
    from pyhpo import Ontology
    result = None
    for g in Ontology.genes:
        if g.name == gene:
            result = g.hpo_set()
            break
    _gene_cache[gene] = result
    return result


def phenotype_verdict(patient_terms: list, gene: str) -> dict:
    base = {"source": "HPO (pyhpo)", "link": "https://hpo.jax.org/"}
    if not patient_terms:
        return {**base, "score": 0.0, "verdict": "Фенотип не указан", "kind": "missing",
                "detail": "Нет HPO-терминов."}

    gs = gene_hpo_set(gene)
    if gs is None or len(gs) == 0:
        return {**base, "score": 0.0, "verdict": "Связь гена с HPO неизвестна", "kind": "missing",
                "detail": f"Нет аннотаций для {gene}."}

    from pyhpo import HPOSet
    try:
        ps = HPOSet.from_queries(patient_terms)
        sim = float(ps.similarity(gs))
    except Exception:
        sim = 0.0

    if sim >= 0.6:
        v = "Сильное фенотипическое соответствие"
    elif sim >= 0.3:
        v = "Умеренное фенотипическое соответствие"
    else:
        v = "Слабое фенотипическое соответствие"

    return {**base, "score": sim, "verdict": v, "kind": "computed",
            "detail": f"Сходство {sim:.2f} ({len(patient_terms)} терминов пациента vs {len(gs)} гена)."}


def rarity_verdict(af_value, max_af: float) -> dict:
    base = {"source": "gnomAD", "link": "https://gnomad.broadinstitute.org/"}
    try:
        af = float(af_value)
    except (TypeError, ValueError):
        return {**base, "score": 1.0, "verdict": "Частота неизвестна", "kind": "missing",
                "detail": "Считаем редким (AF отсутствует)."}

    if af < 0.001:
        return {**base, "score": 1.0, "verdict": "Очень редкий", "kind": "confirmed",
                "detail": f"AF = {af:.5f} (< 0.1%)."}
    if af < max_af:
        return {**base, "score": 0.7, "verdict": "Редкий", "kind": "confirmed",
                "detail": f"AF = {af:.5f} (< {max_af})."}
    return {**base, "score": 0.0, "verdict": "Частый", "kind": "confirmed",
            "detail": f"AF = {af:.5f} (>= {max_af})."}


DEFAULT_WEIGHTS = {"phenotype": 0.40, "clinvar": 0.30, "impact": 0.20, "rare": 0.10}


def score_case(case: dict, weights: dict = None, max_af: float = 0.01) -> pd.DataFrame:
    weights = weights or DEFAULT_WEIGHTS
    df = case["variants"].copy()
    if df.empty:
        return pd.DataFrame()

    patient_terms = patient_hpo_terms(case["phenotype"])
    df = filter_rare(df, max_af)

    rows = []
    for _, r in df.iterrows():
        gene = str(r.get("gene", "NA")).strip()
        if gene in ("NA", "nan", ""):
            continue

        cv = clinvar_verdict(r.get("clinical_significance", "NA"))
        imp = impact_verdict(r.get("consequence", "NA"))
        ph = phenotype_verdict(patient_terms, gene)
        ra = rarity_verdict(r.get("population_AF"), max_af)

        total = (weights["phenotype"] * ph["score"] +
                 weights["clinvar"] * cv["score"] +
                 weights["impact"] * imp["score"] +
                 weights["rare"] * ra["score"])

        rows.append({
            "gene": gene,
            "variant": f"{r.get('chrom')}:{r.get('pos')}:{r.get('ref')}>{r.get('alt')}",
            "genotype": r.get("genotype", "NA"),
            "consequence": r.get("consequence", "NA"),
            "clinical_significance": r.get("clinical_significance", "NA"),
            "population_AF": r.get("population_AF", "NA"),
            "priority_score": round(total, 3),
            "phenotype_score": round(ph["score"], 3),
            "clinvar_score": round(cv["score"], 3),
            "impact_score": round(imp["score"], 3),
            "rare_score": round(ra["score"], 3),
            "phenotype_verdict": ph,
            "clinvar_verdict": cv,
            "impact_verdict": imp,
            "rare_verdict": ra,
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out = (out.sort_values("priority_score", ascending=False)
              .groupby("gene", as_index=False).first())

    def get_conf(score, cv_kind):
        if score >= 0.65 and cv_kind == "confirmed": return "Высокая"
        if score >= 0.45: return "Средняя"
        return "Низкая"

    out["confidence"] = out.apply(lambda x: get_conf(x["priority_score"], x["clinvar_verdict"]["kind"]), axis=1)
    return out.sort_values("priority_score", ascending=False).reset_index(drop=True)


def no_candidate_reason(result: pd.DataFrame) -> str:
    if result.empty:
        return "Не осталось аннотированных вариантов после фильтрации."
    top = result.iloc[0]
    problems = []
    if top["phenotype_score"] < 0.3: problems.append("низкое фенотипическое соответствие")
    if top["clinvar_verdict"]["kind"] == "missing": problems.append("нет данных ClinVar")
    if top["impact_score"] < 0.3: problems.append("слабое воздействие")
    if not problems: problems.append("все признаки слабые")
    return "Убедительный кандидат не найден. Причины: " + ", ".join(problems) + "."
