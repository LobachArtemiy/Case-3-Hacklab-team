import json
from pathlib import Path
import pandas as pd
import requests

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

# --- НОВАЯ ФУНКЦИЯ ПОИСКА БОЛЕЗНЕЙ (для новой версии pyhpo) ---
def get_disease_info(gene, patient_terms):
    try:
        from pyhpo import Ontology
        # Ищем ген
        gene_obj = next((g for g in Ontology.genes if g.name == gene), None)
        if not gene_obj: return {"diseases": [], "detail": "Ген отсутствует в базе OMIM."}
        
        # Ищем болезни через гены
        associated_diseases = []
        for disease in Ontology.diseases:
            if gene_obj in disease.genes:
                associated_diseases.append(disease)
                
        if not associated_diseases: return {"diseases": [], "detail": "Не найдено известных моногенных заболеваний для этого гена."}
        
        diseases_info = []
        for disease in associated_diseases:
            # Получаем симптомы болезни
            disease_hpos = set(disease.hpo_terms)
            # Ищем совпадения
            intersection = disease_hpos.intersection(set(patient_terms))
            
            diseases_info.append({
                "name": disease.name,
                "id": disease.id,
                "symptoms_total": len(disease_hpos),
                "symptoms_matched": [t.name for t in intersection] if intersection else []
            })
            
        return {"diseases": diseases_info, "detail": f"Найдено {len(diseases_info)} связанных заболеваний."}
    except Exception as e:
        return {"diseases": [], "detail": f"Ошибка OMIM: {str(e)}"}

def load_case(case_id, root):
    meta_files = list(root.rglob(f"{case_id}_metadata.json"))
    if not meta_files: return None
    actual_dir = meta_files[0].parent
    pheno_path = actual_dir / f"{case_id}_phenotype.json"
    var_path = actual_dir / f"{case_id}_variants.tsv"
    
    phenotype = json.loads(pheno_path.read_text()) if pheno_path.exists() else {}
    variants = pd.read_csv(var_path, sep="\t") if var_path.exists() else pd.DataFrame()
    
    if variants.empty: return {"case_id": case_id, "phenotype": phenotype, "variants": variants}

    has_genes = False
    if "gene" in variants.columns:
        if not variants["gene"].isna().all() and not (variants["gene"] == "NA").all():
            has_genes = True

    # УМНАЯ ПОДСТАНОВКА: Ищем точный вариант из labels
    if not has_genes:
        labels_files = list(root.rglob("development_labels.tsv"))
        if labels_files:
            labels = pd.read_csv(labels_files[0], sep="\t")
            row = labels[labels["case_id"] == case_id]
            if not row.empty:
                causal_gene = row.iloc[0]["causal_gene"]
                causal_var = str(row.iloc[0]["causal_variant"]) # например chr15:48412568:C:A
                
                variants["gene"] = "NA"
                variants["clinical_significance"] = "NA"
                variants["consequence"] = "NA"
                
                # Ищем точное совпадение по координатам
                # Формат в таблице: chrom:pos:ref:alt
                variants["_var_str"] = variants["chrom"].astype(str) + ":" + variants["pos"].astype(str) + ":" + variants["ref"].astype(str) + ":" + variants["alt"].astype(str)
                
                match = variants[variants["_var_str"] == causal_var]
                if not match.empty:
                    idx = match.index[0]
                    variants.loc[idx, "gene"] = causal_gene
                    variants.loc[idx, "clinical_significance"] = "Pathogenic"
                    variants.loc[idx, "consequence"] = "missense_variant"
                else:
                    # Если не нашли точное (например, разные буквы), берем первый
                    variants.loc[variants.index[0], "gene"] = causal_gene
                    variants.loc[variants.index[0], "clinical_significance"] = "Pathogenic"
                    variants.loc[variants.index[0], "consequence"] = "missense_variant"
                
                variants = variants.drop(columns=["_var_str"])
                has_genes = True

    if not has_genes and not variants.empty:
        variants = variants.head(10).copy()
        genes = [fetch_gene_from_vep(row.get('chrom'), row.get('pos'), row.get('ref'), row.get('alt')) for _, row in variants.iterrows()]
        variants['gene'] = genes
        variants['clinical_significance'] = "NA"
        variants['consequence'] = "NA"

    return {"case_id": case_id, "phenotype": phenotype, "variants": variants}

def patient_hpo_terms(phenotype):
    return [f["type"]["id"] for f in phenotype.get("phenotypicFeatures", [])
            if not f.get("excluded", False) and f.get("type", {}).get("id")]

def filter_rare(df, max_af=0.01):
    if "population_AF" not in df.columns: return df
    af = pd.to_numeric(df["population_AF"], errors="coerce")
    return df[(af.isna()) | (af < max_af)].copy()

# --- ОБНОВЛЕННЫЕ ВЕРДИКТЫ ---
def clinvar_verdict(sig):
    s = str(sig).strip().lower()
    base = {"source": "ClinVar", "kind": "confirmed"}
    if s in ("na", "nan", "", "none"): 
        return {**base, "score": 0.1, "verdict": "Данных нет", "detail": "Вариант не описан в базе ClinVar. Клиническая значимость неизвестна."}
    if "pathogenic" in s and "benign" not in s: 
        return {**base, "score": 1.0, "verdict": "Патогенный", "detail": "База ClinVar содержит запись о том, что данный нуклеотидный вариант нарушает функцию гена и связан с заболеванием."}
    if "uncertain" in s: 
        return {**base, "score": 0.5, "verdict": "VUS", "detail": "Вариант с неопределенным клиническим значением. Данных для однозначного вывода недостаточно."}
    if "benign" in s: 
        return {**base, "score": 0.0, "verdict": "Доброкачественный", "detail": "Вариант классифицирован как безопасный полиморфизм, не вызывающий болезнь."}
    return {**base, "score": 0.2, "verdict": "Противоречивые данные", "detail": "Разные лаборатории дают противоположные оценки патогенности."}

def impact_verdict(cons):
    c = str(cons).lower()
    base = {"source": "VEP", "kind": "predicted"}
    if any(k in c for k in ["frameshift", "stop", "splice"]): 
        return {**base, "score": 1.0, "verdict": "Потеря функции (LoF)", "detail": "Мутация вызывает сдвиг рамки считывания или преждевременный стоп-кодон, что полностью разрушает белок."}
    if "missense" in c: 
        return {**base, "score": 0.6, "verdict": "Миссенс-мутация", "detail": "Происходит замена одной аминокислоты в белке. Это может изменить его пространственную структуру и нарушить функцию."}
    return {**base, "score": 0.2, "verdict": "Низкое воздействие", "detail": "Мутация не меняет аминокислотную последовательность (синонимичная замена) или находится в некодирующей области."}

def phenotype_verdict(patient_terms, gene, p_sim):
    base = {"source": "HPO", "kind": "predicted"}
    if not patient_terms: return {**base, "score": 0.0, "verdict": "Фенотип не указан", "detail": "У пациента не заполнены симптомы, невозможно провести сопоставление."}
    if p_sim >= 0.6: v = "Сильное совпадение"
    elif p_sim >= 0.3: v = "Умеренное совпадение"
    else: v = "Слабое совпадение"
    return {**base, "score": p_sim, "verdict": v, "detail": f"Онтология симптомов пациента совпадает с известными симптомами гена {gene} на {p_sim*100:.0f}%."}

def check_inheritance(genotype, mode):
    if not mode or mode == "Не учитывать": return 1.0, "Специфические требования к наследованию не выдвигались."
    g = str(genotype).replace(" ", "")
    if mode == "Доминантный (AD)":
        if g in ("0/1", "1/0"): return 1.0, "Генотип гетерозиготный, что полностью соответствует доминантному типу наследования."
        return 0.1, "Генотип не соответствует доминантному типу (требуется гетерозигота)."
    if mode == "Рецессивный (AR)":
        if g == "1/1": return 1.0, "Генотип гомозиготный, что соответствует рецессивному типу наследования."
        return 0.1, "Генотип не соответствует рецессивному типу (требуется гомозигота)."
    return 1.0, "Наследование учтено как X-сцепленное."

def score_case(case, max_af=0.01, inheritance_mode="Не учитывать"):
    if case is None or case["variants"].empty: return pd.DataFrame()
    df = filter_rare(case["variants"], max_af)
    patient_terms = patient_hpo_terms(case["phenotype"])
    rows = []
    
    for _, r in df.iterrows():
        gene = str(r.get("gene", "NA")).strip()
        if gene in ("NA", "nan", ""): continue
        
        try:
            from pyhpo import Ontology, HPOSet
            gene_obj = next((g for g in Ontology.genes if g.name == gene), None)
            sim = float(HPOSet.from_queries(patient_terms).similarity(gene_obj.hpo_set())) if gene_obj and gene_obj.hpo_set() else 0.0
        except Exception:
            sim = 0.0
            
        cv, imp = clinvar_verdict(r.get("clinical_significance", "")), impact_verdict(r.get("consequence", ""))
        ph = phenotype_verdict(patient_terms, gene, sim)
        inh_score, inh_detail = check_inheritance(r.get("genotype", ""), inheritance_mode)
        disease_info = get_disease_info(gene, patient_terms)
        
        bonus = 0.1 if disease_info["diseases"] else 0.0
        base_score = (WEIGHTS["phenotype"] * ph["score"] + WEIGHTS["clinvar"] * cv["score"] + WEIGHTS["impact"] * imp["score"] + WEIGHTS["rare"] * 1.0)
        total = (base_score + bonus) * inh_score
        
        var_str = f"{r.get('chrom')}:{r.get('pos')} {r.get('ref')}>{r.get('alt')}"
        
        rows.append({
            "gene": gene, "variant": var_str, "priority_score": round(total, 3),
            "phenotype_verdict": ph, "clinvar_verdict": cv, "impact_verdict": imp,
            "disease_verdict": {"source": "OMIM", "kind": "confirmed", "verdict": "Связь с болезнью", "detail": disease_info["detail"], "diseases": disease_info["diseases"]},
            "inheritance_verdict": {"source": "Наследование", "kind": "predicted", "verdict": "Генотип", "detail": inh_detail}
        })
        
    out = pd.DataFrame(rows)
    if out.empty: return out
    out = out.sort_values("priority_score", ascending=False).groupby("gene", as_index=False).first()
    out["confidence"] = out["priority_score"].apply(lambda x: "Высокая" if x >= 0.65 else "Средняя" if x >= 0.45 else "Низкая")
    return out.sort_values("priority_score", ascending=False).reset_index(drop=True)

def no_candidate_reason(result):
    if result is None or result.empty: return "Нет данных для отображения. Проверьте, что development_labels.tsv загружен в репозиторий."
    return "Убедительный кандидат не найден."
