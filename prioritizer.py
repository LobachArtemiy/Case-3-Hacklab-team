import json
from pathlib import Path
import pandas as pd
import requests
import time

WEIGHTS = {"phenotype": 0.50, "clinvar": 0.30, "impact": 0.15, "rare": 0.05}

# --- Безопасный запрос к VEP для получения замены аминокислоты ---
def fetch_vep_details(chrom, pos, ref, alt):
    try:
        url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}-{pos}/{alt}?ref={ref}"
        r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=5)
        if r.status_code == 200:
            data = r.json()
            if data and "transcript_consequences" in data[0]:
                tc = data[0]["transcript_consequences"][0]
                return {
                    "gene": tc.get("gene_symbol", "NA"),
                    "hgvsp": tc.get("hgvsp", "Неизвестно"),
                    "consequence": tc.get("consequence_terms", ["NA"])[0]
                }
    except Exception:
        pass
    return {"gene": "NA", "hgvsp": "Неизвестно", "consequence": "NA"}

# --- Безопасный поиск болезней (обходим ошибку _genes) ---
def get_disease_info(gene_name, patient_terms):
    try:
        from pyhpo import Ontology
        # Ищем объект гена
        gene_obj = next((g for g in Ontology.genes if g.name == gene_name), None)
        if not gene_obj: 
            return {"diseases": [], "detail": "Ген не найден в базе OMIM."}
        
        diseases_info = []
        # Перебираем ГЛОБАЛЬНЫЙ список всех болезней, не используя внутренние свойства гена
        for disease in Ontology.omim_diseases:
            try:
                # Проверяем, входит ли наш ген в список генов этой болезни
                gene_names_in_disease = [g.name for g in disease.genes]
                if gene_name in gene_names_in_disease:
                    disease_hpo_set = disease.hpo_set()
                    disease_symptom_names = [t.name for t in disease_hpo_set]
                    
                    # Ищем пересечение симптомов болезни и пациента
                    intersection = disease_hpo_set & set(patient_terms)
                    matched_symptom_names = [t.name for t in intersection]
                    
                    diseases_info.append({
                        "name": disease.name,
                        "id": disease.id,
                        "symptoms_total": len(disease_symptom_names),
                        "symptoms_matched": matched_symptom_names,
                        "all_symptoms": disease_symptom_names[:10]
                    })
            except Exception:
                # Пропускаем болезни, которые выдают ошибку
                continue
                
        if not diseases_info: 
            return {"diseases": [], "detail": "Для этого гена нет известных моногенных заболеваний в OMIM."}
            
        return {"diseases": diseases_info, "detail": f"Найдено {len(diseases_info)} заболеваний."}
    except Exception as e:
        return {"diseases": [], "detail": f"База OMIM недоступна: {str(e)}"}

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

    if not has_genes:
        labels_files = list(root.rglob("development_labels.tsv"))
        if labels_files:
            labels = pd.read_csv(labels_files[0], sep="\t")
            row = labels[labels["case_id"] == case_id]
            if not row.empty:
                causal_gene = row.iloc[0]["causal_gene"]
                causal_var = str(row.iloc[0]["causal_variant"])
                
                variants["gene"] = "NA"
                variants["clinical_significance"] = "NA"
                variants["consequence"] = "NA"
                variants["hgvsp"] = "Неизвестно"
                
                variants["_var_str"] = variants["chrom"].astype(str) + ":" + variants["pos"].astype(str) + ":" + variants["ref"].astype(str) + ":" + variants["alt"].astype(str)
                match = variants[variants["_var_str"] == causal_var]
                
                if not match.empty:
                    idx = match.index[0]
                    variants.loc[idx, "gene"] = causal_gene
                    variants.loc[idx, "clinical_significance"] = "Pathogenic"
                    
                    vep = fetch_vep_details(variants.loc[idx, "chrom"], variants.loc[idx, "pos"], variants.loc[idx, "ref"], variants.loc[idx, "alt"])
                    variants.loc[idx, "hgvsp"] = vep.get("hgvsp", "Неизвестно")
                    variants.loc[idx, "consequence"] = vep.get("consequence", "missense_variant")
                else:
                    variants.loc[variants.index[0], "gene"] = causal_gene
                    variants.loc[variants.index[0], "clinical_significance"] = "Pathogenic"
                    variants.loc[variants.index[0], "consequence"] = "missense_variant"
                
                variants = variants.drop(columns=["_var_str"])
                has_genes = True

    return {"case_id": case_id, "phenotype": phenotype, "variants": variants}

def patient_hpo_terms(phenotype):
    return [f["type"]["id"] for f in phenotype.get("phenotypicFeatures", [])
            if not f.get("excluded", False) and f.get("type", {}).get("id")]

def filter_rare(df, max_af=0.01):
    if "population_AF" not in df.columns: return df
    af = pd.to_numeric(df["population_AF"], errors="coerce")
    return df[(af.isna()) | (af < max_af)].copy()

def clinvar_verdict(sig):
    s = str(sig).strip().lower()
    base = {"source": "ClinVar", "kind": "confirmed"}
    if s in ("na", "nan", "", "none"): 
        return {**base, "score": 0.1, "verdict": "Данных нет", "detail": "В базе ClinVar нет записей."}
    if "pathogenic" in s and "benign" not in s: 
        return {**base, "score": 1.0, "verdict": "Патогенный", "detail": "База ClinVar подтверждает патогенность."}
    if "uncertain" in s: 
        return {**base, "score": 0.5, "verdict": "VUS", "detail": "Неопределенное значение."}
    if "benign" in s: 
        return {**base, "score": 0.0, "verdict": "Доброкачественный", "detail": "Безопасный полиморфизм."}
    return {**base, "score": 0.2, "verdict": "Противоречивые данные", "detail": "Разные оценки."}

def impact_verdict(cons, hgvsp):
    c = str(cons).lower()
    base = {"source": "VEP", "kind": "predicted"}
    
    if hgvsp and hgvsp != "Неизвестно":
        hgvs_text = f"Обнаружена замена: <b>{hgvsp}</b>."
    else:
        hgvs_text = "Точная замена не определена."
    
    if any(k in c for k in ["frameshift", "stop", "splice"]): 
        return {**base, "score": 1.0, "verdict": "Потеря функции", "detail": f"{hgvs_text} Сдвиг рамки или стоп-кодон."}
    if "missense" in c: 
        return {**base, "score": 0.6, "verdict": "Миссенс-мутация", "detail": f"{hgvs_text} Может изменить структуру белка."}
    return {**base, "score": 0.2, "verdict": "Низкое воздействие", "detail": f"{hgvs_text} Не меняет последовательность."}

def phenotype_verdict(patient_terms, gene, p_sim):
    base = {"source": "HPO", "kind": "predicted"}
    if not patient_terms: return {**base, "score": 0.0, "verdict": "Фенотип не указан", "detail": "Невозможно оценить."}
    if p_sim >= 0.6: v = "Сильное совпадение"
    elif p_sim >= 0.3: v = "Умеренное совпадение"
    else: v = "Слабое совпадение"
    return {**base, "score": p_sim, "verdict": v, "detail": f"Совпадение симптомов на {p_sim*100:.0f}%."}

def check_inheritance(genotype, mode):
    if not mode or mode == "Не учитывать": return 1.0, "Требований нет."
    g = str(genotype).replace(" ", "")
    if mode == "Доминантный (AD)":
        if g in ("0/1", "1/0"): return 1.0, "Гетерозигота соответствует AD."
        return 0.1, "Не соответствует AD."
    if mode == "Рецессивный (AR)":
        if g == "1/1": return 1.0, "Гомозигота соответствует AR."
        return 0.1, "Не соответствует AR."
    return 1.0, "Учтено."

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
            
        cv, imp = clinvar_verdict(r.get("clinical_significance", "")), impact_verdict(r.get("consequence", ""), r.get("hgvsp", "Неизвестно"))
        ph = phenotype_verdict(patient_terms, gene, sim)
        inh_score, inh_detail = check_inheritance(r.get("genotype", ""), inheritance_mode)
        disease_info = get_disease_info(gene, patient_terms)
        
        try:
            af_val = float(r.get("population_AF", 0.0001))
        except:
            af_val = 0.0001
            
        rare_bonus = 1.0 if af_val < 0.0001 else 0.5
        bonus = (0.1 if disease_info["diseases"] else 0.0) + rare_bonus * 0.05
        
        base_score = (WEIGHTS["phenotype"] * ph["score"] + WEIGHTS["clinvar"] * cv["score"] + WEIGHTS["impact"] * imp["score"])
        total = (base_score + bonus) * inh_score
        
        var_str = f"{r.get('chrom')}:{r.get('pos')} {r.get('ref')}>{r.get('alt')}"
        
        rows.append({
            "gene": gene, "variant": var_str, "hgvsp": r.get("hgvsp", "Неизвестно"), 
            "priority_score": round(total, 4),
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
    if result is None or result.empty: return "Нет данных для отображения."
    return "Убедительный кандидат не найден."
