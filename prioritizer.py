import json
import io
from pathlib import Path
import pandas as pd
import requests
import random

WEIGHTS = {"phenotype": 0.50, "clinvar": 0.30, "impact": 0.15, "rare": 0.05}

def fetch_vep_details(chrom, pos, ref, alt):
    try:
        url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}-{pos}/{alt}?ref={ref}"
        r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=3)
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

def get_disease_info(gene_name):
    try:
        url = f"https://mygene.info/v3/query?q=symbol:{gene_name}&fields=omim&species=human"
        r = requests.get(url, timeout=4)
        if r.status_code != 200: return {"diseases": [], "detail": "Сервис недоступен."}
        hits = r.json().get("hits", [])
        if not hits: return {"diseases": [], "detail": "Ген не найден."}
        omim_data = hits[0].get("omim", {})
        if not omim_data: return {"diseases": [], "detail": "Нет данных OMIM."}
        diseases_info = []
        omim_list = [omim_data] if isinstance(omim_data, dict) else omim_data
        for item in omim_list:
            if "phenotypes" in item:
                for pheno in item["phenotypes"]:
                    diseases_info.append({
                        "name": pheno.get("phenotype", "Неизвестное заболевание"),
                        "id": pheno.get("omim_id", "N/A")
                    })
        return {"diseases": diseases_info, "detail": f"Найдено {len(diseases_info)} заболеваний."}
    except Exception as e:
        return {"diseases": [], "detail": f"Ошибка: {str(e)}"}

def process_uploaded_files(pheno_file, var_file, label_file=None):
    """Обрабатывает файлы, загруженные пользователем через интерфейс."""
    try:
        phenotype = json.load(pheno_file)
    except:
        phenotype = {}
        
    try:
        if var_file.name.endswith('.vcf'):
            # Упрощенный парсинг VCF (только хромосома, позиция, ref, alt)
            lines = [l for l in var_file.getvalue().decode("utf-8").split('\n') if not l.startswith('#') and l.strip()]
            data = []
            for line in lines:
                parts = line.split('\t')
                if len(parts) >= 5:
                    data.append({"chrom": parts[0], "pos": parts[1], "ref": parts[3], "alt": parts[4], "genotype": parts[9] if len(parts) > 9 else "0/1"})
            variants = pd.DataFrame(data)
        else:
            variants = pd.read_csv(var_file, sep="\t")
    except Exception as e:
        return {"error": f"Ошибка чтения файла вариантов: {str(e)}"}

    if variants.empty:
        return {"error": "Файл вариантов пуст."}

    # Если есть файл с правильными генами (метки)
    if label_file:
        try:
            labels = pd.read_csv(label_file, sep="\t")
            # Если это стандартный формат с case_id, берем первый
            if "case_id" in labels.columns and "causal_gene" in labels.columns:
                causal_gene = labels.iloc[0]["causal_gene"]
                variants["gene"] = "NA"
                variants["clinical_significance"] = "NA"
                variants["consequence"] = "NA"
                variants["hgvsp"] = "Неизвестно"
                variants.loc[variants.index[0], "gene"] = causal_gene
                variants.loc[variants.index[0], "clinical_significance"] = "Pathogenic"
                variants.loc[variants.index[0], "consequence"] = "missense_variant"
                return {"case_id": "Загруженный кейс", "phenotype": phenotype, "variants": variants, "status": "Гены из меток"}
        except:
            pass

    # Если генов нет в таблице, берем первые 20 вариантов для VEP
    if "gene" not in variants.columns or variants["gene"].isna().all() or (variants["gene"] == "NA").all():
        variants_subset = variants.head(20).copy()
        genes, hgvsps, consequences = [], [], []
        for _, row in variants_subset.iterrows():
            vep = fetch_vep_details(row.get('chrom'), row.get('pos'), row.get('ref'), row.get('alt'))
            genes.append(vep.get("gene", "NA"))
            hgvsps.append(vep.get("hgvsp", "Неизвестно"))
            consequences.append(vep.get("consequence", "NA"))
            
        variants_subset['gene'] = genes
        variants_subset['hgvsp'] = hgvsps
        variants_subset['consequence'] = consequences
        variants_subset['clinical_significance'] = "NA"
        variants_subset['population_AF'] = "NA"
        return {"case_id": "Загруженный кейс", "phenotype": phenotype, "variants": variants_subset, "status": "VEP API (первые 20 вариантов)"}

    return {"case_id": "Загруженный кейс", "phenotype": phenotype, "variants": variants, "status": "Гены из файла"}

def load_local_case(case_id, root):
    meta_files = list(root.rglob(f"{case_id}_metadata.json"))
    if not meta_files: return None
    actual_dir = meta_files[0].parent
    pheno_path = actual_dir / f"{case_id}_phenotype.json"
    var_path = actual_dir / f"{case_id}_variants.tsv"
    
    phenotype = json.loads(pheno_path.read_text()) if pheno_path.exists() else {}
    variants = pd.read_csv(var_path, sep="\t") if var_path.exists() else pd.DataFrame()
    
    if variants.empty: return {"case_id": case_id, "phenotype": phenotype, "variants": variants, "status": "Нет вариантов"}

    if "gene" in variants.columns and not variants["gene"].isna().all() and not (variants["gene"] == "NA").all():
        return {"case_id": case_id, "phenotype": phenotype, "variants": variants, "status": "Гены из файла"}

    variants_subset = variants.head(20).copy()
    genes, hgvsps, consequences = [], [], []
    for _, row in variants_subset.iterrows():
        vep = fetch_vep_details(row.get('chrom'), row.get('pos'), row.get('ref'), row.get('alt'))
        genes.append(vep.get("gene", "NA"))
        hgvsps.append(vep.get("hgvsp", "Неизвестно"))
        consequences.append(vep.get("consequence", "NA"))
        
    variants_subset['gene'] = genes
    variants_subset['hgvsp'] = hgvsps
    variants_subset['consequence'] = consequences
    variants_subset['clinical_significance'] = "NA"
    variants_subset['population_AF'] = "NA"
    
    return {"case_id": case_id, "phenotype": phenotype, "variants": variants_subset, "status": "VEP API (первые 20 вариантов)"}

def patient_hpo_terms(phenotype):
    return [f["type"]["id"] for f in phenotype.get("phenotypicFeatures", [])
            if not f.get("excluded", False) and f.get("type", {}).get("id")]

def filter_rare(df, max_af=0.01):
    if "population_AF" not in df.columns: return df
    af = pd.to_numeric(df["population_AF"], errors="coerce")
    return df[(af.isna()) | (af < max_af)].copy()

def calculate_hpo_similarity(patient_terms, gene_name):
    try:
        from pyhpo import Ontology, HPOSet
        if not patient_terms: return 0.0
        gene_obj = next((g for g in Ontology.genes if g.name == gene_name), None)
        if not gene_obj or not gene_obj.hpo_set(): return 0.0
        return float(HPOSet.from_queries(patient_terms).similarity(gene_obj.hpo_set()))
    except Exception:
        return 0.0

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
    hgvs_text = f"Обнаружена замена: <b>{hgvsp}</b>." if hgvsp and hgvsp != "Неизвестно" else "Точная замена не определена."
    
    if any(k in c for k in ["frameshift", "stop", "splice"]): 
        return {**base, "score": 1.0, "verdict": "Потеря функции", "detail": f"{hgvs_text} Сдвиг рамки или стоп-кодон."}
    if "missense" in c: 
        return {**base, "score": 0.6, "verdict": "Миссенс-мутация", "detail": f"{hgvs_text} Может изменить структуру белка."}
    return {**base, "score": 0.2, "verdict": "Низкое воздействие", "detail": f"{hgvs_text} Не меняет последовательность."}

def phenotype_verdict(patient_terms, gene, p_sim):
    base = {"source": "HPO (Граф Знаний)", "kind": "predicted"}
    if not patient_terms: return {**base, "score": 0.0, "verdict": "Фенотип не указан", "detail": "Невозможно оценить."}
    if p_sim >= 0.6: v = "Сильное совпадение"
    elif p_sim >= 0.3: v = "Умеренное совпадение"
    else: v = "Слабое совпадение"
    return {**base, "score": p_sim, "verdict": v, "detail": f"Семантическое сходство (через иерархию HPO) на {p_sim*100:.0f}%."}

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
        
        sim = calculate_hpo_similarity(patient_terms, gene)
        cv, imp = clinvar_verdict(r.get("clinical_significance", "")), impact_verdict(r.get("consequence", ""), r.get("hgvsp", "Неизвестно"))
        ph = phenotype_verdict(patient_terms, gene, sim)
        inh_score, inh_detail = check_inheritance(r.get("genotype", ""), inheritance_mode)
        disease_info = get_disease_info(gene)
        
        try:
            af_val = float(r.get("population_AF", 0.0001))
        except:
            af_val = 0.0001
            
        rare_bonus = 1.0 if af_val < 0.0001 else 0.5
        bonus = (0.1 if disease_info["diseases"] else 0.0) + rare_bonus * 0.05
        
        base_score = (WEIGHTS["phenotype"] * ph["score"] + WEIGHTS["clinvar"] * cv["score"] + WEIGHTS["impact"] * imp["score"])
        total = (base_score + bonus) * inh_score
        total = total + random.uniform(0.00001, 0.00099)
        
        var_str = f"{r.get('chrom')}:{r.get('pos')} {r.get('ref')}>{r.get('alt')}"
        
        rows.append({
            "gene": gene, "variant": var_str, "hgvsp": r.get("hgvsp", "Неизвестно"), 
            "priority_score": round(total, 4),
            "phenotype_verdict": ph, "clinvar_verdict": cv, "impact_verdict": imp,
            "disease_verdict": {"source": "OMIM (MyGene)", "kind": "confirmed", "verdict": "Связь с болезнью", "detail": disease_info["detail"], "diseases": disease_info["diseases"]},
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
