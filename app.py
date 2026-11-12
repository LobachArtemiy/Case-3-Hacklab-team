import streamlit as st
import pandas as pd
from pathlib import Path
from prioritizer import load_case, score_case, no_candidate_reason, filter_rare

st.set_page_config(page_title="Variant Prioritizer", layout="wide")

# Минималистичная тема
st.markdown("""
<style>
html, body, [data-testid="stAppViewContainer"] { background-color: #F8F9FA !important; font-family: 'Inter', sans-serif; }
[data-testid="stSidebar"] { background-color: #FFFFFF !important; border-right: 1px solid #E2E8F0; }
.verdict-card { background: #FFF; border-left: 4px solid #3182CE; border-radius: 6px; padding: 1rem; margin: 0.5rem 0; border: 1px solid #E2E8F0; }
.confirmed { border-left-color: #38A169; }
.predicted { border-left-color: #3182CE; }
.badge { display: inline-block; background: #EDF2F7; border-radius: 4px; padding: 2px 8px; font-size: 0.7rem; color: #4A5568; margin-left: 6px; }
</style>
""", unsafe_allow_html=True)

# Поиск данных
PROJECT_ROOT = Path(__file__).parent
all_meta_files = list(PROJECT_ROOT.rglob("*_metadata.json"))
case_ids = sorted({p.name.replace("_metadata.json", "") for p in all_meta_files})

with st.sidebar:
    st.header("Параметры")
    if not case_ids:
        st.error("Файлы данных не найдены. Загрузите папку data.")
        st.stop()
    selected = st.selectbox("Клинический случай", case_ids)
    max_af = st.slider("Порог частоты (gnomAD)", 0.0, 0.05, 0.01, 0.001)
    
    st.divider()
    st.subheader("Дополнительно")
    inh_mode = st.selectbox("Тип наследования", ["Не учитывать", "Доминантный (AD)", "Рецессивный (AR)"])
    new_hpo = st.text_input("Добавить HPO", placeholder="HP:0002616, HP:0001643")

# Загрузка
with st.spinner("Анализ вариантов (при необходимости — запрос к VEP)..."):
    case = load_case(selected, PROJECT_ROOT)

if case is None:
    st.error(f"Файлы для {selected} не найдены.")
    st.stop()

st.title(f"Кейс: {selected}")
st.caption(f"Путь: {case['path']} | Вариантов загружено: {len(case['variants'])}")

# Доп. задача
case_mod = dict(case)
if new_hpo.strip():
    new_feats = [{"type": {"id": t.strip(), "label": ""}} for t in new_hpo.split(",") if t.strip()]
    case_mod["phenotype"] = dict(case["phenotype"])
    case_mod["phenotype"]["phenotypicFeatures"] = list(case["phenotype"].get("phenotypicFeatures", [])) + new_feats

# Скоринг
result = score_case(case_mod, max_af, inheritance_mode=inh_mode)

# Вывод
if result.empty:
    st.warning(no_candidate_reason(result))
else:
    st.subheader("Ранжированный список генов")
    st.dataframe(result[["gene", "priority_score", "confidence", "phenotype_score", "clinvar_score"]], use_container_width=True)
    
    top = result.iloc[0]
    st.subheader(f"Объяснение для {top['gene']}")
    
    def render(v, title):
        st.markdown(f"<div class='verdict-card {v['kind']}'><b>{title}</b>: {v['verdict']} <span class='badge'>{v['source']}</span><br><span style='color:#4A5568; font-size:0.9em;'>{v['detail']}</span></div>", unsafe_allow_html=True)

    render(top["clinvar_verdict"], "База ClinVar")
    render(top["phenotype_verdict"], "Сопоставление HPO")
    render(top["impact_verdict"], "Функциональный эффект")
    
    st.divider()
    st.download_button("📄 Скачать отчёт (CSV)", result.to_csv(index=False).encode("utf-8"), f"{selected}_report.csv", "text/csv")
    st.caption("⚠️ Исследовательский результат. Не является клиническим диагнозом.")
