import streamlit as st
import pandas as pd
from pathlib import Path
from prioritizer import load_case, score_case, no_candidate_reason

st.set_page_config(page_title="Variant Prioritizer", layout="wide")

st.markdown("""
<style>
html, body, [data-testid="stAppViewContainer"] { background-color: #F8F9FA !important; font-family: sans-serif; }
[data-testid="stSidebar"] { background-color: #FFFFFF !important; border-right: 1px solid #E2E8F0; }
.verdict-box { background: #FFF; border-left: 4px solid #3182CE; border-radius: 6px; padding: 1rem; margin-bottom: 1rem; border: 1px solid #E2E8F0; box-shadow: 0 1px 3px rgba(0,0,0,0.05); }
.verdict-title { font-size: 1.2rem; font-weight: 700; margin-bottom: 0.5rem; color: #1A202C; }
.verdict-badge { display: inline-block; background: #EDF2F7; border-radius: 4px; padding: 2px 8px; font-size: 0.75rem; color: #4A5568; margin-left: 6px; font-weight: normal; }
.verdict-row { margin-top: 0.5rem; padding: 0.5rem; border-radius: 4px; font-size: 0.9rem; }
.db-confirmed { background: #F0FFF4; border-left: 3px solid #38A169; }
.db-predicted { background: #EBF8FF; border-left: 3px solid #3182CE; }
.db-missing { background: #EDF2F7; border-left: 3px solid #A0AEC0; }
</style>
""", unsafe_allow_html=True)

PROJECT_ROOT = Path(__file__).parent
objects_files = list(PROJECT_ROOT.rglob("objects.tsv"))
if objects_files:
    df_objects = pd.read_csv(objects_files[0], sep="\t")
    case_ids = df_objects["case_id"].tolist()
else:
    all_meta_files = list(PROJECT_ROOT.rglob("*_metadata.json"))
    case_ids = sorted({p.name.replace("_metadata.json", "") for p in all_meta_files})

if "result" not in st.session_state:
    st.session_state["result"] = None
    st.session_state["case_id"] = None
    st.session_state["case_raw"] = None

with st.sidebar:
    st.header("Параметры")
    if not case_ids:
        st.error("Файлы данных не найдены.")
        st.stop()
    selected = st.selectbox("Клинический случай", case_ids)
    max_af = st.slider("Порог частоты (gnomAD)", 0.0, 0.05, 0.01, 0.001)

st.title(f"Кейс: {selected}")

if st.button("🔬 Проанализировать", use_container_width=True, type="primary"):
    with st.spinner("Анализ баз данных (ClinVar, OMIM, VEP)..."):
        case = load_case(selected, PROJECT_ROOT)
        if case is None:
            st.error("Файлы не найдены.")
            st.stop()
        st.session_state["case_raw"] = case
        st.session_state["case_id"] = selected
        st.session_state["result"] = score_case(case, max_af, inheritance_mode="Не учитывать")

if st.session_state["result"] is not None and st.session_state["case_id"] == selected:
    result = st.session_state["result"]
    case = st.session_state["case_raw"]
    
    st.markdown("---")
    st.subheader("➕ Дополнительная задача: уточнить клиническую картину")
    st.markdown("Измените симптомы или тип наследования, чтобы увидеть, как это повлияет на рейтинг.")
    
    col1, col2 = st.columns(2)
    with col1:
        inh_mode = st.selectbox("Ожидаемый тип наследования", ["Не учитывать", "Доминантный (AD)", "Рецессивный (AR)", "X-сцепленный (XL)"])
    with col2:
        new_hpo = st.text_input("Добавить HPO-термины (через запятую)", placeholder="HP:0002616, HP:0001643")
        
    if st.button("🔄 Пересчитать рейтинг", type="secondary"):
        with st.spinner("Пересчет..."):
            case_mod = dict(case)
            if new_hpo.strip():
                new_feats = [{"type": {"id": t.strip(), "label": ""}} for t in new_hpo.split(",") if t.strip()]
                case_mod["phenotype"] = dict(case["phenotype"])
                case_mod["phenotype"]["phenotypicFeatures"] = list(case["phenotype"].get("phenotypicFeatures", [])) + new_feats
            st.session_state["result"] = score_case(case_mod, max_af, inheritance_mode=inh_mode)
            st.rerun()
            
    st.markdown("---")
    
    result = st.session_state["result"]
    if result.empty:
        st.warning(no_candidate_reason(result))
    else:
        top5 = result.head(5)
        st.subheader("Топ-5 генов-кандидатов")
        
        for i, (_, row) in enumerate(top5.iterrows(), 1):
            st.markdown(f"""
            <div class="verdict-box">
                <div class="verdict-title">{i}. {row['gene']} 
                    <span class="verdict-badge">Скор: {row['priority_score']}</span>
                    <span class="verdict-badge">Уверенность: {row['confidence']}</span>
                </div>
            """, unsafe_allow_html=True)
            
            def render_verdict(v, title):
                css_class = "db-confirmed" if v['kind'] == "confirmed" else ("db-predicted" if v['kind'] == "predicted" else "db-missing")
                st.markdown(f"<div class='verdict-row {css_class}'><b>{title}</b> ({v['source']}): {v['verdict']}. <span style='color:#4A5568;font-size:0.85em;'>{v['detail']}</span></div>", unsafe_allow_html=True)

            render_verdict(row["clinvar_verdict"], "Клиническая значимость")
            render_verdict(row["omim_verdict"], "Связь с болезнью")
            render_verdict(row["clingen_verdict"], "Валидация гена")
            render_verdict(row["phenotype_verdict"], "Совпадение симптомов")
            st.markdown("</div>", unsafe_allow_html=True)
