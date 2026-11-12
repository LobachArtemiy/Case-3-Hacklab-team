import streamlit as st
import pandas as pd
from pathlib import Path
from prioritizer import load_local_case, process_uploaded_files, score_case, no_candidate_reason

st.set_page_config(page_title="Variant Prioritizer", layout="wide")

st.markdown("""
<style>
/* ... стили ... */
html, body, [data-testid="stAppViewContainer"] { background-color: #F8F9FA !important; font-family: sans-serif; }
[data-testid="stSidebar"] { background-color: #FFFFFF !important; border-right: 1px solid #E2E8F0; }
.verdict-box { background: #FFF; border-left: 4px solid #3182CE; border-radius: 6px; padding: 1.25rem; margin-bottom: 1.5rem; border: 1px solid #E2E8F0; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }
.verdict-title { font-size: 1.3rem; font-weight: 700; margin-bottom: 0.25rem; color: #1A202C; }
.verdict-var { font-size: 0.85rem; color: #718096; margin-bottom: 1rem; font-family: monospace; }
.verdict-badge { display: inline-block; background: #EDF2F7; border-radius: 4px; padding: 2px 8px; font-size: 0.75rem; color: #4A5568; margin-left: 8px; font-weight: normal; }
.verdict-row { margin-top: 0.5rem; padding: 0.75rem; border-radius: 4px; font-size: 0.9rem; line-height: 1.5; }
.db-confirmed { background: #F0FFF4; border-left: 3px solid #38A169; }
.db-predicted { background: #EBF8FF; border-left: 3px solid #3182CE; }
.db-missing { background: #EDF2F7; border-left: 3px solid #A0AEC0; }
.disease-list { margin-top: 0.5rem; padding-left: 1rem; }
.disease-item { margin-bottom: 1rem; padding: 0.5rem; background: #F7FAFC; border-radius: 4px; border: 1px solid #E2E8F0; }
.disease-name { font-weight: 700; color: #2D3748; font-size: 1rem; }
.aa-change { color: #D69E2E; font-weight: bold; }
</style>
""", unsafe_allow_html=True)

PROJECT_ROOT = Path(__file__).parent
objects_files = list(PROJECT_ROOT.rglob("objects.tsv"))
case_ids = []
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
    
    # --- Блок загрузки файлов ---
    st.subheader("Загрузить свой кейс")
    st.markdown("Поддерживается VCF или TSV.")
    uploaded_pheno = st.file_uploader("Фенотип (JSON)", type=["json"])
    uploaded_var = st.file_uploader("Варианты (TSV или VCF)", type=["tsv", "vcf"])
    uploaded_label = st.file_uploader("Метки (TSV, опционально)", type=["tsv"])
    
    st.divider()
    
    st.subheader("Или выберите встроенный")
    selected = st.selectbox("Клинический случай", case_ids if case_ids else ["Нет"])
    max_af = st.slider("Порог частоты (gnomAD)", 0.0, 0.05, 0.01, 0.001)

st.title("Гибридная система приоритизации")
st.caption("Экспертные правила + Граф Знаний (HPO) + Статистика")

if st.button("Проанализировать", use_container_width=True, type="primary"):
    with st.spinner("Анализ..."):
        # Приоритет у загруженных файлов
        if uploaded_pheno and uploaded_var:
            case = process_uploaded_files(uploaded_pheno, uploaded_var, uploaded_label)
        else:
            case = load_local_case(selected, PROJECT_ROOT)
            
        if case is None or "error" in case:
            st.error(case.get("error", "Ошибка загрузки данных."))
        else:
            st.session_state["case_raw"] = case
            st.session_state["case_id"] = case["case_id"]
            st.session_state["result"] = score_case(case, max_af, inheritance_mode="Не учитывать")

if st.session_state["result"] is not None:
    result = st.session_state["result"]
    case = st.session_state["case_raw"]
    
    st.caption(f"Источник данных: {case.get('status', 'Неизвестно')}")
    
    if result.empty:
        st.warning(no_candidate_reason(result))
    else:
        num_found = len(result)
        if num_found == 1: st.subheader("Найден 1 ген-кандидат")
        elif num_found <= 5: st.subheader(f"Найдено {num_found} гена-кандидата")
        else: st.subheader("Топ-5 генов-кандидатов")
            
        top5 = result.head(5)
        
        for i, (_, row) in enumerate(top5.iterrows(), 1):
            hgvsp = row.get('hgvsp', 'Неизвестно')
            hgvsp_display = f" <span class='aa-change'>({hgvsp})</span>" if hgvsp != "Неизвестно" else ""
            
            st.markdown(f"""
            <div class="verdict-box">
                <div class="verdict-title">{i}. {row['gene']} {hgvsp_display}
                    <span class="verdict-badge">Score: {row['priority_score']}</span>
                    <span class="verdict-badge">Confidence: {row['confidence']}</span>
                </div>
                <div class="verdict-var">Вариант: {row['variant']}</div>
            """, unsafe_allow_html=True)
            
            def render_verdict(v, title):
                css_class = "db-confirmed" if v['kind'] == "confirmed" else ("db-predicted" if v['kind'] == "predicted" else "db-missing")
                st.markdown(f"<div class='verdict-row {css_class}'><b>{title}</b> ({v['source']}): <b>{v['verdict']}</b>.<br><span style='color:#4A5568;'>{v['detail']}</span></div>", unsafe_allow_html=True)

            render_verdict(row["clinvar_verdict"], "Клиническая значимость")
            render_verdict(row["impact_verdict"], "Функциональный эффект")
            render_verdict(row["phenotype_verdict"], "Совпадение симптомов")
            render_verdict(row["inheritance_verdict"], "Соответствие наследованию")
            
            dv = row["disease_verdict"]
            st.markdown(f"<div class='verdict-row db-confirmed'><b>Связь с болезнью</b> ({dv['source']}): <b>{dv['verdict']}</b>.<br><span style='color:#4A5568;'>{dv['detail']}</span></div>", unsafe_allow_html=True)
            
            if dv["diseases"]:
                st.markdown("<div class='disease-list'>", unsafe_allow_html=True)
                for dis in dv["diseases"]:
                    st.markdown(f"<div class='disease-item'><div class='disease-name'>{dis['name']} (OMIM:{dis['id']})</div></div>", unsafe_allow_html=True)
                st.markdown("</div>", unsafe_allow_html=True)
                
            st.markdown("</div>", unsafe_allow_html=True)
            
    st.markdown("---")
    st.subheader("Дополнительная задача")
    col1, col2 = st.columns(2)
    with col1:
        inh_mode = st.selectbox("Тип наследования", ["Не учитывать", "Доминантный (AD)", "Рецессивный (AR)", "X-сцепленный (XL)"])
    with col2:
        new_hpo = st.text_input("Добавить HPO (через запятую)", placeholder="HP:0002616")
        
    if st.button("Пересчитать рейтинг", type="secondary"):
        with st.spinner("Пересчет..."):
            case_mod = dict(case)
            if new_hpo.strip():
                new_feats = [{"type": {"id": t.strip(), "label": ""}} for t in new_hpo.split(",") if t.strip()]
                case_mod["phenotype"] = dict(case["phenotype"])
                case_mod["phenotype"]["phenotypicFeatures"] = list(case["phenotype"].get("phenotypicFeatures", [])) + new_feats
            st.session_state["result"] = score_case(case_mod, max_af, inheritance_mode=inh_mode)
            st.rerun()
