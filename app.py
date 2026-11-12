import streamlit as st
import pandas as pd
from pathlib import Path
from prioritizer import load_case, patient_hpo_terms, score_case, no_candidate_reason, filter_rare

st.set_page_config(page_title="Приоритизация вариантов", layout="wide")

# --- Тёплая пастельная тема (без белого) ---
st.markdown("""
<style>
:root {
    --bg: #F4EBE1;
    --sidebar: #EADBC8;
    --card-bg: #FCF8F2;
    --text: #3E2723;
    --accent: #D84315;
    --success: #558B2F;
    --info: #E65100;
    --muted: #A1887F;
}
.stApp { background: var(--bg); color: var(--text); }
h1, h2, h3, h4 { color: var(--text) !important; }
[data-testid="stSidebar"] { background: var(--sidebar); }
.stButton>button { background: var(--accent); color: #FFF; border-radius: 8px; border: none; }
.stButton>button:hover { background: #BF360C; }
.stDataFrame { background: var(--card-bg); border-radius: 8px; }
.verdict-card { background: var(--card-bg); border-left: 5px solid var(--muted); border-radius: 8px; padding: 1rem; margin: 0.5rem 0; box-shadow: 0 2px 5px rgba(62,39,35,0.1); }
.confirmed { border-left-color: var(--success); }
.predicted { border-left-color: var(--info); }
.missing { border-left-color: var(--muted); }
.badge { display: inline-block; font-size: 0.75rem; padding: 2px 8px; border-radius: 10px; background: var(--sidebar); color: var(--text); margin-left: 6px; }
</style>
""", unsafe_allow_html=True)

st.title("Приоритизация генетических вариантов")

# --- Поиск файлов (максимально надежный) ---
project_root = Path(__file__).parent
all_meta_files = list(project_root.rglob("*_metadata.json"))
case_ids = sorted({p.name.replace("_metadata.json", "") for p in all_meta_files})

# --- Sidebar ---
with st.sidebar:
    st.header("Настройки")
    if not case_ids:
        st.error("Файлы данных не найдены. Убедитесь, что папка data с файлами case_* загружена в репозиторий.")
        st.stop()
        
    selected = st.selectbox("Клинический случай", case_ids)
    max_af = st.slider("Порог частоты (gnomAD)", 0.0, 0.05, 0.01, 0.001)
    st.info("Веса фиксированы по ACMG: Фенотип (0.5), ClinVar (0.3), VEP (0.15), Редкость (0.05).")

# Загружаем данные
case = load_case(selected, project_root / "data")

# Проверка на отсутствие аннотаций
if case["variants"].empty:
    st.error(f"Файлы вариантов для {selected} пусты.")
    st.stop()
    
if "gene" not in case["variants"].columns or case["variants"]["gene"].dropna().eq("NA").all():
    st.warning("В файлах вариантов нет аннотаций генов (колонка gene = NA). Скоринг не сможет работать. Убедитесь, что вы загрузили правильные файлы `_annotated.tsv`.")

# --- Дополнительная задача ---
st.header("➕ Дополнительная задача: Изменение параметров")
col1, col2 = st.columns(2)
with col1:
    inh_mode = st.selectbox("Ожидаемый тип наследования:", ["Не учитывать", "Доминантный (AD)", "Рецессивный (AR)", "X-сцепленный (XL)"])
with col2:
    new_hpo = st.text_input("Добавить HPO-термины (через запятую):")

case_mod = dict(case)
if new_hpo.strip():
    new_feats = [{"type": {"id": t.strip(), "label": ""}} for t in new_hpo.split(",") if t.strip()]
    case_mod["phenotype"] = dict(case["phenotype"])
    case_mod["phenotype"]["phenotypicFeatures"] = list(case["phenotype"].get("phenotypicFeatures", [])) + new_feats

result = score_case(case_mod, max_af, inheritance_mode=inh_mode)

# --- Отображение ---
st.header("01. Логика анализа")
raw_count = len(case["variants"])
rare_count = len(filter_rare(case["variants"], max_af))
candidate_count = len(result)
st.dataframe(pd.DataFrame({"Этап": ["Сырые", "После фильтра", "Кандидаты"], "Количество": [raw_count, rare_count, candidate_count]}), hide_index=True)

st.header("02. Точность вариантов")
st.markdown(f"Отсеяно по частоте: **{raw_count - rare_count}** вариантов. Ранжирование по фиксированным весам.")

st.header("03. Связь с болезнью")
if result.empty:
    st.warning(no_candidate_reason(result))
else:
    if result["priority_score"].max() < 0.3:
        st.warning(no_candidate_reason(result))

    st.dataframe(result[["gene", "priority_score", "confidence", "phenotype_score", "clinvar_score", "inheritance_match"]], use_container_width=True)
    
    top = result.iloc[0]
    st.subheader(f"Топ-кандидат: {top['gene']} (скор {top['priority_score']})")
    
    def render(v, title):
        st.markdown(f"<div class='verdict-card {v['kind']}'><b>{title}</b> — {v['verdict']} <span class='badge'>{v['source']}</span><br><span style='color:#5D4037'>{v['detail']}</span></div>", unsafe_allow_html=True)

    render(top["clinvar_verdict"], "ClinVar")
    render(top["phenotype_verdict"], "Фенотип")
    render(top["impact_verdict"], "VEP")

st.header("04. Польза для врача")
if not result.empty:
    st.download_button("Скачать отчёт (CSV)", result.to_csv(index=False).encode("utf-8"), f"{selected}_report.csv")
st.caption("⚠️ Исследовательский результат. Не является диагнозом.")
