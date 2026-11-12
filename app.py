import streamlit as st
import pandas as pd
from pathlib import Path
from prioritizer import load_case, score_case, no_candidate_reason, filter_rare

st.set_page_config(page_title="Variant Prioritizer", layout="wide", initial_sidebar_state="expanded")

# --- Светлая тема с акцентами как на слайде ---
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
#MainMenu, footer, header { visibility: hidden; }
:root {
    --bg: #F8F9FA; --card: #FFFFFF; --border: #E0E0E0;
    --primary: #5C83EA; --text: #1A202C; --text-muted: #718096;
}
html, body, [data-testid="stAppViewContainer"], [data-testid="stMain"] { 
    background-color: var(--bg) !important; color: var(--text) !important; font-family: 'Inter', sans-serif !important;
}
[data-testid="stSidebar"] { background-color: #FFFFFF !important; border-right: 1px solid var(--border) !important; }
[data-testid="stSidebar"] * { color: var(--text) !important; }

/* Заголовки */
.page-header { padding: 1.5rem; background: var(--primary); color: white; border-radius: 8px; margin-bottom: 1.5rem; }
.page-header h1 { font-size: 1.5rem; margin: 0; color: white !important; text-transform: uppercase; font-weight: 700; }
.page-header p { margin: 0; font-size: 0.9rem; opacity: 0.9; }

/* Блоки критериев */
.crit-card { background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 1.25rem; height: 100%; margin-bottom: 1rem; box-shadow: 0 2px 4px rgba(0,0,0,0.02); }
.crit-title { font-size: 1.1rem; font-weight: 700; color: var(--text); margin-bottom: 0.25rem; display: flex; align-items: center; gap: 8px; text-transform: uppercase; }
.crit-num { color: var(--primary); font-size: 0.9rem; font-weight: 700; }
.crit-desc { font-size: 0.85rem; color: var(--text-muted); margin-bottom: 1rem; line-height: 1.4; }
.crit-content { font-size: 0.9rem; color: var(--text); }

/* Вердикты */
.verdict-card { background: #F7FAFC; border-left: 4px solid var(--primary); border-radius: 4px; padding: 0.75rem 1rem; margin: 0.5rem 0; }
.confirmed { border-left-color: #38A169; background: #F0FFF4; }
.predicted { border-left-color: #3182CE; background: #EBF8FF; }
.badge { display: inline-block; background: #EDF2F7; border-radius: 4px; padding: 2px 8px; font-size: 0.7rem; font-weight: 600; color: #4A5568; margin-left: 6px; }
[data-testid="stMetric"] { background: #FFF !important; border: 1px solid var(--border) !important; border-radius: 8px !important; padding: 1rem !important; }
[data-testid="stMetricValue"] { color: var(--primary) !important; }
</style>
""", unsafe_allow_html=True)

# --- Поиск файлов ---
PROJECT_ROOT = Path(__file__).parent
all_meta_files = list(PROJECT_ROOT.rglob("*_metadata.json"))
case_ids = sorted({p.name.replace("_metadata.json", "") for p in all_meta_files})

# --- Sidebar ---
with st.sidebar:
    st.header("⚙️ Настройки")
    if not case_ids:
        st.error("❌ Файлы данных не найдены!")
        st.write("Убедитесь, что папки `development` и `hidden` лежат внутри папки `data` в репозитории.")
        st.stop()
        
    selected = st.selectbox("Клинический случай", case_ids)
    max_af = st.slider("Порог частоты (gnomAD)", 0.0, 0.05, 0.01, 0.001)
    
    st.divider()
    st.subheader("Доп. задача")
    inh_mode = st.selectbox("Тип наследования", ["Не учитывать", "Доминантный (AD)", "Рецессивный (AR)"])
    new_hpo = st.text_input("Добавить HPO-термины", placeholder="HP:0002616, HP:0001643")

# Загрузка
case = load_case(selected, PROJECT_ROOT)
if case is None:
    st.error(f"Не удалось загрузить файлы для {selected}.")
    st.stop()

# Применяем доп. задачу
case_mod = dict(case)
if new_hpo.strip():
    new_feats = [{"type": {"id": t.strip(), "label": ""}} for t in new_hpo.split(",") if t.strip()]
    case_mod["phenotype"] = dict(case["phenotype"])
    case_mod["phenotype"]["phenotypicFeatures"] = list(case["phenotype"].get("phenotypicFeatures", [])) + new_feats

# Скоринг
result = score_case(case_mod, max_af, inheritance_mode=inh_mode)
top = result.iloc[0] if not result.empty else None

# --- Header ---
st.markdown("""
<div class="page-header">
  <p>КЕЙС 3 • ПИШ ИТМО</p>
  <h1>Поиск генетических причин редких наследственных заболеваний</h1>
</div>
""", unsafe_allow_html=True)

# --- Строгое соответствие 4 критериям ---
col1, col2 = st.columns(2)

with col1:
    # 01 Логика анализа
    raw = len(case["variants"])
    rare = len(filter_rare(case["variants"], max_af))
    cand = len(result)
    st.markdown(f"""
    <div class="crit-card">
        <div class="crit-title"><span class="crit-num">01</span> Логика анализа</div>
        <div class="crit-desc">Путь от сырых геномных данных к кандидатам: фильтрация, аннотация и приоритизация вариантов.</div>
        <div class="crit-content">
            Воронка: <b>{raw}</b> вариантов → фильтр по частоте (<b>{rare}</b>) → аннотация → <b>{cand}</b> кандидатов.
            <br><small style="color:#E53E3E;">Отсеяно: {raw - rare} частых полиморфизмов.</small>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # 03 Связь с болезнью
    if top:
        cv_badge = f"<span class='badge'>{top['clinvar_verdict']['source']}</span>"
        ph_badge = f"<span class='badge'>{top['phenotype_verdict']['source']}</span>"
        st.markdown(f"""
        <div class="crit-card">
            <div class="crit-title"><span class="crit-num">03</span> Связь с болезнью</div>
            <div class="crit-desc">Чем подтверждена связь варианта с заболеванием: базы данных, литература, критерии патогенности.</div>
            <div class="crit-content">
                Топ-ген: <b>{top['gene']}</b> (скор {top['priority_score']})<br>
                <div style="margin-top:8px;" class="confirmed">
                    <b>ClinVar:</b> {top['clinvar_verdict']['verdict']} {cv_badge}
                </div>
                <div style="margin-top:8px;" class="predicted">
                    <b>Фенотип (HPO):</b> {top['phenotype_verdict']['verdict']} {ph_badge}
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)

with col2:
    # 02 Точность вариантов
    st.markdown(f"""
    <div class="crit-card">
        <div class="crit-title"><span class="crit-num">02</span> Точность вариантов</div>
        <div class="crit-desc">Насколько найденные варианты действительно могут быть причиной болезни и как отсеяны ложные находки.</div>
        <div class="crit-content">
            Отсеяны варианты с частотой > {max_af} (gnomAD) и доброкачественные (ClinVar).
            <br>Использованы фиксированные веса ACMG: Фенотип (0.5), База (0.3), VEP (0.15).
        </div>
    </div>
    """, unsafe_allow_html=True)

    # 04 Польза для врача
    st.markdown(f"""
    <div class="crit-card">
        <div class="crit-title"><span class="crit-num">04</span> Польза для врача-генетика</div>
        <div class="crit-desc">Насколько результат понятен специалисту и пригоден для принятия клинического решения.</div>
        <div class="crit-content">
            Разделение подтвержденных данных БД (зеленый) и прогнозов (синий).
            Реализовано доп. задание: учет наследования ({inh_mode}) и добавление симптомов.
        </div>
    </div>
    """, unsafe_allow_html=True)

st.divider()

# --- Таблица кандидатов ---
st.subheader("Ранжированные кандидаты")
if result.empty:
    st.warning(no_candidate_reason(result))
else:
    st.dataframe(result[["gene", "priority_score", "confidence", "phenotype_score", "clinvar_score"]], use_container_width=True)
    
    st.download_button("📄 Скачать отчёт (CSV)", result.to_csv(index=False).encode("utf-8"), f"{selected}_report.csv", "text/csv")
    st.caption("⚠️ Исследовательский результат. Не является клиническим диагнозом.")
