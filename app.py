import streamlit as st
import pandas as pd
from pathlib import Path
from prioritizer import load_case, score_case, no_candidate_reason, filter_rare

st.set_page_config(page_title="Variant Prioritizer", layout="wide", initial_sidebar_state="expanded")

# --- Стили из ChemX ---
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap');
#MainMenu, footer, header { visibility: hidden; }
:root {
    --bg: #0f1117; --surface: #161b27; --border: #1e2535; --accent: #4f8ef7;
    --text: #e2e8f0; --text-muted: #64748b; --success: #10b981;
    --font: 'Inter', sans-serif; --mono: 'IBM Plex Mono', monospace;
}
html, body, [data-testid="stAppViewContainer"], [data-testid="stMain"] { background-color: var(--bg) !important; color: var(--text) !important; font-family: var(--font) !important; }
[data-testid="stSidebar"] { background-color: var(--surface) !important; border-right: 1px solid var(--border) !important; }
[data-testid="stSidebar"] * { color: var(--text) !important; }
.sidebar-logo { display: flex; align-items: center; gap: 10px; margin-bottom: 1.5rem; padding-bottom: 1.25rem; border-bottom: 1px solid var(--border); }
.sidebar-logo-mark { width: 32px; height: 32px; background: var(--accent); border-radius: 6px; display: flex; align-items: center; justify-content: center; }
.sidebar-logo-text { font-size: 0.9rem; font-weight: 600; }
.sidebar-logo-sub { font-size: 0.7rem; color: var(--text-muted); text-transform: uppercase; }
.page-header { padding: 2rem 0 1.5rem; border-bottom: 1px solid var(--border); margin-bottom: 1.75rem; }
.page-header h1 { font-size: 1.5rem; color: var(--text); margin: 0; }
.page-header p { font-size: 0.85rem; color: var(--text-muted); margin: 0; }
.section-label { font-size: 0.65rem; font-weight: 600; letter-spacing: 0.1em; text-transform: uppercase; color: var(--text-muted); margin: 1.5rem 0 0.75rem; }
.verdict-card { background: var(--surface); border-left: 4px solid var(--accent); border-radius: 6px; padding: 1rem; margin: 0.5rem 0; border: 1px solid var(--border); }
.confirmed { border-left-color: var(--success); }
.predicted { border-left-color: var(--accent); }
.badge { display: inline-block; background: var(--bg); border: 1px solid var(--border); border-radius: 4px; padding: 2px 8px; font-size: 0.72rem; font-family: var(--mono); color: var(--accent); margin-left: 6px; }
[data-testid="stMetric"] { background: var(--surface) !important; border: 1px solid var(--border) !important; border-radius: 6px !important; padding: 1rem !important; }
[data-testid="stMetricValue"] { color: var(--accent) !important; font-family: var(--mono) !important; }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="page-header">
  <h1>Variant Prioritizer</h1>
  <p>Prioritization of genetic variants for rare disease diagnostics</p>
</div>
""", unsafe_allow_html=True)

# --- Архитектура: поиск данных ---
PROJECT_ROOT = Path(__file__).parent

with st.sidebar:
    st.markdown("""
    <div class="sidebar-logo">
        <div class="sidebar-logo-mark"><svg width="16" height="16" viewBox="0 0 16 16" fill="none"><path d="M8 2L14 8L8 14L2 8L8 2Z" stroke="white" stroke-width="1.5"/></svg></div>
        <div><div class="sidebar-logo-text">GenoX</div><div class="sidebar-logo-sub">Prioritizer v1</div></div>
    </div>
    """, unsafe_allow_html=True)

    st.markdown('<div class="section-label">Case Selection</div>', unsafe_allow_html=True)
    all_meta_files = list(PROJECT_ROOT.rglob("*_metadata.json"))
    case_ids = sorted({p.name.replace("_metadata.json", "") for p in all_meta_files})
    
    if not case_ids:
        st.error("❌ Файлы данных не найдены. Убедитесь, что папка data загружена в репозиторий.")
        st.stop()
        
    selected = st.selectbox("Clinical Case", case_ids, label_visibility="collapsed")
    max_af = st.slider("gnomAD AF Threshold", 0.0, 0.05, 0.01, 0.001)
    inh_mode = st.selectbox("Inheritance Mode", ["Не учитывать", "Доминантный (AD)", "Рецессивный (AR)", "X-сцепленный (XL)"], label_visibility="collapsed")
    new_hpo = st.text_input("Add HPO Terms", placeholder="HP:0002616, HP:0001643", label_visibility="collapsed")

# Загрузка
case = load_case(selected, PROJECT_ROOT)
if case is None:
    st.error(f"Не удалось загрузить файлы для {selected}.")
    st.stop()

# Доп. задача
case_mod = dict(case)
if new_hpo.strip():
    new_feats = [{"type": {"id": t.strip(), "label": ""}} for t in new_hpo.split(",") if t.strip()]
    case_mod["phenotype"] = dict(case["phenotype"])
    case_mod["phenotype"]["phenotypicFeatures"] = list(case["phenotype"].get("phenotypicFeatures", [])) + new_feats

# Скоринг
result = score_case(case_mod, max_af, inheritance_mode=inh_mode)

# --- Вкладки ---
tab1, tab2 = st.tabs(["Analysis", "Verdicts"])

with tab1:
    st.markdown('<div class="section-label">Pipeline Funnel</div>', unsafe_allow_html=True)
    raw_count = len(case["variants"])
    rare_count = len(filter_rare(case["variants"], max_af))
    candidate_count = len(result)
    
    c1, c2, c3 = st.columns(3)
    c1.metric("Raw Variants", raw_count)
    c2.metric("After AF Filter", rare_count)
    c3.metric("Gene Candidates", candidate_count)

with tab2:
    if result.empty:
        st.warning(no_candidate_reason(result))
    else:
        st.dataframe(result[["gene", "priority_score", "confidence", "phenotype_score", "clinvar_score"]], use_container_width=True)
        
        top = result.iloc[0]
        st.markdown(f"### Top: {top['gene']} <span class='badge'>Score: {top['priority_score']}</span>", unsafe_allow_html=True)
        
        def render(v, title):
            st.markdown(f"<div class='verdict-card {v['kind']}'><b>{title}</b> — {v['verdict']} <span class='badge'>{v['source']}</span><br><span style='color:#94a3b8'>{v['detail']}</span></div>", unsafe_allow_html=True)

        render(top["clinvar_verdict"], "ClinVar")
        render(top["phenotype_verdict"], "Phenotype (HPO)")
        
    st.download_button("Download CSV", result.to_csv(index=False).encode("utf-8"), f"{selected}_report.csv", "text/csv")
