
import streamlit as st
import pandas as pd
from pathlib import Path
from prioritizer import (load_case, patient_hpo_terms, excluded_hpo_terms,
                         score_case, no_candidate_reason, filter_rare, DEFAULT_WEIGHTS)

st.set_page_config(page_title="Приоритизация вариантов", layout="wide")

# ---------- Тема: синий, зелёный, оранжевый ----------
st.markdown("""
<style>
:root {
    --bg: #F8F9FA;
    --sidebar: #E9ECEF;
    --text: #1A365D;
    --accent: #F57C00;
    --success: #2E7D32;
    --warning: #F9A825;
    --danger: #C62828;
    --info: #1976D2;
    --muted: #78909C;
}
.stApp { background: var(--bg); color: var(--text); }
h1, h2, h3, h4 { color: var(--text) !important; }
[data-testid="stSidebar"] { background: var(--sidebar); }
.stButton>button { 
    background: var(--accent); color: #fff; border-radius: 8px; border: none; 
}
.stButton>button:hover { background: #E65100; }
.verdict-card {
    background: #fff; border-left: 4px solid var(--muted);
    border-radius: 8px; padding: 1rem; margin: 0.5rem 0;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08);
}
.confirmed { border-left-color: var(--success); }
.predicted { border-left-color: var(--info); }
.missing   { border-left-color: var(--muted); }
.uncertain { border-left-color: var(--warning); }
.conflicting { border-left-color: var(--danger); }
.badge {
    display: inline-block; font-size: 0.75rem; padding: 2px 8px;
    border-radius: 10px; background: var(--sidebar); color: var(--text); margin-left: 6px;
}
</style>
""", unsafe_allow_html=True)

DATA_DIR = Path("/content")
st.title("Приоритизация генетических вариантов")
st.caption("Прототип. Исследовательское использование, не клинический диагноз.")

with st.sidebar:
    st.header("Настройки")
    case_ids = sorted({p.name.replace("_metadata.json", "") for p in DATA_DIR.glob("case_*_metadata.json")})
    selected = st.selectbox("Клинический случай", case_ids or ["нет данных"])
    max_af = st.slider("Порог частоты (gnomAD)", 0.0, 0.05, 0.01, 0.001)

    st.subheader("Веса признаков")
    st.caption("Ползунки задают относительный вклад. Система нормирует их так, чтобы сумма была 1.0.")
    
    w_ph = st.slider("Фенотип (HPO)", 0.0, 1.0, 0.40, 0.05)
    w_cv = st.slider("ClinVar", 0.0, 1.0, 0.30, 0.05)
    w_im = st.slider("VEP impact", 0.0, 1.0, 0.20, 0.05)
    w_ra = st.slider("Редкость (gnomAD)", 0.0, 1.0, 0.10, 0.05)
    
    total_w = w_ph + w_cv + w_im + w_ra
    if total_w > 0:
        weights = {"phenotype": w_ph/total_w, "clinvar": w_cv/total_w, 
                   "impact": w_im/total_w, "rare": w_ra/total_w}
    else:
        weights = DEFAULT_WEIGHTS
        
    st.markdown(f"**Итоговые веса (сумма = 1.0):**\n"
                f"• Фенотип: `{weights['phenotype']:.2f}`\n"
                f"• ClinVar: `{weights['clinvar']:.2f}`\n"
                f"• VEP impact: `{weights['impact']:.2f}`\n"
                f"• Редкость: `{weights['rare']:.2f}`")

if selected == "нет данных":
    st.warning("Нет данных в /content")
    st.stop()

case = load_case(selected, DATA_DIR)
result = score_case(case, weights, max_af)

# 01 Логика анализа
st.header("01. Логика анализа")
raw_count = len(case["variants"])
rare_df = filter_rare(case["variants"], max_af)
rare_count = len(rare_df)
annotated_count = len(rare_df[rare_df["gene"].notna() & (rare_df["gene"] != "NA")]) if "gene" in rare_df else 0
candidate_count = len(result)

funnel = pd.DataFrame({
    "Этап": ["1. Сырые варианты", "2. После фильтра по частоте", "3. С аннотацией гена", "4. Генов-кандидатов"],
    "Количество": [raw_count, rare_count, annotated_count, candidate_count]
})
st.dataframe(funnel, use_container_width=True, hide_index=True)
st.progress(min(1.0, candidate_count / max(1, raw_count)))
st.caption(f"От {raw_count} вариантов до {candidate_count} генов.")

# 02 Точность вариантов
st.header("02. Точность вариантов")
af_cut = raw_count - rare_count
benign_cut = len(rare_df[rare_df["clinical_significance"].astype(str).str.lower().str.contains("benign", na=False)]) if "clinical_significance" in rare_df else 0

st.markdown(f"""
- **{af_cut}** вариантов отсеяно по частоте.
- **{benign_cut}** помечены как доброкачественные в ClinVar.
""")

labels_path = DATA_DIR / "development_labels.tsv"
if labels_path.exists():
    labels = pd.read_csv(labels_path, sep="\t")
    dev = labels[labels["split"] == "development"]
    hits = []
    for _, r in dev.iterrows():
        c = load_case(r["case_id"], DATA_DIR)
        if c["variants"].empty: continue
        ranked = score_case(c, weights, max_af)
        genes = ranked["gene"].tolist() if not ranked.empty else []
        true = r["causal_gene"]
        hits.append({"case": r["case_id"], "true": true, "top1": genes[0] if genes else None,
                     "hit@1": int(bool(genes) and genes[0] == true),
                     "hit@3": int(true in genes[:3])})
    ev = pd.DataFrame(hits)
    if not ev.empty:
        st.dataframe(ev, use_container_width=True, hide_index=True)
        st.metric("Top-1 accuracy", f"{ev['hit@1'].mean():.0%}")
        st.metric("Top-3 accuracy", f"{ev['hit@3'].mean():.0%}")

# 03 Связь с болезнью
st.header("03. Связь с болезнью")
if result.empty:
    st.warning(no_candidate_reason(result))
else:
    if result["priority_score"].max() < 0.3:
        st.warning(no_candidate_reason(result))

    st.dataframe(result[["gene", "priority_score", "confidence", "phenotype_score", "clinvar_score", "impact_score"]],
                 use_container_width=True, height=300)

    top = result.iloc[0]
    st.subheader("Подробный вердикт топ-кандидата")
    st.markdown(f"### Ген **{top['gene']}** <span class='badge'>скор {top['priority_score']}</span> "
                f"<span class='badge'>уверенность: {top['confidence']}</span>", unsafe_allow_html=True)
    st.markdown(f"**Вариант:** `{top['variant']}` • **Генотип:** `{top['genotype']}`")

    def render_verdict(v, title):
        css = v.get("kind", "unknown")
        link = v.get("link", "")
        st.markdown(
            f"<div class='verdict-card {css}'>"
            f"<b>{title}</b> — {v['verdict']} <span class='badge'>{v['source']}</span><br>"
            f"<span style='color:#455A64'>{v['detail']}</span><br>"
            f"<small>Вес: {v['score']:.2f} | <a href='{link}' target='_blank'>Источник</a></small>"
            f"</div>", unsafe_allow_html=True)

    st.markdown("#### Подтверждено базами данных")
    render_verdict(top["clinvar_verdict"], "ClinVar")
    render_verdict(top["rare_verdict"], "gnomAD")
    st.markdown("#### Вычислительные прогнозы")
    render_verdict(top["impact_verdict"], "VEP")
    render_verdict(top["phenotype_verdict"], "HPO")

    st.markdown("#### Как получен итоговый скор")
    formula = (f"`{weights['phenotype']:.2f}×{top['phenotype_score']} (фенотип) + "
               f"{weights['clinvar']:.2f}×{top['clinvar_score']} (ClinVar) + "
               f"{weights['impact']:.2f}×{top['impact_score']} (VEP) + "
               f"{weights['rare']:.2f}×{top['rare_score']} (gnomAD) = {top['priority_score']}`")
    st.code(formula, language="text")

# 04 Польза для врача
st.header("04. Польза для врача-генетика")
if not result.empty:
    st.download_button(
        label="Скачать отчёт (CSV)",
        data=result.to_csv(index=False).encode("utf-8"),
        file_name=f"{case['case_id']}_report.csv",
        mime="text/csv"
    )
    if result["priority_score"].max() < 0.3:
        st.info("Рекомендации при низкой уверенности: расширить фенотип, проверить аннотацию, рассмотреть трио-секвенирование.")
    st.caption("Исследовательский результат, не клинический диагноз.")

# Доп. задача
st.divider()
st.header("Дополнительная задача: изменить фенотип")
pos = patient_hpo_terms(case["phenotype"])
new_hpo = st.text_input("Новые HPO-термины (через запятую):")
to_exclude = st.multiselect("Исключить из фенотипа:", options=pos)

if st.button("Пересчитать рейтинг"):
    pheno2 = {"phenotypicFeatures": []}
    for feat in case["phenotype"].get("phenotypicFeatures", []):
        tid = feat.get("type", {}).get("id")
        if tid in to_exclude:
            continue
        pheno2["phenotypicFeatures"].append(feat)
    for t in [x.strip() for x in new_hpo.split(",") if x.strip()]:
        pheno2["phenotypicFeatures"].append({"type": {"id": t, "label": ""}})

    case2 = dict(case)
    case2["phenotype"] = pheno2
    res2 = score_case(case2, weights, max_af)

    if res2.empty:
        st.warning(no_candidate_reason(res2))
    else:
        st.success(f"Новый топ-ген: **{res2.iloc[0]['gene']}** (скор {res2.iloc[0]['priority_score']})")
        st.dataframe(res2[["gene", "priority_score", "confidence", "phenotype_score", "clinvar_score"]], use_container_width=True)
