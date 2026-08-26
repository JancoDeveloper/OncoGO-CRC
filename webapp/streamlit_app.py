"""OncoGO-CRC — App web Nivel 2 (dashboard interactivo).

Explora la **tabla de hipótesis priorizadas** (Capa C), el **volcano** diferencial y un **explorador de
red**. Lee artefactos reales del pipeline si existen; si no, muestra un conjunto de demostración sintético.

Ejecutar:
    pip install -r webapp/requirements.txt
    streamlit run webapp/streamlit_app.py

Fuente de datos: variable de entorno ``ONCOGO_CRC_RESULTS`` o el campo de la barra lateral.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from lib import data_access as da
from lib import plots

st.set_page_config(page_title="OncoGO-CRC · Nivel 2", page_icon="🧬", layout="wide")


@st.cache_data(show_spinner=False)
def _load(results_dir: str) -> da.Dataset:
    return da.load_dataset(Path(results_dir) if results_dir else None)


# ----------------------------------------------------------------------------- Sidebar
st.sidebar.title("🧬 OncoGO-CRC")
st.sidebar.caption("Dashboard interactivo de priorización de hipótesis (Nivel 2)")
default_dir = str(da.default_results_dir())
results_dir = st.sidebar.text_input("Directorio de resultados", value=default_dir)
if st.sidebar.button("↻ Recargar datos"):
    _load.clear()

ds = _load(results_dir)

if ds.is_demo:
    st.sidebar.warning("Modo DEMOSTRACIÓN (datos sintéticos). Ejecuta el pipeline para ver resultados reales.")
else:
    st.sidebar.success(f"Datos reales: {ds.results_dir}")

st.sidebar.markdown("---")
st.sidebar.markdown(
    "**Aviso.** Los resultados son **hipótesis exploratorias** con evidencia convergente, "
    "no dianas terapéuticas validadas."
)

# ----------------------------------------------------------------------------- Header
st.title("OncoGO-CRC — priorización de hipótesis (proteómica de CCR)")
st.markdown(
    "Herramienta metodológica reproducible: integra proteómica pública de cáncer colorrectal, "
    "identifica una firma reproducible y **prioriza candidatos por evidencia convergente** "
    "(diferencial, reproducibilidad, esencialidad, tractabilidad, especificidad, seguridad y genómica)."
)

hyp = ds.hypotheses.copy()
tab_overview, tab_table, tab_volcano, tab_network, tab_prov = st.tabs(
    ["📊 Resumen", "🧾 Hipótesis", "🌋 Volcano", "🕸️ Red", "🔎 Procedencia"]
)

# ============================================================== Resumen
with tab_overview:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Candidatos priorizados", len(hyp))
    if "convergence_score" in hyp.columns:
        c2.metric("Score máx.", f"{hyp['convergence_score'].max():.3f}")
    if "evidence_count" in hyp.columns:
        c3.metric("Evidencias medias/gen", f"{hyp['evidence_count'].mean():.1f}")
    n_prov = sum(1 for v in ds.provenance.get("providers", {}).values() if v.get("available"))
    c4.metric("Fuentes externas activas", n_prov)

    if "evidence_tier" in hyp.columns:
        st.subheader("Distribución por tier de evidencia")
        counts = hyp["evidence_tier"].value_counts().sort_index()
        st.bar_chart(counts)

    st.subheader("Top 15 por score de convergencia")
    top_cols = [c for c in ["rank", "gene_symbol", "convergence_score", "evidence_tier",
                            "evidence_count", "log2fc", "padj"] if c in hyp.columns]
    st.dataframe(hyp[top_cols].head(15), use_container_width=True, hide_index=True)

# ============================================================== Tabla de hipótesis
with tab_table:
    st.subheader("Tabla de hipótesis priorizadas")
    fcol1, fcol2, fcol3 = st.columns([2, 1, 1])
    query = fcol1.text_input("Filtrar por gen (contiene)", "")
    tiers = sorted(hyp["evidence_tier"].dropna().unique()) if "evidence_tier" in hyp.columns else []
    sel_tiers = fcol2.multiselect("Tier", tiers, default=tiers)
    min_ev = fcol3.slider("Mín. evidencias", 0,
                          int(hyp["evidence_count"].max()) if "evidence_count" in hyp.columns else 0, 0)

    view = hyp
    if query:
        view = view[view["gene_symbol"].astype(str).str.contains(query, case=False, na=False)]
    if sel_tiers and "evidence_tier" in view.columns:
        view = view[view["evidence_tier"].isin(sel_tiers)]
    if "evidence_count" in view.columns:
        view = view[view["evidence_count"] >= min_ev]

    st.caption(f"{len(view)} de {len(hyp)} candidatos")
    st.dataframe(view, use_container_width=True, hide_index=True)
    st.download_button(
        "⬇️ Descargar CSV filtrado", view.to_csv(index=False).encode("utf-8"),
        "hypothesis_priority_filtrado.csv", "text/csv",
    )

    if not ds.evidence_matrix.empty:
        st.subheader("Perfil de evidencia por gen")
        gene = st.selectbox("Gen", view["gene_symbol"].tolist() if not view.empty else hyp["gene_symbol"].tolist())
        prof = ds.evidence_matrix[ds.evidence_matrix["gene_symbol"].astype(str).str.upper() == str(gene).upper()]
        if not prof.empty and "favorability" in prof.columns:
            chart = prof.set_index("signal")["favorability"].dropna()
            st.bar_chart(chart)
            st.caption("Favorabilidad (0–1) por señal disponible. Señales ausentes no se muestran ni penalizan.")

# ============================================================== Volcano
with tab_volcano:
    st.subheader("Volcano diferencial")
    if not ds.differentials:
        st.info("No hay tablas diferenciales disponibles.")
    else:
        comp = st.selectbox("Comparación", list(ds.differentials.keys()))
        cc1, cc2 = st.columns(2)
        alpha = cc1.slider("Umbral padj", 0.001, 0.2, 0.05, 0.001)
        min_fc = cc2.slider("Mín. |log2FC|", 0.0, 3.0, 0.0, 0.1)
        top_genes = set(hyp.get("gene_symbol", pd.Series(dtype=str)).head(10).astype(str))
        hl = st.multiselect("Resaltar genes", sorted(hyp["gene_symbol"].astype(str).unique()),
                            default=sorted(top_genes))
        fig = plots.volcano_figure(ds.differentials[comp], alpha=alpha, min_abs_log2fc=min_fc,
                                   highlight=set(hl))
        st.plotly_chart(fig, use_container_width=True)

# ============================================================== Red
with tab_network:
    st.subheader("Explorador de red (coexpresión)")
    if ds.processed_matrix is None:
        st.info("No hay matriz procesada disponible para calcular correlaciones.")
    else:
        n_top = st.slider("Top-N candidatos", 5, 60, 25, 5)
        thr = st.slider("Umbral |correlación|", 0.3, 0.95, 0.6, 0.05)
        top = hyp.head(n_top)["gene_symbol"].astype(str).tolist()
        G = plots.correlation_network(ds.processed_matrix, top, threshold=thr)
        scores = dict(zip(hyp["gene_symbol"].astype(str),
                          hyp.get("convergence_score", pd.Series(dtype=float))))
        st.caption(f"{G.number_of_nodes()} nodos · {G.number_of_edges()} aristas (|r| ≥ {thr})")
        st.plotly_chart(plots.network_figure(G, scores=scores), use_container_width=True)
        st.caption("Color/tamaño de nodo = score de convergencia / grado. "
                   "Aristas = |correlación| entre perfiles proteicos por encima del umbral.")

# ============================================================== Procedencia
with tab_prov:
    st.subheader("Procedencia y reproducibilidad")
    prov = ds.provenance
    st.markdown(f"**Comparación primaria:** `{prov.get('primary_comparison', 'n/d')}`  ·  "
                f"**Candidatos:** {prov.get('n_candidates', len(hyp))}  ·  "
                f"**Generado:** {prov.get('generated_at', 'n/d')}")
    st.markdown(f"> {prov.get('note', '')}")

    providers = prov.get("providers", {})
    if providers:
        rows = [{"fuente": k, "disponible": v.get("available"),
                 "versión": v.get("version"), "n_genes": v.get("n_genes"),
                 "estado": v.get("status")} for k, v in providers.items()]
        st.markdown("**Proveedores de evidencia externa**")
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    weights = prov.get("weights", {})
    if weights:
        st.markdown("**Pesos del score de convergencia**")
        st.dataframe(pd.DataFrame([{"señal": k, "peso": v} for k, v in weights.items()]),
                     use_container_width=True, hide_index=True)
