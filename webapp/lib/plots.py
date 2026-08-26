"""Gráficas interactivas (Plotly) para la app web: volcano y explorador de red."""

from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objects as go


def volcano_figure(
    diff: pd.DataFrame,
    *,
    alpha: float = 0.05,
    min_abs_log2fc: float = 0.0,
    highlight: set[str] | None = None,
) -> go.Figure:
    """Volcano interactivo: log2FC vs -log10(padj), color por estado, hover por gen."""
    d = diff.copy()
    gene_col = "gene_symbol" if "gene_symbol" in d.columns else "protein_id"
    d["_gene"] = d[gene_col].astype(str)
    d["log2fc"] = pd.to_numeric(d.get("log2fc"), errors="coerce")
    d["padj"] = pd.to_numeric(d.get("padj"), errors="coerce").clip(lower=1e-300)
    d = d.dropna(subset=["log2fc", "padj"])
    d["neglog10_padj"] = -np.log10(d["padj"])

    up = (d["padj"] < alpha) & (d["log2fc"] >= max(min_abs_log2fc, 1e-9))
    down = (d["padj"] < alpha) & (d["log2fc"] <= -max(min_abs_log2fc, 1e-9))
    d["estado"] = np.select([up, down], ["Sube (tumor)", "Baja (tumor)"], default="No significativo")

    color_map = {"Sube (tumor)": "#c0392b", "Baja (tumor)": "#2471a3", "No significativo": "#b2babb"}
    fig = go.Figure()
    for estado, sub in d.groupby("estado"):
        fig.add_trace(
            go.Scattergl(
                x=sub["log2fc"],
                y=sub["neglog10_padj"],
                mode="markers",
                name=estado,
                marker={"size": 6, "color": color_map.get(estado, "#888"), "opacity": 0.75},
                text=sub["_gene"],
                hovertemplate="<b>%{text}</b><br>log2FC=%{x:.2f}<br>-log10(padj)=%{y:.2f}<extra></extra>",
            )
        )
    if highlight:
        h = d[d["_gene"].str.upper().isin({g.upper() for g in highlight})]
        if not h.empty:
            fig.add_trace(
                go.Scattergl(
                    x=h["log2fc"], y=h["neglog10_padj"], mode="markers+text",
                    name="Seleccionados", text=h["_gene"], textposition="top center",
                    marker={"size": 12, "color": "#f1c40f", "line": {"width": 1.5, "color": "#000"}},
                    hovertemplate="<b>%{text}</b><br>log2FC=%{x:.2f}<br>-log10(padj)=%{y:.2f}<extra></extra>",
                )
            )
    fig.add_hline(y=-np.log10(alpha), line_dash="dash", line_color="#7f8c8d",
                  annotation_text=f"padj={alpha}")
    if min_abs_log2fc > 0:
        fig.add_vline(x=min_abs_log2fc, line_dash="dot", line_color="#7f8c8d")
        fig.add_vline(x=-min_abs_log2fc, line_dash="dot", line_color="#7f8c8d")
    fig.update_layout(
        xaxis_title="log2 fold-change (tumor vs normal)",
        yaxis_title="-log10(padj)",
        legend_title="Estado",
        height=560,
        margin={"l": 40, "r": 20, "t": 30, "b": 40},
    )
    return fig


def correlation_network(
    processed: pd.DataFrame,
    genes: list[str],
    *,
    threshold: float = 0.6,
) -> nx.Graph:
    """Red proteína–proteína por correlación de Pearson entre las filas seleccionadas."""
    d = processed.copy()
    d["_g"] = d["protein_id"].astype(str).str.upper()
    sel = d[d["_g"].isin({g.upper() for g in genes})]
    if sel.empty:
        return nx.Graph()
    sample_cols = [c for c in sel.columns if c not in ("protein_id", "_g", "gene_symbol", "uniprot_accession")]
    mat = sel.set_index("_g")[sample_cols].astype(float)
    labels = list(mat.index)
    C = np.corrcoef(mat.to_numpy())
    G = nx.Graph()
    for g in labels:
        G.add_node(g)
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            if abs(C[i, j]) >= threshold:
                G.add_edge(labels[i], labels[j], weight=float(C[i, j]))
    return G


def network_figure(G: nx.Graph, scores: dict[str, float] | None = None) -> go.Figure:
    """Dibuja un grafo networkx con layout de resorte en Plotly (interactivo)."""
    if G.number_of_nodes() == 0:
        return go.Figure().update_layout(title="Sin nodos para mostrar")
    pos = nx.spring_layout(G, seed=42, k=0.6)
    edge_x, edge_y = [], []
    for u, v in G.edges():
        edge_x += [pos[u][0], pos[v][0], None]
        edge_y += [pos[u][1], pos[v][1], None]
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y, mode="lines",
        line={"width": 0.6, "color": "#bbb"}, hoverinfo="none", showlegend=False,
    )
    nx_, ny_, ntext, ncolor, nsize = [], [], [], [], []
    for node in G.nodes():
        nx_.append(pos[node][0])
        ny_.append(pos[node][1])
        deg = G.degree(node)
        sc = float(scores.get(node, np.nan)) if scores else np.nan
        ntext.append(f"{node}<br>grado={deg}" + (f"<br>score={sc:.2f}" if scores and not np.isnan(sc) else ""))
        ncolor.append(sc if scores else deg)
        nsize.append(8 + 2.2 * deg)
    node_trace = go.Scatter(
        x=nx_, y=ny_, mode="markers+text", text=list(G.nodes()), textposition="top center",
        hovertext=ntext, hoverinfo="text",
        marker={
            "size": nsize, "color": ncolor, "colorscale": "Viridis", "showscale": True,
            "colorbar": {"title": "score" if scores else "grado"},
            "line": {"width": 1, "color": "#333"},
        },
        showlegend=False,
    )
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        height=600, margin={"l": 10, "r": 10, "t": 30, "b": 10},
        xaxis={"visible": False}, yaxis={"visible": False},
    )
    return fig
