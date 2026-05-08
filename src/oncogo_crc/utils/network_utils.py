from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd


def correlation_network(
    mat: pd.DataFrame,
    *,
    threshold: float = 0.6,
) -> nx.Graph:
    """
    Grafo no dirigido proteína–proteína por correlación de Pearson entre filas (proteínas).
    `mat`: filas = proteínas, columnas = muestras.
    """
    idx = list(mat.index)
    C = np.corrcoef(mat.values.astype(float))
    G = nx.Graph()
    for p in idx:
        G.add_node(p)
    n = len(idx)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(C[i, j]) >= threshold:
                G.add_edge(idx[i], idx[j], weight=float(C[i, j]))
    return G


def hub_proteins(G: nx.Graph, *, min_degree: int = 5) -> pd.DataFrame:
    deg = dict(G.degree())
    rows = [{"protein_id": n, "degree": d, "is_hub": d >= min_degree} for n, d in deg.items()]
    df = pd.DataFrame(rows).sort_values("degree", ascending=False)
    return df
