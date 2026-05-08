from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def differential_ttest(
    mat: pd.DataFrame,
    labels: pd.Series,
    *,
    group_a: str,
    group_b: str,
    alpha: float = 0.05,
    fdr_method: str = "fdr_bh",
) -> pd.DataFrame:
    """t-test por fila entre dos grupos; ajuste FDR."""
    cols = mat.columns.intersection(labels.index)
    mat = mat[cols]
    lab = labels.loc[cols]
    idx_a = lab[lab == group_a].index
    idx_b = lab[lab == group_b].index
    rows = []
    for name, row in mat.iterrows():
        a = row[idx_a].dropna().astype(float)
        b = row[idx_b].dropna().astype(float)
        if len(a) < 2 or len(b) < 2:
            rows.append((name, np.nan, np.nan, np.nan))
            continue
        t_stat, p_val = stats.ttest_ind(a, b, equal_var=False)
        fc = float(a.mean() - b.mean())
        rows.append((name, fc, float(t_stat), float(p_val)))
    res = pd.DataFrame(rows, columns=["protein_id", "log2fc", "t_stat", "p_value"])
    ok = res["p_value"].notna()
    adj = np.full(len(res), np.nan)
    if ok.any():
        _, q, _, _ = multipletests(res.loc[ok, "p_value"], alpha=alpha, method=fdr_method)
        adj[ok.to_numpy()] = q
    res["padj"] = adj
    res["significant"] = res["padj"] < alpha
    return res
