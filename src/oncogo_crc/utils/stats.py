from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def build_paired_sample_order(
    sm: pd.DataFrame,
    *,
    condition_column: str,
    group_a: str,
    group_b: str,
    participant_col: str,
) -> list[tuple[str, str]] | None:
    """Pares (sample_id grupo A, sample_id grupo B) con exactamente un par por participante.

    Devuelve None si faltan columnas o no hay al menos 2 pares completos.
    """
    need = {"sample_id", condition_column, participant_col}
    if not need.issubset(sm.columns):
        return None
    pairs: list[tuple[str, str]] = []
    for _, g in sm.groupby(participant_col):
        sa = g.loc[g[condition_column].astype(str) == str(group_a), "sample_id"].astype(str).tolist()
        sb = g.loc[g[condition_column].astype(str) == str(group_b), "sample_id"].astype(str).tolist()
        if len(sa) == 1 and len(sb) == 1:
            pairs.append((sa[0], sb[0]))
    return pairs if len(pairs) >= 2 else None


def differential_ttest(
    mat: pd.DataFrame,
    labels: pd.Series,
    *,
    group_a: str,
    group_b: str,
    alpha: float = 0.05,
    fdr_method: str = "fdr_bh",
    min_abs_log2fc: float = 0.0,
    paired_sample_order: list[tuple[str, str]] | None = None,
) -> pd.DataFrame:
    """Contraste por fila entre dos grupos sobre matriz **log2** (Δ de medias ≈ log-FC).

    - Sin pares: t-test de Welch entre grupos (``equal_var=False``).
    - Con ``paired_sample_order``: lista de ``(id_muestra_A, id_muestra_B)``; estadístico = t
      sobre las diferencias emparejadas :math:`x_B - x_A` (``ttest_1samp`` vs 0).
    """
    cols = mat.columns.intersection(labels.index)
    mat = mat[cols]
    lab = labels.loc[cols]
    idx_a = lab[lab == group_a].index
    idx_b = lab[lab == group_b].index
    rows = []
    method = "paired_ttest_mean_diff" if paired_sample_order else "welch_ttest"

    if paired_sample_order:
        for name, row in mat.iterrows():
            diffs: list[float] = []
            for sid_a, sid_b in paired_sample_order:
                if sid_a not in row.index or sid_b not in row.index:
                    continue
                va, vb = row[sid_a], row[sid_b]
                if pd.isna(va) or pd.isna(vb):
                    continue
                diffs.append(float(va) - float(vb))
            if len(diffs) < 2:
                rows.append((name, np.nan, np.nan, np.nan, method))
                continue
            arr = np.asarray(diffs, dtype=float)
            t_stat, p_val = stats.ttest_1samp(arr, popmean=0.0, nan_policy="omit")
            fc = float(np.mean(arr))
            rows.append((name, fc, float(t_stat), float(p_val), method))
    else:
        for name, row in mat.iterrows():
            a = row[idx_a].dropna().astype(float)
            b = row[idx_b].dropna().astype(float)
            if len(a) < 2 or len(b) < 2:
                rows.append((name, np.nan, np.nan, np.nan, method))
                continue
            t_stat, p_val = stats.ttest_ind(a, b, equal_var=False)
            fc = float(a.mean() - b.mean())
            rows.append((name, fc, float(t_stat), float(p_val), method))
    res = pd.DataFrame(rows, columns=["protein_id", "log2fc", "t_stat", "p_value", "stat_method"])
    ok = res["p_value"].notna()
    adj = np.full(len(res), np.nan)
    if ok.any():
        _, q, _, _ = multipletests(res.loc[ok, "p_value"], alpha=alpha, method=fdr_method)
        adj[ok.to_numpy()] = q
    res["padj"] = adj
    sig = res["padj"] < alpha
    mfc = float(min_abs_log2fc or 0.0)
    if mfc > 0 and "log2fc" in res.columns:
        fc_ok = res["log2fc"].notna() & (res["log2fc"].abs() >= mfc)
        sig = sig & fc_ok
    res["significant"] = sig
    return res
