from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer
from sklearn.preprocessing import StandardScaler


def filter_missingness(
    expr: pd.DataFrame,
    meta_cols: list[str],
    *,
    protein_missing_max: float,
    sample_missing_max: float,
) -> pd.DataFrame:
    """Elimina proteínas y muestras con exceso de NA."""
    sample_cols = [c for c in expr.columns if c not in meta_cols]
    mat = expr[sample_cols]
    prot_frac = mat.isna().mean(axis=1)
    keep_prot = prot_frac <= protein_missing_max
    expr2 = expr.loc[keep_prot].copy()
    mat2 = expr2[sample_cols]
    samp_frac = mat2.isna().mean(axis=0)
    keep_samp = samp_frac <= sample_missing_max
    kept_samples = [c for c, k in zip(sample_cols, keep_samp) if k]
    return expr2[meta_cols + kept_samples]


def log2_transform(mat: pd.DataFrame, *, offset: float = 1e-6) -> pd.DataFrame:
    """log2 con offset (legado). Preferir :func:`log2_transform_masked` para intensidades no positivas."""
    x = mat.astype(float)
    return np.log2(x + offset)


def log2_transform_masked(mat: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """``log2(x)`` solo donde ``x > 0``; valores no positivos → NaN (antes de imputación KNN).

    Devuelve ``(matriz_log2, n_celdas_enmascaradas)``.
    """
    x = mat.astype(float)
    n_bad = int((x <= 0).sum().sum())
    safe = x.where(x > 0)
    out = np.log2(safe)
    return out, n_bad


def knn_impute(mat: pd.DataFrame, n_neighbors: int = 5) -> pd.DataFrame:
    imp = KNNImputer(n_neighbors=n_neighbors)
    arr = imp.fit_transform(mat.values)
    return pd.DataFrame(arr, index=mat.index, columns=mat.columns)


def zscore_by_protein(mat: pd.DataFrame) -> pd.DataFrame:
    scaler = StandardScaler(with_mean=True, with_std=True)
    arr = scaler.fit_transform(mat.T.values).T
    return pd.DataFrame(arr, index=mat.index, columns=mat.columns)
