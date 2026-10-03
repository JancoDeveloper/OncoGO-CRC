"""t moderada y contraste treat sobre valores observados.

Sigue el estimador de momentos de la varianza a priori de limma (Smyth, 2004)
y el p-valor de ``treat`` (McCarthy y Smyth, 2009): se contrasta
|log2FC| > lfc, no un filtro aplicado después de un test contra cero.
No imputa ausencias. Si el prior no se puede estimar, d0 = 0 y el
resultado coincide con una t ordinaria con el mismo desplazamiento.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import digamma, polygamma
from scipy.optimize import brentq
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests


def trigamma_inverse(y: float) -> float:
    """z tal que trigamma(z) = y. trigamma decrece de +inf a 0."""
    if not np.isfinite(y) or y <= 0:
        raise ValueError("trigamma inversa solo está definida para y > 0")
    lower, upper = 1e-8, 1e8
    if polygamma(1, lower) < y:
        raise ValueError("la varianza empírica no admite un prior positivo")
    if polygamma(1, upper) > y:
        return upper
    return float(brentq(lambda z: float(polygamma(1, z) - y), lower, upper))


def fit_f_prior(s2: np.ndarray, df: np.ndarray) -> tuple[float, float]:
    """Grados de libertad y escala a priori (d0, s0²). d0 = 0 si no hay prior."""
    ok = np.isfinite(s2) & (s2 > 0) & np.isfinite(df) & (df > 0)
    if int(ok.sum()) < 2:
        return 0.0, float("nan")
    x = np.asarray(s2, dtype=float)[ok]
    d = np.asarray(df, dtype=float)[ok]
    e = np.log(x) - digamma(d / 2.0) + np.log(d / 2.0)
    n = int(e.size)
    emean = float(np.mean(e))
    evar = float(np.sum((e - emean) ** 2) / (n - 1))
    if not np.isfinite(evar) or evar <= 0:
        return 0.0, float("nan")
    try:
        z = trigamma_inverse(evar)
    except ValueError:
        return 0.0, float("nan")
    d0 = float(2.0 * z)
    if not np.isfinite(d0) or d0 <= 0:
        return 0.0, float("nan")
    s02 = float(np.exp(emean + digamma(d0 / 2.0) - np.log(d0 / 2.0)))
    if not np.isfinite(s02) or s02 <= 0:
        return 0.0, float("nan")
    return d0, s02


def posterior_scale(s2: np.ndarray, df: np.ndarray, d0: float, s02: float) -> tuple[np.ndarray, np.ndarray]:
    """s² posterior y grados de libertad totales."""
    s2 = np.asarray(s2, dtype=float)
    df = np.asarray(df, dtype=float)
    if not np.isfinite(d0) or d0 <= 0 or not np.isfinite(s02):
        return s2.copy(), df.copy()
    post = (d0 * s02 + df * s2) / (d0 + df)
    return post, df + d0


def treat_pvalues(coef: np.ndarray, se: np.ndarray, df_total: np.ndarray, lfc: float) -> np.ndarray:
    """p de |coef| > lfc. El mínimo de las dos colas desplazadas es el p-valor válido en el borde |coef| = lfc."""
    coef = np.asarray(coef, dtype=float)
    se = np.asarray(se, dtype=float)
    df_total = np.asarray(df_total, dtype=float)
    threshold = abs(float(lfc))
    with np.errstate(divide="ignore", invalid="ignore"):
        t_right = (coef - threshold) / se
        t_left = (coef + threshold) / se
        p_right = student_t.sf(t_right, df_total)
        p_left = student_t.cdf(t_left, df_total)
        p = np.minimum(p_right, p_left)
    bad = ~(np.isfinite(coef) & np.isfinite(se) & (se > 0) & np.isfinite(df_total) & (df_total > 0))
    p = np.where(bad, np.nan, p)
    return p


def _bh(p: np.ndarray, alpha: float) -> tuple[np.ndarray, np.ndarray]:
    padj = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        _, q, _, _ = multipletests(p[ok], alpha=alpha, method="fdr_bh")
        padj[ok] = q
    sig = padj < alpha
    return padj, sig


def moderate_table(
    *,
    protein_ids: np.ndarray,
    coef: np.ndarray,
    s2: np.ndarray,
    df: np.ndarray,
    stdev_unscaled: np.ndarray,
    n_obs: np.ndarray,
    lfc: float,
    alpha: float,
    stat_method: str,
) -> tuple[pd.DataFrame, float, float]:
    """Tabla por proteína. t_stat es la t moderada contra 0; p_value es el de treat."""
    d0, s02 = fit_f_prior(s2, df)
    s2_post, df_total = posterior_scale(s2, df, d0, s02)
    with np.errstate(divide="ignore", invalid="ignore"):
        se = stdev_unscaled * np.sqrt(s2_post)
        t_stat = coef / se
    p = treat_pvalues(coef, se, df_total, lfc)
    padj, sig = _bh(p, alpha)
    out = pd.DataFrame(
        {
            "protein_id": protein_ids,
            "log2fc": coef,
            "t_stat": t_stat,
            "se": se,
            "df_total": df_total,
            "n_obs": n_obs,
            "p_value": p,
            "padj": padj,
            "significant": sig,
            "stat_method": stat_method,
        }
    )
    return out, d0, s02


def paired_observed_stats(diffs: np.ndarray, *, min_pairs: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """diffs: proteínas × pares, ya en escala log2. Un par con ausencia no se imputa."""
    n = np.sum(np.isfinite(diffs), axis=1).astype(float)
    coef = np.nanmean(diffs, axis=1)
    s2 = np.nanvar(diffs, axis=1, ddof=1)
    df = n - 1.0
    stdev = 1.0 / np.sqrt(n)
    too_few = n < float(min_pairs)
    coef = np.where(too_few, np.nan, coef)
    s2 = np.where(too_few, np.nan, s2)
    df = np.where(too_few, np.nan, df)
    stdev = np.where(too_few, np.nan, stdev)
    return coef, s2, df, stdev, n


def paired_covariate_stats(
    diffs: np.ndarray,
    covariate: np.ndarray,
    *,
    min_pairs: int,
    at: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Efecto emparejado evaluado en ``covariate = at``.

    Modelo por proteína: ``diff = a + b * x + e``. El efecto reportado es la
    predicción en ``x = at``, no la media de las diferencias: con ``at = 0`` es
    el efecto que quedaría si el par no tuviera diferencia de composición.
    Evaluar en la media de ``x`` devolvería la media de ``y`` sin cambio alguno,
    que no responde a si el efecto es composición de tejido.

    Devuelve además la pendiente y su t para saber qué proteínas dependen de
    la covariable.
    """
    n_prot = diffs.shape[0]
    coef = np.full(n_prot, np.nan)
    s2 = np.full(n_prot, np.nan)
    df = np.full(n_prot, np.nan)
    stdev = np.full(n_prot, np.nan)
    n_used = np.zeros(n_prot)
    slope = np.full(n_prot, np.nan)
    slope_t = np.full(n_prot, np.nan)
    x_all = np.asarray(covariate, dtype=float)
    for i in range(n_prot):
        row = diffs[i]
        ok = np.isfinite(row) & np.isfinite(x_all)
        n = int(ok.sum())
        n_used[i] = n
        # Dos parámetros: hacen falta al menos min_pairs y un grado de libertad residual.
        if n < max(int(min_pairs), 3):
            continue
        y = row[ok]
        xr = x_all[ok]
        xbar = float(xr.mean())
        x = xr - xbar
        sxx = float(np.dot(x, x))
        if sxx <= 0:
            continue
        b = float(np.dot(x, y) / sxx)
        ybar = float(y.mean())
        resid = y - (ybar + b * x)
        dof = n - 2
        var = float(np.dot(resid, resid) / dof)
        coef[i] = ybar + b * (float(at) - xbar)
        s2[i] = var
        df[i] = dof
        stdev[i] = float(np.sqrt(1.0 / n + (float(at) - xbar) ** 2 / sxx))
        slope[i] = b
        if var > 0:
            slope_t[i] = b / float(np.sqrt(var / sxx))
    return coef, s2, df, stdev, n_used, slope, slope_t


def unpaired_observed_stats(
    group_a: np.ndarray,
    group_b: np.ndarray,
    *,
    min_n: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Dos grupos, varianza común (el modelo de limma), solo celdas observadas."""
    na = np.sum(np.isfinite(group_a), axis=1).astype(float)
    nb = np.sum(np.isfinite(group_b), axis=1).astype(float)
    ma = np.nanmean(group_a, axis=1)
    mb = np.nanmean(group_b, axis=1)
    va = np.nanvar(group_a, axis=1, ddof=1)
    vb = np.nanvar(group_b, axis=1, ddof=1)
    df = na + nb - 2.0
    with np.errstate(divide="ignore", invalid="ignore"):
        s2 = ((na - 1.0) * va + (nb - 1.0) * vb) / df
        stdev = np.sqrt(1.0 / na + 1.0 / nb)
    coef = ma - mb
    n = na + nb
    too_few = (na < float(min_n)) | (nb < float(min_n))
    coef = np.where(too_few, np.nan, coef)
    s2 = np.where(too_few, np.nan, s2)
    df = np.where(too_few, np.nan, df)
    stdev = np.where(too_few, np.nan, stdev)
    return coef, s2, df, stdev, n


def n_significant_after_signs(
    diffs: np.ndarray,
    signs: np.ndarray,
    *,
    s2_post: np.ndarray,
    df_total: np.ndarray,
    stdev_unscaled: np.ndarray,
    lfc: float,
    alpha: float,
) -> int:
    """Cuántas proteínas pasan treat+BH si se invierte el signo de algunos pares.

    La varianza de una diferencia no cambia al cambiarle el signo, así que el
    prior se calcula una vez sobre los datos observados.
    """
    coef = np.nanmean(diffs * signs.reshape(1, -1), axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        se = stdev_unscaled * np.sqrt(s2_post)
    p = treat_pvalues(coef, se, df_total, lfc)
    _, sig = _bh(p, alpha)
    return int(sig.sum())
