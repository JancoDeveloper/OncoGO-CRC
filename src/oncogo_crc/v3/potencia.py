"""TOST, alfa efectivo de BH por simulación y efecto mínimo detectable.

«No hay señal» solo corresponde a la clase equivalente. Una proteína no
significativa con un efecto mínimo detectable grande es indeterminada.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import t as student_t

from oncogo_crc.utils.moderated import treat_pvalues


def p_tost(coef: np.ndarray, se: np.ndarray, df: np.ndarray, margen: float) -> np.ndarray:
    """p de Schuirmann: el máximo de las dos colas con márgenes ±margen."""
    coef = np.asarray(coef, dtype=float)
    se = np.asarray(se, dtype=float)
    df = np.asarray(df, dtype=float)
    limite = abs(float(margen))
    with np.errstate(divide="ignore", invalid="ignore"):
        p_inferior = student_t.sf((coef + limite) / se, df)
        p_superior = student_t.cdf((coef - limite) / se, df)
        p = np.maximum(p_inferior, p_superior)
    mal = ~(np.isfinite(coef) & np.isfinite(se) & (se > 0) & np.isfinite(df) & (df > 0))
    return np.where(mal, np.nan, p)


def umbral_bh(p: np.ndarray, alpha: float) -> float:
    """Mayor p rechazado por BH. Si no hay rechazos, alpha / m."""
    valores = np.asarray(p, dtype=float)
    finitos = valores[np.isfinite(valores)]
    m = int(finitos.size)
    if m == 0:
        return float(alpha)
    ordenados = np.sort(finitos)
    limite = alpha * (np.arange(1, m + 1) / m)
    debajo = np.flatnonzero(ordenados <= limite)
    if debajo.size == 0:
        return float(alpha / m)
    return float(ordenados[int(debajo[-1])])


def alfa_efectivo(
    se: np.ndarray,
    df: np.ndarray,
    *,
    alpha: float,
    lfc: float,
    n_sim: int,
    rng: np.random.Generator,
) -> float:
    """Mediana del umbral BH al simular el catálogo con efecto cero."""
    se = np.asarray(se, dtype=float)
    df = np.asarray(df, dtype=float)
    ok = np.isfinite(se) & (se > 0) & np.isfinite(df) & (df > 0)
    se_ok = se[ok]
    df_ok = df[ok]
    if se_ok.size == 0:
        return float(alpha)
    t = rng.standard_t(df_ok, size=(n_sim, se_ok.size))
    coef = t * se_ok
    umbrales = np.empty(n_sim, dtype=float)
    for i in range(n_sim):
        p = treat_pvalues(coef[i], se_ok, df_ok, lfc)
        umbrales[i] = umbral_bh(p, alpha)
    return float(np.median(umbrales))


def efecto_minimo_detectable(
    se: np.ndarray,
    df: np.ndarray,
    *,
    alpha_eff: float,
    lfc: float,
    potencia: float,
    delta_max: float,
    n_sim: int,
    n_delta: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Menor |log2FC| de la rejilla con potencia al menos la pedida.

    NaN si la proteína no tiene error estándar. Inf si ni el extremo de la
    rejilla llega a la potencia.
    """
    se = np.asarray(se, dtype=float)
    df = np.asarray(df, dtype=float)
    salida = np.full(se.shape, np.nan)
    ok = np.isfinite(se) & (se > 0) & np.isfinite(df) & (df > 0)
    if not np.any(ok):
        return salida
    deltas = np.linspace(float(lfc) + 1e-3, float(delta_max), int(n_delta))
    se_ok = se[ok]
    df_ok = df[ok]
    alcanza = np.zeros(se_ok.size, dtype=bool)
    mejor = np.full(se_ok.size, np.inf)
    for delta in deltas:
        t = rng.standard_t(df_ok, size=(n_sim, se_ok.size))
        p = treat_pvalues(delta + se_ok * t, np.broadcast_to(se_ok, t.shape), np.broadcast_to(df_ok, t.shape), lfc)
        poder = np.mean(p < alpha_eff, axis=0)
        nuevos = (~alcanza) & (poder >= potencia)
        mejor[nuevos] = delta
        alcanza |= poder >= potencia
    salida[ok] = mejor
    return salida


def clasificar(treat_significativo: np.ndarray, tost_significativo: np.ndarray) -> np.ndarray:
    """diferente, equivalente, conflicto o indeterminada. Una clase por proteína."""
    treat = np.asarray(treat_significativo, dtype=bool)
    tost = np.asarray(tost_significativo, dtype=bool)
    salida = np.full(treat.shape, "indeterminada", dtype=object)
    salida[treat & tost] = "conflicto"
    salida[treat & ~tost] = "diferente"
    salida[~treat & tost] = "equivalente"
    return salida
