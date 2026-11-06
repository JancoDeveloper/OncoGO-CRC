"""Concordancia de signo, Spearman y concordancia de Lin.

La de signo usa el mismo intervalo de Clopper–Pearson que v1. Spearman y Lin
se calculan sobre los mismos pares: un Pearson solo, y solo en las proteínas
ya elegidas, comprime la lectura.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from oncogo_crc.v2.metricas import Fraccion


def pares_evaluables(efecto_a: np.ndarray, efecto_b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    a = np.asarray(efecto_a, dtype=float)
    b = np.asarray(efecto_b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b) & (a != 0) & (b != 0)
    return a[ok], b[ok]


def concordancia_de_signo(
    efecto_a: np.ndarray,
    efecto_b: np.ndarray,
    *,
    nombre: str,
    universo: int,
    universo_nombre: str,
) -> Fraccion:
    a, b = pares_evaluables(efecto_a, efecto_b)
    mismo = int(np.sum(np.sign(a) == np.sign(b))) if a.size else 0
    return Fraccion(
        nombre=nombre,
        numerador=mismo,
        denominador=int(a.size),
        unidad_denominador="proteinas_unicas",
        universo=universo,
        universo_nombre=universo_nombre,
    )


def spearman(efecto_a: np.ndarray, efecto_b: np.ndarray) -> dict[str, float | int | None]:
    a, b = pares_evaluables(efecto_a, efecto_b)
    if a.size < 3:
        return {"n": int(a.size), "rho": None, "p": None}
    rho, p = spearmanr(a, b)
    return {"n": int(a.size), "rho": float(rho), "p": float(p)}


def concordancia_de_lin(efecto_a: np.ndarray, efecto_b: np.ndarray) -> dict[str, float | int | None]:
    """Coeficiente de concordancia de Lin (1989)."""
    a, b = pares_evaluables(efecto_a, efecto_b)
    n = int(a.size)
    if n < 3:
        return {"n": n, "ccc": None}
    media_a = float(a.mean())
    media_b = float(b.mean())
    var_a = float(a.var(ddof=1))
    var_b = float(b.var(ddof=1))
    cov = float(np.cov(a, b, ddof=1)[0, 1])
    denominador = var_a + var_b + (media_a - media_b) ** 2
    if denominador <= 0:
        return {"n": n, "ccc": None}
    return {"n": n, "ccc": float((2.0 * cov) / denominador)}


def alinear_por_clave(
    descubrimiento: pd.DataFrame,
    replica: pd.DataFrame,
    *,
    clave_descubrimiento: pd.Series,
    clave_replica: pd.Series,
) -> pd.DataFrame:
    """Una fila por clave compartida. Las claves vacías no se cruzan."""
    izq = descubrimiento.copy()
    der = replica.copy()
    izq["_clave"] = clave_descubrimiento.to_numpy()
    der["_clave"] = clave_replica.to_numpy()
    izq = izq.loc[izq["_clave"].notna() & (izq["_clave"].astype(str) != "")]
    der = der.loc[der["_clave"].notna() & (der["_clave"].astype(str) != "")]
    izq = izq.drop_duplicates("_clave", keep="first")
    der = der.drop_duplicates("_clave", keep="first")
    return izq.merge(der, on="_clave", how="inner", suffixes=("_cptac", "_replica"))
