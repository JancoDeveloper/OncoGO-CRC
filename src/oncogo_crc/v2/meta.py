"""Meta-análisis de efectos aleatorios por proteína (corrección C10).

Sustituye a ``metafor::rma`` porque no hay R instalado (enmienda E2). Implementa
el modelo de efectos aleatorios de un solo nivel:

    y_i = mu + u_i + e_i,    u_i ~ N(0, tau2),    e_i ~ N(0, v_i)

con ``v_i`` el error estándar al cuadrado del efecto en la cohorte ``i``, que
viene del contraste dentro de esa cohorte. Las matrices nunca se fusionan: solo
se combinan los efectos ya estimados.

``tau2`` se estima por REML con iteración de punto fijo, que es el mismo
estimador que ``metafor::rma(method = "REML")``. Se incluye DerSimonian–Laird
porque tiene forma cerrada y sirve de arranque y de comprobación.

Referencias: DerSimonian y Laird 1986; Viechtbauer 2005 (Journal of Educational
and Behavioral Statistics 30:261) para la forma de la ecuación REML; Higgins y
Thompson 2002 para I².
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.stats import chi2, norm

MAX_ITER = 200
TOL = 1e-10


@dataclass(frozen=True)
class ResultadoMeta:
    k: int
    efecto: float
    se: float
    ci_low: float
    ci_high: float
    z: float
    p: float
    tau2: float
    i2: float
    h2: float
    q: float
    q_df: int
    q_p: float
    metodo: str

    def to_dict(self) -> dict[str, float | int | str]:
        return asdict(self)


def _q_de_cochran(y: np.ndarray, v: np.ndarray) -> tuple[float, float]:
    """Q de Cochran con pesos de efectos fijos. Devuelve Q y la media fija."""
    w = 1.0 / v
    mu_fijo = float(np.sum(w * y) / np.sum(w))
    q = float(np.sum(w * (y - mu_fijo) ** 2))
    return q, mu_fijo


def tau2_dersimonian_laird(y: np.ndarray, v: np.ndarray) -> float:
    """Estimador de momentos, en forma cerrada. Nunca negativo."""
    k = y.size
    if k < 2:
        return 0.0
    w = 1.0 / v
    q, _ = _q_de_cochran(y, v)
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    if c <= 0:
        return 0.0
    return max(0.0, (q - (k - 1)) / c)


def tau2_reml(y: np.ndarray, v: np.ndarray, *, inicio: float | None = None) -> float:
    """REML por punto fijo.

    La ecuación de estimación de REML para este modelo es

        tau2 = sum(w_i^2 * ((y_i - mu)^2 + 1/sum(w) - v_i)) / sum(w_i^2)

    con ``w_i = 1 / (v_i + tau2)`` y ``mu`` la media ponderada con esos pesos.
    Iterar esa igualdad converge de forma monótona desde DerSimonian–Laird en los
    casos con los que se trabaja aquí (k pequeño, varianzas positivas).
    """
    k = y.size
    if k < 2:
        return 0.0
    tau2 = tau2_dersimonian_laird(y, v) if inicio is None else float(inicio)
    for _ in range(MAX_ITER):
        w = 1.0 / (v + tau2)
        sw = float(np.sum(w))
        mu = float(np.sum(w * y) / sw)
        w2 = w**2
        nuevo = float(np.sum(w2 * ((y - mu) ** 2 + 1.0 / sw - v)) / np.sum(w2))
        nuevo = max(0.0, nuevo)
        if abs(nuevo - tau2) < TOL:
            tau2 = nuevo
            break
        tau2 = nuevo
    return tau2


def meta_efectos_aleatorios(
    y: np.ndarray | list[float],
    se: np.ndarray | list[float],
    *,
    metodo: str = "reml",
) -> ResultadoMeta | None:
    """Combina efectos de varias cohortes para una sola proteína.

    ``y`` son los log2FC por cohorte y ``se`` sus errores estándar. Devuelve
    ``None`` si quedan menos de dos cohortes con error estándar positivo y
    finito: una sola cohorte no es un meta-análisis.
    """
    y = np.asarray(y, dtype=float)
    se = np.asarray(se, dtype=float)
    ok = np.isfinite(y) & np.isfinite(se) & (se > 0)
    y, se = y[ok], se[ok]
    k = int(y.size)
    if k < 2:
        return None

    v = se**2
    if metodo == "reml":
        tau2 = tau2_reml(y, v)
    elif metodo == "dl":
        tau2 = tau2_dersimonian_laird(y, v)
    elif metodo == "fijo":
        tau2 = 0.0
    else:
        raise ValueError(f"método de meta-análisis desconocido: {metodo!r}")

    w = 1.0 / (v + tau2)
    sw = float(np.sum(w))
    efecto = float(np.sum(w * y) / sw)
    se_efecto = float(np.sqrt(1.0 / sw))
    z = efecto / se_efecto if se_efecto > 0 else np.nan
    p = float(2.0 * norm.sf(abs(z))) if np.isfinite(z) else np.nan

    q, _ = _q_de_cochran(y, v)
    q_df = k - 1
    q_p = float(chi2.sf(q, q_df)) if q_df > 0 else np.nan
    # Higgins y Thompson: I² = (Q - df) / Q, acotado en cero.
    i2 = float(max(0.0, (q - q_df) / q)) if q > 0 else 0.0
    h2 = float(q / q_df) if q_df > 0 and q > 0 else 1.0

    crit = float(norm.ppf(0.975))
    return ResultadoMeta(
        k=k,
        efecto=efecto,
        se=se_efecto,
        ci_low=efecto - crit * se_efecto,
        ci_high=efecto + crit * se_efecto,
        z=float(z),
        p=p,
        tau2=float(tau2),
        i2=i2,
        h2=h2,
        q=q,
        q_df=q_df,
        q_p=q_p,
        metodo=metodo,
    )
