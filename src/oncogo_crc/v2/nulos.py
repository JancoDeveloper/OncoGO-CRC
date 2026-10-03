"""Nulos de permutación con el número de permutaciones que el p necesita (C1).

v1 usó 100 permutaciones y reportó p = 0,0099, que es exactamente 1/101: el
mínimo alcanzable. Un p que coincide con el piso no es una estimación, es un
tope. Aquí hay dos remedios.

En CPTAC, 10.000 permutaciones del signo de los 94 pares. Eso baja el piso a
1/10001 y permite un intervalo de Clopper–Pearson con sentido.

En la réplica, 16 pares significan 2^16 = 65.536 inversiones de signo posibles.
Ese espacio se enumera completo: el p deja de ser una estimación y pasa a ser
exacto.

Para que 65.536 permutaciones sobre miles de proteínas sean viables, el cálculo
se hace en forma cerrada en vez de recomputar la tabla. Invertir el signo de un
par no cambia qué pares están observados ni la suma de cuadrados, así que para
la permutación con vector de signos ``s``:

    media_i = (sum_j s_j * d_ij) / n_i
    s2_i    = (sum_j d_ij^2 - n_i * media_i^2) / (n_i - 1)

porque ``s_j^2 = 1``. Las dos salen de productos de matrices, y el prior de la
varianza se estima una sola vez sobre los datos observados, igual que en v1.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from statsmodels.stats.multitest import multipletests

from oncogo_crc.utils.moderated import (
    fit_f_prior,
    paired_observed_stats,
    posterior_scale,
    treat_pvalues,
)

# Permutaciones por bloque. Acota la memoria del producto de matrices.
BLOQUE = 512


@dataclass
class Preparado:
    """Lo que no cambia al invertir signos, precalculado una vez."""

    d0_llenado: np.ndarray      # diferencias con las ausencias en cero
    n_obs: np.ndarray           # pares observados por proteína
    suma_cuadrados: np.ndarray  # sum_j d_ij^2 sobre los observados
    evaluable: np.ndarray       # proteínas con pares suficientes
    prior_df: float
    prior_s02: float
    n_pares: int
    min_pares: int


def preparar(diffs: np.ndarray, *, min_pares: int) -> Preparado:
    """Precalcula el prior de la varianza sobre los datos observados.

    El prior se estima una vez y se reutiliza en todas las permutaciones. Es lo
    correcto: invertir el signo de un par no cambia la varianza de la proteína,
    así que la distribución de varianzas del nulo es la observada.
    """
    diffs = np.asarray(diffs, dtype=float)
    finito = np.isfinite(diffs)
    llenado = np.where(finito, diffs, 0.0)
    n_obs = finito.sum(axis=1).astype(float)
    suma_cuadrados = (llenado**2).sum(axis=1)
    _coef, s2, df, _stdev, _n = paired_observed_stats(diffs, min_pairs=min_pares)
    prior_df, prior_s02 = fit_f_prior(s2, df)
    return Preparado(
        d0_llenado=llenado,
        n_obs=n_obs,
        suma_cuadrados=suma_cuadrados,
        evaluable=n_obs >= float(min_pares),
        prior_df=float(prior_df),
        prior_s02=float(prior_s02),
        n_pares=int(diffs.shape[1]),
        min_pares=int(min_pares),
    )


def _conteos_del_bloque(
    prep: Preparado,
    signos: np.ndarray,
    *,
    lfc: float,
    alpha: float,
) -> np.ndarray:
    """Número de proteínas significativas para cada vector de signos del bloque.

    ``signos`` tiene forma (b, n_pares) con valores en {-1, +1}.
    """
    n = prep.n_obs
    evaluable = prep.evaluable
    # medias: proteínas x b
    medias = (prep.d0_llenado @ signos.T) / np.where(n[:, None] > 0, n[:, None], np.nan)
    df = n - 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        s2 = (prep.suma_cuadrados[:, None] - n[:, None] * medias**2) / df[:, None]
    # Una varianza negativa por error de redondeo se recorta a cero.
    s2 = np.where(s2 < 0, 0.0, s2)
    s2_post, df_total = posterior_scale(
        s2, np.repeat(df[:, None], signos.shape[0], axis=1), prep.prior_df, prep.prior_s02
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.sqrt(s2_post) / np.sqrt(n[:, None])
    p = treat_pvalues(medias, se, df_total, lfc)
    p = np.where(evaluable[:, None], p, np.nan)

    conteos = np.zeros(signos.shape[0], dtype=int)
    for j in range(signos.shape[0]):
        col = p[:, j]
        ok = np.isfinite(col)
        if not ok.any():
            continue
        _rechaza, q, _a, _b = multipletests(col[ok], alpha=alpha, method="fdr_bh")
        conteos[j] = int((q < alpha).sum())
    return conteos


def _observado(prep: Preparado, *, lfc: float, alpha: float) -> int:
    unos = np.ones((1, prep.n_pares))
    return int(_conteos_del_bloque(prep, unos, lfc=lfc, alpha=alpha)[0])


def nulo_aleatorio(
    diffs: np.ndarray,
    *,
    lfc: float,
    alpha: float,
    min_pares: int,
    n_permutaciones: int,
    semilla: int,
) -> dict:
    """Invierte al azar el signo de cada par, ``n_permutaciones`` veces."""
    prep = preparar(diffs, min_pares=min_pares)
    observado = _observado(prep, lfc=lfc, alpha=alpha)
    rng = np.random.default_rng(semilla)
    conteos: list[int] = []
    restantes = int(n_permutaciones)
    while restantes > 0:
        b = min(BLOQUE, restantes)
        signos = rng.choice(np.array([-1.0, 1.0]), size=(b, prep.n_pares))
        conteos.extend(_conteos_del_bloque(prep, signos, lfc=lfc, alpha=alpha).tolist())
        restantes -= b
    nulo = np.asarray(conteos, dtype=int)
    return {
        "tipo": "aleatorio",
        "observado": observado,
        "n_permutaciones": int(n_permutaciones),
        "b_al_menos_observado": int((nulo >= observado).sum()),
        "nulo": nulo,
        "nulo_mediana": float(np.median(nulo)),
        "nulo_p95": float(np.percentile(nulo, 95)),
        "nulo_max": int(nulo.max()) if nulo.size else 0,
        "semilla": int(semilla),
        "prior_df": prep.prior_df,
    }


def nulo_exacto(
    diffs: np.ndarray,
    *,
    lfc: float,
    alpha: float,
    min_pares: int,
    max_pares: int = 20,
) -> dict:
    """Enumera las 2^k inversiones de signo. Con 16 pares son 65.536.

    El límite de ``max_pares`` evita que alguien lance 2^94 por accidente.
    """
    k = int(np.asarray(diffs).shape[1])
    if k > max_pares:
        raise ValueError(
            f"{k} pares son 2^{k} inversiones: demasiadas para enumerar. Use nulo_aleatorio."
        )
    prep = preparar(diffs, min_pares=min_pares)
    observado = _observado(prep, lfc=lfc, alpha=alpha)
    total = 2**k
    conteos = np.zeros(total, dtype=int)
    idx = 0
    buffer: list[tuple[float, ...]] = []
    for combo in product((-1.0, 1.0), repeat=k):
        buffer.append(combo)
        if len(buffer) == BLOQUE:
            signos = np.asarray(buffer)
            conteos[idx : idx + len(buffer)] = _conteos_del_bloque(
                prep, signos, lfc=lfc, alpha=alpha
            )
            idx += len(buffer)
            buffer = []
    if buffer:
        signos = np.asarray(buffer)
        conteos[idx : idx + len(buffer)] = _conteos_del_bloque(prep, signos, lfc=lfc, alpha=alpha)
    b = int((conteos >= observado).sum())
    # En una enumeración exhaustiva, el vector de signos sin invertir está dentro
    # del espacio y empata con lo observado por definición. Su complemento también,
    # porque invertir todos los signos no cambia el valor absoluto del estadístico.
    # Así que b nunca baja de 2 y el p mínimo alcanzable es 3/(B+1), no 1/(B+1).
    b_triviales = 2
    return {
        "tipo": "exacto",
        "observado": observado,
        "n_permutaciones": total,
        "b_al_menos_observado": b,
        "b_triviales": b_triviales,
        "b_no_triviales": max(0, b - b_triviales),
        "p_minimo_alcanzable_exacto": (b_triviales + 1) / (total + 1),
        "es_minimo_alcanzable": b <= b_triviales,
        "nulo": conteos,
        "nulo_mediana": float(np.median(conteos)),
        "nulo_p95": float(np.percentile(conteos, 95)),
        "nulo_max": int(conteos.max()),
        "semilla": None,
        "prior_df": prep.prior_df,
        "nota": (
            f"Enumeración completa de 2^{k} = {total} inversiones de signo. El p es "
            "exacto, no una estimación. La identidad y su complemento siempre empatan "
            f"con lo observado, así que el p no puede bajar de {(b_triviales + 1) / (total + 1):.3e}."
        ),
    }
