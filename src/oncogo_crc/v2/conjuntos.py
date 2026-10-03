"""Prueba competitiva de conjuntos de proteínas (correcciones C4 y C13).

Sustituye a ``limma::camera`` porque no hay R instalado (enmienda E2).

Por qué competitiva y no autocontenida. Una prueba autocontenida pregunta si el
conjunto cambia; con 688 proteínas significativas de 7394, casi cualquier
conjunto grande cambia. Una prueba competitiva pregunta si el conjunto cambia
**más que el resto de proteínas medidas**, que es la pregunta del enriquecimiento.

Por qué hay que corregir por correlación. Las proteínas de una misma vía no son
independientes: suben y bajan juntas. Si se ignora, el error estándar del
estadístico del conjunto se subestima y salen vías significativas por montones.
``camera`` corrige con un factor de inflación de la varianza que depende de la
correlación intergénica media:

    VIF = 1 + (m - 1) * rho

con ``m`` el tamaño del conjunto. Aquí ``rho`` se estima de los residuos de la
matriz, no se fija a mano.

Referencia: Wu y Smyth 2012, Nucleic Acids Research 40:e133.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import rankdata, t as t_dist


@dataclass(frozen=True)
class ResultadoConjunto:
    termino: str
    n_en_conjunto: int
    n_fuera: int
    estadistico: float
    direccion: int
    p_value: float
    vif: float
    rho: float
    metodo: str


def correlacion_intergenica(residuos: np.ndarray, *, max_proteinas: int = 2000) -> float:
    """Correlación media entre proteínas, estimada de los residuos del modelo.

    Toma la media de las correlaciones por pares fuera de la diagonal. Con muchas
    proteínas se submuestrea de forma determinista (las primeras ``max_proteinas``
    por varianza) para no construir una matriz enorme. Un valor negativo se
    recorta a cero: un VIF menor que uno haría la prueba anticonservadora.
    """
    x = np.asarray(residuos, dtype=float)
    if x.ndim != 2 or x.shape[0] < 2 or x.shape[1] < 3:
        return 0.0
    # Solo filas con varianza; una fila constante no aporta correlación.
    var = np.nanvar(x, axis=1)
    keep = np.isfinite(var) & (var > 0)
    x = x[keep]
    if x.shape[0] < 2:
        return 0.0
    if x.shape[0] > max_proteinas:
        orden = np.argsort(-np.nanvar(x, axis=1), kind="stable")[:max_proteinas]
        x = x[np.sort(orden)]
    x = x - np.nanmean(x, axis=1, keepdims=True)
    sd = np.sqrt(np.nanmean(x**2, axis=1, keepdims=True))
    sd[sd == 0] = np.nan
    z = x / sd
    z = np.nan_to_num(z, nan=0.0)
    n = z.shape[0]
    corr = (z @ z.T) / z.shape[1]
    suma = float(corr.sum() - np.trace(corr))
    pares = n * (n - 1)
    rho = suma / pares if pares > 0 else 0.0
    return float(max(0.0, rho))


def camera_un_conjunto(
    estadisticos: np.ndarray,
    en_conjunto: np.ndarray,
    *,
    termino: str,
    rho: float,
    df: float,
) -> ResultadoConjunto | None:
    """Compara el conjunto contra el resto con el VIF de camera.

    ``estadisticos`` son los t moderados por proteína. El contraste es la
    diferencia de medias entre el conjunto y el resto, con la varianza del
    conjunto inflada por ``1 + (m - 1) * rho``.
    """
    stats = np.asarray(estadisticos, dtype=float)
    sel = np.asarray(en_conjunto, dtype=bool)
    ok = np.isfinite(stats)
    stats, sel = stats[ok], sel[ok]
    m = int(sel.sum())
    n_fuera = int((~sel).sum())
    if m < 2 or n_fuera < 2:
        return None

    vif = 1.0 + (m - 1) * rho
    media_dentro = float(stats[sel].mean())
    media_fuera = float(stats[~sel].mean())
    n = stats.size
    # Varianza común, como en camera: una sola estimación para los dos grupos.
    var_comun = float(((stats - stats.mean()) ** 2).sum() / (n - 1))
    if var_comun <= 0:
        return None
    se = float(np.sqrt(var_comun * (vif / m + 1.0 / n_fuera)))
    if se <= 0:
        return None
    estadistico = (media_dentro - media_fuera) / se
    df_total = max(1.0, float(df))
    p = float(2.0 * t_dist.sf(abs(estadistico), df_total))
    return ResultadoConjunto(
        termino=termino,
        n_en_conjunto=m,
        n_fuera=n_fuera,
        estadistico=float(estadistico),
        direccion=int(np.sign(estadistico)),
        p_value=p,
        vif=vif,
        rho=float(rho),
        metodo="competitiva_camera",
    )


def camera(
    estadisticos: np.ndarray,
    pertenencia: dict[str, np.ndarray],
    *,
    rho: float,
    df: float,
    min_proteinas: int,
) -> list[ResultadoConjunto]:
    """Aplica la prueba a cada término. Descarta los que no llegan al mínimo."""
    salida = []
    for termino, sel in pertenencia.items():
        if int(np.asarray(sel, dtype=bool).sum()) < min_proteinas:
            continue
        res = camera_un_conjunto(estadisticos, sel, termino=termino, rho=rho, df=df)
        if res is not None:
            salida.append(res)
    return salida


def nulo_etiquetas_de_gen(
    estadisticos: np.ndarray,
    pertenencia: dict[str, np.ndarray],
    *,
    rho: float,
    df: float,
    min_proteinas: int,
    n_permutaciones: int,
    semilla: int,
) -> dict[str, float]:
    """Nulo de C13: permuta qué proteína pertenece a qué término.

    Rompe la asociación entre el estadístico y el conjunto conservando el tamaño
    de cada conjunto. Es el nulo que faltaba: R5 de v1 solo permutaba los pares,
    así que no cubría la etapa de vías.

    Devuelve, por término, la fracción de permutaciones en las que el
    estadístico absoluto fue al menos el observado.
    """
    rng = np.random.default_rng(semilla)
    stats = np.asarray(estadisticos, dtype=float)
    ok = np.isfinite(stats)
    stats = stats[ok]
    observado: dict[str, float] = {}
    tamanos: dict[str, int] = {}
    for termino, sel in pertenencia.items():
        sel = np.asarray(sel, dtype=bool)[ok]
        m = int(sel.sum())
        if m < min_proteinas:
            continue
        res = camera_un_conjunto(stats, sel, termino=termino, rho=rho, df=df)
        if res is None:
            continue
        observado[termino] = abs(res.estadistico)
        tamanos[termino] = m

    conteo = {k: 0 for k in observado}
    n = stats.size
    for _ in range(n_permutaciones):
        perm = rng.permutation(n)
        for termino, m in tamanos.items():
            sel = np.zeros(n, dtype=bool)
            sel[perm[:m]] = True
            res = camera_un_conjunto(stats, sel, termino=termino, rho=rho, df=df)
            if res is not None and abs(res.estadistico) >= observado[termino]:
                conteo[termino] += 1
    return {k: (conteo[k] + 1) / (n_permutaciones + 1) for k in observado}


def ranks_normalizados(valores: np.ndarray) -> np.ndarray:
    """Rangos en (0, 1). Útil para una prueba de conjuntos libre de distribución."""
    v = np.asarray(valores, dtype=float)
    r = rankdata(v, nan_policy="omit")
    return r / (np.nansum(np.isfinite(v)) + 1.0)
