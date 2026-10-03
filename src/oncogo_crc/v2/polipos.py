"""Secuencia adenoma–carcinoma sin maldición del ganador (corrección C5).

El problema de v1. Las proteínas se eligieron por ser significativas en el brazo
tumor de PXD034575 y después se midió su magnitud en ese mismo brazo. Una
proteína entra en la lista en parte porque el ruido de esa cohorte la empujó
hacia arriba, así que su magnitud en tumor está inflada. La conclusión «el
pólipo tiene menor magnitud que el tumor en el 91,9 % de los casos» sale en
parte por construcción: el tumor es el brazo donde se seleccionó.

Tres remedios, los tres preregistrados.

1. **Selección independiente**: las proteínas se eligen en CPTAC, que es otra
   cohorte y otra tecnología, y se miden tumor y pólipo en PXD034575. El ruido
   que las seleccionó no es el ruido que las mide.
2. **Mitades divididas**: se parten los 16 pacientes en dos, se selecciona en una
   mitad y se mide en la otra. Mil repeticiones con semilla derivada.
3. **Nulo de intercambio**: para cada proteína se sortea cuál de los dos brazos
   hace de tumor y se vuelve a aplicar la regla de selección sobre ese brazo. Si
   la fracción «menor magnitud en el pólipo» se mantiene, no dice nada del
   adenoma: es regresión a la media. El preregistro especificaba un nulo que
   intercambiaba etiquetas dentro de cada paciente; ese nulo estaba mal y la
   desviación se declara en ``docs/enmiendas_v2.md`` como enmienda E7.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from oncogo_crc.utils.moderated import moderate_table, paired_observed_stats
from oncogo_crc.v2.metricas import Fraccion


def _concordancia(
    lfc_tumor: np.ndarray,
    lfc_polipo: np.ndarray,
    *,
    etiqueta: str,
) -> dict:
    """Mismo signo, y entre las de mismo signo, menor magnitud en el pólipo."""
    ok = (
        np.isfinite(lfc_tumor)
        & np.isfinite(lfc_polipo)
        & (lfc_tumor != 0)
        & (lfc_polipo != 0)
    )
    t, p = lfc_tumor[ok], lfc_polipo[ok]
    mismo = np.sign(t) == np.sign(p)
    direccion = Fraccion(
        nombre=f"{etiqueta}_direccion",
        numerador=int(mismo.sum()),
        denominador=int(t.size),
        unidad_denominador="proteinas_unicas",
    )
    intermedio = Fraccion(
        nombre=f"{etiqueta}_intermedio",
        numerador=int((np.abs(p[mismo]) < np.abs(t[mismo])).sum()),
        denominador=int(mismo.sum()),
        unidad_denominador="proteinas_unicas",
        nota="Condicionado a que el signo coincida.",
    )
    return {
        "direccion": direccion.to_dict(),
        "intermedio_dado_mismo_signo": intermedio.to_dict(),
        "mediana_abs_lfc_tumor": float(np.median(np.abs(t))) if t.size else None,
        "mediana_abs_lfc_polipo": float(np.median(np.abs(p))) if p.size else None,
    }


def seleccion_independiente(
    tumor: pd.DataFrame,
    polipo: pd.DataFrame,
    seleccionadas: set[str],
    *,
    mapear: dict[str, str] | None = None,
) -> dict:
    """Las proteínas se eligen fuera y se miden en los dos brazos de la réplica.

    ``seleccionadas`` son identificadores de la cohorte de selección (CPTAC, que
    usa símbolo de gen). ``mapear`` traduce el identificador de la réplica
    (UniProt) a ese espacio.
    """
    t = tumor[["protein_id", "log2fc"]].rename(columns={"log2fc": "lfc_tumor"})
    p = polipo[["protein_id", "log2fc"]].rename(columns={"log2fc": "lfc_polipo"})
    m = t.merge(p, on="protein_id", how="inner")
    if mapear is not None:
        clave = m["protein_id"].astype(str).map(lambda v: mapear.get(v, ""))
    else:
        clave = m["protein_id"].astype(str).str.upper()
    m = m.loc[clave.isin({s.strip().upper() for s in seleccionadas})]
    res = _concordancia(
        m["lfc_tumor"].to_numpy(dtype=float),
        m["lfc_polipo"].to_numpy(dtype=float),
        etiqueta="seleccion_independiente",
    )
    res["n_seleccionadas_en_cohorte_externa"] = len(seleccionadas)
    res["n_medidas_en_la_replica"] = int(len(m))
    res["nota"] = (
        "Las proteínas se eligieron en CPTAC y se midieron en PXD034575. El ruido "
        "que las seleccionó no es el que las mide, así que la magnitud en tumor no "
        "está inflada por la selección."
    )
    return res


def mitades_divididas(
    diffs_tumor: np.ndarray,
    diffs_polipo: np.ndarray,
    protein_ids: np.ndarray,
    *,
    lfc: float,
    alpha: float,
    min_pares: int,
    n_repeticiones: int,
    semilla: int,
) -> dict:
    """Selecciona en una mitad de los pacientes y mide en la otra.

    Las dos matrices tienen que venir con los pacientes en el mismo orden de
    columnas: la columna ``j`` de las dos es el mismo paciente.
    """
    if diffs_tumor.shape[1] != diffs_polipo.shape[1]:
        raise ValueError("tumor y pólipo deben tener los mismos pacientes en el mismo orden")
    n_pac = diffs_tumor.shape[1]
    if n_pac < 2 * min_pares:
        return {
            "estado": "no_calculable",
            "motivo": f"{n_pac} pacientes no alcanzan para dos mitades de {min_pares}",
        }
    rng = np.random.default_rng(semilla)
    direcciones: list[float] = []
    intermedios: list[float] = []
    n_evaluables: list[int] = []
    for _ in range(n_repeticiones):
        orden = rng.permutation(n_pac)
        a, b = orden[: n_pac // 2], orden[n_pac // 2 :]
        # Selección en la mitad A, usando solo el brazo tumor.
        coef, s2, df, stdev, n_obs = paired_observed_stats(
            diffs_tumor[:, a], min_pairs=min_pares
        )
        sel_tabla, _d0, _s02 = moderate_table(
            protein_ids=protein_ids,
            coef=coef,
            s2=s2,
            df=df,
            stdev_unscaled=stdev,
            n_obs=n_obs,
            lfc=lfc,
            alpha=alpha,
            stat_method="mitad_seleccion",
        )
        sel = sel_tabla["significant"].to_numpy()
        if sel.sum() < 5:
            continue
        # Medición en la mitad B, en los dos brazos.
        lfc_t = np.nanmean(diffs_tumor[:, b], axis=1)
        lfc_p = np.nanmean(diffs_polipo[:, b], axis=1)
        res = _concordancia(lfc_t[sel], lfc_p[sel], etiqueta="mitad")
        d = res["direccion"]["valor"]
        i = res["intermedio_dado_mismo_signo"]["valor"]
        if d is not None:
            direcciones.append(d)
            n_evaluables.append(res["direccion"]["denominador"])
        if i is not None:
            intermedios.append(i)

    def _pct(arr: list[float]) -> dict:
        if not arr:
            return {"estado": "no_calculable", "n_repeticiones_utiles": 0}
        a = np.asarray(arr)
        return {
            "mediana": float(np.median(a)),
            "p2_5": float(np.percentile(a, 2.5)),
            "p97_5": float(np.percentile(a, 97.5)),
            "n_repeticiones_utiles": int(a.size),
        }

    return {
        "n_repeticiones_pedidas": int(n_repeticiones),
        "semilla": int(semilla),
        "n_pacientes": int(n_pac),
        "direccion": _pct(direcciones),
        "intermedio_dado_mismo_signo": _pct(intermedios),
        "n_proteinas_evaluables_mediana": (
            float(np.median(n_evaluables)) if n_evaluables else None
        ),
        "nota": (
            "Selección en una mitad de los pacientes y medición en la otra. Ninguna "
            "proteína se mide en los pacientes que la seleccionaron."
        ),
    }


def _estadistico_pareado(diffs: np.ndarray, min_pares: int) -> tuple[np.ndarray, np.ndarray]:
    """Media por par y su t ordinario. Las proteínas con pocos pares salen NaN.

    El filtro por número de pares no es cosmético. Sin él, la proteína medida en
    un solo paciente puede tener la media más grande de la matriz y se cuela en
    cualquier selección por magnitud, aunque no tenga error estándar estimable.
    """
    limpio = np.where(np.isfinite(diffs), diffs, np.nan)
    n = np.isfinite(limpio).sum(axis=1)
    ok = n >= min_pares
    media = np.full(n.shape, np.nan)
    sd = np.full(n.shape, np.nan)
    if ok.any():
        media[ok] = np.nanmean(limpio[ok], axis=1)
        sd[ok] = np.nanstd(limpio[ok], axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = media / (sd / np.sqrt(np.maximum(n, 1)))
    return np.where(ok, media, np.nan), np.where(ok & np.isfinite(t), t, np.nan)


def _seleccion_por_estadistico(
    estadistico: np.ndarray, elegibles: np.ndarray, n_elegir: int
) -> np.ndarray:
    """Las ``n_elegir`` proteínas elegibles con mayor |estadístico|."""
    puntaje = np.where(elegibles & np.isfinite(estadistico), np.abs(estadistico), -np.inf)
    orden = np.argsort(-puntaje, kind="stable")
    sel = np.zeros(estadistico.size, dtype=bool)
    tope = min(n_elegir, int(np.isfinite(puntaje[puntaje > -np.inf]).sum()))
    sel[orden[:tope]] = True
    return sel


def nulo_intercambio(
    diffs_tumor: np.ndarray,
    diffs_polipo: np.ndarray,
    *,
    n_elegir: int,
    n_permutaciones: int,
    semilla: int,
    min_pares: int = 3,
) -> dict:
    """Nulo de regresión a la media para la asimetría pólipo–tumor.

    La hipótesis nula es que los dos brazos son intercambiables. Para cada
    proteína se decide al azar cuál brazo hace de «tumor», y después se vuelve a
    seleccionar con la misma regla sobre ese brazo. Eso es lo importante: el nulo
    tiene que repetir **todo** el procedimiento, selección incluida. Si solo se
    intercambian las etiquetas sin reseleccionar, se destruye la ventaja de la
    selección y el nulo deja de medir el sesgo que se quería medir.

    Un nulo centrado en el valor observado dice que la asimetría es regresión a
    la media. Un nulo por debajo dice que el adenoma es de verdad intermedio.
    """
    if diffs_tumor.shape != diffs_polipo.shape:
        raise ValueError("las dos matrices deben tener la misma forma")
    rng = np.random.default_rng(semilla)
    n_prot = diffs_tumor.shape[0]
    media_t, t_t = _estadistico_pareado(diffs_tumor, min_pares)
    media_p, t_p = _estadistico_pareado(diffs_polipo, min_pares)
    # Solo proteínas medibles en los dos brazos: si no, intercambiarlos no es
    # una permutación válida porque los dos lados no son comparables.
    elegibles = np.isfinite(media_t) & np.isfinite(media_p) & np.isfinite(t_t) & np.isfinite(t_p)

    sel_obs = _seleccion_por_estadistico(t_t, elegibles, n_elegir)
    observado = _concordancia(media_t[sel_obs], media_p[sel_obs], etiqueta="observado")
    obs_inter = observado["intermedio_dado_mismo_signo"]["valor"]

    nulos: list[float] = []
    for _ in range(n_permutaciones):
        intercambia = rng.random(n_prot) < 0.5
        primero = np.where(intercambia, media_p, media_t)
        segundo = np.where(intercambia, media_t, media_p)
        est = np.where(intercambia, t_p, t_t)
        sel = _seleccion_por_estadistico(est, elegibles, n_elegir)
        res = _concordancia(primero[sel], segundo[sel], etiqueta="nulo")
        v = res["intermedio_dado_mismo_signo"]["valor"]
        if v is not None:
            nulos.append(v)
    arr = np.asarray(nulos)
    if arr.size == 0 or obs_inter is None:
        return {"estado": "no_calculable"}
    b = int((arr >= obs_inter).sum())
    return {
        "observado_intermedio": float(obs_inter),
        "observado_direccion": observado["direccion"]["valor"],
        "n_elegir": int(n_elegir),
        "n_elegidas": int(sel_obs.sum()),
        "n_elegibles_en_los_dos_brazos": int(elegibles.sum()),
        "min_pares": int(min_pares),
        "nulo_mediana": float(np.median(arr)),
        "nulo_p2_5": float(np.percentile(arr, 2.5)),
        "nulo_p97_5": float(np.percentile(arr, 97.5)),
        "n_permutaciones": int(arr.size),
        "b_al_menos_observado": b,
        "p_empirico": (b + 1) / (arr.size + 1),
        "semilla": int(semilla),
        "nota": (
            "Para cada proteína se sortea qué brazo hace de tumor y se vuelve a "
            "seleccionar con la misma regla. El nulo incluye la selección, así que "
            "lo que mide es cuánta de la asimetría observada es regresión a la media."
        ),
    }
