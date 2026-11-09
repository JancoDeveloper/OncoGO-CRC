"""Etiquetas de PXD034575 a partir del sufijo de la columna.

P es el tumor primario, M la metástasis hepática y N la mucosa del mismo
participante. El contraste principal es el primario frente a la mucosa.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from oncogo_crc.utils.moderated import (
    moderate_table,
    paired_observed_stats,
    unpaired_observed_stats,
)

SUFIJO_A_TEJIDO = {"P": "primary", "M": "liver_metastasis", "N": "mucosa"}
CONTRASTE_PRINCIPAL = ("primary", "mucosa")
TEJIDOS = frozenset(SUFIJO_A_TEJIDO.values())
_COLUMNA = re.compile(r"(CRC\d+)([PMN])\b", re.IGNORECASE)


def filas_desde_columnas(columnas: list[str]) -> pd.DataFrame:
    """Una fila por columna cuyo nombre trae CRC<n> y un sufijo P, M o N."""
    filas = []
    for col in columnas:
        hallado = _COLUMNA.search(str(col))
        if not hallado:
            continue
        sufijo = hallado.group(2).upper()
        filas.append(
            {
                "sample_id": str(col),
                "participant_id": hallado.group(1).upper(),
                "suffix": sufijo,
                "tissue": SUFIJO_A_TEJIDO[sufijo],
                "accession": "PXD034575",
            }
        )
    if not filas:
        raise RuntimeError("PXD034575: ninguna columna trae sufijo P, M o N")
    tabla = pd.DataFrame(filas)
    ajenos = sorted(set(tabla["tissue"]) - TEJIDOS)
    if ajenos:
        raise RuntimeError(f"PXD034575: tejido fuera de la leyenda: {ajenos}")
    return tabla


def matriz_log2_lfq(cruda: pd.DataFrame) -> pd.DataFrame:
    """LFQ en log2. El cero y los negativos son ausencia. No se imputa."""
    if "protein_id" not in cruda.columns:
        raise RuntimeError("la matriz de PXD034575 no trae protein_id")
    mat = cruda.drop_duplicates("protein_id", keep="first").set_index("protein_id")
    valores = mat.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        logged = np.where(valores > 0, np.log2(valores), np.nan)
    return pd.DataFrame(logged, index=mat.index, columns=mat.columns)


def _pares(meta: pd.DataFrame, grupo_a: str, grupo_b: str) -> list[tuple[str, str]]:
    pares: list[tuple[str, str]] = []
    for _, grupo in meta.groupby("participant_id"):
        a = grupo.loc[grupo["tissue"] == grupo_a, "sample_id"].astype(str).tolist()
        b = grupo.loc[grupo["tissue"] == grupo_b, "sample_id"].astype(str).tolist()
        if len(a) == 1 and len(b) == 1:
            pares.append((a[0], b[0]))
    return pares


def contraste_pareado(
    log2mat: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    grupo_a: str,
    grupo_b: str,
    nombre: str,
    lfc: float,
    alpha: float,
    min_pares: int,
) -> tuple[pd.DataFrame, int]:
    pares = _pares(meta, grupo_a, grupo_b)
    if len(pares) < 2:
        raise RuntimeError(f"{nombre}: no hay pares {grupo_a} frente a {grupo_b}")
    a_cols = [a for a, _ in pares]
    b_cols = [b for _, b in pares]
    diffs = log2mat[a_cols].to_numpy(dtype=float) - log2mat[b_cols].to_numpy(dtype=float)
    coef, s2, df, stdev, n_obs = paired_observed_stats(diffs, min_pairs=min_pares)
    tabla, _d0, _s02 = moderate_table(
        protein_ids=log2mat.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method=f"moderated_t_treat_paired_{nombre}",
    )
    tabla["comparison_id"] = nombre
    tabla["n_pares_diseno"] = len(pares)
    return tabla, len(pares)


def contraste_no_pareado(
    log2mat: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    grupo_a: str,
    grupo_b: str,
    nombre: str,
    lfc: float,
    alpha: float,
    min_n: int,
) -> pd.DataFrame:
    """El mismo modelo no pareado de v1, con las columnas de la leyenda nueva."""
    cols_a = meta.loc[meta["tissue"] == grupo_a, "sample_id"].astype(str).tolist()
    cols_b = meta.loc[meta["tissue"] == grupo_b, "sample_id"].astype(str).tolist()
    grupo_x = log2mat[cols_a].to_numpy(dtype=float)
    grupo_y = log2mat[cols_b].to_numpy(dtype=float)
    coef, s2, df, stdev, n_obs = unpaired_observed_stats(grupo_x, grupo_y, min_n=min_n)
    tabla, _d0, _s02 = moderate_table(
        protein_ids=log2mat.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method=f"moderated_t_treat_unpaired_{nombre}",
    )
    tabla["comparison_id"] = nombre
    return tabla


def centrar_por_mediana(log2mat: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series, str]:
    """Resta a cada muestra la mediana de las proteínas de referencia. Devuelve el desplazamiento."""
    presencia = log2mat.notna().mean(axis=1)
    referencia, regla = presencia >= 1.0, "100 %"
    if int(referencia.sum()) < 100:
        referencia, regla = presencia >= 0.9, "90 %"
    medianas = log2mat.loc[referencia].median(axis=0)
    desplazamiento = medianas - medianas.median()
    return log2mat - desplazamiento, desplazamiento, regla
