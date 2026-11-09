"""Controles de mucosa v3.4. No detectado no es falla.

La regla del bloque se escribió después de ver los resultados v3.3.
El veredicto v3.3 se conserva en results/v3_3/ y se informa al lado.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd


def filas_fisicas(grupos: pd.DataFrame, gen: str, mapa_uni: dict[str, str]) -> dict:
    """Cuenta filas de proteinGroups que son del gen, antes y después de los filtros de MaxQuant."""
    buscado = gen.upper()
    gn = grupos["Gene names"].fillna("").str.upper().str.split(";")
    fasta = grupos["Fasta headers"].fillna("").str.upper().str.contains(
        rf"GN={re.escape(buscado)}(?:\s|;|$)", regex=True
    )
    ids = grupos["Protein IDs"].fillna("").str.split(";").apply(
        lambda xs: any(mapa_uni.get(str(x).split("-")[0].replace("CON__", "").replace("REV__", "")) == buscado for x in xs)
    )
    es = gn.apply(lambda xs: buscado in [p.strip() for p in xs]) | fasta | ids
    malos = (
        grupos["Protein IDs"].fillna("").str.startswith(("CON__", "REV__"))
        | grupos[["Reverse", "Potential contaminant", "Only identified by site"]].eq("+").any(axis=1)
    )
    return {"antes_de_filtrar": int(es.sum()), "despues_de_filtrar": int((es & ~malos).sum())}


def estado_control(fisico: dict, n_obs: float | None, log2fc: float, min_pares: int) -> str:
    if fisico["antes_de_filtrar"] == 0:
        return "no_detectado"
    if fisico["despues_de_filtrar"] == 0:
        return "filtrado"
    if n_obs is None or not np.isfinite(n_obs) or n_obs < min_pares or not np.isfinite(log2fc):
        return "no_evaluable_por_pares"
    return "pasa" if log2fc < 0 else "falla"


def estado_bloque(estados: list[str], qc_mucosa_todos_pasan: bool) -> str:
    evaluables = [e for e in estados if e in ("pasa", "falla")]
    if not evaluables:
        return "no_evaluable"
    if "falla" in evaluables or not qc_mucosa_todos_pasan:
        return "falla"
    return "pasa"
