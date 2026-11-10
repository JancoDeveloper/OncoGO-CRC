"""Rejilla v3.4. 288 celdas de suero, 72 de ML y las dos normalizaciones de R1."""

from __future__ import annotations

import itertools

CUANTIFICACION = ("PG.MaxLFQ", "Genes.MaxLFQ")
CONTAMINANTES = ("libreria_v3_3", "sin_quitar", "con_queratinas_epiteliales", "sin_grupos_mixtos")
PANELES = ("excluidos", "covariable_coagulacion", "incluidos")
PRESENCIA = ("min_n_3", "0.5", "0.7")
LFC = (0.585, 1.0)
ALFA = (0.05, 0.10)


def celdas_suero() -> list[dict]:
    salida = []
    for cuant, cont, panel, presencia, lfc, alfa in itertools.product(
        CUANTIFICACION, CONTAMINANTES, PANELES, PRESENCIA, LFC, ALFA
    ):
        salida.append({
            "cuantificacion": cuant,
            "contaminantes": cont,
            "paneles": panel,
            "presencia": presencia,
            "lfc": lfc,
            "alfa": alfa,
            "es_la_celda_preregistrada": (
                cuant == "PG.MaxLFQ" and cont == "libreria_v3_3" and panel == "excluidos"
                and presencia == "min_n_3" and lfc == 0.585 and alfa == 0.05
            ),
        })
    return salida


def celdas_ml() -> list[dict]:
    salida = []
    for cuant, cont, panel, presencia in itertools.product(CUANTIFICACION, CONTAMINANTES, PANELES, PRESENCIA):
        salida.append({
            "cuantificacion": cuant,
            "contaminantes": cont,
            "paneles": panel,
            "presencia": presencia,
            "es_la_celda_preregistrada": (
                cuant == "PG.MaxLFQ" and cont == "libreria_v3_3" and panel == "excluidos" and presencia == "min_n_3"
            ),
        })
    return salida
