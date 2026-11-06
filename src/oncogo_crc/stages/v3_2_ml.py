"""Prueba global del suero. Tres versiones y el nulo de permutación."""

from __future__ import annotations

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.stages.v3_2_controles import _abrir
from oncogo_crc.v3.ml_suero import evaluar, matrices_de_suero
from oncogo_crc.v3.suero_v32 import semilla_v32
from oncogo_crc.v3.salida_v32 import volcar_json


def run(settings: Settings) -> None:
    prereg, salida, sello = _abrir(settings, "v3_2_ml")
    principal = salida / "suero_matriz_principal.tsv"
    humana = salida / "suero_matriz_humana_fdr.tsv"
    calidad = salida / "suero_calidad_muestra.tsv"
    if not principal.is_file() or not humana.is_file() or not calidad.is_file():
        raise FileNotFoundError("faltan las matrices del suero; corre v3-2-suero antes")
    bloques = matrices_de_suero(
        pd.read_csv(principal, sep="\t"),
        pd.read_csv(humana, sep="\t"),
        pd.read_csv(calidad, sep="\t"),
        prereg,
    )
    semilla_cv = semilla_v32(prereg, "ml_cv")
    semilla_nulo = semilla_v32(prereg, "ml_permutacion")
    versiones = {}
    for nombre, (x, y) in bloques.items():
        print(f"[v3.2] ML {nombre}: validación cruzada y {1000} permutaciones")
        versiones[nombre] = evaluar(x, y, semilla_cv, semilla_nulo)
    volcar_json({
        "estado": "calculado",
        "pregunta": "cáncer frente a colonoscopia negativa",
        "n_por_clase": "sale de las columnas CP_ y CN_",
        "versiones": versiones,
        "no_es_biomarcador": True,
        "fase_pepe": 1,
    }, salida / "ml_suero_permutacion.json", sello)
    print("[v3.2] ml_suero_permutacion.json escrito")
