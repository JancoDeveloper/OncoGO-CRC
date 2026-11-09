"""Prueba global del suero. Seis modelos y el titular de la regla escrita antes de correr."""

from __future__ import annotations

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.stages.v3_3_controles import _abrir
from oncogo_crc.v3.ml_suero import _columnas_cp_cn, evaluar, matrices_de_suero, titular
from oncogo_crc.v3.salida import volcar_json, volcar_tsv
from oncogo_crc.v3.suero_v32 import semilla_v33


def run(settings: Settings) -> None:
    prereg, salida, sello = _abrir(settings, "v3_3_ml")
    principal = salida / "suero_matriz_principal.tsv"
    humana = salida / "suero_matriz_humana_fdr.tsv"
    calidad = salida / "suero_calidad_muestra.tsv"
    if not principal.is_file() or not humana.is_file() or not calidad.is_file():
        raise FileNotFoundError("faltan las matrices del suero; corre v3-3-suero antes")
    principal_df = pd.read_csv(principal, sep="\t")
    humana_df = pd.read_csv(humana, sep="\t")
    calidad_df = pd.read_csv(calidad, sep="\t")
    bloques = matrices_de_suero(principal_df, humana_df, calidad_df, prereg)
    semilla_cv = semilla_v33(prereg, "ml_cv")
    semilla_nulo = semilla_v33(prereg, "ml_permutacion")
    versiones = {}
    for nombre, (x, y, tipo) in bloques.items():
        if x.shape[0] != 40 or int(np.sum(y)) != 20:
            raise RuntimeError(f"{nombre} no tiene 40 filas (20 CP y 20 CN)")
        print(f"[v3.3] ML {nombre}: validación cruzada y 1000 permutaciones")
        versiones[nombre] = evaluar(x, y, semilla_cv, semilla_nulo, tipo)
    textos = prereg.get("enmienda_v3_3.ml_suero.titular.textos")
    decision = titular(
        versiones,
        float(prereg.get("enmienda_v3_3.ml_suero.titular.alfa")),
        float(prereg.get("enmienda_v3_3.ml_suero.titular.auc_minima")),
        textos,
    )
    crc, cn = _columnas_cp_cn(principal_df, humana_df)
    volcar_json({
        "estado": "calculado",
        "pregunta": "cáncer frente a colonoscopia negativa",
        "n_por_clase": {"crc": 20, "colonoscopy_negative": 20},
        "columnas_usadas": crc + cn,
        "escala": prereg.get("enmienda_v3_3.ml_suero.escala"),
        "versiones": versiones,
        "titular": decision,
        "no_es_biomarcador": True,
        "fase_pepe": 1,
    }, salida / "ml_suero_permutacion.json", sello)
    filas = []
    for nombre, bloque in versiones.items():
        filas.append({"modelo": nombre, **{k: v for k, v in bloque.items() if not isinstance(v, (dict, list))}})
    volcar_tsv(pd.DataFrame(filas), salida / "ml_suero_modelos.tsv", sello)
    print(f"[v3.3] titular: {decision['clave']}")
