"""H5 v3.4: el submuestreo completa las repeticiones y el universo es el mismo."""

from __future__ import annotations

import numpy as np
import pandas as pd

from oncogo_crc.v3.h5 import concordancias_mismo_universo, submuestreo


def _matriz(n_genes: int, n_pares: int, semilla: int, desplazamiento: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(semilla)
    genes = [f"G{i}" for i in range(n_genes)]
    columnas = [f"P{i}" for i in range(n_pares)]
    tumor = pd.DataFrame(rng.normal(desplazamiento, 0.2, size=(n_genes, n_pares)), index=genes, columns=columnas)
    normal = pd.DataFrame(rng.normal(0.0, 0.2, size=(n_genes, n_pares)), index=genes, columns=columnas)
    return tumor, normal


def test_el_bucle_de_v3_3_no_completa_ninguna_repeticion():
    """Con n_colon = 96 y un carcinoma de 80 pares, el break de v3.3 deja el conteo en 0."""
    n_colon = 96
    pares = {"CCRCC": 80, "OV": 10}
    iguales = 0
    for _ in range(200):
        completos = 1
        for n in pares.values():
            if n < n_colon:
                break
            completos += 1
        if completos == 1 + len(pares):
            iguales += 1
    assert iguales == 0


def test_submuestreo_completa_200_y_n_sub_es_10():
    matrices = {
        "COAD": _matriz(8, 96, 1, 2.0),
        "CCRCC": _matriz(8, 80, 2, 1.0),
        "OV": _matriz(8, 10, 3, 0.5),
    }
    genes = [f"G{i}" for i in range(8)]
    mapa = {gen: gen for gen in genes}
    resultado = submuestreo(
        matrices, ["COAD", "CCRCC", "OV"], genes, mapa, np.random.default_rng(930),
        n_rep=200, min_pares=8, lfc=0.585, alpha=0.05, umbral_pan=0.20, fraccion_pan=0.10,
    )
    assert resultado["estado"] == "calculado"
    assert resultado["n_sub"] == 10
    assert resultado["n_repeticiones"] == 200


def test_concordancias_comparten_el_denominador():
    tanaka = pd.Series({"G1": 1.2, "G2": -0.8, "G3": 0.4, "SOLO": 1.0})
    efectos = {
        "COAD": pd.DataFrame({"log2fc": [1.0, -0.5, np.nan]}, index=["G1", "G2", "G3"]),
        "OV": pd.DataFrame({"log2fc": [0.7, np.nan, 0.2]}, index=["G1", "G2", "G3"]),
        "PDAC": pd.DataFrame({"log2fc": [-0.2, -0.4, 0.9]}, index=["G1", "G2", "G3"]),
    }
    salida, universo = concordancias_mismo_universo(tanaka, efectos, ["COAD", "OV", "PDAC"])
    denominadores = {bloque["denominador"] for bloque in salida.values()}
    assert len(denominadores) == 1
    assert "SOLO" not in universo
    assert "G3" not in universo


def test_submuestreo_no_calculado_si_hay_menos_de_8_pares():
    matrices = {"COAD": _matriz(4, 7, 4, 1.0), "OV": _matriz(4, 12, 5, 1.0)}
    genes = [f"G{i}" for i in range(4)]
    resultado = submuestreo(
        matrices, ["COAD", "OV"], genes, {g: g for g in genes}, np.random.default_rng(1),
        n_rep=3, min_pares=8, lfc=0.585, alpha=0.05, umbral_pan=0.20, fraccion_pan=0.10,
    )
    assert resultado["estado"] == "no_calculado"
