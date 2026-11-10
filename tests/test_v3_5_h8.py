"""H8 v3.5: mapa de Windows, estadio, directLFQ y Fisher. Fallan antes del parche."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from scipy.stats import fisher_exact

from oncogo_crc.v3.h8 import (
    agregar_directlfq,
    clasificar_corridas,
    fisher_carga,
    limpiar_estadio,
    unir_mapa,
)


def test_el_mapa_une_rutas_de_windows_y_excluye_estandares():
    columnas = [
        "Protein.Group",
        r"D:\X\SH_01.mzML",
        "SH_02.mzML",
        r"D:\X\ST_01_Pre_Batch_01.mzML",
    ]
    muestras, estandares = clasificar_corridas(columnas)
    assert muestras == {"SH_01": r"D:\X\SH_01.mzML", "SH_02": "SH_02.mzML"}
    assert estandares == ["ST_01_Pre_Batch_01.mzML"]
    meta = pd.DataFrame({"source name": ["SH_01", "SH_02"], "characteristic [stage]": ["4", "1"]})
    mapa = unir_mapa(muestras, meta)
    assert mapa["uno_a_uno"] is True
    assert mapa["n_unidas"] == 2


def test_el_estadio_limpia_el_retorno_de_carro_y_la_columna_en_singular():
    meta = pd.DataFrame({
        "source name": ["SH_01"],
        "characteristic [stage]\r": [" 4\r "],
    })
    limpio = limpiar_estadio(meta)
    assert list(limpio.columns) == ["muestra", "estadio"]
    assert limpio.iloc[0]["estadio"] == "4"


def test_directlfq_agrega_dos_proteinas_en_cuatro_corridas(tmp_path: Path):
    filas = []
    for proteina, iones in (("P1", ("i1", "i2")), ("P2", ("i3",))):
        for ion in iones:
            filas.append({
                "protein": proteina,
                "ion": ion,
                "SH_01": 100.0,
                "SH_02": 200.0,
                "SH_03": 50.0,
                "SH_04": 80.0 if ion == "i1" else float("nan"),
            })
    entrada = tmp_path / "sintetica.aq_reformat.tsv"
    pd.DataFrame(filas).to_csv(entrada, sep="\t", index=False)
    salida, usado = agregar_directlfq(entrada)
    tabla = pd.read_csv(salida, sep="\t")
    assert len(tabla) == 2
    muestras = [c for c in tabla.columns if c.startswith("SH_")]
    assert len(muestras) == 4
    assert usado in {"nan", "cero"}


def test_fisher_de_una_tabla_conocida():
    tabla = [[15, 5], [10, 70]]
    razon, p = fisher_exact(tabla, alternative="greater")
    resultado = fisher_carga(tabla)
    assert resultado["odds_ratio"] == pytest.approx(float(razon))
    assert resultado["p"] == pytest.approx(float(p))
    assert resultado["estado"] == "sostiene"
    assert resultado["n_tumorales_epiteliales"] == 20
