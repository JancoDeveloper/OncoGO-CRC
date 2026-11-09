"""Controles de mucosa: no detectado no es falla."""

from __future__ import annotations

import pandas as pd

from oncogo_crc.stages.v3_4_controles import estado_bloque, estado_control, filas_fisicas


def _grupos(gen: str | None) -> pd.DataFrame:
    fila = {
        "Gene names": gen or "",
        "Fasta headers": f"sp|P00915|{gen}_HUMAN GN={gen}" if gen else "",
        "Protein IDs": "P00915" if gen else "P99999",
        "Reverse": "",
        "Potential contaminant": "",
        "Only identified by site": "",
    }
    return pd.DataFrame([fila])


def test_gen_ausente_es_no_detectado():
    fisico = filas_fisicas(_grupos(None), "AQP8", {})
    assert fisico["antes_de_filtrar"] == 0
    assert estado_control(fisico, None, float("nan"), 8) == "no_detectado"


def test_pocos_pares_no_es_falla():
    fisico = {"antes_de_filtrar": 1, "despues_de_filtrar": 1}
    assert estado_control(fisico, 4, float("nan"), 8) == "no_evaluable_por_pares"


def test_ca1_con_descenso_pasa():
    fisico = {"antes_de_filtrar": 1, "despues_de_filtrar": 1}
    assert estado_control(fisico, 13, -2.31, 8) == "pasa"


def test_bloque_pasa_si_el_qc_de_mucosa_pasa_en_todos():
    estados = ["pasa", "no_detectado", "no_evaluable_por_pares", "no_detectado"]
    assert estado_bloque(estados, True) == "pasa"


def test_bloque_falla_si_el_qc_no_pasa_en_todos():
    estados = ["pasa", "no_detectado", "no_evaluable_por_pares", "no_detectado"]
    assert estado_bloque(estados, False) == "falla"
