"""Tests de C9: R1 recalculado con el efecto ya ajustado por composición.

Lo que estos tests fijan es que la selección y el signo salgan los dos de la
tabla ajustada. Si uno de los dos se quedara en la tabla cruda, el cálculo
volvería a ser el estrato de v1 y no respondería a la pregunta de C9.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.stages.v2_certified import _r1_tras_el_ajuste

R1_BASE = {
    "global": {
        "nombre": "R1_direccion",
        "valor": 0.5,
        "numerador": 2,
        "denominador": 4,
        "ci_low": 0.0676,
        "ci_high": 0.9324,
    }
}


def _overlap(genes, a_crudo, b, padj_crudo):
    return pd.DataFrame(
        {
            "protein_id_cptac": genes,
            "log2fc_cptac": a_crudo,
            "log2fc_pride": b,
            "padj_cptac": padj_crudo,
        }
    )


def _ajustada(genes, log2fc, padj):
    return pd.DataFrame({"protein_id": genes, "log2fc": log2fc, "padj": padj})


def test_sin_tabla_ajustada_queda_no_calculable():
    ov = _overlap(["A"], [1.0], [1.0], [0.01])
    res = _r1_tras_el_ajuste(ov, None, alpha=0.05, r1_sin_ajustar=R1_BASE)
    assert res["estado"] == "no_calculable"

    vacia = pd.DataFrame(columns=["protein_id", "log2fc", "padj"])
    res2 = _r1_tras_el_ajuste(ov, vacia, alpha=0.05, r1_sin_ajustar=R1_BASE)
    assert res2["estado"] == "no_calculable"


def test_la_seleccion_usa_el_padj_ajustado_y_no_el_crudo():
    """B es significativa sin ajustar y deja de serlo al ajustar: no entra."""
    genes = ["A", "B"]
    ov = _overlap(genes, [1.0, 1.0], [1.0, -1.0], [0.01, 0.01])
    aj = _ajustada(genes, [1.0, 1.0], [0.01, 0.40])
    res = _r1_tras_el_ajuste(ov, aj, alpha=0.05, r1_sin_ajustar=R1_BASE)
    assert res["fraccion"]["denominador"] == 1
    assert res["fraccion"]["numerador"] == 1


def test_el_signo_usa_el_efecto_ajustado_y_no_el_crudo():
    """El ajuste le cambia el signo a A: la concordancia tiene que cambiar con él."""
    genes = ["A"]
    ov = _overlap(genes, [1.0], [1.0], [0.01])
    crudo_concuerda = _r1_tras_el_ajuste(
        ov, _ajustada(genes, [1.0], [0.01]), alpha=0.05, r1_sin_ajustar=R1_BASE
    )
    assert crudo_concuerda["fraccion"]["numerador"] == 1

    volteado = _r1_tras_el_ajuste(
        ov, _ajustada(genes, [-1.0], [0.01]), alpha=0.05, r1_sin_ajustar=R1_BASE
    )
    assert volteado["fraccion"]["numerador"] == 0


def test_un_efecto_ajustado_de_cero_no_es_evaluable():
    """Con efecto cero no hay signo que comparar."""
    genes = ["A", "B"]
    ov = _overlap(genes, [1.0, 1.0], [1.0, 1.0], [0.01, 0.01])
    aj = _ajustada(genes, [0.0, 1.0], [0.01, 0.01])
    res = _r1_tras_el_ajuste(ov, aj, alpha=0.05, r1_sin_ajustar=R1_BASE)
    assert res["fraccion"]["denominador"] == 1


def test_si_nada_sobrevive_el_ajuste_queda_no_calculable():
    genes = ["A", "B"]
    ov = _overlap(genes, [1.0, 1.0], [1.0, 1.0], [0.01, 0.01])
    aj = _ajustada(genes, [1.0, 1.0], [0.9, 0.8])
    res = _r1_tras_el_ajuste(ov, aj, alpha=0.05, r1_sin_ajustar=R1_BASE)
    assert res["estado"] == "no_calculable"
    assert res["n_significativas_tras_el_ajuste"] == 0


def test_una_proteina_que_la_replica_no_mide_no_entra():
    genes = ["A", "B"]
    ov = _overlap(genes, [1.0, 1.0], [1.0, np.nan], [0.01, 0.01])
    aj = _ajustada(genes, [1.0, 1.0], [0.01, 0.01])
    res = _r1_tras_el_ajuste(ov, aj, alpha=0.05, r1_sin_ajustar=R1_BASE)
    assert res["fraccion"]["denominador"] == 1


def test_el_cruce_de_genes_no_depende_de_mayusculas_ni_espacios():
    ov = _overlap([" tp53 "], [1.0], [1.0], [0.01])
    aj = _ajustada(["TP53"], [1.0], [0.01])
    res = _r1_tras_el_ajuste(ov, aj, alpha=0.05, r1_sin_ajustar=R1_BASE)
    assert res["fraccion"]["denominador"] == 1


def test_se_reporta_cuanto_denominador_cuesta_el_ajuste():
    genes = [f"G{i}" for i in range(4)]
    ov = _overlap(genes, [1.0] * 4, [1.0] * 4, [0.01] * 4)
    aj = _ajustada(genes, [1.0] * 4, [0.01, 0.01, 0.9, 0.9])
    res = _r1_tras_el_ajuste(ov, aj, alpha=0.05, r1_sin_ajustar=R1_BASE)
    assert res["fraccion"]["denominador"] == 2
    assert res["n_proteinas_que_pierde_el_denominador"] == 2


def test_el_veredicto_sale_de_los_intervalos_y_no_de_los_puntos():
    """Dos valores distintos con intervalos que se solapan no son un cambio."""
    genes = [f"G{i}" for i in range(4)]
    ov = _overlap(genes, [1.0] * 4, [1.0, 1.0, 1.0, -1.0], [0.01] * 4)
    aj = _ajustada(genes, [1.0] * 4, [0.01] * 4)
    res = _r1_tras_el_ajuste(ov, aj, alpha=0.05, r1_sin_ajustar=R1_BASE)
    assert res["fraccion"]["valor"] == pytest.approx(0.75)
    assert res["cambio_en_la_concordancia"] == pytest.approx(0.25)
    assert res["los_dos_intervalos_se_solapan"] is True
    assert "se solapan" in res["lectura"]


def test_sin_intervalo_en_la_base_no_se_declara_solape():
    base = {"global": {"valor": 0.5, "denominador": 4, "ci_low": None, "ci_high": None}}
    genes = ["A", "B", "C"]
    ov = _overlap(genes, [1.0] * 3, [1.0] * 3, [0.01] * 3)
    aj = _ajustada(genes, [1.0] * 3, [0.01] * 3)
    res = _r1_tras_el_ajuste(ov, aj, alpha=0.05, r1_sin_ajustar=base)
    assert res["los_dos_intervalos_se_solapan"] is False
