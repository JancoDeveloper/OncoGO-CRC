"""Tests de C5: la maldición del ganador en la secuencia adenoma–carcinoma."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.v2.polipos import (
    mitades_divididas,
    nulo_intercambio,
    seleccion_independiente,
)


def test_seleccion_independiente_solo_mide_las_proteinas_elegidas_fuera():
    tumor = pd.DataFrame({"protein_id": ["A", "B", "C", "D"], "log2fc": [2.0, -2.0, 1.0, 3.0]})
    polipo = pd.DataFrame({"protein_id": ["A", "B", "C", "D"], "log2fc": [1.0, -1.0, -0.5, 2.0]})
    res = seleccion_independiente(tumor, polipo, {"A", "B", "D"})
    assert res["n_medidas_en_la_replica"] == 3
    assert res["direccion"]["denominador"] == 3
    assert res["direccion"]["numerador"] == 3  # A, B y D coinciden en signo
    # De las tres, todas tienen menor magnitud en el pólipo.
    assert res["intermedio_dado_mismo_signo"]["numerador"] == 3


def test_seleccion_independiente_traduce_identificadores():
    tumor = pd.DataFrame({"protein_id": ["P1", "P2"], "log2fc": [1.0, -1.0]})
    polipo = pd.DataFrame({"protein_id": ["P1", "P2"], "log2fc": [0.5, -0.5]})
    mapa = {"P1": "GENA", "P2": "GENB"}
    res = seleccion_independiente(tumor, polipo, {"GENA"}, mapear=mapa)
    assert res["n_medidas_en_la_replica"] == 1


def test_el_nulo_de_intercambio_reproduce_la_asimetria_por_seleccion():
    """Sin señal real, seleccionar en un brazo ya hace que el otro parezca menor.

    Es la maldición del ganador en estado puro: los dos brazos son ruido
    independiente, así que el nulo tiene que igualar lo observado.
    """
    rng = np.random.default_rng(20260914)
    n_prot, n_pac = 600, 16
    t = rng.normal(0, 1, size=(n_prot, n_pac))
    p = rng.normal(0, 1, size=(n_prot, n_pac))
    res = nulo_intercambio(t, p, n_elegir=60, n_permutaciones=300, semilla=20260914)
    # El observado está alto por construcción, no por biología.
    assert res["observado_intermedio"] > 0.8
    # Y el nulo lo iguala, porque incluye la misma selección.
    assert res["p_empirico"] > 0.05, (
        "con ruido puro la asimetría debe ser indistinguible del nulo que reselecciona"
    )
    assert res["nulo_mediana"] > 0.8


def test_el_nulo_de_intercambio_separa_una_asimetria_real():
    """Si el pólipo es de verdad intermedio, el nulo no llega a lo observado."""
    rng = np.random.default_rng(7)
    n_prot, n_pac = 600, 16
    efecto = rng.normal(0, 2, size=(n_prot, 1))
    t = efecto + rng.normal(0, 0.3, size=(n_prot, n_pac))
    p = 0.4 * efecto + rng.normal(0, 0.3, size=(n_prot, n_pac))
    res = nulo_intercambio(t, p, n_elegir=200, n_permutaciones=300, semilla=20260914)
    assert res["observado_intermedio"] > 0.95
    assert res["p_empirico"] < 0.05
    assert res["nulo_mediana"] < res["observado_intermedio"]


def test_mitades_divididas_no_mide_en_los_pacientes_que_seleccionaron():
    rng = np.random.default_rng(3)
    n_prot, n_pac = 300, 16
    efecto = rng.normal(0, 2, size=(n_prot, 1))
    t = efecto + rng.normal(0, 0.4, size=(n_prot, n_pac))
    p = 0.5 * efecto + rng.normal(0, 0.4, size=(n_prot, n_pac))
    ids = np.array([f"P{i}" for i in range(n_prot)])
    res = mitades_divididas(
        t,
        p,
        ids,
        lfc=0.585,
        alpha=0.05,
        min_pares=3,
        n_repeticiones=50,
        semilla=20260914,
    )
    assert res["n_pacientes"] == 16
    assert res["direccion"]["n_repeticiones_utiles"] > 0
    assert 0.0 <= res["direccion"]["mediana"] <= 1.0
    assert res["direccion"]["p2_5"] <= res["direccion"]["mediana"] <= res["direccion"]["p97_5"]


def test_mitades_divididas_avisa_cuando_no_hay_pacientes():
    t = np.zeros((10, 4))
    res = mitades_divididas(
        t, t, np.arange(10).astype(str), lfc=0.585, alpha=0.05, min_pares=3,
        n_repeticiones=10, semilla=1,
    )
    assert res["estado"] == "no_calculable"
    assert "pacientes" in res["motivo"]


def test_las_matrices_tienen_que_tener_los_mismos_pacientes():
    with pytest.raises(ValueError, match="mismos pacientes"):
        mitades_divididas(
            np.zeros((5, 16)), np.zeros((5, 8)), np.arange(5).astype(str),
            lfc=0.585, alpha=0.05, min_pares=3, n_repeticiones=5, semilla=1,
        )
    with pytest.raises(ValueError, match="misma forma"):
        nulo_intercambio(
            np.zeros((5, 16)), np.zeros((5, 8)), n_elegir=3, n_permutaciones=5, semilla=1
        )


def test_el_nulo_descarta_proteinas_con_pocos_pares():
    """El fallo que apareció con los datos reales: una proteína de un solo par
    puede tener la media más grande de la matriz y colarse en la selección."""
    rng = np.random.default_rng(11)
    n_prot, n_pac = 200, 16
    t = rng.normal(0, 1, size=(n_prot, n_pac))
    p = rng.normal(0, 1, size=(n_prot, n_pac))
    # Las primeras 20 solo están medidas en un paciente, y con un valor enorme.
    t[:20, :] = np.nan
    t[:20, 0] = 50.0
    p[:20, :] = np.nan
    res = nulo_intercambio(
        t, p, n_elegir=30, n_permutaciones=100, semilla=20260914, min_pares=3
    )
    assert res["n_elegibles_en_los_dos_brazos"] == 180
    assert res["n_elegidas"] == 30
    # Sin el filtro, el resultado sería no_calculable porque el pólipo es NaN.
    assert res["observado_intermedio"] is not None


def test_el_nulo_de_intercambio_es_reproducible():
    rng = np.random.default_rng(5)
    t = rng.normal(size=(100, 16))
    p = rng.normal(size=(100, 16))
    kw = dict(n_elegir=20, n_permutaciones=50, semilla=20260914)
    a = nulo_intercambio(t, p, **kw)
    b = nulo_intercambio(t, p, **kw)
    assert a["p_empirico"] == b["p_empirico"]
    assert a["nulo_mediana"] == b["nulo_mediana"]
