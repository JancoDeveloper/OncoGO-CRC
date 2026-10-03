"""Tests de C8: la ventana de detección separa «no se mide» de «no está»."""

from __future__ import annotations

import numpy as np
import pandas as pd

from oncogo_crc.v2 import deteccion


def _universo(n=400, semilla=20260914):
    """Concentraciones repartidas en seis órdenes de magnitud."""
    rng = np.random.default_rng(semilla)
    ids = [f"P{i:04d}" for i in range(n)]
    log10 = rng.uniform(3.0, 9.0, size=n)
    return ids, pd.Series(10.0**log10, index=ids), log10


def test_la_curva_recupera_una_dependencia_de_la_concentracion():
    ids, conc, log10 = _universo()
    # Detección determinista por encima de 10^6 pg/L, con un poco de ruido.
    rng = np.random.default_rng(1)
    p = 1.0 / (1.0 + np.exp(-(log10 - 6.0) * 2.0))
    detectada = pd.Series(rng.random(len(ids)) < p, index=ids)
    res = deteccion.curva_de_deteccion(conc, detectada, n_min=50)
    assert res["estado"] == "calculado"
    assert res["pendiente_por_orden_de_magnitud"] > 1.0
    assert res["z_pendiente"] > 4.0
    # El umbral tiene que caer cerca del punto medio real.
    assert 5.5 < res["umbral_log10_pg_l"] < 6.5


def test_sin_dependencia_la_pendiente_no_se_separa_de_cero():
    ids, conc, _ = _universo()
    rng = np.random.default_rng(2)
    detectada = pd.Series(rng.random(len(ids)) < 0.5, index=ids)
    res = deteccion.curva_de_deteccion(conc, detectada, n_min=50)
    assert res["estado"] == "calculado"
    assert abs(res["z_pendiente"]) < 2.0


def test_pocas_proteinas_dejan_la_curva_no_evaluable():
    ids = ["A", "B", "C"]
    conc = pd.Series([1e4, 1e6, 1e8], index=ids)
    det = pd.Series([False, True, True], index=ids)
    res = deteccion.curva_de_deteccion(conc, det, n_min=50)
    assert res["estado"] == "no_evaluable"
    assert res["n"] == 3


def test_todo_detectado_deja_la_curva_no_evaluable():
    ids, conc, _ = _universo()
    det = pd.Series(True, index=ids)
    res = deteccion.curva_de_deteccion(conc, det, n_min=50)
    assert res["estado"] == "no_evaluable"
    assert "del mismo lado" in res["motivo"]


def test_la_ventana_rescata_una_concordancia_escondida_por_la_deteccion():
    """Caso construido: concuerdan las abundantes y las escasas son ruido."""
    rng = np.random.default_rng(3)
    n_alto, n_bajo = 60, 120
    ids = [f"A{i}" for i in range(n_alto)] + [f"B{i}" for i in range(n_bajo)]
    conc = pd.Series(
        list(10.0 ** rng.uniform(8.0, 9.0, n_alto)) + list(10.0 ** rng.uniform(3.0, 4.0, n_bajo)),
        index=ids,
    )
    tejido = pd.Series(rng.normal(size=len(ids)), index=ids)
    plasma = tejido.copy()
    # Las abundantes concuerdan; las escasas van al azar.
    plasma.iloc[n_alto:] = rng.normal(size=n_bajo)
    res = deteccion.concordancia_dentro_de_la_ventana(
        tejido, plasma, conc, umbral_log10=7.0, n_min=10
    )
    assert res["estado"] == "calculado"
    assert res["dentro_de_la_ventana"]["valor"] == 1.0
    assert res["fuera_de_la_ventana"]["valor"] < 0.7
    assert res["todas"]["valor"] < 1.0
    assert res["veredicto"]["la_ventana_explica_la_brecha"] is True


def test_si_dentro_de_la_ventana_sigue_siendo_azar_la_brecha_es_real():
    """El caso de los datos reales: restringir no cambia nada.

    Se repite con varias semillas porque con ruido puro un intervalo del 95 %
    excluye el 0,5 una de cada veinte veces. La propiedad que se comprueba es
    que eso sea la excepción, no que no ocurra nunca.
    """
    azar = 0
    for semilla in range(20):
        rng = np.random.default_rng(semilla)
        ids = [f"P{i}" for i in range(200)]
        conc = pd.Series(10.0 ** rng.uniform(3.0, 9.0, 200), index=ids)
        tejido = pd.Series(rng.normal(size=200), index=ids)
        plasma = pd.Series(rng.normal(size=200), index=ids)
        res = deteccion.concordancia_dentro_de_la_ventana(
            tejido, plasma, conc, umbral_log10=6.0, n_min=10
        )
        azar += int(res["veredicto"]["dentro_de_la_ventana_es_compatible_con_azar"])
        if res["veredicto"]["dentro_de_la_ventana_es_compatible_con_azar"]:
            assert res["veredicto"]["la_ventana_explica_la_brecha"] is False
    assert azar >= 17, f"solo {azar} de 20 repeticiones de ruido salieron compatibles con azar"


def test_sin_umbral_la_ventana_queda_no_evaluable():
    ids = [f"P{i}" for i in range(50)]
    conc = pd.Series(np.full(50, 1e6), index=ids)
    t = pd.Series(np.ones(50), index=ids)
    res = deteccion.concordancia_dentro_de_la_ventana(
        t, t, conc, umbral_log10=None, n_min=10
    )
    assert res["estado"] == "no_evaluable"
    assert "umbral" in res["motivo"]


def test_la_tabla_por_proteina_trae_una_fila_por_identificador():
    ids = ["A", "B", "C"]
    t = pd.Series([1.0, -1.0, 2.0], index=ids)
    p = pd.Series([1.0, 1.0, np.nan], index=ids)
    conc = pd.Series([1e8, 1e4, np.nan], index=ids)
    det = pd.Series([True, False, False], index=ids)
    tabla = deteccion.tabla_por_proteina(t, p, conc, det, umbral_log10=6.0)
    assert len(tabla) == 3
    assert list(tabla["dentro_de_la_ventana"]) == [True, False, False]
    assert list(tabla["mismo_signo"]) == [True, False, None]


def test_las_proteinas_sin_concentracion_no_entran_en_la_curva():
    ids = [f"P{i}" for i in range(100)]
    conc = pd.Series([1e6] * 50 + [np.nan] * 50, index=ids)
    det = pd.Series([True] * 25 + [False] * 25 + [True] * 50, index=ids)
    res = deteccion.curva_de_deteccion(conc, det, n_min=10)
    assert res["n"] == 50
