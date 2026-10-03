"""Tests de C6: el remuestreo y la validación se agrupan por paciente."""

from __future__ import annotations

import numpy as np
import pytest
from sklearn.model_selection import GroupKFold

from oncogo_crc.v2.modelos import (
    bootstrap_por_conglomerado,
    veredicto_sin_solape,
)


def _datos_por_paciente(n_pac=16, semilla=20260914):
    """Dos muestras por paciente: una tumor y una normal."""
    rng = np.random.default_rng(semilla)
    grupos = np.repeat(np.arange(n_pac), 2)
    y = np.tile([1, 0], n_pac)
    # Un puntaje informativo con efecto de paciente, que es lo que importa aquí.
    efecto_paciente = rng.normal(0, 1.0, size=n_pac)
    puntaje = y * 2.0 + np.repeat(efecto_paciente, 2) + rng.normal(0, 0.5, size=y.size)
    return y, grupos, puntaje


def test_el_bootstrap_remuestrea_pacientes_no_especimenes():
    y, grupos, s = _datos_por_paciente()
    res = bootstrap_por_conglomerado(
        y, grupos, {"m": s}, n_bootstrap=400, semilla=20260914
    )
    assert res["unidad_de_remuestreo"] == "paciente"
    assert res["n_pacientes"] == 16
    assert res["n_muestras"] == 32


def test_agrupar_por_paciente_cambia_el_intervalo():
    """La unidad de remuestreo importa, y en un diseño emparejado lo estrecha.

    El efecto de paciente se cancela al comparar su tumor con su normal.
    Remuestrear especímenes rompe esa cancelación y añade varianza que el diseño
    no tiene. No se afirma una dirección universal: se comprueba que cambia y se
    documenta la que corresponde a este diseño.
    """
    y, grupos, s = _datos_por_paciente()
    agrupado = bootstrap_por_conglomerado(
        y, grupos, {"m": s}, n_bootstrap=600, semilla=20260914
    )
    # El mismo cálculo tratando cada espécimen como su propio grupo.
    suelto = bootstrap_por_conglomerado(
        y, np.arange(y.size), {"m": s}, n_bootstrap=600, semilla=20260914
    )
    ancho_agrupado = agrupado["por_modelo"]["m"]["ci_high"] - agrupado["por_modelo"]["m"]["ci_low"]
    ancho_suelto = suelto["por_modelo"]["m"]["ci_high"] - suelto["por_modelo"]["m"]["ci_low"]
    assert ancho_agrupado != pytest.approx(ancho_suelto, abs=1e-6)
    assert ancho_agrupado < ancho_suelto


def test_la_diferencia_por_pares_se_calcula_en_la_misma_remuestra():
    y, grupos, s = _datos_por_paciente()
    rng = np.random.default_rng(1)
    peor = y * 0.5 + rng.normal(0, 1.5, size=y.size)
    res = bootstrap_por_conglomerado(
        y, grupos, {"bueno": s, "malo": peor}, n_bootstrap=400, semilla=20260914
    )
    dif = res["diferencias_por_pares"]["bueno_menos_malo"]
    assert dif["mediana"] > 0
    assert "incluye_cero" in dif
    assert dif["ci_low"] <= dif["mediana"] <= dif["ci_high"]


def test_dos_modelos_identicos_dan_diferencia_cero_sin_incertidumbre():
    y, grupos, s = _datos_por_paciente()
    res = bootstrap_por_conglomerado(
        y, grupos, {"a": s, "b": s.copy()}, n_bootstrap=200, semilla=20260914
    )
    dif = res["diferencias_por_pares"]["a_menos_b"]
    assert dif["mediana"] == pytest.approx(0.0, abs=1e-12)
    assert dif["incluye_cero"] is True


def test_el_veredicto_exige_intervalos_disjuntos():
    por_modelo = {
        "a": {"ci_low": 0.80, "ci_high": 0.99},
        "b": {"ci_low": 0.60, "ci_high": 0.85},
    }
    res = veredicto_sin_solape(por_modelo, {"a": 0.92, "b": 0.74})
    assert res["intervalos_disjuntos"] is False
    assert res["veredicto"].startswith("sin ganador")

    por_modelo_disjunto = {
        "a": {"ci_low": 0.90, "ci_high": 0.99},
        "b": {"ci_low": 0.50, "ci_high": 0.70},
    }
    res2 = veredicto_sin_solape(por_modelo_disjunto, {"a": 0.95, "b": 0.60})
    assert res2["intervalos_disjuntos"] is True
    assert res2["veredicto"] == "gana a"


def test_ningun_paciente_cae_en_entrenamiento_y_validacion():
    """La garantía de GroupKFold, comprobada sobre los grupos reales."""
    _y, grupos, _s = _datos_por_paciente()
    gkf = GroupKFold(n_splits=4)
    vistos = []
    for tr, va in gkf.split(np.zeros((grupos.size, 3)), np.tile([1, 0], 16), groups=grupos):
        assert set(grupos[tr]).isdisjoint(set(grupos[va]))
        vistos.append(set(grupos[va]))
    # Los pliegues cubren a todos los pacientes y no se repiten.
    union = set().union(*vistos)
    assert union == set(np.unique(grupos))
    assert sum(len(v) for v in vistos) == len(union)


def test_tamanos_desiguales_son_un_error():
    with pytest.raises(ValueError, match="mismo tamaño"):
        bootstrap_por_conglomerado(
            np.array([1, 0, 1]), np.array([0, 0]), {"m": np.array([1.0, 0.0, 1.0])},
            n_bootstrap=5, semilla=1,
        )


def test_el_bootstrap_es_reproducible():
    y, grupos, s = _datos_por_paciente()
    kw = dict(n_bootstrap=100, semilla=20260914)
    a = bootstrap_por_conglomerado(y, grupos, {"m": s}, **kw)
    b = bootstrap_por_conglomerado(y, grupos, {"m": s}, **kw)
    assert a["por_modelo"]["m"] == b["por_modelo"]["m"]
