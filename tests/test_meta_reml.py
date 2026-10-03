"""Valida el meta-análisis de Python contra casos con resultado conocido (C10, E2).

No hay R instalado, así que la equivalencia con ``metafor::rma`` se comprueba
contra propiedades y contra casos con forma cerrada, no contra su salida.
"""

from __future__ import annotations

import numpy as np
import pytest

from oncogo_crc.v2.meta import (
    meta_efectos_aleatorios,
    tau2_dersimonian_laird,
    tau2_reml,
)


def test_sin_heterogeneidad_el_efecto_es_la_media_ponderada_fija():
    """Si los efectos son iguales, tau2 = 0 y sale la media de efectos fijos."""
    y = np.array([0.5, 0.5, 0.5])
    se = np.array([0.1, 0.2, 0.4])
    res = meta_efectos_aleatorios(y, se)
    assert res is not None
    assert res.tau2 == pytest.approx(0.0, abs=1e-9)
    assert res.i2 == pytest.approx(0.0, abs=1e-9)
    w = 1.0 / se**2
    assert res.efecto == pytest.approx(float(np.sum(w * y) / np.sum(w)), abs=1e-12)
    assert res.se == pytest.approx(float(np.sqrt(1.0 / np.sum(w))), abs=1e-12)


def test_dos_cohortes_identicas_dan_se_de_un_sobre_raiz_dos():
    """Caso con forma cerrada: dos efectos iguales con el mismo error estándar."""
    res = meta_efectos_aleatorios([1.0, 1.0], [0.2, 0.2])
    assert res is not None
    assert res.efecto == pytest.approx(1.0)
    assert res.se == pytest.approx(0.2 / np.sqrt(2.0), abs=1e-12)
    assert res.q == pytest.approx(0.0, abs=1e-12)


def test_heterogeneidad_grande_infla_tau2_y_el_intervalo():
    """Efectos opuestos con error pequeño: tau2 > 0, I² alto y el IC se ensancha."""
    homog = meta_efectos_aleatorios([0.80, 0.82, 0.78], [0.05, 0.05, 0.05])
    hetero = meta_efectos_aleatorios([-1.0, 0.0, 1.0], [0.05, 0.05, 0.05])
    assert homog is not None and hetero is not None
    assert hetero.tau2 > homog.tau2
    assert hetero.i2 > 0.9
    assert (hetero.ci_high - hetero.ci_low) > (homog.ci_high - homog.ci_low)


def test_reml_no_es_menor_que_cero_y_arranca_de_dersimonian_laird():
    rng = np.random.default_rng(20260914)
    for _ in range(50):
        y = rng.normal(0, 1, size=4)
        se = rng.uniform(0.05, 0.5, size=4)
        dl = tau2_dersimonian_laird(y, se**2)
        reml = tau2_reml(y, se**2)
        assert dl >= 0.0
        assert reml >= 0.0


def test_reml_converge_al_punto_fijo_de_su_ecuacion():
    """Metido en su propia ecuación, el tau2 de REML se reproduce."""
    y = np.array([0.4, 0.9, -0.2, 0.6])
    v = np.array([0.02, 0.05, 0.03, 0.04])
    tau2 = tau2_reml(y, v)
    w = 1.0 / (v + tau2)
    sw = w.sum()
    mu = float((w * y).sum() / sw)
    w2 = w**2
    recalculado = max(0.0, float((w2 * ((y - mu) ** 2 + 1.0 / sw - v)).sum() / w2.sum()))
    assert recalculado == pytest.approx(tau2, abs=1e-8)


def test_i2_coincide_con_higgins_thompson():
    y = np.array([0.2, 1.4, -0.3])
    se = np.array([0.1, 0.1, 0.1])
    res = meta_efectos_aleatorios(y, se)
    assert res is not None
    esperado = max(0.0, (res.q - res.q_df) / res.q)
    assert res.i2 == pytest.approx(esperado, abs=1e-12)
    assert res.q_df == 2


def test_una_sola_cohorte_no_es_meta_analisis():
    assert meta_efectos_aleatorios([0.5], [0.1]) is None
    # Un error estándar no positivo no entra: no aporta peso finito.
    assert meta_efectos_aleatorios([0.5, 0.7], [0.1, 0.0]) is None
    assert meta_efectos_aleatorios([0.5, np.nan], [0.1, 0.2]) is None


def test_el_modelo_fijo_ignora_la_heterogeneidad():
    y = np.array([-1.0, 1.0])
    se = np.array([0.1, 0.1])
    fijo = meta_efectos_aleatorios(y, se, metodo="fijo")
    aleatorio = meta_efectos_aleatorios(y, se, metodo="reml")
    assert fijo is not None and aleatorio is not None
    assert fijo.tau2 == 0.0
    # Con heterogeneidad, el modelo de efectos aleatorios reconoce más incertidumbre.
    assert aleatorio.se > fijo.se
