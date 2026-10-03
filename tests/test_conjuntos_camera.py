"""Valida la prueba competitiva de conjuntos y su nulo (C4, C13, E2).

Las matrices de estos tests son sintéticas y viven solo aquí, como exige el
preregistro: ningún dato simulado entra en results/.
"""

from __future__ import annotations

import numpy as np
import pytest

from oncogo_crc.v2.conjuntos import (
    camera,
    camera_un_conjunto,
    correlacion_intergenica,
    nulo_etiquetas_de_gen,
)


def test_conjunto_desplazado_sale_en_la_direccion_correcta():
    rng = np.random.default_rng(20260914)
    stats = rng.normal(0, 1, size=400)
    sel = np.zeros(400, dtype=bool)
    sel[:40] = True
    stats[sel] += 2.0
    res = camera_un_conjunto(stats, sel, termino="T", rho=0.0, df=90)
    assert res is not None
    assert res.direccion == 1
    assert res.p_value < 1e-6
    assert res.n_en_conjunto == 40
    assert res.n_fuera == 360


def test_conjunto_sin_senal_no_sale_significativo():
    rng = np.random.default_rng(1)
    stats = rng.normal(0, 1, size=400)
    sel = np.zeros(400, dtype=bool)
    sel[:40] = True
    res = camera_un_conjunto(stats, sel, termino="T", rho=0.0, df=90)
    assert res is not None
    assert res.p_value > 0.05


def test_la_correlacion_vuelve_la_prueba_mas_conservadora():
    """Es el punto de camera: con rho > 0 el mismo efecto da un p mayor."""
    rng = np.random.default_rng(7)
    stats = rng.normal(0, 1, size=300)
    sel = np.zeros(300, dtype=bool)
    sel[:30] = True
    stats[sel] += 0.8
    sin_corr = camera_un_conjunto(stats, sel, termino="T", rho=0.0, df=90)
    con_corr = camera_un_conjunto(stats, sel, termino="T", rho=0.05, df=90)
    assert sin_corr is not None and con_corr is not None
    assert con_corr.vif > sin_corr.vif == 1.0
    assert con_corr.p_value > sin_corr.p_value
    assert abs(con_corr.estadistico) < abs(sin_corr.estadistico)


def test_el_vif_es_uno_mas_m_menos_uno_por_rho():
    stats = np.arange(100, dtype=float)
    sel = np.zeros(100, dtype=bool)
    sel[:10] = True
    res = camera_un_conjunto(stats, sel, termino="T", rho=0.1, df=50)
    assert res is not None
    assert res.vif == pytest.approx(1.0 + (10 - 1) * 0.1)


def test_conjuntos_demasiado_pequenos_se_descartan():
    stats = np.arange(100, dtype=float)
    pert = {"grande": np.zeros(100, dtype=bool), "chico": np.zeros(100, dtype=bool)}
    pert["grande"][:20] = True
    pert["chico"][:4] = True
    res = camera(stats, pert, rho=0.0, df=50, min_proteinas=10)
    assert [r.termino for r in res] == ["grande"]


def test_correlacion_intergenica_detecta_bloques_y_no_baja_de_cero():
    rng = np.random.default_rng(3)
    # Diez proteínas que comparten un factor latente: correlación positiva.
    factor = rng.normal(0, 1, size=30)
    bloque = np.array([factor + rng.normal(0, 0.3, size=30) for _ in range(10)])
    rho_bloque = correlacion_intergenica(bloque)
    independiente = rng.normal(0, 1, size=(10, 30))
    rho_indep = correlacion_intergenica(independiente)
    assert rho_bloque > 0.5
    assert rho_indep >= 0.0
    assert rho_bloque > rho_indep


def test_nulo_de_etiquetas_respeta_el_piso_del_p():
    rng = np.random.default_rng(11)
    stats = rng.normal(0, 1, size=200)
    sel = np.zeros(200, dtype=bool)
    sel[:20] = True
    stats[sel] += 3.0
    n_perm = 200
    ps = nulo_etiquetas_de_gen(
        stats,
        {"T": sel},
        rho=0.0,
        df=90,
        min_proteinas=10,
        n_permutaciones=n_perm,
        semilla=20260914,
    )
    assert set(ps) == {"T"}
    # C1: ningún p de permutación puede bajar del piso.
    assert ps["T"] >= 1.0 / (n_perm + 1)
    # Un desplazamiento de 3 sigmas tiene que quedar en el extremo del nulo.
    assert ps["T"] <= 2.0 / (n_perm + 1)


def test_nulo_de_etiquetas_no_marca_un_conjunto_sin_senal():
    rng = np.random.default_rng(13)
    stats = rng.normal(0, 1, size=200)
    sel = np.zeros(200, dtype=bool)
    sel[:20] = True
    ps = nulo_etiquetas_de_gen(
        stats,
        {"T": sel},
        rho=0.0,
        df=90,
        min_proteinas=10,
        n_permutaciones=200,
        semilla=20260914,
    )
    assert ps["T"] > 0.05


def test_el_nulo_es_reproducible_con_la_misma_semilla():
    rng = np.random.default_rng(17)
    stats = rng.normal(0, 1, size=150)
    sel = np.zeros(150, dtype=bool)
    sel[:15] = True
    kw = dict(rho=0.0, df=60, min_proteinas=10, n_permutaciones=100, semilla=20260914)
    assert nulo_etiquetas_de_gen(stats, {"T": sel}, **kw) == nulo_etiquetas_de_gen(
        stats, {"T": sel}, **kw
    )
