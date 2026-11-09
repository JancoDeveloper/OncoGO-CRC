import numpy as np
import pandas as pd

import oncogo_crc.v3.ml_suero as ml


def _tablas():
    cols = ["protein_id"] + [f"CP_{i}" for i in range(20)] + [f"CN_{i}" for i in range(20)] + [f"Normal_{i}" for i in range(20)]
    rng = np.random.default_rng(0)
    tabla = pd.DataFrame([[f"P{j}"] + list(rng.lognormal(20, 1, 60)) for j in range(30)], columns=cols)
    calidad = pd.DataFrame({
        "sample_id": cols[1:],
        "eritrocitos": 0.0,
        "plaquetas": 0.0,
        "coagulacion": rng.normal(size=60),
        "piel": 0.0,
        "n_precursores": 1000,
    })
    return tabla, calidad


def test_cada_muestra_entra_una_vez():
    tabla, calidad = _tablas()
    for nombre, (X, y, _tipo) in ml.matrices_de_suero(tabla, tabla.copy(), calidad).items():
        assert X.shape[0] == 40, nombre
        assert int(y.sum()) == 20, nombre


def test_la_matriz_entra_en_log2():
    tabla, calidad = _tablas()
    X, _, _ = ml.matrices_de_suero(tabla, tabla.copy(), calidad)["M2_proteinas_principal"]
    assert np.nanmax(X) < 40


def test_ruido_no_separa(monkeypatch):
    monkeypatch.setattr(ml, "N_REPETICIONES", 2)
    monkeypatch.setattr(ml, "N_PERMUTACIONES", 49)
    rng = np.random.default_rng(1)
    X = rng.normal(size=(40, 50))
    y = np.array([1] * 20 + [0] * 20)
    r = ml.evaluar(X, y, 1, 2)
    assert 0.2 < r["auc_oof_media"] < 0.8
    assert r["p_permutacion"] > 0.01


def test_ledell_contiene_la_auc():
    y = np.array([1] * 20 + [0] * 20)
    s = np.r_[np.linspace(0.3, 1, 20), np.linspace(0, 0.7, 20)]
    auc, bajo, alto = ml.ic_ledell(y, s)
    assert bajo <= auc <= alto
