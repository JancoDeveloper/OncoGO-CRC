import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer

import oncogo_crc.v3.h7 as h7


def _tabla(y, semilla=0):
    rng = np.random.default_rng(semilla)
    n = len(y)
    tabla = pd.DataFrame({
        "v": y.astype(float) + rng.normal(scale=0.01, size=n),
        "y_a": y,
        "y_b": y,
    }, index=[f"P{i}" for i in range(n)])
    return tabla


def test_etiquetas_aleatorias_no_sostienen(monkeypatch):
    monkeypatch.setattr(h7, "N_BOOTSTRAP", 200)
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, size=80)
    p = rng.random(size=80)
    por = {"a": h7.auc_bootstrap(y, p, 1), "b": h7.auc_bootstrap(y, rng.random(size=80), 2)}
    assert h7.veredicto(por) != "sostiene"


def test_una_variable_que_determina_la_etiqueta_sostiene(monkeypatch):
    monkeypatch.setattr(h7, "N_BOOTSTRAP", 200)
    y = np.array([1] * 40 + [0] * 40)
    tabla = _tabla(y, semilla=1)
    por = {}
    for recurso, otro in (("a", "b"), ("b", "a")):
        pred = h7.fuera_de_recurso(tabla, ["a", "b"], recurso, ["v"], [], 3)
        por[recurso] = h7.auc_bootstrap(tabla[f"y_{recurso}"].to_numpy(), pred.to_numpy(), 4)
    assert h7.veredicto(por) == "sostiene"


def test_la_pca_se_ajusta_dentro_del_tubo(monkeypatch):
    vistos = {}

    class Espia(ColumnTransformer):
        def fit_transform(self, X, y=None, **kwargs):
            vistos["filas"] = list(getattr(X, "index", range(len(X))))
            return super().fit_transform(X, y, **kwargs)

    monkeypatch.setattr(h7, "ColumnTransformer", Espia)
    y = np.array([1] * 20 + [0] * 20)
    X = pd.DataFrame({"v": y.astype(float), **{f"e{i}": np.random.default_rng(i).normal(size=40) for i in range(6)}})
    modelo = h7.tubo(["v"], [f"e{i}" for i in range(6)], 0)
    modelo.fit(X.iloc[:30], y[:30])
    assert len(vistos["filas"]) == 30
    assert 39 not in vistos["filas"]


def test_auc_bootstrap_pide_diez_por_clase():
    y = np.array([1] * 5 + [0] * 5)
    p = np.linspace(0, 1, 10)
    assert h7.auc_bootstrap(y, p, 0)["estado"] == "no_evaluable"
