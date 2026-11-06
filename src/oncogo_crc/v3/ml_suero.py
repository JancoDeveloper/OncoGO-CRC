"""Prueba global de señal en el suero. No es un biomarcador.

La imputación y el escalado se ajustan dentro de cada pliegue. La permutación
repite ese mismo ajuste. Con 20 muestras por grupo esto pone a prueba el nulo,
no produce una cifra clínica.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, permutation_test_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from oncogo_crc.v3.prereg import PreregistroV3

N_SPLITS = 5
N_REPETICIONES = 20
N_PERMUTACIONES = 1000
PRESENCIA_MINIMA = 0.70


class ImputarDentroDelPliegue(BaseEstimator, TransformerMixin):
    """Mediana del entrenamiento. Una columna vista en menos del 70 % queda en cero."""

    def fit(self, X, y=None):
        valores = np.asarray(X, dtype=float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            self.mediana_ = np.nanmedian(valores, axis=0)
        self.fraccion_ = np.isfinite(valores).mean(axis=0)
        self.n_features_in_ = valores.shape[1]
        return self

    def transform(self, X):
        valores = np.array(np.asarray(X, dtype=float), copy=True)
        for j in range(valores.shape[1]):
            if self.fraccion_[j] < PRESENCIA_MINIMA or not np.isfinite(self.mediana_[j]):
                valores[:, j] = 0.0
            else:
                hueco = ~np.isfinite(valores[:, j])
                valores[hueco, j] = self.mediana_[j]
        return valores


def _tubo(semilla: int) -> Pipeline:
    return Pipeline([
        ("imputar", ImputarDentroDelPliegue()),
        ("escalar", StandardScaler()),
        ("modelo", LogisticRegressionCV(
            Cs=(0.1, 1.0, 10.0),
            cv=N_SPLITS,
            penalty="l2",
            solver="lbfgs",
            max_iter=500,
            scoring="roc_auc",
            random_state=semilla,
        )),
    ])


def _auc_repetida(X: np.ndarray, y: np.ndarray, semilla: int) -> dict:
    from sklearn.model_selection import cross_val_score
    tubo = _tubo(semilla)
    partes = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPETICIONES, random_state=semilla)
    puntajes = cross_val_score(tubo, X, y, cv=partes, scoring="roc_auc")
    return {
        "auc_media": float(np.mean(puntajes)),
        "ic_bajo": float(np.quantile(puntajes, 0.025)),
        "ic_alto": float(np.quantile(puntajes, 0.975)),
        "n_pliegues": int(puntajes.size),
        "intervalo": "cuantiles de la validación cruzada repetida, no un intervalo de DeLong",
    }


def _nulo(X: np.ndarray, y: np.ndarray, semilla: int) -> dict:
    tubo = _tubo(semilla)
    partes = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=semilla)
    estimacion, permutaciones, p = permutation_test_score(
        tubo,
        X,
        y,
        cv=partes,
        n_permutations=N_PERMUTACIONES,
        scoring="roc_auc",
        random_state=semilla,
        n_jobs=1,
    )
    return {
        "auc_de_la_permutacion": float(estimacion),
        "p": float(p),
        "media_del_nulo": float(np.mean(permutaciones)),
        "n_permutaciones": N_PERMUTACIONES,
    }


def evaluar(X: np.ndarray, y: np.ndarray, semilla_cv: int, semilla_nulo: int) -> dict:
    if len(np.unique(y)) < 2 or min(np.bincount(y.astype(int))) < N_SPLITS:
        return {"estado": "no_evaluable", "motivo": "un pliegue estratificado de 5 necesita al menos 5 por clase"}
    return {"estado": "calculado", **_auc_repetida(X, y, semilla_cv), **_nulo(X, y, semilla_nulo)}


def matrices_de_suero(principal: pd.DataFrame, humana: pd.DataFrame, calidad: pd.DataFrame, prereg: PreregistroV3) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Cáncer frente a colonoscopia negativa. Los sanos no entran en esta prueba."""
    grupos = {}
    for columnas in (principal.columns, humana.columns):
        for col in columnas:
            if col == "protein_id":
                continue
            if str(col).startswith("CP_"):
                grupos.setdefault("crc", []).append(col)
            elif str(col).startswith("CN_"):
                grupos.setdefault("colonoscopy_negative", []).append(col)
    columnas = grupos.get("crc", []) + grupos.get("colonoscopy_negative", [])
    y = np.array([1] * len(grupos.get("crc", [])) + [0] * len(grupos.get("colonoscopy_negative", [])))
    salida = {}
    for nombre, tabla in (("filtrada_humana", humana), ("sin_contaminantes_ni_paneles", principal)):
        mat = tabla.set_index("protein_id")[columnas].apply(pd.to_numeric, errors="coerce")
        salida[nombre] = (mat.to_numpy(dtype=float).T, y)
    coag = calidad.set_index("sample_id")["coagulacion"].reindex(columnas)
    salida["solo_coagulacion"] = (coag.to_numpy(dtype=float).reshape(-1, 1), y)
    return salida
