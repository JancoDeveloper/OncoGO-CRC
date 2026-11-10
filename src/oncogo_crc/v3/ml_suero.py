"""Prueba global de señal en el suero. No es un biomarcador.

La imputación, el residual y el escalado se ajustan dentro de cada pliegue.
La permutación repite ese mismo ajuste. Con 20 muestras por grupo esto pone
a prueba el nulo, no produce una cifra clínica.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, permutation_test_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

N_SPLITS = 5
N_REPETICIONES = 20
N_PERMUTACIONES = 1000
N_POR_CLASE = 20
PRESENCIA_MINIMA = 0.70
QC_COLUMNAS = ("eritrocitos", "plaquetas", "coagulacion", "piel", "log2_n_precursores")


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


class ResidualizarQC(BaseEstimator, TransformerMixin):
    """Las últimas k columnas son las QC. Ajusta, con el entrenamiento, la regresión de cada proteína sobre las QC y la quita."""

    def __init__(self, k_qc: int = len(QC_COLUMNAS)):
        self.k_qc = k_qc

    def _partir(self, X):
        valores = np.asarray(X, dtype=float)
        return valores[:, : -self.k_qc], valores[:, -self.k_qc :]

    def fit(self, X, y=None):
        prot, qc = self._partir(X)
        self.qc_media_ = np.nanmean(qc, axis=0)
        Q = np.column_stack([np.ones(len(qc)), np.nan_to_num(qc - self.qc_media_)])
        self.beta_ = np.zeros((Q.shape[1], prot.shape[1]))
        for j in range(prot.shape[1]):
            ok = np.isfinite(prot[:, j])
            if ok.sum() > Q.shape[1] + 2:
                self.beta_[:, j] = np.linalg.lstsq(Q[ok], prot[ok, j], rcond=None)[0]
        self.n_features_in_ = np.asarray(X).shape[1]
        return self

    def transform(self, X):
        prot, qc = self._partir(X)
        Q = np.nan_to_num(qc - self.qc_media_)
        return prot - Q @ self.beta_[1:]


def _tubo(semilla: int, tipo: str = "proteinas") -> Pipeline:
    pasos = []
    if tipo == "residualizado":
        pasos.append(("residualizar", ResidualizarQC()))
    pasos += [("imputar", ImputarDentroDelPliegue()), ("escalar", StandardScaler())]
    if tipo == "pca5":
        pasos.append(("pca", PCA(n_components=5, random_state=semilla)))
    pasos.append(("modelo", LogisticRegressionCV(
        Cs=(0.1, 1.0, 10.0), cv=N_SPLITS, penalty="l2", solver="lbfgs",
        max_iter=500, scoring="roc_auc", random_state=semilla,
    )))
    return Pipeline(pasos)


def ic_ledell(y: np.ndarray, s: np.ndarray, alfa: float = 0.05) -> tuple[float, float, float]:
    """AUC e intervalo por la curva de influencia (LeDell et al. 2015)."""
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    pos, neg = s[y == 1], s[y == 0]
    n1, n0, n = len(pos), len(neg), len(s)
    comp = (pos[:, None] > neg[None, :]) + 0.5 * (pos[:, None] == neg[None, :])
    auc = float(comp.mean())
    influencia = np.empty(n)
    influencia[y == 1] = (comp.mean(axis=1) - auc) * n / n1
    influencia[y == 0] = (comp.mean(axis=0) - auc) * n / n0
    se = float(np.sqrt(np.var(influencia) / n))
    z = float(norm.ppf(1 - alfa / 2))
    return auc, max(0.0, auc - z * se), min(1.0, auc + z * se)


def _auc_repetida(X, y, semilla, tipo) -> dict:
    partes = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPETICIONES, random_state=semilla)
    puntajes, oof = [], np.full((N_REPETICIONES, len(y)), np.nan)
    for k, (tr, te) in enumerate(partes.split(X, y)):
        tubo = _tubo(semilla, tipo).fit(X[tr], y[tr])
        p = tubo.predict_proba(X[te])[:, 1]
        puntajes.append(roc_auc_score(y[te], p))
        oof[k // N_SPLITS, te] = p
    auc_oof = [roc_auc_score(y, fila) for fila in oof]
    _, bajo, alto = ic_ledell(y, oof[0])
    return {
        "auc_oof_media": float(np.mean(auc_oof)),
        "ic95_ledell_bajo": bajo,
        "ic95_ledell_alto": alto,
        "ic95_nota": "LeDell 2015 sobre las predicciones fuera del pliegue de la primera repetición",
        "auc_media_por_pliegue": float(np.mean(puntajes)),
        "dispersion_pliegues_p2_5": float(np.quantile(puntajes, 0.025)),
        "dispersion_pliegues_p97_5": float(np.quantile(puntajes, 0.975)),
        "n_pliegues": len(puntajes),
    }


def _nulo(X, y, semilla, tipo, n_permutaciones: int = N_PERMUTACIONES) -> dict:
    partes = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=semilla)
    observada, permutaciones, p = permutation_test_score(
        _tubo(semilla, tipo), X, y, cv=partes, n_permutations=n_permutaciones,
        scoring="roc_auc", random_state=semilla, n_jobs=-1,
    )
    return {
        "auc_observada_5_pliegues": float(observada),
        "p_permutacion": float(p),
        "media_del_nulo": float(np.mean(permutaciones)),
        "percentil_95_del_nulo": float(np.quantile(permutaciones, 0.95)),
        "n_permutaciones": int(n_permutaciones),
    }


def evaluar(X, y, semilla_cv, semilla_nulo, tipo="proteinas", n_permutaciones: int = N_PERMUTACIONES) -> dict:
    if X.shape[0] != len(y):
        raise RuntimeError("X e y no tienen el mismo número de filas")
    if len(np.unique(y)) < 2 or min(np.bincount(y.astype(int))) < N_SPLITS:
        return {"estado": "no_evaluable", "motivo": "un pliegue estratificado de 5 necesita al menos 5 por clase"}
    return {
        "estado": "calculado",
        "tipo": tipo,
        "n": int(len(y)),
        "n_variables": int(X.shape[1]),
        **_auc_repetida(X, y, semilla_cv, tipo),
        **_nulo(X, y, semilla_nulo, tipo, n_permutaciones),
    }


def _columnas_cp_cn(*tablas: pd.DataFrame) -> tuple[list[str], list[str]]:
    todas = list(dict.fromkeys(str(c) for t in tablas for c in t.columns if c != "protein_id"))
    crc = [c for c in todas if c.startswith("CP_")]
    cn = [c for c in todas if c.startswith("CN_")]
    for tabla in tablas:
        faltan = sorted(set(crc + cn) - set(map(str, tabla.columns)))
        if faltan:
            raise RuntimeError(f"una matriz no trae las columnas {faltan[:5]}")
    if len(crc) != N_POR_CLASE or len(cn) != N_POR_CLASE:
        raise RuntimeError(f"se esperaban {N_POR_CLASE} CP y {N_POR_CLASE} CN; hay {len(crc)} y {len(cn)}")
    return crc, cn


def _log2_si_hace_falta(mat: np.ndarray) -> np.ndarray:
    finitos = mat[np.isfinite(mat)]
    if finitos.size and np.nanmax(finitos) > 100:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(mat > 0, np.log2(mat), np.nan)
    return mat


def matrices_de_suero(principal, humana, calidad, prereg=None) -> dict[str, tuple[np.ndarray, np.ndarray, str]]:
    """Cáncer frente a colonoscopia negativa. Una fila por muestra, sin repetir."""
    del prereg
    crc, cn = _columnas_cp_cn(principal, humana)
    columnas = crc + cn
    if len(set(columnas)) != len(columnas) or len(columnas) != 2 * N_POR_CLASE:
        raise RuntimeError("hay muestras repetidas en el ML")
    y = np.array([1] * len(crc) + [0] * len(cn))

    def _x(tabla):
        mat = tabla.assign(protein_id=tabla["protein_id"].astype(str)).set_index("protein_id")[columnas]
        return _log2_si_hace_falta(mat.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float).T)

    qc_tabla = calidad.assign(sample_id=calidad["sample_id"].astype(str)).set_index("sample_id")
    qc_tabla = qc_tabla.assign(log2_n_precursores=np.log2(pd.to_numeric(qc_tabla["n_precursores"], errors="coerce")))
    qc = qc_tabla.reindex(columnas)[list(QC_COLUMNAS)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if np.isnan(qc).all(axis=1).any():
        raise RuntimeError("hay muestras sin puntuación de calidad")
    xh, xp = _x(humana), _x(principal)
    coag = qc[:, [QC_COLUMNAS.index("coagulacion")]]
    return {
        "M1_proteinas_filtrada_humana": (xh, y, "proteinas"),
        "M2_proteinas_principal": (xp, y, "proteinas"),
        "M3_solo_qc": (qc, y, "proteinas"),
        "M4_solo_coagulacion": (coag, y, "proteinas"),
        "M5_principal_residualizada_qc": (np.hstack([xp, qc]), y, "residualizado"),
        "M6_principal_pca5": (xp, y, "pca5"),
    }


def titular(versiones: dict, alfa: float, auc_minima: float, textos: dict) -> dict:
    m2, m5 = versiones["M2_proteinas_principal"], versiones["M5_principal_residualizada_qc"]
    if m2.get("estado") != "calculado" or m2["p_permutacion"] >= alfa:
        clave = "sin_separacion_detectable"
    elif m5.get("estado") == "calculado" and m5["p_permutacion"] < alfa and m5["auc_oof_media"] >= auc_minima:
        clave = "senal_mas_alla_de_la_calidad"
    else:
        clave = "separacion_explicada_por_calidad"
    return {"clave": clave, "texto": textos[clave], "regla": "enmienda_v3_3.ml_suero.titular"}
