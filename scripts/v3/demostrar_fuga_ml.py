"""Ruido con y sin filas duplicadas. Evidencia de la fuga; no es un test.

Semilla 822, la de enmienda_v3_3.semillas.desplazamientos.ml_demostracion_fuga.
"""

from __future__ import annotations

import json

import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from oncogo_crc.settings import load_settings
from oncogo_crc.v3.ml_suero import _tubo
from oncogo_crc.v3.salida import carpeta

SEMILLA = 822


def _auc(X: np.ndarray, y: np.ndarray) -> float:
    partes = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEMILLA)
    puntajes = []
    for tr, te in partes.split(X, y):
        tubo = _tubo(SEMILLA, "proteinas").fit(X[tr], y[tr])
        p = tubo.predict_proba(X[te])[:, 1]
        puntajes.append(float(roc_auc_score(y[te], p)))
    return float(np.mean(puntajes))


def main() -> None:
    rng = np.random.default_rng(SEMILLA)
    X = rng.normal(size=(40, 50))
    y = np.array([1] * 20 + [0] * 20)
    payload = {
        "semilla": SEMILLA,
        "forma": [40, 50],
        "auc_sin_duplicar": _auc(X, y),
        "auc_con_filas_duplicadas": _auc(np.vstack([X, X]), np.concatenate([y, y])),
        "nota": "Ruido normal. Duplicar las filas puede meter la misma muestra en el entrenamiento y en la prueba.",
    }
    destino = carpeta(load_settings()) / "ml_demostracion_fuga.json"
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(destino)
    print(payload["auc_sin_duplicar"], payload["auc_con_filas_duplicadas"])


if __name__ == "__main__":
    main()
