"""H7 no entrena con un recurso que no tiene 10 observaciones por clase."""

from __future__ import annotations

import numpy as np
import pandas as pd

from oncogo_crc.v3.h7 import partir_recursos


def test_un_recurso_de_1_sobre_125_no_entra_en_el_entrenamiento():
    tabla = pd.DataFrame({
        "y_alvez": [1] * 53 + [0] * 73,
        "y_hpa": [1] * 53 + [0] * 63 + [np.nan] * 10,
        "y_papier": [1] + [0] * 125,
    })
    entrenables, no = partir_recursos(tabla, ["alvez", "hpa", "papier"])
    assert entrenables == ["alvez", "hpa"]
    assert no[0]["recurso"] == "papier"
    assert no[0]["n_pos"] == 1
    assert no[0]["n_neg"] == 125
    assert "papier" not in entrenables
