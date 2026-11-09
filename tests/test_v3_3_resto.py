"""Reglas fijadas antes de la corrida. No leen el report ni las matrices reales."""

import numpy as np
import pandas as pd

from oncogo_crc.v3.h5 import clase_de_proteina
from oncogo_crc.v3.h10 import clasificar_gen
from oncogo_crc.v3.k2k3 import filas_de_independencia, pendiente_deming
from oncogo_crc.v3.ml_suero import ImputarDentroDelPliegue


def test_deming_recupera_una_recta():
    x = np.arange(10, dtype=float)
    y = 2.0 * x + 1.0
    ajuste = pendiente_deming(x, y)
    assert ajuste["n"] == 10
    assert abs(ajuste["pendiente"] - 2.0) < 1e-8
    assert abs(ajuste["intercepto"] - 1.0) < 1e-8


def test_un_hash_compartido_impide_fusionar():
    filas = filas_de_independencia({
        "nivel_a_hashes_de_pride": [{"par": "A_vs_B", "n_hashes_compartidos": 15}],
        "nivel_b_correlacion_entre_muestras": {},
    })
    decisiones = {fila["nivel"]: fila["decision"] for fila in filas}
    assert decisiones["A_hash"] == "una_sola_cohorte"
    assert decisiones["B_correlacion"] == "una_sola_cohorte"


def test_pan_carcinoma_exige_cuatro_del_mismo_signo():
    otros = [(1.0, True)] * 4 + [(-1.0, True)] * 3
    assert clase_de_proteina(True, 1.0, otros) == "pan_carcinoma"
    assert clase_de_proteina(True, 1.0, [(1.0, False)] * 7) == "especifica_colon"
    assert clase_de_proteina(False, 1.0, [(1.0, True)] * 2 + [(1.0, False)] * 5) == "intermedia"


def test_h10_solo_es_especifica_si_las_demas_no_pasan():
    # Enmienda v3.3: la primera fila de colorrectal ya no decide. Este bloque,
    # sin las tres filas de Control, deja de ser "especifica".
    alpha = 0.05
    especifica = pd.DataFrame({
        "Disease": ["Colorectal cancer", "Lung cancer", "Myositis"],
        "Class": ["Cancer", "Cancer", "Autoimmune"],
        "padj": [0.01, 0.2, 0.4],
        "logFC": [1.0, 0.2, 0.2],
        "tipo_control": ["sanos", "sanos", "sanos"],
    })
    assert clasificar_gen(especifica, alpha, 0.585)["clase"] == "no_medida"


def test_la_imputacion_no_ve_el_valor_de_prueba():
    modelo = ImputarDentroDelPliegue()
    modelo.fit(np.array([[1.0, np.nan], [3.0, np.nan]]))
    salida = modelo.transform(np.array([[np.nan, 999.0]]))
    assert salida[0, 0] == 2.0
    assert salida[0, 1] == 0.0
