import pandas as pd
import pytest

from oncogo_crc.v3.h5 import clase_de_proteina
from oncogo_crc.v3.h10 import clasificar_gen, tipo_de_control


def _fila(control, logfc, padj):
    return {
        "Disease": "Colorectal cancer",
        "Control": control,
        "tipo_control": tipo_de_control(control),
        "logFC": logfc,
        "padj": padj,
    }


def test_siete_carcinomas_dan_lo_mismo_que_en_v3():
    assert clase_de_proteina(True, 1.0, [(1.0, True)] * 7) == "pan_carcinoma"
    assert clase_de_proteina(True, 1.0, [(1.0, True)] * 4 + [(-1.0, True)] * 3) == "pan_carcinoma"


def test_seis_con_cuatro_concordantes_es_pan_carcinoma():
    otros = [(1.0, True)] * 4 + [(-1.0, True)] * 2
    assert clase_de_proteina(True, 1.0, otros) == "pan_carcinoma"


def test_cinco_carcinomas_no_se_clasifican():
    assert clase_de_proteina(True, 1.0, [(1.0, True)] * 5) == "no_clasificable"


def test_h10_clases_con_los_tres_controles():
    alpha, lfc = 0.05, 0.585
    base = [
        _fila("Healthy", 1.0, 0.01),
        _fila("Class", 1.0, 0.01),
        _fila("All other diseases", 1.0, 0.01),
    ]
    assert clasificar_gen(pd.DataFrame(base), alpha, lfc)["clase"] == "especifica_crc"
    solo_enfermedades = base.copy()
    solo_enfermedades[1] = _fila("Class", 0.1, 0.2)
    assert clasificar_gen(pd.DataFrame(solo_enfermedades), alpha, lfc)["clase"] == "especifica_frente_a_enfermedades"
    solo_canceres = base.copy()
    solo_canceres[2] = _fila("All other diseases", 0.1, 0.2)
    assert clasificar_gen(pd.DataFrame(solo_canceres), alpha, lfc)["clase"] == "especifica_frente_a_canceres"
    ninguna = base.copy()
    ninguna[1] = _fila("Class", -1.0, 0.01)
    ninguna[2] = _fila("All other diseases", -1.0, 0.01)
    assert clasificar_gen(pd.DataFrame(ninguna), alpha, lfc)["clase"] == "no_especifica"
    sanos = base.copy()
    sanos[0] = _fila("Healthy", 0.1, 0.2)
    assert clasificar_gen(pd.DataFrame(sanos), alpha, lfc)["clase"] == "no_elevada_frente_a_sanos"


def test_un_control_desconocido_detiene_la_etapa():
    with pytest.raises(ValueError):
        tipo_de_control("adjacent tissue")
