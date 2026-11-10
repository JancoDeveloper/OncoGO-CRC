"""Orden de adquisición: los AppleDouble no cuentan y 54 crudos no cierran la prueba."""

from __future__ import annotations

import json

from oncogo_crc.stages.v3_5_preanalitica import contar_crudos, fecha_de_metadato, resumir_legibles


def test_el_conteo_ignora_appledouble():
    nombres = [
        "._CP_02_DIA.raw",
        "CP_02_DIA.raw",
        "CN_07_DIA.raw",
        "HC_15_DIA.raw",
        "nota.txt",
    ]
    validos, apple = contar_crudos(nombres)
    assert validos == ["CN_07_DIA.raw", "CP_02_DIA.raw", "HC_15_DIA.raw"]
    assert apple == 1


def test_el_analizador_lee_el_campo_fijado_en_la_enmienda():
    meta = {"FileProperties": [
        {"name": "Content Creation Date", "value": "10/02/2023 04:14:10"},
    ]}
    assert fecha_de_metadato(meta, "Content Creation Date") == "2023-10-02T04:14:10"


def test_con_54_crudos_el_estado_es_no_evaluable_y_las_cifras_existen():
    resumen = resumir_legibles(n_validos=60, n_con_fecha=54, auc=0.62, spearman={"piel": {"rho": 0.1, "p": 0.4}})
    assert resumen["estado"] == "no_evaluable"
    assert resumen["motivo"] == "54 de 60 crudos legibles; la prueba preregistrada pide 60"
    assert resumen["calculado_con_legibles"]["auc"] == 0.62
    assert resumen["calculado_con_legibles"]["n_con_fecha"] == 54
    json.dumps(resumen)
