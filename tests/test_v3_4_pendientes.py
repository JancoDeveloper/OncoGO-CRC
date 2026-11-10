"""Piezas de v3.4 que no leen la cohorte."""

from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.v3.h6 import clase_de_expresion
from oncogo_crc.v3.h9 import _alinear, wald
from oncogo_crc.v3.rejilla_v34 import celdas_ml, celdas_suero


def test_la_rejilla_tiene_288_y_72():
    suero = celdas_suero()
    ml = celdas_ml()
    assert len(suero) == 288
    assert len(ml) == 72
    assert sum(celda["es_la_celda_preregistrada"] for celda in suero) == 1
    assert sum(celda["es_la_celda_preregistrada"] for celda in ml) == 1


def test_wald_de_un_instrumento():
    salida = wald(2.0, 0.1, 0.4, 0.2)
    assert salida["estado"] == "exploratorio"
    assert abs(salida["beta"] - 0.2) < 1e-12
    assert salida["se"] > 0


def test_wald_no_divide_por_cero():
    assert wald(0.0, 0.1, 0.4, 0.2)["estado"] == "no_evaluable"


def test_alinea_el_alelo_del_gwas_con_a1():
    assert _alinear("A", "G", "A", "G", 0.5, 0.1) == (0.5, 0.1)
    assert _alinear("G", "A", "A", "G", 0.5, 0.1) == (-0.5, 0.1)
    assert _alinear("C", "T", "A", "G", 0.5, 0.1) is None


def test_vesicula_exige_el_umbral_del_archivo_y_tumor():
    from oncogo_crc.v3.vesiculas import etiqueta, por_gen
    assert etiqueta(0.01, 0.01, True) == 1.0
    assert etiqueta(0.02, 0.01, True) == 0.0
    assert etiqueta(0.001, 0.01, False) == 0.0
    assert np.isnan(etiqueta(np.nan, 0.01, True))
    genes = por_gen(pd.DataFrame({"gen": ["A", "A", "B"], "etiqueta": [0.0, 1.0, 0.0]}))
    assert genes["A"] == 1.0 and genes["B"] == 0.0


def test_tmbed_cuenta_helice_y_senal():
    from oncogo_crc.v3.secuencia_features import helices_y_senal
    assert helices_y_senal("SS..HHH...hh..") == (2.0, 1.0)
    assert helices_y_senal("....") == (0.0, 0.0)


def test_la_pagina_v34_explica_la_corrida_sin_imprimir_none():
    import importlib.util
    base = Path(__file__).resolve().parents[1] / "results" / "v3_4"
    if not (base / "certificado_v3_4.json").is_file():
        return
    ruta = Path(__file__).resolve().parents[1] / "scripts" / "v3" / "exportar_sitio_v3_4.py"
    spec = importlib.util.spec_from_file_location("exportar_sitio_v3_4", ruta)
    modulo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(modulo)
    html = modulo.pagina(modulo.datos_publicos(base))
    assert "None" not in html
    assert "decidir sobre un paciente" in html
    assert 'id="recorrido"' in html
    assert 'id="conclusion"' in html


def test_el_sitio_v3_4_exige_certificado_publicable():
    from oncogo_crc.stages.v3_4_report import puede_publicar
    assert puede_publicar({"estado": "publicable", "faltan": []})
    assert not puede_publicar({"estado": "incompleto", "faltan": [{"id": "h7_variables"}]})
    assert not puede_publicar({"estado": "publicable", "faltan": [{"id": "h7_variables"}]})


def test_clase_h6_usa_el_doble():
    assert clase_de_expresion(4.0, 1.0) == "epitelial"
    assert clase_de_expresion(1.0, 4.0) == "no_epitelial"
    assert clase_de_expresion(3.0, 2.0) == "mixta"
