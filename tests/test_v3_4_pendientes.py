"""Piezas de v3.4 que no leen la cohorte."""

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


def test_clase_h6_usa_el_doble():
    assert clase_de_expresion(4.0, 1.0) == "epitelial"
    assert clase_de_expresion(1.0, 4.0) == "no_epitelial"
    assert clase_de_expresion(3.0, 2.0) == "mixta"
