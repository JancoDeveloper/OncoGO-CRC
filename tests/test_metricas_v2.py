"""Tests de las correcciones C1, C3 y C15 sobre fracciones y p de permutación."""

from __future__ import annotations

import pytest

from oncogo_crc.v2.metricas import Fraccion, p_empirico


# ------------------------------------------------------------------- C3
def test_una_fraccion_no_existe_sin_numerador_y_denominador():
    with pytest.raises(TypeError):
        Fraccion(nombre="r2")  # type: ignore[call-arg]


def test_el_valor_se_deriva_y_no_se_pasa_a_mano():
    f = Fraccion(
        nombre="R2_detectabilidad",
        numerador=226,
        denominador=688,
        unidad_denominador="proteinas_unicas",
    )
    assert f.valor == pytest.approx(226 / 688, abs=1e-12)
    d = f.to_dict()
    assert d["numerador"] == 226
    assert d["denominador"] == 688
    assert d["valor"] == pytest.approx(d["numerador"] / d["denominador"], abs=1e-9)


def test_la_unidad_del_denominador_es_obligatoria_y_restringida():
    """El bug de v1 fue contar filas del cruce como si fueran proteínas."""
    with pytest.raises(ValueError, match="no válida"):
        Fraccion(
            nombre="mala",
            numerador=119,
            denominador=702,
            unidad_denominador="filas_del_cruce",
        )


def test_numerador_mayor_que_denominador_es_un_error():
    with pytest.raises(ValueError, match="numerador"):
        Fraccion(nombre="mala", numerador=10, denominador=5, unidad_denominador="proteinas_unicas")


def test_denominador_cero_no_inventa_un_valor():
    f = Fraccion(nombre="vacia", numerador=0, denominador=0, unidad_denominador="terminos")
    assert f.valor is None
    assert f.intervalo() == (None, None)
    assert f.p_vs_azar() is None


# ------------------------------------------------------------------ C15
def test_r1_reporta_la_cobertura_sobre_los_hallazgos():
    """R1 de v1: 104/119, y 119 es el 17,3 % de las 688 certificadas."""
    r1 = Fraccion(
        nombre="R1_direccion",
        numerador=104,
        denominador=119,
        unidad_denominador="proteinas_unicas",
        universo=688,
        universo_nombre="n_certificadas",
    )
    assert r1.valor == pytest.approx(104 / 119, abs=1e-12)
    assert r1.cobertura == pytest.approx(119 / 688, abs=1e-12)
    assert r1.cobertura < 0.50  # umbral del preregistro para no generalizar
    d = r1.to_dict()
    assert d["universo"] == 688
    assert d["cobertura"] == pytest.approx(119 / 688, abs=1e-12)


def test_sin_universo_no_se_inventa_cobertura():
    f = Fraccion(nombre="x", numerador=1, denominador=2, unidad_denominador="proteinas_unicas")
    assert f.cobertura is None
    assert "cobertura" not in f.to_dict()


def test_el_intervalo_de_r1_reproduce_el_de_v1():
    r1 = Fraccion(
        nombre="R1_direccion",
        numerador=104,
        denominador=119,
        unidad_denominador="proteinas_unicas",
    )
    low, high = r1.intervalo()
    assert low == pytest.approx(0.8006330017494038, abs=1e-9)
    assert high == pytest.approx(0.9277026405016786, abs=1e-9)


# ------------------------------------------------------------------- C1
def test_ningun_p_de_permutacion_baja_del_piso():
    for b in range(0, 6):
        for n in (100, 1000, 10000, 65536):
            res = p_empirico(b, n)
            assert res["p"] >= 1.0 / (n + 1)
            assert res["p"] == pytest.approx((b + 1) / (n + 1), abs=1e-15)


def test_cero_aciertos_se_declara_como_tope_no_como_estimacion():
    """El p de v1 era 1/101 = 0,0099: el mínimo alcanzable con 100 permutaciones."""
    v1 = p_empirico(0, 100)
    assert v1["p"] == pytest.approx(0.009900990099009901, abs=1e-15)
    assert v1["es_tope"] is True
    assert v1["p"] == v1["p_minimo_alcanzable"]

    v2 = p_empirico(0, 10000)
    assert v2["es_tope"] is True
    assert v2["p"] < v1["p"]  # más permutaciones bajan el piso


def test_el_p_de_permutacion_trae_intervalo():
    res = p_empirico(3, 10000)
    assert res["es_tope"] is False
    assert res["ci_low"] is not None and res["ci_high"] is not None
    assert res["ci_low"] <= res["p"] <= res["ci_high"]


def test_permutaciones_no_positivas_son_un_error():
    with pytest.raises(ValueError):
        p_empirico(0, 0)
