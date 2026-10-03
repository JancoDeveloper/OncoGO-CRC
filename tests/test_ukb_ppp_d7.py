"""Tests de D7: pQTL del plasma en UK Biobank.

Esta fuente requiere cuenta y queda marcada ``requiere_acceso``. Estos tests son
la razón por la que «script preparado» significa algo: el análisis se ejecuta
entero sobre estadísticas de resumen sintéticas, así que el día que alguien con
acceso descargue el fichero real no hay nada que programar ni que depurar.

Los datos son sintéticos y viven solo aquí.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pandas as pd
import pytest

RAIZ = Path(__file__).resolve().parents[1]


def _modulo():
    ruta = RAIZ / "scripts" / "download" / "ukb_ppp_pqtl.py"
    spec = importlib.util.spec_from_file_location("ukb_ppp_pqtl", ruta)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ukb_ppp_pqtl"] = mod
    spec.loader.exec_module(mod)
    return mod


U = _modulo()


def _stats(filas):
    return pd.DataFrame(
        filas,
        columns=["protein", "variant", "cis_trans", "beta", "se", "log10p", "maf"],
    )


# --------------------------------------------------- varianza explicada
def test_la_varianza_explicada_sigue_la_formula_aditiva():
    assert U.varianza_explicada(0.5, 0.3) == pytest.approx(2 * 0.3 * 0.7 * 0.25)
    # Un efecto nulo no explica nada, por frecuente que sea la variante.
    assert U.varianza_explicada(0.0, 0.5) == 0.0


def test_el_signo_del_efecto_no_cambia_la_varianza():
    assert U.varianza_explicada(0.4, 0.2) == U.varianza_explicada(-0.4, 0.2)


def test_una_frecuencia_imposible_no_devuelve_un_numero():
    for maf in (0.0, 1.0, -0.1, 1.5, float("nan")):
        assert math.isnan(U.varianza_explicada(0.5, maf))
    assert math.isnan(U.varianza_explicada(float("nan"), 0.3))


def test_la_varianza_es_maxima_con_la_variante_mas_equilibrada():
    """A igual efecto, una variante al 50 % mueve más el nivel que una rara."""
    assert U.varianza_explicada(0.5, 0.5) > U.varianza_explicada(0.5, 0.05)


# --------------------------------------------------- resumen por proteína
def test_se_queda_con_el_cis_pqtl_que_mas_varianza_explica_no_con_el_de_menor_p():
    """Un p minúsculo en 54 000 personas puede explicar casi nada del nivel."""
    d = _stats(
        [
            # p mucho menor, pero efecto pequeño: explica poco.
            ["SAA1", "rs_p_minusculo", "cis", 0.05, 0.001, 300.0, 0.5],
            # p mayor, pero efecto grande: explica mucho.
            ["SAA1", "rs_efecto_grande", "cis", 0.90, 0.05, 50.0, 0.4],
        ]
    )
    r = U.resumir_por_proteina(d).iloc[0]
    assert r["variante"] == "rs_efecto_grande"
    assert bool(r["determinada_por_genetica"]) is True


def test_las_variantes_que_no_pasan_el_umbral_del_articulo_no_entran():
    d = _stats(
        [
            ["SAA1", "rs_flojo", "cis", 0.9, 0.05, 3.0, 0.4],
            ["LBP", "rs_fuerte", "cis", 0.9, 0.05, 50.0, 0.4],
        ]
    )
    res = U.resumir_por_proteina(d)
    assert list(res["proteina"]) == ["LBP"]


def test_un_trans_pqtl_no_se_usa_como_el_mejor_pero_se_cuenta():
    """El cis es el que habla del nivel propio de la proteína; el trans se anota."""
    d = _stats(
        [
            ["SAA1", "rs_cis", "cis", 0.30, 0.02, 50.0, 0.4],
            ["SAA1", "rs_trans", "trans", 0.95, 0.02, 60.0, 0.4],
        ]
    )
    r = U.resumir_por_proteina(d).iloc[0]
    assert r["variante"] == "rs_cis"
    assert r["n_cis_pqtl"] == 1
    assert r["n_trans_pqtl"] == 1


def test_una_proteina_solo_con_trans_no_aparece_en_el_resumen():
    d = _stats([["SAA1", "rs_trans", "trans", 0.95, 0.02, 60.0, 0.4]])
    assert U.resumir_por_proteina(d).empty


def test_el_umbral_de_varianza_decide_el_veredicto_y_no_el_p():
    """Justo por debajo del umbral no se declara determinada por genética."""
    # 2·0,4·0,6·beta² = 0,05  =>  beta ≈ 0,3227
    d = _stats(
        [
            ["A", "r1", "cis", 0.30, 0.01, 99.0, 0.4],
            ["B", "r2", "cis", 0.35, 0.01, 99.0, 0.4],
        ]
    )
    res = U.resumir_por_proteina(d).set_index("proteina")
    assert res.loc["A", "varianza_explicada"] < U.UMBRAL_VARIANZA
    assert bool(res.loc["A", "determinada_por_genetica"]) is False
    assert bool(res.loc["B", "determinada_por_genetica"]) is True


def test_un_formato_distinto_detiene_el_analisis_en_lugar_de_leer_mal():
    d = _stats([["SAA1", "r", "cis", 0.3, 0.01, 99.0, 0.4]]).rename(
        columns={"beta": "effect"}
    )
    with pytest.raises(KeyError, match="beta"):
        U.resumir_por_proteina(d)


def test_el_resumen_sale_ordenado_por_varianza_explicada():
    d = _stats(
        [
            ["A", "r1", "cis", 0.10, 0.01, 99.0, 0.4],
            ["B", "r2", "cis", 0.90, 0.01, 99.0, 0.4],
            ["C", "r3", "cis", 0.50, 0.01, 99.0, 0.4],
        ]
    )
    assert list(U.resumir_por_proteina(d)["proteina"]) == ["B", "C", "A"]


# --------------------------------------------------- lectura
def test_una_proteina_ausente_no_se_declara_sin_pqtl_a_secas():
    """Igual que en D4b: no medido y no encontrado no son lo mismo."""
    d = _stats([["LBP", "r", "cis", 0.90, 0.01, 99.0, 0.4]])
    res = U.lectura(U.resumir_por_proteina(d))
    assert res["por_proteina"]["SAA1"]["estado"] == "no_en_el_panel_o_sin_cis_pqtl"
    assert "2 923" in res["por_proteina"]["SAA1"]["aviso"]
    assert res["por_proteina"]["LBP"]["estado"] == "con_cis_pqtl"


def test_la_lectura_lista_las_proteinas_determinadas_por_genetica():
    d = _stats(
        [
            ["SAA1", "r1", "cis", 0.90, 0.01, 99.0, 0.4],
            ["LBP", "r2", "cis", 0.05, 0.01, 99.0, 0.4],
        ]
    )
    res = U.lectura(U.resumir_por_proteina(d))
    assert res["determinadas_por_genetica"] == ["SAA1"]
    assert res["n_determinadas_por_genetica"] == 1


def test_la_lectura_declara_lo_que_un_pqtl_no_dice():
    res = U.lectura(U.resumir_por_proteina(_stats([])))
    assert "no invalida un marcador" in res["lo_que_no_dice"]


def test_cada_proteina_objetivo_dice_por_que_se_consulta():
    res = U.lectura(U.resumir_por_proteina(_stats([])))
    assert set(res["por_proteina"]) == set(U.OBJETIVO)
    for v in res["por_proteina"].values():
        assert v["motivo_de_interes"]


def test_las_tres_proteinas_de_plasma_de_v1_estan_en_la_lista_objetivo():
    """La lista se fija antes de ver los datos. Si cambia, es una enmienda."""
    assert {"SAA1", "LBP", "CHMP4A"} <= set(U.OBJETIVO)
