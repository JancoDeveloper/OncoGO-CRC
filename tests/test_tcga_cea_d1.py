"""Tests de D1: determinantes del CEA sérico en TCGA.

Las matrices son sintéticas y viven solo aquí. Lo que estos tests fijan son las
decisiones de codificación que, si se equivocan, cambian la conclusión sin que
nadie lo note: qué cuenta como colon derecho, que MX no es M0, y que el modelo
recupera un efecto que se ha puesto a mano.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.v2 import tcga_cea as T


# ------------------------------------------------------------ codificación
def test_el_lado_sale_de_los_codigos_icd_y_no_de_colon_contra_recto():
    icd = pd.Series(["C18.0", "C18.2", "C18.7", "C20.9", "C18.9", "C49.4", None])
    lado = T.lado_del_tumor(icd)
    assert list(lado) == [
        "derecho", "derecho", "izquierdo", "izquierdo",
        "sin_especificar", "otro_o_ausente", "otro_o_ausente",
    ]


def test_el_angulo_esplenico_es_izquierdo_y_el_hepatico_derecho():
    """La frontera embriológica cae entre los dos ángulos, no en el transverso."""
    lado = T.lado_del_tumor(pd.Series(["C18.3", "C18.4", "C18.5", "C18.6"]))
    assert list(lado) == ["derecho", "derecho", "izquierdo", "izquierdo"]


def test_mx_no_se_cuenta_como_m0():
    """MX es que no se buscó metástasis. Meterlo en M0 inventaría información."""
    m = T.metastasis(pd.Series(["M0", "M1", "M1A", "M1B", "MX", None]))
    assert list(m) == ["M0", "M1", "M1", "M1", "desconocido", "desconocido"]


def test_los_subestadios_se_agrupan_sin_confundir_ii_con_iii():
    e = T.estadio_agrupado(
        pd.Series(["STAGE I", "STAGE IIA", "STAGE IIIB", "STAGE IVA", "STAGE IA", None])
    )
    assert list(e) == ["I", "II", "III", "IV", "I", "desconocido"]


def test_el_msi_usa_el_umbral_publicado_de_mantis():
    s = T.msi(pd.Series([0.2, 0.39, 0.4, 0.8, np.nan]))
    assert list(s) == ["estable", "estable", "inestable", "inestable", "desconocido"]


def test_preparar_no_logaritma_un_cea_de_cero():
    d = T.preparar(
        pd.DataFrame(
            {
                "cea_pretratamiento_ng_ml": [0.0, 1.0, 100.0],
                "rnaseq_CEACAM5": [0.0, 10.0, 1000.0],
                "icd_o_3_sitio": ["C18.0"] * 3,
                "estadio_ajcc": ["STAGE I"] * 3,
                "estadio_m": ["M0"] * 3,
                "msi_mantis": [0.3] * 3,
            }
        )
    )
    assert np.isnan(d.loc[0, "log10_cea"])
    assert d.loc[1, "log10_cea"] == pytest.approx(0.0)
    assert d.loc[2, "log10_cea"] == pytest.approx(2.0)
    # log2(x+1) para que una expresión de cero no se vaya a menos infinito.
    assert d.loc[0, "log2_CEACAM5"] == pytest.approx(0.0)


# ------------------------------------------------------------ censura
def test_sin_acumulacion_en_un_suelo_no_hace_falta_tobit():
    cea = pd.Series([0.0] + list(np.linspace(0.9, 50, 200)))
    res = T.censura_del_ensayo(cea)
    assert res["n_en_el_valor_minimo"] == 1
    assert res["hace_falta_tobit"] is False
    assert "mínimos cuadrados" in res["decision"]


def test_con_muchos_valores_pegados_al_minimo_si_hace_falta():
    cea = pd.Series([0.5] * 40 + list(np.linspace(1, 50, 100)))
    res = T.censura_del_ensayo(cea)
    assert res["hace_falta_tobit"] is True
    assert "Tobit" in res["decision"]


def test_sin_valores_de_cea_la_censura_no_es_evaluable():
    assert T.censura_del_ensayo(pd.Series([np.nan, np.nan]))["estado"] == "no_evaluable"


# ------------------------------------------------------------ modelo
def _cohorte_sintetica(n=300, efecto_gen=0.30, efecto_m1=1.0, semilla=20260914):
    rng = np.random.default_rng(semilla)
    gen = rng.normal(10, 2, n)
    m1 = rng.random(n) < 0.2
    pureza = rng.uniform(0.3, 0.9, n)
    log_cea = (
        0.5 + efecto_gen * (gen - 10) + efecto_m1 * m1 + rng.normal(0, 0.4, n)
    )
    return pd.DataFrame(
        {
            "cea_pretratamiento_ng_ml": 10**log_cea,
            "rnaseq_CEACAM5": 2**gen - 1,
            "pureza_absolute": pureza,
            "edad": rng.integers(40, 85, n),
            "sexo": rng.choice(["Male", "Female"], n),
            "icd_o_3_sitio": rng.choice(["C18.0", "C18.7"], n),
            "estadio_ajcc": np.where(m1, "STAGE IV", "STAGE IIA"),
            "estadio_m": np.where(m1, "M1", "M0"),
            "msi_mantis": rng.uniform(0.1, 0.6, n),
            "subtipo_pancanatlas": rng.choice(["COAD_CIN", "COAD_MSI"], n),
        }
    )


def test_el_modelo_recupera_un_efecto_puesto_a_mano():
    d = T.preparar(_cohorte_sintetica(efecto_gen=0.30))
    res = T.modelo_del_cea(d)
    assert res["estado"] == "calculado"
    gen = next(c for c in res["coeficientes"] if c["termino"] == "log2_CEACAM5")
    assert gen["coeficiente"] == pytest.approx(0.30, abs=0.05)
    assert gen["ci_low"] < 0.30 < gen["ci_high"]
    assert gen["p"] < 1e-6


def test_el_modelo_no_inventa_un_efecto_que_no_existe():
    d = T.preparar(_cohorte_sintetica(efecto_gen=0.0))
    res = T.modelo_del_cea(d)
    gen = next(c for c in res["coeficientes"] if c["termino"] == "log2_CEACAM5")
    assert gen["ci_low"] < 0 < gen["ci_high"]
    assert res["aporte_de_ceacam5_al_r2"] < 0.02


def test_se_reporta_cuanto_aporta_el_gen_por_encima_del_resto():
    """Un R² alto del modelo completo no significa que el gen explique nada."""
    d = T.preparar(_cohorte_sintetica(efecto_gen=0.0, efecto_m1=2.0))
    res = T.modelo_del_cea(d)
    assert res["r2"] > 0.5
    assert res["r2_solo_ceacam5"] < 0.05
    assert res["aporte_de_ceacam5_al_r2"] < 0.02


def test_con_pocos_casos_el_modelo_se_declara_no_calculable():
    d = T.preparar(_cohorte_sintetica(n=30))
    res = T.modelo_del_cea(d)
    assert res["estado"] == "no_calculable"
    assert "mínimo" in res["motivo"]


def test_la_categoria_de_referencia_es_la_mas_frecuente_y_se_declara():
    d = T.preparar(_cohorte_sintetica())
    res = T.modelo_del_cea(d)
    refs = res["categorias_de_referencia"]
    assert refs["metastasis"] == "M0"
    assert refs["estadio_grupo"] == "II"
    # La referencia no aparece como término: si apareciera, el diseño sería singular.
    terminos = {c["termino"] for c in res["coeficientes"]}
    assert "metastasis:M0" not in terminos


def test_una_categoria_ausente_no_rompe_el_diseno():
    d = T.preparar(_cohorte_sintetica())
    d["subtipo_pancanatlas"] = np.nan
    res = T.modelo_del_cea(d)
    assert res["estado"] == "calculado"


# ------------------------------------------------------------ no secretor
def test_el_fenotipo_no_secretor_se_define_con_el_umbral_clinico():
    rng = np.random.default_rng(7)
    n = 200
    gen = rng.normal(10, 2, n)
    # La mitad de los que expresan mucho no elevan el CEA.
    cea = np.where((gen > 11) & (rng.random(n) < 0.5), 2.0, 20.0)
    d = T.preparar(
        pd.DataFrame(
            {
                "cea_pretratamiento_ng_ml": cea,
                "rnaseq_CEACAM5": 2**gen - 1,
                "icd_o_3_sitio": ["C18.0"] * n,
                "estadio_ajcc": ["STAGE IIA"] * n,
                "estadio_m": ["M0"] * n,
                "msi_mantis": rng.uniform(0.1, 0.6, n),
                "pureza_absolute": rng.uniform(0.3, 0.9, n),
                "edad": rng.integers(40, 85, n),
            }
        )
    )
    res = T.fenotipo_no_secretor(d)
    assert res["estado"] == "calculado"
    assert res["umbral_de_cea_ng_ml"] == 5.0
    assert res["grupos"]["no_secretor"]["n"] > 0
    # Por construcción, los no secretores tienen el CEA por debajo del umbral.
    assert res["grupos"]["no_secretor"]["mediana_cea"] < 5.0
    assert res["grupos"]["secretor"]["mediana_cea"] >= 5.0


def test_la_multiplicidad_del_fenotipo_queda_declarada():
    d = T.preparar(_cohorte_sintetica())
    res = T.fenotipo_no_secretor(d)
    assert "exploratorio" in res["nota_de_multiplicidad"]
    assert str(len(res["contrastes_no_secretor_vs_secretor"])) in res["nota_de_multiplicidad"]


def test_con_pocos_casos_el_fenotipo_no_se_calcula():
    d = T.preparar(_cohorte_sintetica(n=20))
    assert T.fenotipo_no_secretor(d)["estado"] == "no_calculable"


def test_sin_ceros_la_sensibilidad_no_aplica():
    d = T.preparar(_cohorte_sintetica())
    assert T.sensibilidad_con_el_cero(d, suelo=0.1)["estado"] == "no_aplica"


def test_con_un_cero_la_sensibilidad_lo_imputa_y_lo_dice():
    d = _cohorte_sintetica()
    d.loc[0, "cea_pretratamiento_ng_ml"] = 0.0
    res = T.sensibilidad_con_el_cero(T.preparar(d), suelo=0.9)
    assert res["estado"] == "calculado"
    assert res["n_ceros_imputados"] == 1
    assert res["suelo_usado_ng_ml"] == 0.9
