"""La tabla de hipótesis y el certificado leen results/v3_3. El changelog copia los JSON."""

import json

from oncogo_crc.stages.v3_3_hipotesis import FUENTE
from oncogo_crc.stages.v3_3_report import armar_certificado
from oncogo_crc.v3.changelog_v33 import filas


def test_fuentes_de_hipotesis_apuntan_a_v3_3():
    assert FUENTE == {
        "H1": "results/v3_3/k1_etiquetas.json",
        "H5": "results/v3_3/h5_especificidad_organo.json",
        "H7": "results/v3_3/h7_visibilidad_en_sangre.json",
        "H10": "results/v3_3/h10_especificidad.json",
    }


def test_certificado_copia_el_titular_y_no_lo_elige():
    certificado = armar_certificado({
        "suero": {"frase": "No se detecta diferencia; 166 proteínas son equivalentes"},
        "tejido": {"h1_veredicto": "rechaza_h1", "h1_fragil_a_normalizacion": False},
        "k2": {},
        "k3": {},
        "k5": {},
        "controles": {},
        "multiverso": {},
        "h5": {"estado": "refuta", "fraccion_pan_carcinoma": 0.12},
        "h7": {"estado": "refuta"},
        "h10": {"estado": "descriptiva", "n_especificas_crc": 5, "fase_pepe": 1},
        "ml": {"titular": {"clave": "separacion_explicada_por_calidad", "texto": "sale de la regla"}},
    })
    assert certificado["estado"] == "calculado"
    assert certificado["titular_ml"]["clave"] == "separacion_explicada_por_calidad"
    assert certificado["h1_fragil_a_normalizacion"] is False
    assert certificado["h7"]["estado"] == "refuta"
    assert certificado["faltan"] == []


def test_changelog_solo_anota_lo_que_cambia(tmp_path):
    viejo = tmp_path / "v3_2"
    nuevo = tmp_path / "v3_3"
    viejo.mkdir()
    nuevo.mkdir()
    (viejo / "suero_resumen.json").write_text(json.dumps({
        "n_principal": 205,
        "frase": "No se detecta diferencia; 105 proteínas son equivalentes",
        "clases": {"principal": {"crc_frente_a_colonoscopia_negativa": {"n": 205, "diferente": 0, "equivalente": 105, "indeterminada": 100, "conflicto": 0}}},
    }), encoding="utf-8")
    (nuevo / "suero_resumen.json").write_text(json.dumps({
        "n_principal": 284,
        "frase": "No se detecta diferencia; 166 proteínas son equivalentes",
        "clases": {"principal": {"crc_frente_a_colonoscopia_negativa": {"n": 284, "diferente": 0, "equivalente": 166, "indeterminada": 118, "conflicto": 0}}},
    }), encoding="utf-8")
    (viejo / "k1_etiquetas.json").write_text(json.dumps({
        "h1_veredicto": "rechaza_h1",
        "contrastes": {"primario_frente_a_mucosa": {"numerador": 89, "denominador": 90}},
    }), encoding="utf-8")
    (nuevo / "k1_etiquetas.json").write_text(json.dumps({
        "h1_veredicto": "rechaza_h1",
        "h1_fragil_a_normalizacion": False,
        "contrastes": {"primario_frente_a_mucosa": {"numerador": 147, "denominador": 147}},
    }), encoding="utf-8")
    tabla = filas(viejo, nuevo)
    por_metrica = dict(zip(tabla["metrica"], tabla["causa"], strict=False))
    assert por_metrica["n_principal"] == "2.1"
    assert por_metrica["primario_frente_a_mucosa"] == "2.3"
    assert por_metrica["h1_fragil_a_normalizacion"] == "2.4"
    assert "h1_veredicto" not in por_metrica
