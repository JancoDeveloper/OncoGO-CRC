"""El certificado v3.5 acepta una limitación solo si la pieza y el estado coinciden."""

from __future__ import annotations

from oncogo_crc.stages.v3_5_report import armar_certificado


def _piezas() -> dict:
    return {
        "h5": {"submuestreo": {"estado": "calculado", "n_repeticiones": 200, "n_sub": 10},
               "concordancias_universo_principal": {"COAD": {"denominador": 29}},
               "potencia_por_carcinoma": {}, "sensibilidad_sin_ov": {}},
        "controles": {"controles_fisicos": True, "veredicto_v3_3": "falla", "veredicto_v3_4": "pasa"},
        "h6": {"estado": "no_sostiene", "motivo": "Fisher", "tabla_degenerada": True,
               "motivo_degenerada": "ninguna proteína tiene signo distinto; la prueba no puede detectar una diferencia"},
        "h7": {"constructo_por_recurso": {}, "independientes": {}, "concordancia_de_etiquetas": {},
               "modelo_completo": True, "variables_faltantes": [],
               "por_recurso": {"papier_2024": {"estado": "no_evaluable", "n_pos": 1, "n_neg": 125}}},
        "h8": {"estado": "no_evaluable", "motivo": "directLFQ no se pudo evaluar"},
        "h9": {"estado": "exploratorio"},
        "multiverso": {"n_celdas_suero": 288, "n_celdas_ml": 72, "n_r1": 2, "resumen_presente": True},
        "preanalitica": {"estado": "no_evaluable", "n_con_fecha": 54, "n_crudos_validos": 60},
        "triangulacion": {"n_filas": 14, "celdas_vacias": 0},
        "ml_reproduccion": {"auc": 0.88025},
        "regresion_v1": {"estado": "pasa"},
        "manifiestos": {"estado": "pasa"},
    }


def test_h8_no_evaluable_con_solo_la_regla_de_h6_queda_en_rojo():
    prereg = {"enmienda_v3_5": {"certificado": {"limitaciones_aceptadas": [
        {"pieza": "h6", "estados": ["no_evaluable"], "motivo": "sin leyenda verificable de Pelka ni del HPA"},
    ]}}}
    certificado = armar_certificado(_piezas(), prereg)
    h8 = next(c for c in certificado["casillas"] if c["id"] == "h8")
    assert h8["estado"] == "rojo"
    assert certificado["estado"] == "incompleto"
