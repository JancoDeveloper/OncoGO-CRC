"""El certificado nombra las casillas en rojo. No escribe faltan vacío."""

from oncogo_crc.stages.v3_4_report import armar_certificado

PREREG = {
    "enmienda_v3_4": {
        "certificado": {
            "limitaciones_aceptadas": [
                "H9 requiere_acceso",
                "ESM-2 o TMbed no ejecutados si el autor dijo que no",
                "un recurso de H7 no_disponible con motivo",
                "H6 no_evaluable si ni Pelka ni el HPA tienen leyenda verificable",
            ]
        }
    }
}


def test_h7_de_v3_3_deja_el_certificado_incompleto():
    piezas = {
        "h7": {
            "estado": "refuta",
            "variables_faltantes": [
                "especificidad_tisular",
                "origen_celular_h6",
                "presencia_en_vesiculas_pxd033236",
                "papier_2024",
            ],
            "modelo_completo": False,
        }
    }
    certificado = armar_certificado(piezas, PREREG)
    assert certificado["estado"] == "incompleto"
    assert certificado["faltan"]
    nombres = [c["id"] for c in certificado["faltan"]]
    assert "h7_variables" in nombres
    assert "h7_constructos" in nombres
    assert all(c["estado"] == "rojo" for c in certificado["faltan"])
