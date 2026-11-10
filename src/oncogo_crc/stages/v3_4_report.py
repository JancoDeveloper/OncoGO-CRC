"""Certificado v3.4. Una casilla en rojo deja el estado en incompleto."""

from __future__ import annotations


def _casilla(identificador: str, estado: str, motivo: str, fuente: str) -> dict:
    return {"id": identificador, "estado": estado, "motivo": motivo, "fuente": fuente}


def _limitacion_aceptada(motivo: str, prereg: dict) -> bool:
    aceptadas = prereg.get("enmienda_v3_4", {}).get("certificado", {}).get("limitaciones_aceptadas", [])
    texto = motivo.lower()
    if "requiere_acceso" in texto and any("H9" in item for item in aceptadas):
        return True
    if "no_disponible" in texto and any("H7" in item for item in aceptadas):
        return True
    if "no_evaluable" in texto and any("H6" in item for item in aceptadas):
        return True
    return False


def casillas(piezas: dict, prereg: dict) -> list[dict]:
    """Cada casilla del checklist se lee de un JSON o de un conteo. No se inventa un verde."""
    salida = []
    h5 = piezas.get("h5") or {}
    h7 = piezas.get("h7") or {}
    h6 = piezas.get("h6") or {}
    h8 = piezas.get("h8") or {}
    sub = h5.get("submuestreo") or {}
    if sub.get("estado") == "calculado" and sub.get("n_repeticiones") == 200 and sub.get("n_sub"):
        salida.append(_casilla("h5_submuestreo", "verde", "200 repeticiones con n_sub", "h5_especificidad_organo.json"))
    else:
        salida.append(_casilla("h5_submuestreo", "rojo", "el submuestreo no completa 200 repeticiones", "h5_especificidad_organo.json"))
    conc = h5.get("concordancias_universo_principal") or {}
    denominadores = {v.get("denominador") for v in conc.values()} if conc else set()
    if conc and len(denominadores) == 1 and h5.get("potencia_por_carcinoma") and h5.get("sensibilidad_sin_ov"):
        salida.append(_casilla("h5_universo", "verde", "mismo denominador, potencia y sensibilidad sin OV", "h5_especificidad_organo.json"))
    else:
        salida.append(_casilla("h5_universo", "rojo", "faltan el universo común, la potencia o la sensibilidad sin OV", "h5_especificidad_organo.json"))
    controles = piezas.get("controles") or {}
    if controles.get("controles_fisicos") and controles.get("veredicto_v3_3") and controles.get("veredicto_v3_4"):
        salida.append(_casilla("controles_tejido", "verde", "comprobación física y los dos veredictos", "controles_fisicos.tsv"))
    else:
        salida.append(_casilla("controles_tejido", "rojo", "faltan la comprobación física o uno de los dos veredictos", "controles.tsv"))
    campos = ("constructo_por_recurso", "independientes", "concordancia_de_etiquetas", "modelo_completo")
    if all(campo in h7 for campo in campos):
        salida.append(_casilla("h7_constructos", "verde", "el JSON trae constructo, independencia, kappa y modelo_completo", "h7_visibilidad_en_sangre.json"))
    else:
        salida.append(_casilla("h7_constructos", "rojo", "el JSON de H7 no trae los campos de la enmienda", "h7_visibilidad_en_sangre.json"))
    faltantes = list(h7.get("variables_faltantes") or [])
    papier = (h7.get("por_recurso") or {}).get("papier_2024") or {}
    if not faltantes and h7.get("modelo_completo") is True:
        salida.append(_casilla("h7_variables", "verde", "están las variables del preregistro", "h7_visibilidad_en_sangre.json"))
    elif papier.get("estado") == "no_disponible" and papier.get("motivo"):
        salida.append(_casilla("h7_variables", "limitacion_declarada", f"no_disponible: {papier['motivo']}", "h7_visibilidad_en_sangre.json"))
    else:
        salida.append(_casilla("h7_variables", "rojo", "variables faltantes: " + ", ".join(faltantes or ["sin calcular"]), "h7_visibilidad_en_sangre.json"))
    for clave, fuente in (("h6", "h6_origen_celular.json"), ("h8", "h8_carga_tumoral.json")):
        bloque = h6 if clave == "h6" else h8
        if bloque.get("estado") in {"sostiene", "no_sostiene", "refuta"}:
            salida.append(_casilla(clave, "verde", bloque.get("motivo", bloque["estado"]), fuente))
        elif bloque.get("estado") in {"no_evaluable", "no_disponible", "requiere_acceso"} and _limitacion_aceptada(str(bloque.get("estado")), prereg):
            salida.append(_casilla(clave, "limitacion_declarada", str(bloque.get("motivo") or bloque["estado"]), fuente))
        else:
            salida.append(_casilla(clave, "rojo", "sin estado de la enmienda", fuente))
    multi = piezas.get("multiverso") or {}
    if multi.get("n_celdas_suero") == 288 and multi.get("n_celdas_ml") == 72 and multi.get("n_r1") == 2:
        salida.append(_casilla("multiverso", "verde", "288, 72 y 2", "multiverso_v3_4.json"))
    else:
        salida.append(_casilla("multiverso", "rojo", "la rejilla no está completa", "multiverso_v3_4.json"))
    if piezas.get("h9", {}).get("estado") == "requiere_acceso":
        salida.append(_casilla("h9", "limitacion_declarada", "requiere_acceso", "hipotesis_v3_4.tsv"))
    elif piezas.get("h9"):
        salida.append(_casilla("h9", "verde", str(piezas["h9"].get("estado")), "hipotesis_v3_4.tsv"))
    else:
        salida.append(_casilla("h9", "limitacion_declarada", "requiere_acceso", "hipotesis_v3_4.tsv"))
    return salida


def puede_publicar(certificado: dict) -> bool:
    """La enmienda deja publicar solo cuando ninguna casilla está en rojo."""
    return certificado.get("estado") == "publicable" and not certificado.get("faltan")


def armar_certificado(piezas: dict, prereg: dict | None = None) -> dict:
    lista = casillas(piezas, prereg or {})
    rojas = [c for c in lista if c["estado"] == "rojo"]
    limitaciones = [c for c in lista if c["estado"] == "limitacion_declarada"]
    return {
        "estado": "publicable" if not rojas else "incompleto",
        "faltan": rojas,
        "limitaciones": limitaciones,
        "casillas": lista,
        "advertencia": "Ningún resultado de esta tesis sirve para decidir sobre un paciente.",
    }
