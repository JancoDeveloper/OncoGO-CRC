"""Certificado v3.5. Una limitación se acepta solo si la pieza y el estado están en la lista."""

from __future__ import annotations

from pathlib import Path


def _casilla(identificador: str, estado: str, motivo: str, fuente: str, pieza: str | None = None, estado_pieza: str | None = None) -> dict:
    return {
        "id": identificador,
        "estado": estado,
        "motivo": motivo,
        "fuente": fuente,
        "pieza": pieza or identificador,
        "estado_pieza": estado_pieza,
    }


def _limitacion_aceptada(pieza_id: str, estado: str, motivo: str, aceptadas: list) -> bool:
    del motivo
    for item in aceptadas:
        if not isinstance(item, dict):
            continue
        if item.get("pieza") == pieza_id and estado in item.get("estados", []):
            return True
    return False


def _aceptadas(prereg: dict) -> list:
    bloque = prereg.get("enmienda_v3_5", {})
    if not bloque and hasattr(prereg, "get"):
        crudo = prereg.get("enmienda_v3_5.certificado.limitaciones_aceptadas")
        if isinstance(crudo, list):
            return crudo
    return bloque.get("certificado", {}).get("limitaciones_aceptadas", [])


def _por_estado(clave: str, bloque: dict, fuente: str, aceptadas: list, estados_verdes: set[str]) -> dict:
    estado = str(bloque.get("estado") or "")
    motivo = str(bloque.get("motivo") or estado or "sin estado")
    if clave == "h6" and bloque.get("tabla_degenerada"):
        motivo = f"{motivo}. {bloque.get('motivo_degenerada', 'tabla degenerada')}"
    if estado in estados_verdes:
        return _casilla(clave, "verde", motivo, fuente, clave, estado)
    if _limitacion_aceptada(clave, estado, motivo, aceptadas):
        return _casilla(clave, "limitacion_declarada", motivo, fuente, clave, estado)
    return _casilla(clave, "rojo", motivo or "sin estado de la enmienda", fuente, clave, estado)


def casillas(piezas: dict, prereg: dict) -> list[dict]:
    aceptadas = _aceptadas(prereg)
    salida = []
    h5 = piezas.get("h5") or {}
    h7 = piezas.get("h7") or {}
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
    if controles.get("controles_fisicos") and controles.get("veredicto_v3_3") == "falla" and controles.get("veredicto_v3_4"):
        salida.append(_casilla("controles_tejido", "verde", "comprobación física; veredicto v3.3 = falla", "controles_v3_4.json"))
    else:
        salida.append(_casilla("controles_tejido", "rojo", "faltan la comprobación física o el veredicto v3.3 no es falla", "controles_v3_4.json"))
    campos = ("constructo_por_recurso", "independientes", "concordancia_de_etiquetas", "modelo_completo")
    if all(campo in h7 for campo in campos):
        salida.append(_casilla("h7_constructos", "verde", "el JSON trae constructo, independencia, kappa y modelo_completo", "h7_visibilidad_en_sangre.json"))
    else:
        salida.append(_casilla("h7_constructos", "rojo", "el JSON de H7 no trae los campos de la enmienda", "h7_visibilidad_en_sangre.json"))
    faltantes = list(h7.get("variables_faltantes") or [])
    if not faltantes and h7.get("modelo_completo") is True:
        salida.append(_casilla("h7_variables", "verde", "están las variables del preregistro", "h7_visibilidad_en_sangre.json"))
    else:
        salida.append(_casilla("h7_variables", "rojo", "variables faltantes: " + ", ".join(faltantes or ["sin calcular"]), "h7_visibilidad_en_sangre.json"))
    recursos = h7.get("por_recurso") or {}
    no_entrenable = next((r for r, bloque in recursos.items() if bloque.get("estado") in {"no_evaluable", "no_disponible"}), None)
    if no_entrenable:
        estado_recurso = recursos[no_entrenable]["estado"]
        motivo = f"{no_entrenable} {estado_recurso}"
        if _limitacion_aceptada("h7_recurso", estado_recurso, motivo, aceptadas):
            salida.append(_casilla("h7_recurso", "limitacion_declarada", motivo, "h7_visibilidad_en_sangre.json", "h7_recurso", estado_recurso))
        else:
            salida.append(_casilla("h7_recurso", "rojo", motivo, "h7_visibilidad_en_sangre.json", "h7_recurso", estado_recurso))
    else:
        salida.append(_casilla("h7_recurso", "verde", "cada recurso usado tiene al menos 10 por clase", "h7_visibilidad_en_sangre.json"))
    salida.append(_por_estado("h6", piezas.get("h6") or {}, "h6_origen_celular.json", aceptadas, {"sostiene", "no_sostiene", "refuta"}))
    salida.append(_por_estado("h8", piezas.get("h8") or {}, "h8_carga_tumoral.json", aceptadas, {"sostiene", "no_sostiene"}))
    pre = piezas.get("preanalitica") or {}
    if pre.get("estado") == "calculado" and pre.get("n_con_fecha") == 60:
        salida.append(_casilla("preanalitica", "verde", "calculado con 60", "preanalitica_v3_5.json", "preanalitica", "calculado"))
    elif pre.get("estado") == "no_evaluable" and pre.get("n_con_fecha") and _limitacion_aceptada("preanalitica", "no_evaluable", "", aceptadas):
        salida.append(_casilla("preanalitica", "limitacion_declarada", str(pre.get("motivo") or "calculado con los legibles"), "preanalitica_v3_5.json", "preanalitica", "no_evaluable"))
    else:
        salida.append(_casilla("preanalitica", "rojo", str(pre.get("motivo") or "sin la prueba de orden"), "preanalitica_v3_5.json", "preanalitica", pre.get("estado")))
    tri = piezas.get("triangulacion") or {}
    if tri.get("n_filas") == 14 and tri.get("celdas_vacias") == 0:
        salida.append(_casilla("triangulacion", "verde", "14 filas y ninguna celda vacía", "triangulacion_v3_5.tsv"))
    else:
        salida.append(_casilla("triangulacion", "rojo", "la triangulación no tiene 14 filas completas", "triangulacion_v3_5.tsv"))
    ml = piezas.get("ml_reproduccion") or {}
    auc = ml.get("auc")
    if isinstance(auc, (int, float)) and abs(float(auc) - 0.88025) <= 1e-9:
        salida.append(_casilla("ml_reproduccion", "verde", "la celda preregistrada reproduce 0.88025", "ml_celda_preregistrada.json"))
    else:
        salida.append(_casilla("ml_reproduccion", "rojo", "la celda preregistrada no reproduce 0.88025", "ml_celda_preregistrada.json"))
    multi = piezas.get("multiverso") or {}
    if multi.get("n_celdas_suero") == 288 and multi.get("n_celdas_ml") == 72 and multi.get("n_r1") == 2 and multi.get("resumen_presente"):
        salida.append(_casilla("multiverso", "verde", "288, 72 y 2, con resumen", "multiverso_resumen_v3_5.json"))
    else:
        salida.append(_casilla("multiverso", "rojo", "la rejilla o el resumen no está completo", "multiverso_resumen_v3_5.json"))
    h9 = piezas.get("h9") or {}
    if h9.get("estado") == "exploratorio":
        salida.append(_casilla("h9", "verde", "exploratorio", "h9_wald.json", "h9", "exploratorio"))
        salida.append(_casilla("h9_fuente", "verde", "la fuente es h9_wald.json", "h9_wald.json"))
    elif h9.get("estado") == "requiere_acceso" and _limitacion_aceptada("h9", "requiere_acceso", "", aceptadas):
        salida.append(_casilla("h9", "limitacion_declarada", "requiere_acceso", "h9_wald.json", "h9", "requiere_acceso"))
        salida.append(_casilla("h9_fuente", "verde", "la fuente citada existe como nombre", "h9_wald.json"))
    else:
        salida.append(_casilla("h9", "rojo", "H9 sin estado exploratorio ni limitación aceptada", "h9_wald.json", "h9", h9.get("estado")))
        salida.append(_casilla("h9_fuente", "rojo", "la fuente de H9 no es h9_wald.json", "h9_wald.json"))
    for clave, fuente in (("regresion_v1", "regresion_v1.json"), ("manifiestos", "manifiestos_v3_5.json")):
        bloque = piezas.get(clave) or {}
        if bloque.get("estado") == "pasa":
            salida.append(_casilla(clave, "verde", "pasa", fuente))
        else:
            salida.append(_casilla(clave, "rojo", str(bloque.get("motivo") or "no pasa"), fuente))
    return salida


def puede_publicar(certificado: dict, aceptadas: list | None = None, directorio: Path | None = None) -> bool:
    if certificado.get("estado") != "publicable" or certificado.get("faltan"):
        return False
    lista = aceptadas if aceptadas is not None else []
    for casilla in certificado.get("casillas", []):
        if casilla.get("estado") == "limitacion_declarada" and not _limitacion_aceptada(casilla.get("pieza", casilla["id"]), str(casilla.get("estado_pieza")), "", lista):
            return False
        if directorio is not None and not (directorio / casilla["fuente"]).is_file():
            return False
    return True


def armar_certificado(piezas: dict, prereg: dict | None = None, directorio: Path | None = None) -> dict:
    registro = prereg or {}
    lista = casillas(piezas, registro)
    if directorio is not None:
        for casilla in lista:
            if not (directorio / casilla["fuente"]).is_file() and casilla["estado"] != "rojo":
                casilla["estado"] = "rojo"
                casilla["motivo"] = f"no está {casilla['fuente']}"
    rojas = [c for c in lista if c["estado"] == "rojo"]
    limitaciones = [c for c in lista if c["estado"] == "limitacion_declarada"]
    aceptadas = _aceptadas(registro)
    certificado = {
        "estado": "publicable" if not rojas else "incompleto",
        "faltan": rojas,
        "limitaciones": limitaciones,
        "casillas": lista,
        "advertencia": "Ningún resultado de esta tesis sirve para decidir sobre un paciente.",
    }
    certificado["puede_publicar"] = puede_publicar(certificado, aceptadas, None if directorio is None else directorio)
    if directorio is not None and not certificado["puede_publicar"] and certificado["estado"] == "publicable":
        certificado["estado"] = "incompleto"
    return certificado
