"""Certificado de la corrida limpia. Las cifras esperadas son las de control, no una lectura de v3.5."""

from __future__ import annotations

import json
from pathlib import Path

from oncogo_crc.stages.v3_5_report import _limitacion_aceptada
from oncogo_crc.v3.rutas_publicas import texto_con_ruta_local
from oncogo_crc.v3.salida import volcar_json

_LIMITACIONES = [
    {"pieza": "h8", "estados": ["no_evaluable"]},
    {"pieza": "preanalitica", "estados": ["no_evaluable"]},
    {"pieza": "h7_recurso", "estados": ["no_evaluable", "no_disponible"]},
    {"pieza": "h6", "estados": ["no_evaluable"]},
    {"pieza": "h9", "estados": ["requiere_acceso"]},
]


def _cerca(valor, esperado, tolerancia) -> bool:
    try:
        return abs(float(valor) - float(esperado)) <= tolerancia
    except (TypeError, ValueError):
        return valor == esperado


def casillas(piezas: dict) -> dict[str, dict]:
    suero = piezas.get("suero_celda") or {}
    ml = piezas.get("ml") or {}
    multi = piezas.get("multiverso") or {}
    controles = piezas.get("controles") or {}
    h5 = piezas.get("h5") or {}
    h6 = piezas.get("h6") or {}
    h7 = piezas.get("h7") or {}
    h8 = piezas.get("h8") or {}
    h9 = piezas.get("h9") or {}
    h10 = piezas.get("h10") or {}
    pre = piezas.get("pre") or {}
    v1 = piezas.get("v1") or {}
    salida: dict[str, dict] = {}
    salida["regresion_v1"] = _por_bool(all(v1.get(k) == n for k, n in {
        "cptac_n_significant": 688, "pride_n_significant": 140, "n_both_significant_padj_lt_0.05": 28,
        "n_same_sign": 104, "n_sign_evaluable": 119, "peptide_n_measurable": 226,
        "n_certified_surviving": 307,
    }.items()), "results/v3_6/regresion_v1.json")
    salida["sin_mezcla"] = _por_bool(piezas.get("sin_mezcla") is True, "results/v3_6/lecturas_v3_6.json")
    salida["reproduce_v3_5"] = _por_bool(piezas.get("reproduce") is True, "results/v3_6/changelog_numeros.tsv")
    salida["rutas_publicas"] = _por_bool(piezas.get("rutas_publicas") is True, "results/v3_6/")
    salida["suero"] = _por_bool(
        suero.get("n") == 284 and suero.get("diferente") == 0 and suero.get("equivalente") == 166 and suero.get("indeterminada") == 118,
        "results/v3_6/suero_embudo.json",
    )
    m2 = (ml.get("versiones") or {}).get("M2_proteinas_principal") or {}
    salida["ml_reproduccion"] = _por_bool(_cerca(m2.get("auc_oof_media"), 0.88025, 1e-9), "results/v3_6/ml_celda_preregistrada.json")
    salida["multiverso"] = _por_bool(
        multi.get("n_celdas_suero") == 288 and multi.get("n_celdas_ml") == 72 and multi.get("n_r1") == 2,
        "results/v3_6/multiverso_resumen.json",
    )
    salida["controles_tejido"] = _por_bool(
        controles.get("veredicto_v3_3") == "falla" and bool(controles.get("controles_fisicos")) and bool(controles.get("veredicto_v3_4")),
        "results/v3_6/controles_tejido.json",
    )
    salida["h5"] = _por_bool(h5.get("estado") == "sostiene" and h5.get("universo_principal_n") == 29, "results/v3_6/h5_especificidad_organo.json")
    salida["h6"] = _por_bool(h6.get("estado") == "no_sostiene" and h6.get("tabla_degenerada") is True, "results/v3_6/h6_origen_celular.json")
    alvez = (h7.get("por_recurso") or {}).get("alvez_2023") or {}
    salida["h7"] = _por_bool(h7.get("estado") == "refuta" and _cerca(alvez.get("auc"), 0.48281209614887566, 1e-9), "results/v3_6/h7_visibilidad_en_sangre.json")
    salida["h8"] = _estado_h8(h8)
    auc = (pre.get("calculado_con_legibles") or {}).get("auc") or {}
    if pre.get("estado") == "calculado" and pre.get("n_con_fecha") == 60:
        salida["preanalitica"] = _por_bool(True, "results/v3_6/preanalitica.json")
    elif pre.get("estado") == "no_evaluable" and isinstance(pre.get("cota_con_60"), dict):
        salida["preanalitica"] = {"estado": "limitacion_aceptada", "fuente": "results/v3_6/preanalitica.json", "valor": auc.get("valor")}
    else:
        salida["preanalitica"] = _por_bool(False, "results/v3_6/preanalitica.json")
    salida["h9"] = _por_bool(h9.get("estado") in {"exploratorio", "no_evaluable"}, "results/v3_6/h9_wald.json")
    salida["h10"] = _por_bool(h10.get("estado") == "descriptiva", "results/v3_6/h10_especificidad.json")
    salida["triangulacion"] = _por_bool(piezas.get("triangulacion_ok") is True, "results/v3_6/triangulacion.tsv")
    salida["determinismo"] = _por_bool(piezas.get("determinismo") is True, "results/v3_6/determinismo.json")
    return salida


def _estado_h8(h8: dict) -> dict:
    if h8.get("estado") in {"sostiene", "no_sostiene"}:
        return _por_bool(True, "results/v3_6/h8_carga_tumoral.json")
    if h8.get("estado") == "no_evaluable":
        return {"estado": "limitacion_aceptada", "fuente": "results/v3_6/h8_carga_tumoral.json", "valor": h8.get("n_universo")}
    return _por_bool(False, "results/v3_6/h8_carga_tumoral.json")


def _por_bool(ok: bool, fuente: str) -> dict:
    return {"estado": "pasa" if ok else "rojo", "fuente": fuente}


def rutas_limpias(carpeta: Path) -> bool:
    for path in carpeta.rglob("*"):
        if path.suffix.lower() not in {".json", ".tsv", ".txt", ".md"} or not path.is_file():
            continue
        if path.name.startswith("._"):
            continue
        if texto_con_ruta_local(path.read_text(encoding="utf-8", errors="replace")):
            return False
    return True


def emitir(piezas: dict, salida: Path, sello: dict) -> dict:
    cajas = casillas(piezas)
    limitaciones = []
    for nombre, caja in cajas.items():
        if caja["estado"] == "limitacion_aceptada":
            limitaciones.append({"pieza": nombre, "estado": "no_evaluable" if nombre in {"h8", "preanalitica"} else caja["estado"]})
    aceptadas = all(_limitacion_aceptada(item["pieza"], "no_evaluable", "", _LIMITACIONES) for item in limitaciones)
    rojas = [nombre for nombre, caja in cajas.items() if caja["estado"] == "rojo"]
    payload = {
        "estado": "publicable" if not rojas and aceptadas else "incompleto",
        "puede_publicar": not rojas and aceptadas,
        "casillas": cajas,
        "limitaciones": limitaciones,
        "rojas": rojas,
        "lectura": "Los grupos se midieron en bloques: primero el cáncer, después la colonoscopia negativa y al final los sanos. La piel y la coagulación cambian con el orden. Con estos datos no se puede separar el efecto del grupo del efecto del lote de adquisición.",
    }
    volcar_json(payload, salida / "certificado.json", sello)
    return payload
