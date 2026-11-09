"""Informe que lee los JSON. No escribe una cifra que no esté en un fichero."""

from __future__ import annotations

import json
from pathlib import Path

from oncogo_crc.settings import Settings
from oncogo_crc.stages.v3_3_controles import _abrir
from oncogo_crc.v3.changelog_v33 import filas as filas_changelog
from oncogo_crc.v3.salida import volcar_json, volcar_tsv


def armar_certificado(piezas: dict) -> dict:
    """Copia campos de los JSON. No calcula un veredicto nuevo."""
    tejido = piezas.get("tejido") or {}
    ml = piezas.get("ml") or {}
    h5 = piezas.get("h5") or {}
    h7 = piezas.get("h7") or {}
    h10 = piezas.get("h10") or {}
    suero = piezas.get("suero")
    faltan = [nombre for nombre, valor in piezas.items() if valor is None]
    return {
        "estado": "incompleto" if faltan else "calculado",
        "faltan": faltan,
        "frase_suero": None if suero is None else suero.get("frase"),
        "h1": tejido.get("h1_veredicto"),
        "advertencia": "Ningún resultado de esta tesis sirve para decidir sobre un paciente.",
        "titular_ml": ml.get("titular"),
        "h1_fragil_a_normalizacion": tejido.get("h1_fragil_a_normalizacion"),
        "h5": {
            "estado": h5.get("estado"),
            "fraccion_pan_carcinoma": h5.get("fraccion_pan_carcinoma"),
            "conteo": h5.get("conteo"),
            "carcinomas_disponibles": h5.get("carcinomas_disponibles"),
            "carcinomas_excluidos": h5.get("carcinomas_excluidos"),
        },
        "h7": {
            "estado": h7.get("estado"),
            "por_recurso": h7.get("por_recurso"),
            "variables_faltantes": h7.get("variables_faltantes"),
            "modelo_secundario": h7.get("modelo_secundario"),
        },
        "h10": {
            "estado": h10.get("estado"),
            "fase_pepe": h10.get("fase_pepe"),
            "n_especificas_crc": h10.get("n_especificas_crc"),
            "conteo": h10.get("conteo"),
            "valores_de_control": h10.get("valores_de_control"),
        },
    }


def _leer(path: Path) -> dict | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def run(settings: Settings) -> None:
    _prereg, salida, sello = _abrir(settings, "v3_3_report")
    piezas = {
        "suero": _leer(salida / "suero_resumen.json"),
        "tejido": _leer(salida / "k1_etiquetas.json"),
        "k2": _leer(salida / "k2_independencia.json"),
        "k3": _leer(salida / "k3_deming.json"),
        "k5": _leer(salida / "k5_auc.json"),
        "controles": _leer(salida / "controles.json"),
        "multiverso": _leer(salida / "multiverso.json"),
        "h5": _leer(salida / "h5_especificidad_organo.json"),
        "h7": _leer(salida / "h7_visibilidad_en_sangre.json"),
        "h10": _leer(salida / "h10_especificidad.json"),
        "ml": _leer(salida / "ml_suero_permutacion.json"),
    }
    certificado = armar_certificado(piezas)
    volcar_json(certificado, salida / "certificado_v3_3.json", sello)
    viejo = settings.root / "results" / "v3_2"
    volcar_tsv(filas_changelog(viejo, salida), salida / "changelog_numeros.tsv", sello)
    report = salida / "report"
    report.mkdir(parents=True, exist_ok=True)
    lineas = ["<!DOCTYPE html><meta charset='utf-8'><title>OncoGO-CRC v3.3</title><h1>Corrida v3.3</h1>"]
    lineas.append(f"<p>{certificado['advertencia']}</p>")
    lineas.append(f"<p>Estado: {certificado['estado']}.</p>")
    if certificado["frase_suero"]:
        lineas.append(f"<p>Suero: {certificado['frase_suero']}</p>")
    if certificado["h1"]:
        lineas.append(f"<p>H1: {certificado['h1']}</p>")
    titular = certificado.get("titular_ml") or {}
    if titular.get("texto"):
        lineas.append(f"<p>Modelo: {titular['texto']}</p>")
    lineas.append("<p>Los números están en los TSV y JSON de results/v3_3/. Este informe no los copia a mano.</p>")
    if certificado["faltan"]:
        lineas.append("<p>Faltan: " + ", ".join(certificado["faltan"]) + ".</p>")
    (report / "informe.html").write_text("\n".join(lineas) + "\n", encoding="utf-8")
    print(f"[v3.3] informe en {report / 'informe.html'}; faltan {certificado['faltan']}")
