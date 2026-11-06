"""Informe que lee los JSON. No escribe una cifra que no esté en un fichero."""

from __future__ import annotations

import json
from pathlib import Path

from oncogo_crc.settings import Settings
from oncogo_crc.stages.v3_2_controles import _abrir
from oncogo_crc.v3.salida_v32 import volcar_json


def _leer(path: Path) -> dict | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def run(settings: Settings) -> None:
    _prereg, salida, sello = _abrir(settings, "v3_2_report")
    piezas = {
        "suero": _leer(salida / "suero_resumen.json"),
        "tejido": _leer(salida / "k1_etiquetas.json"),
        "k2": _leer(salida / "k2_independencia.json"),
        "k3": _leer(salida / "k3_deming.json"),
        "k5": _leer(salida / "k5_auc.json"),
        "controles": _leer(salida / "controles.json"),
        "multiverso": _leer(salida / "multiverso.json"),
        "h5": _leer(salida / "h5_especificidad_organo.json"),
        "h10": _leer(salida / "h10_especificidad.json"),
        "ml": _leer(salida / "ml_suero_permutacion.json"),
    }
    faltan = [nombre for nombre, valor in piezas.items() if valor is None]
    certificado = {
        "estado": "incompleto" if faltan else "calculado",
        "faltan": faltan,
        "frase_suero": None if piezas["suero"] is None else piezas["suero"].get("frase"),
        "h1": None if piezas["tejido"] is None else piezas["tejido"].get("h1_veredicto"),
        "advertencia": "Ningún resultado de esta tesis sirve para decidir sobre un paciente.",
    }
    volcar_json(certificado, salida / "certificado_v3_2.json", sello)
    report = salida / "report"
    report.mkdir(parents=True, exist_ok=True)
    lineas = ["<!DOCTYPE html><meta charset='utf-8'><title>OncoGO-CRC v3.2</title><h1>Corrida v3.2</h1>"]
    lineas.append(f"<p>{certificado['advertencia']}</p>")
    lineas.append(f"<p>Estado: {certificado['estado']}.</p>")
    if certificado["frase_suero"]:
        lineas.append(f"<p>Suero: {certificado['frase_suero']}</p>")
    if certificado["h1"]:
        lineas.append(f"<p>H1: {certificado['h1']}</p>")
    lineas.append("<p>Los números están en los TSV y JSON de results/v3_2/. Este informe no los copia a mano.</p>")
    if faltan:
        lineas.append("<p>Faltan: " + ", ".join(faltan) + ".</p>")
    (report / "informe.html").write_text("\n".join(lineas) + "\n", encoding="utf-8")
    print(f"[v3.2] informe en {report / 'informe.html'}; faltan {faltan}")
