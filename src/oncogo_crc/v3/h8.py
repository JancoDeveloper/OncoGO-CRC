"""H8. PXD064655 queda sin Fisher si la matriz sigue siendo de precursores."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.v3.salida import volcar_json


def mapa_de_muestras(columnas: list[str], metadatos: pd.DataFrame) -> dict:
    """Une el nombre SH_ del encabezado con el estadio. Las columnas ST_ no son casos."""
    pedidos = set(metadatos["source name"].astype(str))
    hallados = {}
    for columna in columnas:
        nombre = Path(str(columna)).name
        codigo = re.match(r"(SH_\d+)", nombre)
        if codigo and codigo.group(1) in pedidos:
            hallados[codigo.group(1)] = columna
    return {
        "n_metadatos": len(pedidos),
        "n_unidos": len(hallados),
        "uno_a_uno": len(hallados) == len(pedidos) and len(set(hallados.values())) == len(hallados),
        "columnas": hallados,
    }


def ejecutar(settings: Settings, salida: Path, sello: dict) -> dict:
    matriz = settings.resolve("external/pxd064655/report.pr_matrix.tsv")
    meta = settings.resolve("external/pxd064655/Sample_Metadata.txt")
    if not matriz.is_file() or not meta.is_file():
        payload = {"estado": "no_disponible", "motivo": "no está la matriz o el metadato de PXD064655"}
        volcar_json(payload, salida / "h8_carga_tumoral.json", sello)
        return payload
    cabecera = pd.read_csv(matriz, sep="\t", nrows=0).columns.tolist()
    metadatos = pd.read_csv(meta, sep="\t")
    union = mapa_de_muestras(cabecera, metadatos)
    precursor = "Precursor.Id" in cabecera
    try:
        import directlfq  # noqa: F401
        hay_directlfq = True
    except ImportError:
        hay_directlfq = False
    if precursor and not hay_directlfq:
        payload = {
            "estado": "no_evaluable",
            "motivo": "report.pr_matrix.tsv es una matriz de precursores. La enmienda pide directLFQ y ese paquete no está instalado. No se suma ni se promedia como sustituto.",
            "mapa_uno_a_uno": union["uno_a_uno"],
            "n_muestras_unidas": union["n_unidos"],
        }
        volcar_json(payload, salida / "h8_carga_tumoral.json", sello)
        return payload
    payload = {
        "estado": "no_evaluable",
        "motivo": "la agregación directLFQ no se ejecutó en esta etapa",
        "mapa_uno_a_uno": union["uno_a_uno"],
        "n_muestras_unidas": union["n_unidos"],
    }
    volcar_json(payload, salida / "h8_carga_tumoral.json", sello)
    return payload
