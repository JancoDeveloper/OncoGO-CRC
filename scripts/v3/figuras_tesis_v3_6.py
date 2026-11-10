#!/usr/bin/env python3
"""Figuras de la tesis a partir de results/v3_6/. PNG a 300 ppp y PDF."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
V36 = REPO / "results" / "v3_6"
SALIDA = REPO / "tesis-ucaldas" / "paquete-redaccion-v3_6" / "figuras"
CP, CN, HC = "#D55E00", "#0072B2", "#009E73"
DIF, EQ, IND = "#D55E00", "#56B4E9", "#BFC7D1"


def _sello(path: Path, fuente: str) -> None:
    payload = {
        "fuente": fuente,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    path.with_suffix(path.suffix + ".sello.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _guardar(fig, nombre: str, fuente: str) -> None:
    SALIDA.mkdir(parents=True, exist_ok=True)
    png = SALIDA / f"{nombre}.png"
    pdf = SALIDA / f"{nombre}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    _sello(png, fuente)
    _sello(pdf, fuente)


def embudo() -> None:
    datos = json.loads((V36 / "suero_embudo.json").read_text(encoding="utf-8"))
    celda = datos["celda_preregistrada"]
    controles = datos["controles"]
    etiquetas = ["Grupos", "q global", "Humanas", "Celda"]
    valores = [controles["grupos"], controles["grupos_global_pg_q_0_01"], controles["humanos_pg_y_global"], celda["n"]]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.barh(etiquetas[::-1], valores[::-1], color="#0E7C86")
    ax.set_xlabel("Proteínas")
    ax.set_title("Embudo del suero")
    _guardar(fig, "Fig4-1_embudo_suero", "results/v3_6/suero_embudo.json")


def orden() -> None:
    tabla = pd.read_csv(V36 / "preanalitica_orden.tsv", sep="\t")
    colores = {"CP": CP, "CN": CN, "HC": HC}
    fig, ax = plt.subplots(figsize=(8, 3.2))
    for grupo, bloque in tabla.groupby("grupo"):
        ax.scatter(pd.to_numeric(bloque["orden"]), [grupo] * len(bloque), c=colores.get(grupo, "#475569"), s=28)
    ax.set_xlabel("Orden de adquisición")
    ax.set_title("Las muestras se midieron en tres bloques")
    _guardar(fig, "Fig4-2_orden_adquisicion", "results/v3_6/preanalitica_orden.tsv")


def equivalencia() -> None:
    tabla = pd.read_csv(V36 / "suero_k4_crc_cn.tsv", sep="\t")
    colores = {"diferente": DIF, "equivalente": EQ, "indeterminada": IND}
    fig, ax = plt.subplots(figsize=(7, 4))
    cuentas = tabla["clase"].value_counts()
    nombres = ["diferente", "equivalente", "indeterminada"]
    ax.bar(nombres, [int(cuentas.get(n, 0)) for n in nombres], color=[colores[n] for n in nombres])
    ax.set_ylabel("Proteínas")
    ax.set_title("Cáncer frente a colonoscopia negativa")
    _guardar(fig, "Fig4-3_equivalencia_suero", "results/v3_6/suero_k4_crc_cn.tsv")


def main() -> None:
    embudo()
    orden()
    equivalencia()
    print("figuras en", SALIDA)


if __name__ == "__main__":
    main()
