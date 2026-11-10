"""Escritura bajo results/v3_3/. v3 y v3.2 no se tocan."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.sello import escribir_sello
from oncogo_crc.v3.rutas import rechazar_si_es_pipeline_antiguo

VERSION_VIGENTE = "v3_3"
VERSIONES_DE_ESCRITURA = ("v3_3", "v3_4")


def carpeta(settings) -> Path:
    return settings.resolve("results") / VERSION_VIGENTE


def carpeta_v34(settings) -> Path:
    return settings.resolve("results") / "v3_4"


def exigir_vigente(path: Path) -> None:
    texto = str(path).replace("\\", "/")
    rechazar_si_es_pipeline_antiguo(path)
    if any(f"/results/{version}/" in texto or texto.endswith(f"/results/{version}") for version in VERSIONES_DE_ESCRITURA):
        return
    raise RuntimeError(f"v3 solo escribe bajo results/v3_3 o results/v3_4: {path}")


def entrada_v3(settings, salida: Path, nombre: str) -> Path:
    """La corrida v3.4 lee un resultado congelado de v3.3 si todavía no lo reescribió."""
    local = salida / nombre
    if local.is_file():
        return local
    return settings.resolve("results/v3_3") / nombre


def manifiesto_congelado(settings) -> Path | None:
    """Une los manifiestos ya publicados. El fichero nuevo está ignorado y no ensucia el sello."""
    partes: list[str] = []
    for nombre in ("results/v3_manifest.sha256", "results/v3_2_manifest.sha256"):
        candidato = settings.root / nombre
        if not candidato.is_file():
            candidato = settings.resolve(nombre)
        if candidato.is_file():
            partes.append(candidato.read_text(encoding="utf-8").rstrip("\n"))
    if not partes:
        return None
    destino = settings.resolve("results/congelados.sha256")
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_text("\n".join(partes) + "\n", encoding="utf-8")
    return destino


def volcar_tsv(df: pd.DataFrame, path: Path, sello: dict) -> None:
    exigir_vigente(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(df, path)
    escribir_sello(path.with_name(path.name + ".sello.json"), sello)


def volcar_json(payload: dict, path: Path, sello: dict) -> None:
    exigir_vigente(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cuerpo = dict(payload)
    cuerpo["sello"] = sello
    path.write_text(json.dumps(cuerpo, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    escribir_sello(path.with_name(path.name + ".sello.json"), sello)
