"""Congelado. v3.2 se reproduce desde la etiqueta v3.2-final."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.sello import escribir_sello
from oncogo_crc.v3.rutas import rechazar_si_es_pipeline_antiguo


def exigir_v32(path: Path) -> None:
    texto = str(path).replace("\\", "/")
    rechazar_si_es_pipeline_antiguo(path)
    if "/results/v3_2/" not in texto and not texto.endswith("/results/v3_2"):
        raise RuntimeError(f"v3.2 solo escribe bajo results/v3_2: {path}")


def volcar_tsv(df: pd.DataFrame, path: Path, sello: dict) -> None:
    exigir_v32(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(df, path)
    escribir_sello(path.with_name(path.name + ".sello.json"), sello)


def volcar_json(payload: dict, path: Path, sello: dict) -> None:
    exigir_v32(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cuerpo = dict(payload)
    cuerpo["sello"] = sello
    path.write_text(json.dumps(cuerpo, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    escribir_sello(path.with_name(path.name + ".sello.json"), sello)
