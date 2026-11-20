#!/usr/bin/env python3
"""Compara los ficheros de los manifiestos con su SHA-256. No imprime URLs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    resumen = hashlib.sha256()
    with path.open("rb") as fichero:
        for bloque in iter(lambda: fichero.read(1 << 20), b""):
            resumen.update(bloque)
    return resumen.hexdigest()


def _fila(fuente: str, dataset: str, lado: str, ruta: str, esperado: str, datos: Path) -> dict:
    esperado = (esperado or "").strip().lower()
    if not ruta:
        return {
            "fuente": fuente, "dataset": dataset, "lado": lado, "ruta": "",
            "estado": "sin_ruta", "bytes": "", "sha256_archivo": "", "sha256_esperado": esperado,
        }
    path = Path(ruta) if Path(ruta).is_absolute() else datos / ruta
    if not path.is_file():
        return {
            "fuente": fuente, "dataset": dataset, "lado": lado, "ruta": ruta,
            "estado": "falta", "bytes": "", "sha256_archivo": "", "sha256_esperado": esperado,
        }
    digest = sha256(path)
    if not esperado:
        estado = "presente_sin_huella"
    elif digest == esperado:
        estado = "coincide"
    else:
        estado = "no_coincide"
    return {
        "fuente": fuente, "dataset": dataset, "lado": lado, "ruta": ruta,
        "estado": estado, "bytes": path.stat().st_size,
        "sha256_archivo": digest, "sha256_esperado": esperado,
    }


def filas(datos: Path) -> list[dict]:
    salida: list[dict] = []
    manifiesto = json.loads((datos / "data/raw/download_real_manifest.json").read_text(encoding="utf-8"))
    for item in manifiesto:
        dataset = str(item.get("dataset") or "")
        salida.append(_fila(
            "download_real_manifest.json", dataset, "matriz",
            str(item.get("matrix_path") or ""), str(item.get("matrix_sha256") or ""), datos,
        ))
        salida.append(_fila(
            "download_real_manifest.json", dataset, "metadata",
            str(item.get("metadata_path") or ""), str(item.get("metadata_sha256") or ""), datos,
        ))
    with (datos / "data/raw/manifest.tsv").open(encoding="utf-8") as fichero:
        for row in csv.DictReader(fichero, delimiter="\t"):
            salida.append(_fila(
                "manifest.tsv", str(row.get("dataset") or ""), "matriz",
                str(row.get("matrix_path") or ""), "", datos,
            ))
    return salida


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datos", type=Path, required=True)
    parser.add_argument("--salida", type=Path, required=True)
    args = parser.parse_args()
    columnas = ["fuente", "dataset", "lado", "ruta", "estado", "bytes", "sha256_archivo", "sha256_esperado"]
    tabla = filas(args.datos)
    args.salida.parent.mkdir(parents=True, exist_ok=True)
    with args.salida.open("w", encoding="utf-8") as fichero:
        fichero.write("\t".join(columnas) + "\n")
        for fila in tabla:
            fichero.write("\t".join(str(fila[columna]) for columna in columnas) + "\n")
    conteo: dict[str, int] = {}
    for fila in tabla:
        conteo[fila["estado"]] = conteo.get(fila["estado"], 0) + 1
    print(args.salida)
    for estado, n in sorted(conteo.items()):
        print(f"{estado}\t{n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
