#!/usr/bin/env python3
"""Compara una corrida nueva con results/v3_7.

Las 35 métricas de config/esperado_v3_7.tsv se leen de la corrida nueva.
Los ficheros se comparan por SHA-256. Se admiten como diferencia esperada los
sellos, la fecha, el identificador de proceso y el registro de pruebas, y el
orden de filas de gsea_tejido.tsv. Cualquier otra diferencia sale con código 1.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
ESPERADO = RAIZ / "config" / "esperado_v3_7.tsv"
NOMBRES_ESPERADOS = {
    "corrida.json",
    "tests.txt",
    "procedencia.json",
    "sello_corrida.json",
    "lecturas_v3_7.json",
}


def sha256(path: Path) -> str:
    resumen = hashlib.sha256()
    with path.open("rb") as fichero:
        for bloque in iter(lambda: fichero.read(1 << 20), b""):
            resumen.update(bloque)
    return resumen.hexdigest()


def _bajar(obj, camino: str):
    actual = obj
    for parte in camino.split("."):
        if not isinstance(actual, dict) or parte not in actual:
            raise KeyError(camino)
        actual = actual[parte]
    return actual


def valor_metrica(corrida: Path, fichero: str, campo: str):
    path = corrida / fichero
    if not path.is_file():
        raise FileNotFoundError(path)
    if fichero.endswith(".json"):
        return _bajar(json.loads(path.read_text(encoding="utf-8")), campo)
    if "|" not in campo:
        raise ValueError(campo)
    filtro, columna = campo.split("|", 1)
    clave, valor = filtro.split("=", 1)
    tabla = pd.read_csv(path, sep="\t")
    filas = tabla.loc[tabla[clave].astype(str).eq(valor)]
    if len(filas) != 1:
        raise ValueError(f"{fichero} {campo}: {len(filas)} filas")
    return filas.iloc[0][columna]


def metricas(corrida: Path) -> list[str]:
    fallos = []
    tabla = pd.read_csv(ESPERADO, sep="\t")
    for fila in tabla.itertuples(index=False):
        try:
            observado = valor_metrica(corrida, fila.fichero, fila.campo)
        except (FileNotFoundError, KeyError, ValueError) as exc:
            fallos.append(f"{fila.metrica}: {exc}")
            continue
        tolerancia = float(fila.tolerancia)
        if tolerancia == 0:
            if str(observado) != str(fila.valor):
                fallos.append(f"{fila.metrica}: {observado} != {fila.valor}")
        else:
            if not math.isclose(float(observado), float(fila.valor), abs_tol=tolerancia, rel_tol=0.0):
                fallos.append(f"{fila.metrica}: {observado} != {fila.valor} (tol {tolerancia})")
    return fallos


def _filas_gsea(path: Path) -> list[tuple]:
    tabla = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return sorted(tuple(fila) for fila in tabla.itertuples(index=False, name=None))


def ficheros(nueva: Path, referencia: Path) -> list[str]:
    nombres = sorted(
        {path.name for path in nueva.iterdir() if path.is_file()}
        | {path.name for path in referencia.iterdir() if path.is_file()}
    )
    distintos = []
    for nombre in nombres:
        izquierda, derecha = nueva / nombre, referencia / nombre
        if not izquierda.is_file() or not derecha.is_file():
            distintos.append(nombre)
            continue
        if nombre.endswith(".sello.json") or nombre in NOMBRES_ESPERADOS:
            continue
        if nombre == "gsea_tejido.tsv":
            if _filas_gsea(izquierda) != _filas_gsea(derecha):
                distintos.append(nombre)
            continue
        if sha256(izquierda) != sha256(derecha):
            distintos.append(nombre)
    return distintos


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nueva", type=Path, required=True)
    parser.add_argument("--referencia", type=Path, default=RAIZ / "results" / "v3_7")
    args = parser.parse_args()
    if not args.nueva.is_dir():
        print(f"no está la corrida nueva: {args.nueva}", file=sys.stderr)
        return 1
    fallos = metricas(args.nueva)
    distintos = ficheros(args.nueva, args.referencia)
    for fallo in fallos:
        print(f"metrica {fallo}")
    for nombre in distintos:
        print(f"fichero {nombre}")
    if fallos or distintos:
        print(f"{len(fallos)} métricas y {len(distintos)} ficheros fuera de la referencia")
        return 1
    print("la corrida nueva coincide con la referencia en las 35 métricas y en los ficheros comparables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
