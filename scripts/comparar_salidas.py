#!/usr/bin/env python3
"""Compara una corrida nueva con la referencia sellada, fichero a fichero.

Clasifica cada salida como idéntica en bytes, distinta solo en números
(con la diferencia máxima), distinta en texto, o ausente. No reescribe
la referencia.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

TOLERANCIA = 1e-6
SELLOS = {
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


def _numero(texto: str) -> float | None:
    if texto is None:
        return None
    limpio = str(texto).strip()
    if limpio.lower() in {"", "nan", "none", "null"}:
        return None
    try:
        valor = float(limpio)
    except ValueError:
        return None
    if not math.isfinite(valor):
        return None
    return valor


def _filas(path: Path) -> list[list[str]] | None:
    if path.suffix.lower() not in {".tsv", ".csv"}:
        return None
    texto = path.read_text(encoding="utf-8", errors="replace")
    separador = "\t" if path.suffix.lower() == ".tsv" else ","
    return [linea.split(separador) for linea in texto.splitlines()]


def _max_numerico(izq: list[list[str]], der: list[list[str]]) -> tuple[str, float]:
    if not izq or not der or izq[0] != der[0]:
        return "distinto", math.inf
    if len(izq) != len(der):
        ordenada_izq = [tuple(fila) for fila in izq[1:]]
        ordenada_der = [tuple(fila) for fila in der[1:]]
        if sorted(ordenada_izq) == sorted(ordenada_der):
            return "mismo_contenido_otro_orden", 0.0
        return "distinto", math.inf
    maximo = 0.0
    for fila_izq, fila_der in zip(izq, der, strict=True):
        if len(fila_izq) != len(fila_der):
            return "distinto", math.inf
        for celda_izq, celda_der in zip(fila_izq, fila_der, strict=True):
            if celda_izq == celda_der:
                continue
            num_izq, num_der = _numero(celda_izq), _numero(celda_der)
            if num_izq is None or num_der is None:
                return "distinto", math.inf
            maximo = max(maximo, abs(num_izq - num_der))
    if maximo == 0.0:
        return "mismo_valor_otro_formato", 0.0
    if maximo <= TOLERANCIA:
        return "numerico_dentro_de_tolerancia", maximo
    return "numerico_fuera_de_tolerancia", maximo


def _caminar(izq, der, maximo: list[float]) -> bool:
    if isinstance(izq, dict) and isinstance(der, dict):
        if set(izq) != set(der):
            return False
        return all(_caminar(izq[clave], der[clave], maximo) for clave in izq)
    if isinstance(izq, list) and isinstance(der, list):
        if len(izq) != len(der):
            return False
        return all(_caminar(a, b, maximo) for a, b in zip(izq, der, strict=True))
    if isinstance(izq, bool) or isinstance(der, bool):
        return izq == der
    if isinstance(izq, (int, float)) and isinstance(der, (int, float)):
        if not (math.isfinite(float(izq)) and math.isfinite(float(der))):
            return izq == der
        maximo[0] = max(maximo[0], abs(float(izq) - float(der)))
        return True
    return izq == der


def _json(izq: Path, der: Path) -> tuple[str, float]:
    try:
        izquierda = json.loads(izq.read_text(encoding="utf-8"))
        derecha = json.loads(der.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return "distinto", math.inf
    maximo = [0.0]
    if not _caminar(izquierda, derecha, maximo):
        return "distinto", math.inf
    if maximo[0] == 0.0:
        return "mismo_valor_otro_formato", 0.0
    if maximo[0] <= TOLERANCIA:
        return "numerico_dentro_de_tolerancia", maximo[0]
    return "numerico_fuera_de_tolerancia", maximo[0]


def clasificar(nueva: Path, referencia: Path) -> list[dict]:
    nombres = sorted(
        {path.name for path in nueva.iterdir() if path.is_file() and not path.name.startswith("._")}
        | {path.name for path in referencia.iterdir() if path.is_file() and not path.name.startswith("._")}
    )
    filas = []
    for nombre in nombres:
        izquierda, derecha = nueva / nombre, referencia / nombre
        if not derecha.is_file():
            filas.append({"fichero": nombre, "clase": "extra", "diferencia_maxima": ""})
            continue
        if not izquierda.is_file():
            filas.append({"fichero": nombre, "clase": "falta", "diferencia_maxima": ""})
            continue
        if sha256(izquierda) == sha256(derecha):
            filas.append({"fichero": nombre, "clase": "identico", "diferencia_maxima": "0"})
            continue
        if nombre.endswith(".sello.json") or nombre in SELLOS:
            filas.append({"fichero": nombre, "clase": "sello", "diferencia_maxima": ""})
            continue
        tabla_izq, tabla_der = _filas(izquierda), _filas(derecha)
        if tabla_izq is not None and tabla_der is not None:
            clase, maximo = _max_numerico(tabla_izq, tabla_der)
        elif izquierda.suffix == ".json" and derecha.suffix == ".json":
            clase, maximo = _json(izquierda, derecha)
        else:
            clase, maximo = "distinto", math.inf
        texto = "" if math.isinf(maximo) else f"{maximo:.3e}"
        filas.append({"fichero": nombre, "clase": clase, "diferencia_maxima": texto})
    return filas


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nueva", type=Path, required=True)
    parser.add_argument("--referencia", type=Path, required=True)
    parser.add_argument("--salida", type=Path, required=True)
    args = parser.parse_args()
    if not args.nueva.is_dir() or not args.referencia.is_dir():
        print("falta la corrida nueva o la referencia", file=sys.stderr)
        return 1
    tabla = clasificar(args.nueva, args.referencia)
    args.salida.parent.mkdir(parents=True, exist_ok=True)
    with args.salida.open("w", encoding="utf-8") as fichero:
        fichero.write("fichero\tclase\tdiferencia_maxima\ttolerancia\n")
        for fila in tabla:
            fichero.write(f"{fila['fichero']}\t{fila['clase']}\t{fila['diferencia_maxima']}\t{TOLERANCIA:.0e}\n")
    conteo: dict[str, int] = {}
    for fila in tabla:
        conteo[fila["clase"]] = conteo.get(fila["clase"], 0) + 1
    for clase, n in sorted(conteo.items()):
        print(f"{clase}\t{n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
