#!/usr/bin/env python3
"""Comprueba el SHA-1 de los crudos de PXD055218 sin descargar ni mover nada.

El target `make crudos-v3-7` vuelve a bajar y manda a cuarentena el crudo que
no coincide. Este script solo lee y escribe la tabla en el workspace.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]


def _crudos():
    ruta = RAIZ / "scripts" / "download" / "pxd055218_crudos.py"
    spec = importlib.util.spec_from_file_location("pxd055218_crudos", ruta)
    if spec is None or spec.loader is None:
        raise ImportError(ruta)
    modulo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(modulo)
    return modulo


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crudos", type=Path, required=True)
    parser.add_argument("--salida", type=Path, required=True)
    args = parser.parse_args()
    crudos = _crudos()
    checksum = args.crudos / "checksum.txt"
    if not checksum.is_file():
        print(f"falta {checksum}", file=sys.stderr)
        return 1
    tabla = crudos.parsear_checksum(checksum.read_text(encoding="utf-8", errors="replace"))
    filas = []
    for nombre in sorted(nombre for nombre in tabla if nombre.endswith(".raw")):
        crudo = args.crudos / nombre
        if not crudo.is_file():
            filas.append({
                "crudo": nombre, "bytes_local": 0, "sha1_local": "",
                "sha1_pride": tabla[nombre], "igual": False,
            })
            print(f"falta {nombre}", flush=True)
            continue
        fila = crudos.evaluar_crudo(crudo, tabla[nombre], None)
        filas.append(fila)
        print(f"{'ok' if fila['igual'] else 'distinto'} {nombre}", flush=True)
    args.salida.parent.mkdir(parents=True, exist_ok=True)
    columnas = ["crudo", "bytes_local", "sha1_local", "sha1_pride", "igual"]
    lineas = ["\t".join(columnas)]
    for fila in filas:
        lineas.append("\t".join(
            str(fila[columna]).lower() if isinstance(fila[columna], bool) else str(fila[columna])
            for columna in columnas
        ))
    args.salida.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    distintos = sum(1 for fila in filas if not fila["igual"])
    print(f"{len(filas)} crudos; distintos {distintos}")
    return 1 if distintos else 0


if __name__ == "__main__":
    raise SystemExit(main())
