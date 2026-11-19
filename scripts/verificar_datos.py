#!/usr/bin/env python3
"""Compara el workspace con data/checksums.sha256.

Una línea cuya ruta no existe se informa y cuenta como fallo. No descarga.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    resumen = hashlib.sha256()
    with path.open("rb") as fichero:
        for bloque in iter(lambda: fichero.read(1 << 20), b""):
            resumen.update(bloque)
    return resumen.hexdigest()


def lineas(path: Path) -> list[tuple[str, str]]:
    salida = []
    for linea in path.read_text(encoding="utf-8").splitlines():
        texto = linea.strip()
        if not texto or texto.startswith("#"):
            continue
        digest, rel = texto.split(maxsplit=1)
        salida.append((digest, rel))
    return salida


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--manifiesto", type=Path, default=RAIZ / "data" / "checksums.sha256")
    args = parser.parse_args()
    fallos = 0
    for digest, rel in lineas(args.manifiesto):
        candidatos = [args.workspace / rel, RAIZ / rel]
        path = next((candidato for candidato in candidatos if candidato.is_file()), None)
        if path is None:
            print(f"AUSENTE  {rel}")
            fallos += 1
            continue
        hallado = sha256(path)
        if hallado != digest:
            print(f"DISTINTO  {rel}")
            fallos += 1
        else:
            print(f"ok  {rel}")
    return 1 if fallos else 0


if __name__ == "__main__":
    raise SystemExit(main())
