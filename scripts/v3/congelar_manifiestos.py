#!/usr/bin/env python3
"""Congela v2, el cierre de tesis y las entradas de v3.

No modifica los ficheros que resume. Escribe:

- results/v2_manifest.sha256
- results/tesis_cierre_manifest.sha256
- results/v3/manifest_inputs.sha256

Los dos primeros cubren lo ya calculado. El tercero repite las entradas que
``scripts/v2/manifest_hashes.py`` ya vigila, en la carpeta de v3.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "manifest_hashes_v2", REPO / "scripts/v2/manifest_hashes.py"
)
_modulo = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_modulo)
WORKSPACE = _modulo.WORKSPACE
collect_inputs = _modulo.collect_inputs
is_noise = _modulo.is_noise
render = _modulo.render


def _arbol(relativo: str) -> list[Path]:
    base = WORKSPACE / relativo
    if not base.is_dir():
        return []
    return sorted(p for p in base.rglob("*") if p.is_file() and not is_noise(p))


def _escribir(relativo: str, paths: list[Path]) -> None:
    destino = WORKSPACE / relativo
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_text(render(paths), encoding="utf-8")
    print(f"[v3-manifiestos] {len(paths)} ficheros -> {destino}")


def _congelar_carpeta(carpeta: str, salida: str) -> int:
    """Resume una carpeta ya publicada. No modifica sus ficheros."""
    relativo = Path(carpeta)
    if relativo.is_absolute():
        base = relativo
        raiz = REPO
    else:
        en_repo = REPO / relativo
        base = en_repo if en_repo.is_dir() else WORKSPACE / relativo
        raiz = REPO if en_repo.is_dir() else WORKSPACE
    if not base.is_dir():
        print(f"carpeta no accesible: {base}", file=sys.stderr)
        return 2
    paths = sorted(p for p in base.rglob("*") if p.is_file() and not is_noise(p))
    destino = Path(salida)
    if not destino.is_absolute():
        destino = REPO / destino
    destino.parent.mkdir(parents=True, exist_ok=True)
    lineas = [f"{_modulo.sha256(p)}  {p.relative_to(raiz)}" for p in paths]
    destino.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    print(f"[v3-manifiestos] {len(paths)} ficheros -> {destino}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--carpeta", help="Carpeta a resumir, por ejemplo results/v3_2")
    parser.add_argument("--salida", help="Manifiesto de salida, por ejemplo results/v3_2_manifest.sha256")
    args = parser.parse_args()
    if args.carpeta or args.salida:
        if not args.carpeta or not args.salida:
            print("hacen falta --carpeta y --salida juntas", file=sys.stderr)
            return 2
        return _congelar_carpeta(args.carpeta, args.salida)
    if not WORKSPACE.is_dir():
        print(f"workspace no accesible: {WORKSPACE}", file=sys.stderr)
        return 2
    _escribir("results/v2_manifest.sha256", _arbol("results/v2"))
    _escribir("results/tesis_cierre_manifest.sha256", _arbol("results/tesis_cierre"))
    _escribir("results/v3/manifest_inputs.sha256", collect_inputs())
    _escribir("results/v3_manifest.sha256", _arbol("results/v3"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
