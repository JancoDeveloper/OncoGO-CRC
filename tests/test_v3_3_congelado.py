"""results/v3_2 no cambia después de la etiqueta v3.2-final."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


def _sha256(path: Path) -> str:
    resumen = hashlib.sha256()
    with path.open("rb") as entrada:
        while True:
            bloque = entrada.read(1 << 20)
            if not bloque:
                break
            resumen.update(bloque)
    return resumen.hexdigest()


def test_el_manifiesto_de_v3_2_sigue_igual():
    carpeta = REPO / "results" / "v3_2"
    manifiesto = REPO / "results" / "v3_2_manifest.sha256"
    if not carpeta.is_dir() or not manifiesto.is_file():
        pytest.skip("results/v3_2 no está en este entorno")
    esperado = {}
    for linea in manifiesto.read_text(encoding="utf-8").splitlines():
        if not linea.strip():
            continue
        digest, _, rel = linea.partition("  ")
        esperado[rel] = digest
    actual = {}
    for path in sorted(p for p in carpeta.rglob("*") if p.is_file() and not p.name.startswith("._") and p.name != ".DS_Store"):
        actual[str(path.relative_to(REPO))] = _sha256(path)
    assert actual == esperado
