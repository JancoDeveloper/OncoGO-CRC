"""Manifiestos de hash para congelar v1 y verificar que las entradas no cambian.

Dos usos:

    python scripts/v2/manifest_hashes.py write   # escribe los dos manifiestos
    python scripts/v2/manifest_hashes.py verify  # recalcula y compara

El manifiesto de entradas cubre las matrices crudas y la metadata que el
pipeline lee. El de v1 cubre todo lo que la corrida certificada escribió. Si
``verify`` encuentra una diferencia, sale con código 1.

Los ficheros ``._*`` de AppleDouble que ExFAT crea en el disco externo no son
datos del proyecto y quedan fuera de los dos manifiestos.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WORKSPACE = Path("/Volumes/JANCODATA/OncoGO-CRC")

# Entradas: lo que el pipeline lee y nunca debe editar.
INPUT_GLOBS = (
    "data/raw/cptac/cptac_crc_matrix.tsv",
    "data/raw/cptac/provider_manifest.json",
    "data/raw/pride/pride_crc_matrix.tsv",
    "data/raw/pride/Cohort1_proteinGroups.txt",
    "data/raw/pride/cohort1_unique_peptides.tsv",
    "data/raw/pride/provider_manifest.json",
    "data/raw/plasma/plasma_crc_matrix.tsv",
    "data/raw/plasma/provider_manifest.json",
    "metadata/protein_id_maps/protein_mapping_master.tsv",
    "metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv",
    "metadata/gene_ontology/go_gene_long.tsv",
)

# v1: la corrida certificada del 1 de octubre de 2026. Queda congelada.
V1_DIRS = (
    "results/certified",
    "results/certification",
    "results/secondary",
    "results/axes",
)
V1_FILES = ("OncoGO-CRC_tesis_explicada.pdf",)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def is_noise(path: Path) -> bool:
    """AppleDouble de ExFAT y basura del sistema operativo. No son datos."""
    return path.name.startswith("._") or path.name in {".DS_Store", ".oncogo_write_probe"}


def collect_inputs() -> list[Path]:
    out = []
    for rel in INPUT_GLOBS:
        p = WORKSPACE / rel
        if p.is_file():
            out.append(p)
    return sorted(out)


def collect_v1() -> list[Path]:
    out = []
    for rel in V1_DIRS:
        base = WORKSPACE / rel
        if not base.is_dir():
            continue
        out += [p for p in base.rglob("*") if p.is_file() and not is_noise(p)]
    for rel in V1_FILES:
        p = WORKSPACE / rel
        if p.is_file():
            out.append(p)
    # El manifiesto de v2 no entra en el de v1.
    return sorted(p for p in out if "/v2" not in str(p))


def render(paths: list[Path]) -> str:
    lines = [f"{sha256(p)}  {p.relative_to(WORKSPACE)}" for p in paths]
    return "\n".join(lines) + "\n"


def parse(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        digest, _, rel = line.partition("  ")
        out[rel] = digest
    return out


def compare(name: str, expected: dict[str, str], actual: dict[str, str]) -> int:
    missing = sorted(set(expected) - set(actual))
    added = sorted(set(actual) - set(expected))
    changed = sorted(r for r in set(expected) & set(actual) if expected[r] != actual[r])
    if not (missing or added or changed):
        print(f"[{name}] coincide: {len(expected)} ficheros sin cambios.")
        return 0
    for rel in missing:
        print(f"[{name}] FALTA     {rel}")
    for rel in changed:
        print(f"[{name}] CAMBIADO  {rel}")
    for rel in added:
        print(f"[{name}] NUEVO     {rel}")
    return 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("action", choices=("write", "verify"))
    args = ap.parse_args(argv)

    if not WORKSPACE.is_dir():
        print(f"[manifest] workspace no accesible: {WORKSPACE}", file=sys.stderr)
        return 2

    targets = {
        "entradas": (WORKSPACE / "results/v2/manifest_inputs.sha256", collect_inputs()),
        "v1": (WORKSPACE / "results/v1_manifest.sha256", collect_v1()),
    }
    rc = 0
    for name, (out_path, paths) in targets.items():
        text = render(paths)
        if args.action == "write":
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(text, encoding="utf-8")
            print(f"[{name}] {len(paths)} ficheros -> {out_path}")
        else:
            if not out_path.is_file():
                print(f"[{name}] no existe el manifiesto {out_path}", file=sys.stderr)
                rc = 1
                continue
            rc |= compare(name, parse(out_path.read_text(encoding="utf-8")), parse(text))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
