#!/usr/bin/env python3
"""Exporta un JSON listo para enriquecimiento externo (p. ej. ProteinsLanguage, UniProt, STRING).

Uso (desde la raíz del repo, con el venv activo):
  python scripts/export_protein_context_bundle.py
  python scripts/export_protein_context_bundle.py --comparison cptac_tumor_vs_normal --top 300
"""

from __future__ import annotations

import argparse
from pathlib import Path

from oncogo_crc.settings import load_settings
from oncogo_crc.utils.export_bundle import write_protein_context_bundle


def main() -> None:
    ap = argparse.ArgumentParser(description="Bundle de proteínas priorizadas para contexto GO/estructura.")
    ap.add_argument("--root", type=Path, default=None, help="Raíz del repo OncoGO-CRC")
    ap.add_argument("--comparison", type=str, default=None, help="Comparación (default: export_bundle.primary_comparison)")
    ap.add_argument("--top", type=int, default=None, help="Máximo filas diferencial (default: export_bundle.diff_top)")
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Ruta de salida JSON (por defecto: results/export/protein_context_bundle.json)",
    )
    args = ap.parse_args()

    settings = load_settings(args.root)
    res = settings.resolve(settings.config["paths"]["results_dir"])
    exp = settings.config.get("export_bundle", {}) or {}
    comp = args.comparison or str(exp.get("primary_comparison", "cptac_tumor_vs_normal"))
    top = int(args.top if args.top is not None else exp.get("diff_top", 250))
    out_path = args.out or (res / "export" / "protein_context_bundle.json")
    write_protein_context_bundle(
        settings,
        out_path,
        comparison=comp,
        top=top,
        uniprot_top_n=int(exp.get("uniprot_top_n", 0)),
    )
    print(f"[export] {settings.relative_to_project(out_path)}")


if __name__ == "__main__":
    main()
