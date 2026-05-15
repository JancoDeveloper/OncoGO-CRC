#!/usr/bin/env python3
"""
Intersección de proteínas diferenciales CPTAC vs PRIDE (reproducibilidad multi-cohorte).

Usa el mismo **puente UniProt** que la etapa `cohort_overlap` (mapa gene_symbol ↔ uniprot_accession).

Lee por defecto:
  results/differential/cptac_tumor_vs_normal_diff.tsv
  results/differential/pride_validation_tumor_vs_normal_diff.tsv

Escribe:
  results/reproducibility/cptac_pride_differential_overlap.tsv
  results/reproducibility/cptac_pride_differential_overlap_summary.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import load_settings
from oncogo_crc.utils.reproducibility_overlap import (
    cptac_pride_overlap_tables,
    load_gene_uniprot_from_settings_map,
)


def main() -> None:
    ap = argparse.ArgumentParser(description="Intersección CPTAC ∩ PRIDE en tablas diferenciales.")
    ap.add_argument("--root", type=Path, default=None, help="Raíz del repo (para resolver mapa maestro)")
    ap.add_argument(
        "--cptac",
        type=Path,
        default=Path("results/differential/cptac_tumor_vs_normal_diff.tsv"),
        help="TSV diferencial CPTAC",
    )
    ap.add_argument(
        "--pride",
        type=Path,
        default=Path("results/differential/pride_validation_tumor_vs_normal_diff.tsv"),
        help="TSV diferencial PRIDE",
    )
    ap.add_argument(
        "--map",
        type=Path,
        default=None,
        help="TSV mapa gene_symbol,uniprot_accession (default: desde settings / metadata)",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/reproducibility"),
        help="Directorio de salida",
    )
    ap.add_argument(
        "--lenient",
        action="store_true",
        help="Incluir genes con mapeo UniProt ambiguo (por defecto: mapa estricto como cohort_overlap).",
    )
    args = ap.parse_args()

    if not args.cptac.exists() or not args.pride.exists():
        missing = [str(p) for p in (args.cptac, args.pride) if not p.exists()]
        raise SystemExit(f"Faltan archivos: {missing}")

    c = pd.read_csv(args.cptac, sep="\t")
    p = pd.read_csv(args.pride, sep="\t")
    if "protein_id" not in c.columns or "protein_id" not in p.columns:
        raise SystemExit("Ambos TSV deben tener columna protein_id.")

    if args.map is not None:
        gene_to_uni = load_gene_uniprot_from_settings_map(args.map, exclude_ambiguous=not args.lenient)
    else:
        settings = load_settings(args.root)
        gene_to_uni = load_gene_uniprot_from_settings_map(
            settings.protein_id_map_path(), exclude_ambiguous=not args.lenient
        )

    merged, summary = cptac_pride_overlap_tables(c, p, gene_to_uni=gene_to_uni, alpha=args.alpha)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_merged = args.out_dir / "cptac_pride_differential_overlap.tsv"
    sum_path = args.out_dir / "cptac_pride_differential_overlap_summary.tsv"
    merged.to_csv(out_merged, sep="\t", index=False)
    summary.to_csv(sum_path, sep="\t", index=False)
    print(f"[overlap] detalle -> {out_merged}")
    print(f"[overlap] resumen -> {sum_path}")


if __name__ == "__main__":
    main()
