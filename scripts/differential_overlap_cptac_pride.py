#!/usr/bin/env python3
"""
Intersección de proteínas diferenciales CPTAC vs PRIDE (reproducibilidad multi-cohorte).

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


def main() -> None:
    ap = argparse.ArgumentParser(description="Intersección CPTAC ∩ PRIDE en tablas diferenciales.")
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
        "--out-dir",
        type=Path,
        default=Path("results/reproducibility"),
        help="Directorio de salida",
    )
    ap.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Umbral padj para contar significativos (si existe columna padj)",
    )
    args = ap.parse_args()

    if not args.cptac.exists() or not args.pride.exists():
        missing = [str(p) for p in (args.cptac, args.pride) if not p.exists()]
        raise SystemExit(f"Faltan archivos: {missing}")

    c = pd.read_csv(args.cptac, sep="\t")
    p = pd.read_csv(args.pride, sep="\t")
    if "protein_id" not in c.columns or "protein_id" not in p.columns:
        raise SystemExit("Ambos TSV deben tener columna protein_id.")

    def _suffix(df: pd.DataFrame, tag: str) -> pd.DataFrame:
        skip = {"protein_id", "gene_symbol", "uniprot_accession"}
        cols = [x for x in df.columns if x not in skip]
        return df[["protein_id"] + cols].rename(columns={col: f"{col}_{tag}" for col in cols})

    merged = _suffix(c, "cptac").merge(_suffix(p, "pride"), on="protein_id", how="outer", indicator=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_merged = args.out_dir / "cptac_pride_differential_overlap.tsv"
    merged.to_csv(out_merged, sep="\t", index=False)

    def _sig(df: pd.DataFrame, prefix: str) -> pd.Series:
        col = f"padj_{prefix}"
        if col not in merged.columns:
            return pd.Series(False, index=merged.index)
        return merged[col].astype(float) < args.alpha

    in_both = merged["_merge"] == "both"
    sig_c = _sig(merged, "cptac")
    sig_p = _sig(merged, "pride")
    summary = pd.DataFrame(
        [
            {
                "metric": "n_proteins_cptac",
                "value": int(c["protein_id"].nunique()),
            },
            {
                "metric": "n_proteins_pride",
                "value": int(p["protein_id"].nunique()),
            },
            {"metric": "n_intersection_union", "value": int(in_both.sum())},
            {
                "metric": f"n_both_significant_padj_lt_{args.alpha}",
                "value": int((in_both & sig_c & sig_p).sum()),
            },
            {
                "metric": "n_cptac_sig_only",
                "value": int((in_both & sig_c & ~sig_p).sum()),
            },
            {
                "metric": "n_pride_sig_only",
                "value": int((in_both & ~sig_c & sig_p).sum()),
            },
        ]
    )
    sum_path = args.out_dir / "cptac_pride_differential_overlap_summary.tsv"
    summary.to_csv(sum_path, sep="\t", index=False)
    print(f"[overlap] detalle -> {out_merged}")
    print(f"[overlap] resumen -> {sum_path}")


if __name__ == "__main__":
    main()
