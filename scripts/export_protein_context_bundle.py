#!/usr/bin/env python3
"""Exporta un JSON listo para enriquecimiento externo (p. ej. ProteinsLanguage, UniProt, STRING).

Uso (desde la raíz del repo, con el venv activo):
  python scripts/export_protein_context_bundle.py
  python scripts/export_protein_context_bundle.py --comparison cptac_tumor_vs_normal --top 300
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import load_settings
from oncogo_crc.utils.io import read_tsv


def main() -> None:
    ap = argparse.ArgumentParser(description="Bundle de proteínas priorizadas para contexto GO/estructura.")
    ap.add_argument("--root", type=Path, default=None, help="Raíz del repo OncoGO-CRC")
    ap.add_argument("--comparison", type=str, default="cptac_tumor_vs_normal")
    ap.add_argument("--top", type=int, default=250, help="Máximo de filas del diferencial por padj")
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Ruta de salida JSON (por defecto: results/export/protein_context_bundle.json)",
    )
    args = ap.parse_args()

    settings = load_settings(args.root)
    res = settings.resolve(settings.config["paths"]["results_dir"])
    proc = settings.resolve(settings.config["paths"]["processed_dir"])
    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")

    out_path = args.out or (res / "export" / "protein_context_bundle.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    mapping = read_tsv(settings.protein_id_map_path())
    map_cols = [c for c in ["protein_id", "gene_symbol", "uniprot_accession"] if c in mapping.columns]

    diff_path = res / "differential" / f"{args.comparison}_diff.tsv"
    rows: list[dict] = []
    if diff_path.exists():
        diff = read_tsv(diff_path)
        diff = diff.sort_values("padj", na_position="last").head(args.top)
        merged = diff.merge(mapping[map_cols], on="protein_id", how="left") if "protein_id" in mapping.columns else diff
        for _, r in merged.iterrows():
            rows.append(
                {
                    "source": "differential",
                    "comparison": args.comparison,
                    "protein_id": str(r.get("protein_id", "")),
                    "gene_symbol": str(r.get("gene_symbol", "")) if pd.notna(r.get("gene_symbol")) else "",
                    "uniprot_accession": str(r.get("uniprot_accession", ""))
                    if pd.notna(r.get("uniprot_accession"))
                    else "",
                    "log2fc": float(r["log2fc"]) if pd.notna(r.get("log2fc")) else None,
                    "padj": float(r["padj"]) if pd.notna(r.get("padj")) else None,
                }
            )

    ml_path = res / "ml" / f"{primary}_feature_ranking.tsv"
    if ml_path.exists():
        ml = read_tsv(ml_path).head(80)
        ml = ml.merge(mapping[map_cols], on="protein_id", how="left") if "protein_id" in ml.columns else ml
        for _, r in ml.iterrows():
            rows.append(
                {
                    "source": "ml_ranking",
                    "comparison": None,
                    "protein_id": str(r.get("protein_id", "")),
                    "gene_symbol": str(r.get("gene_symbol", "")) if pd.notna(r.get("gene_symbol")) else "",
                    "uniprot_accession": str(r.get("uniprot_accession", ""))
                    if pd.notna(r.get("uniprot_accession"))
                    else "",
                    "combined_score": float(r["combined_score"]) if "combined_score" in r and pd.notna(r["combined_score"]) else None,
                }
            )

    bundle = {
        "project": settings.config.get("project", {}).get("name", "OncoGO-CRC"),
        "primary_dataset": primary,
        "comparison": args.comparison,
        "n_entries": len(rows),
        "note": "Enriquecer GO/AlphaFold/grafos en ProteinsLanguage u otras herramientas usando gene_symbol o UniProt.",
        "proteins": rows,
    }
    out_path.write_text(json.dumps(bundle, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[export] {settings.relative_to_project(out_path)} ({len(rows)} entradas)")


if __name__ == "__main__":
    main()
