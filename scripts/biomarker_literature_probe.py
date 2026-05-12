#!/usr/bin/env python3
"""
Tabla de concordancia literatura ↔ resultados diferenciales (FDR).

Por defecto busca TFF3, LCN2, CEACAM5 en CPTAC, PRIDE y plasma (config/results).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _rank_row(df: pd.DataFrame, gene: str, alpha: float) -> dict[str, object]:
    if df is None or df.empty or "padj" not in df.columns:
        return {"present_sig": False, "padj": None, "rank_pct": None}
    g = gene.upper()
    if "gene_symbol" in df.columns:
        mask = df["gene_symbol"].astype(str).str.upper() == g
    else:
        mask = df["protein_id"].astype(str).str.upper() == g
    sub = df[mask]
    if sub.empty and "gene_symbol" in df.columns:
        sub = df[df["protein_id"].astype(str).str.upper() == g]
    if sub.empty:
        return {"present_sig": False, "padj": None, "rank_pct": None}
    padj = float(sub.iloc[0]["padj"])
    order = df.sort_values("padj", na_position="last").reset_index(drop=True)
    match_col = "gene_symbol" if "gene_symbol" in order.columns else "protein_id"
    pos = order.index[order[match_col].astype(str).str.upper() == g].tolist()
    rank = int(pos[0]) + 1 if pos else None
    n = len(order)
    pct = round(100.0 * (rank / n), 1) if rank and n else None
    return {"present_sig": bool(padj < alpha), "padj": padj, "rank_pct": pct, "rank": rank}


def main() -> None:
    ap = argparse.ArgumentParser(description="Sondeo TFF3/LCN2/CEACAM5 vs diferenciales.")
    ap.add_argument("--results-dir", type=Path, default=Path("results/differential"))
    ap.add_argument("--out", type=Path, default=Path("results/reproducibility/biomarker_literature_probe.tsv"))
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument(
        "--genes",
        type=str,
        default="TFF3,LCN2,CEACAM5",
        help="Lista separada por comas",
    )
    args = ap.parse_args()
    genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()]
    comps = {
        "CPTAC": args.results_dir / "cptac_tumor_vs_normal_diff.tsv",
        "PRIDE": args.results_dir / "pride_validation_tumor_vs_normal_diff.tsv",
        "plasma": args.results_dir / "plasma_crc_vs_healthy_diff.tsv",
    }
    rows: list[dict[str, object]] = []
    loaded: dict[str, pd.DataFrame] = {}
    for label, p in comps.items():
        loaded[label] = pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()

    for gene in genes:
        for label, df in loaded.items():
            r = _rank_row(df, gene, args.alpha)
            rows.append(
                {
                    "gene_symbol": gene,
                    "cohort": label,
                    "diff_tsv_exists": label in loaded and not loaded[label].empty,
                    "significant_padj_lt_alpha": r["present_sig"],
                    "padj": r["padj"],
                    "rank_by_padj": r.get("rank"),
                    "percentile_rank_approx": r["rank_pct"],
                }
            )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print(f"[biomarker_probe] -> {args.out}")


if __name__ == "__main__":
    main()
