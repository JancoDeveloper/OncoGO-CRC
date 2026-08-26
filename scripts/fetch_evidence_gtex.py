#!/usr/bin/env python3
"""Descarga expresión en tejidos normales (**GTEx**) → snapshot Capa C (seguridad/ventana terapéutica).

Fuente: GTEx v8 medianas de TPM por tejido (GCT). Para cada gen calcula la **máxima mediana de TPM**
entre tejidos normales (menor valor ⇒ mejor ventana terapéutica) y el tejido con el máximo.

Salida (``metadata/evidence/gtex_normal_max.tsv``):
    ``gene_symbol   max_normal_tpm   top_tissue``

Uso::

    python scripts/fetch_evidence_gtex.py --out metadata/evidence/gtex_normal_max.tsv
    # opcional: restringir al universo de genes de la firma
    python scripts/fetch_evidence_gtex.py --differential-tsv results/differential/cptac_tumor_vs_normal_diff.tsv
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _evidence_common as ec  # noqa: E402

GTEX_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=Path("metadata/evidence/gtex_normal_max.tsv"))
    p.add_argument("--cache", type=Path, default=Path("metadata/evidence/.gtex_cache/gene_median_tpm.gct.gz"))
    p.add_argument("--url", default=GTEX_URL)
    p.add_argument("--genes", nargs="*", default=None)
    p.add_argument("--genes-tsv", type=Path, default=None)
    p.add_argument("--differential-tsv", type=Path, default=None)
    p.add_argument("--mapping-tsv", type=Path, default=None)
    args = p.parse_args()

    gz = ec.download_file(args.url, args.cache, min_bytes=1_000_000)

    ec.log("[gtex] leyendo GCT (medianas por tejido) …")
    # GCT: línea 1 = versión, línea 2 = dims, línea 3 = cabecera (Name, Description, tejidos…)
    with gzip.open(gz, "rt") as fh:
        _ = fh.readline()
        _ = fh.readline()
        df = pd.read_csv(fh, sep="\t")

    id_cols = [c for c in ("Name", "Description") if c in df.columns]
    tissue_cols = [c for c in df.columns if c not in id_cols]
    if "Description" not in df.columns or not tissue_cols:
        ec.log(f"[gtex] formato inesperado; columnas: {df.columns[:5].tolist()}")
        return 1

    df["gene_symbol"] = df["Description"].astype(str).str.strip().str.upper()
    tpm = df[tissue_cols].apply(pd.to_numeric, errors="coerce")
    df["max_normal_tpm"] = tpm.max(axis=1).round(4)
    df["top_tissue"] = tpm.idxmax(axis=1)

    genes = ec.load_gene_universe(
        genes_arg=args.genes,
        genes_tsv=args.genes_tsv,
        differential_tsv=args.differential_tsv,
        mapping_tsv=args.mapping_tsv,
    )
    out = df[["gene_symbol", "max_normal_tpm", "top_tissue"]].copy()
    if genes:
        out = out[out["gene_symbol"].isin(set(genes))]
    # dedup: mayor expresión por símbolo
    out = out.sort_values("max_normal_tpm", ascending=False).drop_duplicates("gene_symbol")

    ec.write_snapshot(out, args.out, source="GTEx v8 gene median TPM")
    ec.log(f"[gtex] genes: {len(out)} (mediana max_normal_tpm={out['max_normal_tpm'].median():.2f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
