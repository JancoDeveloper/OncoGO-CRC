#!/usr/bin/env python3
"""Reduce **Human Protein Atlas** a un proxy de especificidad tumoral → snapshot Capa C.

Descarga (o usa en caché) ``proteinatlas.tsv.zip`` desde https://www.proteinatlas.org/about/download
y deriva una **especificidad tumoral (0–1)** por gen a partir de la categoría de especificidad de RNA.

Mapeo por defecto (categoría → score), configurable con ``--specificity-column``:
    Tissue enriched=1.0 · Group enriched=0.8 · Tissue enhanced=0.6 · Low tissue specificity=0.1

> Proxy documentado: HPA no expone en el archivo plano una «especificidad tumoral» numérica única;
> este score usa la especificidad de tejido como aproximación. Puedes apuntar a otra columna con
> ``--specificity-column`` si dispones de una métrica mejor.

Salida (``metadata/evidence/hpa_tumor_specificity.tsv``):
    ``gene_symbol   tumor_specificity   specificity_category``

Uso::

    python scripts/fetch_evidence_hpa.py --out metadata/evidence/hpa_tumor_specificity.tsv
"""

from __future__ import annotations

import argparse
import os
import sys
import zipfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _evidence_common as ec  # noqa: E402

HPA_URL = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
_SCORE_MAP = {
    "tissue enriched": 1.0,
    "group enriched": 0.8,
    "tissue enhanced": 0.6,
    "low tissue specificity": 0.1,
    "not detected": 0.0,
}


def _score(cat: str) -> float:
    return _SCORE_MAP.get(str(cat).strip().lower(), 0.3)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=Path("metadata/evidence/hpa_tumor_specificity.tsv"))
    p.add_argument("--cache", type=Path, default=Path("metadata/evidence/.hpa_cache/proteinatlas.tsv.zip"))
    p.add_argument("--url", default=HPA_URL)
    p.add_argument("--specificity-column", default="RNA tissue specificity")
    p.add_argument("--genes", nargs="*", default=None)
    p.add_argument("--genes-tsv", type=Path, default=None)
    p.add_argument("--differential-tsv", type=Path, default=None)
    args = p.parse_args()

    zpath = ec.download_file(args.url, args.cache, min_bytes=1_000_000)
    ec.log("[hpa] leyendo proteinatlas.tsv (grande) …")
    with zipfile.ZipFile(zpath) as zf:
        inner = next((n for n in zf.namelist() if n.endswith(".tsv")), None)
        if inner is None:
            ec.log("[hpa] el zip no contiene un .tsv.")
            return 1
        with zf.open(inner) as fh:
            df = pd.read_csv(fh, sep="\t", low_memory=False)

    gene_col = next((c for c in ("Gene", "Gene name", "gene") if c in df.columns), None)
    spec_col = args.specificity_column if args.specificity_column in df.columns else None
    if gene_col is None or spec_col is None:
        ec.log(f"[hpa] columnas no encontradas. Disponibles (muestra): {df.columns[:15].tolist()}")
        return 1

    out = pd.DataFrame(
        {
            "gene_symbol": df[gene_col].astype(str).str.strip().str.upper(),
            "specificity_category": df[spec_col].astype(str),
        }
    )
    out["tumor_specificity"] = out["specificity_category"].map(_score)

    genes = ec.load_gene_universe(
        genes_arg=args.genes, genes_tsv=args.genes_tsv, differential_tsv=args.differential_tsv
    )
    if genes:
        out = out[out["gene_symbol"].isin(set(genes))]
    out = out.drop_duplicates("gene_symbol").sort_values("tumor_specificity", ascending=False)
    out = out[["gene_symbol", "tumor_specificity", "specificity_category"]]

    ec.write_snapshot(out, args.out, source=f"Human Protein Atlas ({spec_col})")
    ec.log(f"[hpa] genes: {len(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
