#!/usr/bin/env python3
"""Descarga interacciones gen–fármaco de **DGIdb** (druggability) → snapshot Capa C.

Fuente: DGIdb GraphQL API (``https://dgidb.org/api/graphql``). Consulta por lotes los genes del
universo (firma diferencial por defecto) y agrega el número de fármacos que interactúan por gen.

Salida (``metadata/evidence/dgidb_gene_drugs.tsv``):
    ``gene_symbol   n_drugs   drugs``

Uso::

    python scripts/fetch_evidence_dgidb.py \\
        --differential-tsv results/differential/cptac_tumor_vs_normal_diff.tsv \\
        --out metadata/evidence/dgidb_gene_drugs.tsv
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _evidence_common as ec  # noqa: E402

API = "https://dgidb.org/api/graphql"
QUERY = """
query($names:[String!]!){
  genes(names:$names){
    nodes{ name interactions{ drug{ name } } }
  }
}
"""


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=Path("metadata/evidence/dgidb_gene_drugs.tsv"))
    p.add_argument("--genes", nargs="*", default=None, help="Lista explícita de genes (para pruebas).")
    p.add_argument("--genes-tsv", type=Path, default=None)
    p.add_argument(
        "--differential-tsv",
        type=Path,
        default=Path("results/differential/cptac_tumor_vs_normal_diff.tsv"),
    )
    p.add_argument("--mapping-tsv", type=Path, default=None)
    p.add_argument("--limit", type=int, default=0, help="Máx. genes a consultar (0 = todos).")
    p.add_argument("--batch-size", type=int, default=100)
    p.add_argument("--max-drugs-listed", type=int, default=8)
    args = p.parse_args()

    genes = ec.load_gene_universe(
        genes_arg=args.genes,
        genes_tsv=args.genes_tsv,
        differential_tsv=args.differential_tsv,
        mapping_tsv=args.mapping_tsv,
        limit=args.limit,
    )
    if not genes:
        ec.log("[dgidb] no hay universo de genes; especifica --genes / --genes-tsv / --differential-tsv.")
        return 1
    ec.log(f"[dgidb] consultando {len(genes)} genes en lotes de {args.batch_size} …")

    rows: list[dict] = []
    for batch in ec.chunked(genes, args.batch_size):
        data = ec.http_json(API, payload={"query": QUERY, "variables": {"names": batch}})
        nodes = (((data or {}).get("data") or {}).get("genes") or {}).get("nodes") or []
        for node in nodes:
            name = str(node.get("name", "")).strip().upper()
            drugs = sorted(
                {
                    str((it.get("drug") or {}).get("name", "")).strip().upper()
                    for it in (node.get("interactions") or [])
                    if (it.get("drug") or {}).get("name")
                }
            )
            if name:
                rows.append(
                    {
                        "gene_symbol": name,
                        "n_drugs": len(drugs),
                        "drugs": ";".join(drugs[: args.max_drugs_listed]),
                    }
                )

    if not rows:
        ec.log("[dgidb] sin resultados (¿genes no reconocidos?).")
        return 1
    df = pd.DataFrame(rows).drop_duplicates("gene_symbol").sort_values("n_drugs", ascending=False)
    ec.write_snapshot(df, args.out, source="DGIdb GraphQL (dgidb.org)")
    ec.log(f"[dgidb] genes con ≥1 fármaco: {(df['n_drugs'] > 0).sum()} / {len(df)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
