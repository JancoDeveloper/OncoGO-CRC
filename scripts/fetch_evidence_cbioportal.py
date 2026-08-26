#!/usr/bin/env python3
"""Descarga **frecuencia de mutación** por gen en CCR desde **cBioPortal** → snapshot Capa C.

Fuente: cBioPortal REST API (``https://www.cbioportal.org/api``). Por defecto usa el estudio
TCGA COADREAD PanCancer Atlas. Frecuencia = (muestras con mutación del gen) / (muestras del estudio).

> Nota: es una frecuencia **basada en mutaciones** (no incluye CNA/fusiones); es un proxy razonable
> del contexto genómico. Documentado como tal.

Salida (``metadata/evidence/cbioportal_alteration_freq.tsv``):
    ``gene_symbol   alteration_frequency   n_mutated_samples   n_total_samples   entrez_id``

Uso::

    python scripts/fetch_evidence_cbioportal.py \\
        --differential-tsv results/differential/cptac_tumor_vs_normal_diff.tsv \\
        --out metadata/evidence/cbioportal_alteration_freq.tsv
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _evidence_common as ec  # noqa: E402

BASE = "https://www.cbioportal.org/api"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=Path("metadata/evidence/cbioportal_alteration_freq.tsv"))
    p.add_argument("--genes", nargs="*", default=None)
    p.add_argument("--genes-tsv", type=Path, default=None)
    p.add_argument(
        "--differential-tsv",
        type=Path,
        default=Path("results/differential/cptac_tumor_vs_normal_diff.tsv"),
    )
    p.add_argument("--mapping-tsv", type=Path, default=None)
    p.add_argument("--study", default="coadread_tcga_pan_can_atlas_2018")
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--batch-size", type=int, default=500)
    args = p.parse_args()

    profile = f"{args.study}_mutations"
    sample_list = f"{args.study}_all"

    genes = ec.load_gene_universe(
        genes_arg=args.genes,
        genes_tsv=args.genes_tsv,
        differential_tsv=args.differential_tsv,
        mapping_tsv=args.mapping_tsv,
        limit=args.limit,
    )
    if not genes:
        ec.log("[cbioportal] no hay universo de genes.")
        return 1

    # Tamaño del estudio
    sl = ec.http_json(f"{BASE}/sample-lists/{sample_list}")
    n_total = int(sl.get("sampleCount", 0)) if isinstance(sl, dict) else 0
    if n_total <= 0:
        ec.log(f"[cbioportal] no pude obtener sampleCount de {sample_list}.")
        return 1
    ec.log(f"[cbioportal] estudio {args.study}: {n_total} muestras. Resolviendo {len(genes)} genes a Entrez …")

    # Resolver símbolos → Entrez
    entrez_to_symbol: dict[int, str] = {}
    for batch in ec.chunked(genes, 1000):
        recs = ec.http_json(f"{BASE}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", payload=batch)
        for r in recs or []:
            eid = r.get("entrezGeneId")
            sym = str(r.get("hugoGeneSymbol", "")).strip().upper()
            if eid and sym:
                entrez_to_symbol[int(eid)] = sym
    if not entrez_to_symbol:
        ec.log("[cbioportal] no se resolvieron genes a Entrez.")
        return 1
    ec.log(f"[cbioportal] {len(entrez_to_symbol)} genes resueltos. Consultando mutaciones …")

    # Muestras mutadas por gen
    mutated_samples: dict[int, set[str]] = {}
    entrez_ids = list(entrez_to_symbol.keys())
    for batch in ec.chunked(entrez_ids, args.batch_size):
        recs = ec.http_json(
            f"{BASE}/molecular-profiles/{profile}/mutations/fetch?projection=SUMMARY",
            payload={"sampleListId": sample_list, "entrezGeneIds": batch},
        )
        for r in recs or []:
            eid = r.get("entrezGeneId")
            sid = r.get("sampleId")
            if eid and sid:
                mutated_samples.setdefault(int(eid), set()).add(str(sid))

    rows: list[dict] = []
    for eid, sym in entrez_to_symbol.items():
        n_mut = len(mutated_samples.get(eid, set()))
        rows.append(
            {
                "gene_symbol": sym,
                "alteration_frequency": round(n_mut / n_total, 5),
                "n_mutated_samples": n_mut,
                "n_total_samples": n_total,
                "entrez_id": eid,
            }
        )
    df = pd.DataFrame(rows).drop_duplicates("gene_symbol").sort_values("alteration_frequency", ascending=False)
    ec.write_snapshot(df, args.out, source=f"cBioPortal REST ({args.study}, mutations)")
    ec.log(f"[cbioportal] genes mutados (>0): {(df['alteration_frequency'] > 0).sum()} / {len(df)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
