#!/usr/bin/env python3
"""Descarga **tractabilidad** de dianas de **Open Targets** → snapshot Capa C.

Fuente: Open Targets Platform GraphQL API (``https://api.platform.opentargets.org/api/v4/graphql``).
Resuelve cada gen a su Ensembl ID (búsqueda o columna ``ensembl_gene_id`` del mapa) y consulta los
*buckets* de tractabilidad (SM/AB/…); calcula una tractabilidad agregada 0–1.

Salida (``metadata/evidence/opentargets_tractability.tsv``):
    ``gene_symbol   tractability   tractable_sm   tractable_ab   ensembl_id``

Uso::

    python scripts/fetch_evidence_opentargets.py \\
        --differential-tsv results/differential/cptac_tumor_vs_normal_diff.tsv \\
        --out metadata/evidence/opentargets_tractability.tsv --limit 300
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _evidence_common as ec  # noqa: E402

API = "https://api.platform.opentargets.org/api/v4/graphql"
SEARCH_Q = """
query($q:String!){ search(queryString:$q, entityNames:["target"]){ hits{ id name } } }
"""
TARGET_Q = """
query($id:String!){ target(ensemblId:$id){ approvedSymbol tractability{ modality value } } }
"""


def resolve_ensembl(gene: str, cache: dict[str, str]) -> str | None:
    if gene in cache:
        return cache[gene] or None
    data = ec.http_json(API, payload={"query": SEARCH_Q, "variables": {"q": gene}})
    hits = (((data or {}).get("data") or {}).get("search") or {}).get("hits") or []
    ens = ""
    for h in hits:
        if str(h.get("name", "")).strip().upper() == gene and str(h.get("id", "")).startswith("ENSG"):
            ens = h["id"]
            break
    if not ens and hits and str(hits[0].get("id", "")).startswith("ENSG"):
        ens = hits[0]["id"]
    cache[gene] = ens
    return ens or None


def tractability_for(ensembl_id: str) -> dict | None:
    data = ec.http_json(API, payload={"query": TARGET_Q, "variables": {"id": ensembl_id}})
    tgt = ((data or {}).get("data") or {}).get("target")
    if not tgt:
        return None
    buckets = tgt.get("tractability") or []
    if not buckets:
        return {"tractability": 0.0, "tractable_sm": False, "tractable_ab": False}
    n_true = sum(1 for b in buckets if b.get("value"))
    sm = any(b.get("value") and b.get("modality") == "SM" for b in buckets)
    ab = any(b.get("value") and b.get("modality") == "AB" for b in buckets)
    return {"tractability": round(n_true / len(buckets), 4), "tractable_sm": bool(sm), "tractable_ab": bool(ab)}


def _ensembl_from_mapping(mapping_tsv: Path | None) -> dict[str, str]:
    if not mapping_tsv or not Path(mapping_tsv).exists():
        return {}
    df = pd.read_csv(mapping_tsv, sep="\t")
    if "gene_symbol" not in df.columns or "ensembl_gene_id" not in df.columns:
        return {}
    out: dict[str, str] = {}
    for _, r in df.iterrows():
        g = str(r["gene_symbol"]).strip().upper()
        e = str(r["ensembl_gene_id"]).strip().split(";")[0]
        if g and e.startswith("ENSG"):
            out.setdefault(g, e)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=Path("metadata/evidence/opentargets_tractability.tsv"))
    p.add_argument("--genes", nargs="*", default=None)
    p.add_argument("--genes-tsv", type=Path, default=None)
    p.add_argument(
        "--differential-tsv",
        type=Path,
        default=Path("results/differential/cptac_tumor_vs_normal_diff.tsv"),
    )
    p.add_argument("--mapping-tsv", type=Path, default=None, help="Mapa con ensembl_gene_id (evita búsquedas).")
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--cache", type=Path, default=Path("metadata/evidence/.ot_cache/ensembl.json"))
    p.add_argument("--sleep", type=float, default=0.05, help="Pausa entre peticiones (cortesía).")
    args = p.parse_args()

    genes = ec.load_gene_universe(
        genes_arg=args.genes,
        genes_tsv=args.genes_tsv,
        differential_tsv=args.differential_tsv,
        mapping_tsv=args.mapping_tsv,
        limit=args.limit,
    )
    if not genes:
        ec.log("[opentargets] no hay universo de genes.")
        return 1

    cache: dict[str, str] = {}
    if args.cache.exists():
        try:
            cache = json.loads(args.cache.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            cache = {}
    cache.update(_ensembl_from_mapping(args.mapping_tsv))

    ec.log(f"[opentargets] consultando tractabilidad de {len(genes)} genes …")
    rows: list[dict] = []
    for i, g in enumerate(genes, 1):
        try:
            ens = resolve_ensembl(g, cache)
            if not ens:
                continue
            tr = tractability_for(ens)
            if tr is None:
                continue
            rows.append({"gene_symbol": g, **tr, "ensembl_id": ens})
        except SystemExit:
            ec.log(f"[opentargets] error de red en {g}; continúo.")
        if args.sleep:
            time.sleep(args.sleep)
        if i % 50 == 0:
            ec.log(f"[opentargets] {i}/{len(genes)} …")

    args.cache.parent.mkdir(parents=True, exist_ok=True)
    args.cache.write_text(json.dumps(cache), encoding="utf-8")

    if not rows:
        ec.log("[opentargets] sin resultados.")
        return 1
    df = pd.DataFrame(rows).drop_duplicates("gene_symbol").sort_values("tractability", ascending=False)
    ec.write_snapshot(df, args.out, source="Open Targets Platform GraphQL v4")
    ec.log(f"[opentargets] genes con tractabilidad>0: {(df['tractability'] > 0).sum()} / {len(df)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
