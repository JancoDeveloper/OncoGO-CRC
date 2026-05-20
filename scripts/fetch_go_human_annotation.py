#!/usr/bin/env python3
"""Descarga la anotación GO humana real (GOA) y la convierte al formato ``go_gene_long.tsv``.

Fuente: archivo oficial ``goa_human.gaf.gz`` del Gene Ontology Consortium
(`https://current.geneontology.org/annotations/goa_human.gaf.gz`). El formato GAF 2.x tiene:

    DB  DB_Object_ID  DB_Object_Symbol  Qualifier  GO_ID  DB:Reference  Evidence_Code  ...  Aspect  ...

donde ``DB_Object_Symbol`` es el símbolo HGNC y ``Aspect`` ∈ {P, F, C} (BP, MF, CC).

Filtra evidencia "no IEA" por defecto (más confiable) y opcionalmente la incluye con ``--include-iea``.

Uso típico::

    python scripts/fetch_go_human_annotation.py \\
        --out metadata/gene_ontology/go_gene_long.tsv

Salida (TSV): ``gene_symbol``, ``go_id``, ``namespace`` (``BP|MF|CC``).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import sys
import urllib.error
import urllib.request
from pathlib import Path


USER_AGENT = "OncoGO-CRC/goa-fetch (https://github.com; academic research)"
GOA_URL = "https://current.geneontology.org/annotations/goa_human.gaf.gz"
ASPECT_TO_NS = {"P": "BP", "F": "MF", "C": "CC"}


def _download_bytes(url: str, *, timeout_s: float) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT}, method="GET")
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            return resp.read()
    except urllib.error.HTTPError as e:
        print(f"[goa-fetch] HTTP {e.code}: {e.reason}", file=sys.stderr)
        raise SystemExit(1) from e
    except urllib.error.URLError as e:
        print(f"[goa-fetch] red: {e}", file=sys.stderr)
        raise SystemExit(1) from e


def parse_gaf(
    raw_gz: bytes,
    *,
    include_iea: bool,
    namespaces: set[str],
) -> set[tuple[str, str, str]]:
    seen: set[tuple[str, str, str]] = set()
    with gzip.GzipFile(fileobj=io.BytesIO(raw_gz)) as gf:
        for raw in io.TextIOWrapper(gf, encoding="utf-8", errors="replace"):
            if not raw or raw.startswith("!"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            symbol = cols[2].strip().upper()
            qualifier = cols[3].strip().upper()
            go_id = cols[4].strip()
            evidence = cols[6].strip().upper()
            aspect = cols[8].strip().upper()
            if not symbol or not go_id.startswith("GO:"):
                continue
            ns = ASPECT_TO_NS.get(aspect, "")
            if ns not in namespaces:
                continue
            if "NOT" in qualifier.split("|"):
                continue
            if not include_iea and evidence == "IEA":
                continue
            seen.add((symbol, go_id, ns))
    return seen


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("metadata/gene_ontology/go_gene_long.tsv"),
        help="Ruta TSV de salida (relativa al CWD).",
    )
    p.add_argument("--include-iea", action="store_true", help="Incluye evidencia IEA (inferida electrónicamente).")
    p.add_argument(
        "--namespaces",
        nargs="+",
        default=["BP", "MF", "CC"],
        help="Subontologías a incluir (BP/MF/CC).",
    )
    p.add_argument("--timeout", type=float, default=180.0, help="Timeout HTTP (segundos).")
    p.add_argument(
        "--source-url",
        type=str,
        default=GOA_URL,
        help="URL del archivo GAF gzip (por defecto current GO Consortium).",
    )
    args = p.parse_args()

    namespaces = {ns.upper() for ns in args.namespaces}
    print(f"[goa-fetch] descargando {args.source_url} …", file=sys.stderr)
    raw_gz = _download_bytes(args.source_url, timeout_s=args.timeout)
    print(f"[goa-fetch] {len(raw_gz) / 1_048_576:.1f} MiB recibidos; parseando GAF…", file=sys.stderr)

    rows = parse_gaf(raw_gz, include_iea=args.include_iea, namespaces=namespaces)
    n_genes = len({r[0] for r in rows})
    n_terms = len({r[1] for r in rows})

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_symbol", "go_id", "namespace"])
        for symbol, go_id, ns in sorted(rows):
            w.writerow([symbol, go_id, ns])

    print(
        f"[goa-fetch] escrito: {args.out.resolve()} "
        f"({len(rows)} filas; {n_genes} genes únicos; {n_terms} términos únicos; "
        f"IEA={'incluido' if args.include_iea else 'excluido'})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
