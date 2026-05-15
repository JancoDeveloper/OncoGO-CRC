#!/usr/bin/env python3
"""Descarga paginada UniProt KB (humanos, Swiss-Prot) → TSV compatible con ``mapping.supplementary_map_paths``.

Usa la API de búsqueda ``/uniprotkb/search`` con paginación por ``cursor`` (cabecera HTTP ``Link: …; rel="next"``).

Ejemplo::

    python scripts/fetch_uniprot_human_reviewed_gene_map.py \\
        --out metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv

Opciones útiles: ``--page-size 500``, ``--max-rows 5000`` (prueba), ``--sleep 0.2`` entre páginas.

Referencia: https://rest.uniprot.org/ — usar un ``User-Agent`` identificable (política de uso razonable).
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


USER_AGENT = "OncoGO-CRC/uniprot-fetch (https://github.com; academic research)"
BASE_QUERY = "(organism_id:9606) AND (reviewed:true)"
FIELDS = "accession,gene_primary,gene_synonym"
RE_NEXT = re.compile(r"<([^>]+)>;\s*rel=\"next\"", re.IGNORECASE)

OUT_COLS = [
    "protein_id",
    "gene_symbol",
    "uniprot_accession",
    "uniprot_base",
    "ensembl_gene_id",
    "alias_symbols",
    "reviewed_status",
    "mapping_source",
    "is_ambiguous",
]


def _uniprot_base(accession: str) -> str:
    s = str(accession).strip().split(";")[0].strip()
    if not s:
        return ""
    return s.split("-")[0].strip().upper()


def _parse_next_url(link_header: str | None) -> str | None:
    if not link_header:
        return None
    m = RE_NEXT.search(link_header)
    return m.group(1) if m else None


def _aliases(primary: str, synonym_cell: str) -> str:
    p = primary.strip().upper()
    parts = [p]
    syn = str(synonym_cell or "").strip()
    if syn:
        for tok in syn.replace(",", " ").split():
            t = tok.strip().upper()
            if t and t not in parts and t != "NAN":
                parts.append(t)
    return ";".join(parts)


def fetch_pages(
    *,
    page_size: int,
    max_rows: int,
    sleep_s: float,
    timeout_s: float,
) -> list[dict[str, str]]:
    """Devuelve filas únicas por (gene_symbol, uniprot_base); prioriza la primera aparición."""
    url = (
        "https://rest.uniprot.org/uniprotkb/search?"
        f"query={urllib.parse.quote(BASE_QUERY, safe='')}"
        f"&fields={FIELDS}&format=tsv&size={int(page_size)}"
    )
    seen: set[tuple[str, str]] = set()
    rows_out: list[dict[str, str]] = []
    n_pages = 0

    while url:
        n_pages += 1
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT}, method="GET")
        try:
            with urllib.request.urlopen(req, timeout=timeout_s) as resp:
                link = resp.headers.get("Link")
                raw = resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            print(f"[uniprot-fetch] HTTP {e.code}: {e.reason}\n{e.read()[:500]!r}", file=sys.stderr)
            raise SystemExit(1) from e
        except urllib.error.URLError as e:
            print(f"[uniprot-fetch] red: {e}", file=sys.stderr)
            raise SystemExit(1) from e

        reader = csv.reader(io.StringIO(raw), delimiter="\t")
        header = next(reader, None)
        if not header or header[0] != "Entry":
            print(f"[uniprot-fetch] cabecera TSV inesperada en página {n_pages}: {header!r}", file=sys.stderr)
            raise SystemExit(1)

        for rec in reader:
            if len(rec) < 2:
                continue
            acc, primary = rec[0].strip(), rec[1].strip() if len(rec) > 1 else ""
            syn_cell = rec[2].strip() if len(rec) > 2 else ""
            if not acc or not primary or primary.upper() == "NAN":
                continue
            gene = primary.strip().upper()
            base = _uniprot_base(acc)
            if not base:
                continue
            key = (gene, base)
            if key in seen:
                continue
            seen.add(key)
            rows_out.append(
                {
                    "protein_id": gene,
                    "gene_symbol": gene,
                    "uniprot_accession": base,
                    "uniprot_base": base,
                    "ensembl_gene_id": "",
                    "alias_symbols": _aliases(primary, syn_cell),
                    "reviewed_status": "reviewed",
                    "mapping_source": "UniProt_REST_human_reviewed_paginated",
                    "is_ambiguous": "false",
                }
            )
            if max_rows > 0 and len(rows_out) >= max_rows:
                url = None
                break

        if max_rows > 0 and len(rows_out) >= max_rows:
            break

        url = _parse_next_url(link)
        if url and sleep_s > 0:
            time.sleep(sleep_s)

    print(f"[uniprot-fetch] páginas={n_pages}, filas únicas (gen×UniProt base)={len(rows_out)}", file=sys.stderr)
    return rows_out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--out",
        type=Path,
        default=Path("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv"),
        help="Ruta del TSV de salida (relativa al CWD o absoluta).",
    )
    p.add_argument("--page-size", type=int, default=500, help="Filas por petición (máx. recomendado ~500).")
    p.add_argument("--max-rows", type=int, default=0, help="0 = sin límite; útil para pruebas.")
    p.add_argument("--sleep", type=float, default=0.15, help="Pausa entre páginas (segundos).")
    p.add_argument("--timeout", type=float, default=120.0, help="Timeout HTTP por petición.")
    args = p.parse_args()

    rows = fetch_pages(
        page_size=max(50, min(args.page_size, 500)),
        max_rows=max(0, args.max_rows),
        sleep_s=max(0.0, args.sleep),
        timeout_s=max(30.0, args.timeout),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=OUT_COLS, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)
    print(f"[uniprot-fetch] escrito: {args.out.resolve()} ({len(rows)} filas)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
