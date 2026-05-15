"""Consultas mínimas a UniProt REST con caché local (sin dependencia requests)."""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any


def _parse_entry(data: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {"accession": data.get("primaryAccession", "")}
    desc = data.get("proteinDescription") or {}
    rec = desc.get("recommendedName") or {}
    full = rec.get("fullName") or {}
    out["protein_name"] = str(full.get("value", ""))[:500]
    genes = data.get("genes") or []
    if genes and isinstance(genes[0], dict):
        gn = genes[0].get("geneName") or {}
        out["gene_primary"] = str(gn.get("value", ""))
    else:
        out["gene_primary"] = ""
    out["organism"] = str((data.get("organism") or {}).get("scientificName", ""))
    return out


def enrich_accessions(accessions: list[str], *, cache_path: Path, delay_s: float = 0.35) -> dict[str, Any]:
    cache: dict[str, Any] = {}
    if cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            cache = {}
    out: dict[str, Any] = {}
    for raw in accessions:
        acc = str(raw).strip().split("-")[0]
        if not acc or len(acc) < 4:
            continue
        if acc in cache and isinstance(cache[acc], dict) and cache[acc].get("protein_name"):
            out[acc] = cache[acc]
            continue
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=45) as resp:  # noqa: S310
                data = json.loads(resp.read().decode("utf-8"))
            parsed = _parse_entry(data)
            cache[acc] = parsed
            out[acc] = parsed
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, KeyError) as e:
            cache[acc] = {"error": str(e)[:200]}
            out[acc] = cache[acc]
        time.sleep(delay_s)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
    return out
