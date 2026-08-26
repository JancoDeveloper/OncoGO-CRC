"""Utilidades compartidas por los descargadores de evidencia (Capa C).

- Carga del **universo de genes** a consultar (desde la firma diferencial, un TSV de genes,
  o el mapa maestro proteína→gen).
- Cliente HTTP mínimo (POST/GET JSON con reintentos) sin dependencias externas.
- Escritura del snapshot TSV + sidecar ``.version`` (trazabilidad FAIR).

Cada script de ``fetch_evidence_*`` produce un TSV en el formato que espera
``src/oncogo_crc/utils/evidence_providers.py``.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

USER_AGENT = "OncoGO-CRC/evidence-fetch (academic research)"


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------- HTTP
def http_json(
    url: str,
    *,
    payload: dict | list | None = None,
    method: str | None = None,
    timeout_s: float = 60.0,
    retries: int = 3,
    backoff_s: float = 2.0,
    extra_headers: dict[str, str] | None = None,
):
    """POST (si hay payload) o GET JSON con reintentos. Devuelve el objeto decodificado."""
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    data = None
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    if extra_headers:
        headers.update(extra_headers)
    m = method or ("POST" if data is not None else "GET")
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, data=data, headers=headers, method=m)
            with urllib.request.urlopen(req, timeout=timeout_s) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as e:
            last_err = e
            if attempt < retries:
                wait = backoff_s * attempt
                log(f"[http] intento {attempt}/{retries} falló ({e}); reintento en {wait:.0f}s…")
                time.sleep(wait)
    raise SystemExit(f"[http] fallo tras {retries} intentos: {last_err}")


def download_file(url: str, dst: Path, *, timeout_s: float = 600.0, min_bytes: int = 1024) -> Path:
    """Descarga ``url`` a ``dst`` (con caché por tamaño). Devuelve la ruta destino."""
    import shutil

    if dst.exists() and dst.stat().st_size >= min_bytes:
        log(f"[download] caché: {dst.name} ya existe.")
        return dst
    dst.parent.mkdir(parents=True, exist_ok=True)
    log(f"[download] {url} → {dst.name} …")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout_s) as resp, dst.open("wb") as out:
        shutil.copyfileobj(resp, out, length=1 << 16)
    return dst


# --------------------------------------------------------------------- gene universe
def load_gene_universe(
    *,
    genes_arg: list[str] | None = None,
    genes_tsv: Path | None = None,
    differential_tsv: Path | None = None,
    mapping_tsv: Path | None = None,
    limit: int = 0,
) -> list[str]:
    """Devuelve una lista de símbolos de gen (MAYÚSCULAS, únicos) según la fuente disponible.

    Prioridad: ``genes_arg`` > ``genes_tsv`` > ``differential_tsv`` > ``mapping_tsv``.
    """
    genes: list[str] = []
    if genes_arg:
        genes = [g.strip().upper() for g in genes_arg if g.strip()]
    elif genes_tsv and Path(genes_tsv).exists():
        df = pd.read_csv(genes_tsv, sep="\t")
        col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
        genes = df[col].astype(str).str.strip().str.upper().tolist()
    elif differential_tsv and Path(differential_tsv).exists():
        df = pd.read_csv(differential_tsv, sep="\t")
        col = "protein_id" if "protein_id" in df.columns else df.columns[0]
        genes = df[col].astype(str).str.strip().str.upper().tolist()
    elif mapping_tsv and Path(mapping_tsv).exists():
        df = pd.read_csv(mapping_tsv, sep="\t")
        col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
        genes = df[col].astype(str).str.strip().str.upper().tolist()

    seen: set[str] = set()
    out: list[str] = []
    for g in genes:
        if g and g != "NAN" and g not in seen:
            seen.add(g)
            out.append(g)
    if limit and limit > 0:
        out = out[:limit]
    return out


def chunked(seq: list, size: int):
    for i in range(0, len(seq), size):
        yield seq[i : i + size]


# ------------------------------------------------------------------------ snapshot
def write_snapshot(df: pd.DataFrame, out_path: Path, *, source: str, version: str | None = None) -> None:
    """Escribe el TSV del snapshot y un sidecar ``<out>.version`` con la procedencia."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    stamp = version or datetime.now(timezone.utc).strftime("%Y-%m-%d")
    sidecar = out_path.with_suffix(out_path.suffix + ".version")
    sidecar.write_text(
        f"{source} | {stamp} | n_genes={len(df)} | fetched={datetime.now(timezone.utc).isoformat()}\n",
        encoding="utf-8",
    )
    log(f"[snapshot] {out_path} ({len(df)} genes) + {sidecar.name}")
