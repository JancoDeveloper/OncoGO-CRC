#!/usr/bin/env python3
"""Calcula Cox HR (OS) por gen para TCGA COAD+READ usando UCSC Xena.

Estrategia (consistente con literatura de meta-firmas, p. ej. Gentles et al. 2015):
    1. Para cada gen, dicotomizar muestras tumor en alto vs bajo (mediana de log2 RSEM+1).
    2. Ajustar modelo Cox PH univariado con ``lifelines`` sobre OS = (days_to_event, event).
    3. Reportar HR (hazard ratio para el grupo "alto"), IC95 y p-valor del Wald test.

Salida (TSV): ``gene_symbol  cox_hr  cox_log_hr  cox_ci_low  cox_ci_high  cox_p  cox_n  cox_events``

Notas:
    - Fuente expresión: ``https://tcga.xenahubs.net/download/TCGA.COADREAD.sampleMap/HiSeqV2.gz``
    - Fuente supervivencia: ``https://tcga.xenahubs.net/download/survival/COADREAD_survival.txt.gz``
      (columnas ``sample``, ``OS``, ``OS.time``)
    - Tiempo de cómputo: ~10-15 min para todos los genes (~20 k) en una CPU moderna.
    - ``--max-genes`` y ``--only-genes-tsv`` permiten cómputo parcial (recomendado para
      restringir a los ~6 000 genes de la matriz proteómica curada).
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

USER_AGENT = "OncoGO-CRC/xena-fetch (https://github.com; academic research)"
EXPR_URL = "https://tcga.xenahubs.net/download/TCGA.COADREAD.sampleMap/HiSeqV2.gz"
SURV_URL = "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/survival%2FCOADREAD_survival.txt"


def _download_to(url: str, dst: Path, *, timeout_s: float = 300.0) -> None:
    if dst.exists() and dst.stat().st_size > 1024:
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    print(f"[xena-fetch] descargando {url} → {dst.name} …", file=sys.stderr)
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp, dst.open("wb") as out:
            shutil.copyfileobj(resp, out, length=1 << 16)
    except urllib.error.URLError as e:
        print(f"[xena-fetch] error: {e}", file=sys.stderr)
        raise SystemExit(1) from e


def _maybe_gunzip(src: Path, dst: Path) -> None:
    if dst.exists() and dst.stat().st_size > 1024:
        return
    with gzip.open(src, "rb") as gf, dst.open("wb") as out:
        shutil.copyfileobj(gf, out, length=1 << 16)


def _load_only_genes(p: Path | None) -> set[str] | None:
    if not p:
        return None
    df = pd.read_csv(p, sep="\t")
    col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
    return {str(x).strip().upper() for x in df[col] if str(x).strip()}


def _cox_per_gene(
    expr_row: np.ndarray,
    durations: np.ndarray,
    events: np.ndarray,
) -> tuple[float, float, float, float, float, int, int]:
    """Cox univariado dicotomizado por mediana. Devuelve hr, log_hr, ci_low, ci_high, p, n, events."""
    from lifelines import CoxPHFitter

    mask = (~np.isnan(expr_row)) & (~np.isnan(durations)) & (~np.isnan(events))
    if mask.sum() < 20:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, int(mask.sum()), int(events[mask].sum()))
    x = expr_row[mask]
    med = np.nanmedian(x)
    high = (x > med).astype(int)
    if high.sum() < 5 or (1 - high).sum() < 5:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, int(mask.sum()), int(events[mask].sum()))

    d = pd.DataFrame(
        {
            "duration": durations[mask],
            "event": events[mask].astype(int),
            "high": high,
        }
    )
    try:
        cph = CoxPHFitter(penalizer=0.0)
        cph.fit(d, duration_col="duration", event_col="event", show_progress=False)
        row = cph.summary.loc["high"]
        hr = float(row["exp(coef)"])
        log_hr = float(row["coef"])
        ci_low = float(row["exp(coef) lower 95%"])
        ci_high = float(row["exp(coef) upper 95%"])
        p = float(row["p"])
        return (hr, log_hr, ci_low, ci_high, p, int(mask.sum()), int(d["event"].sum()))
    except Exception:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, int(mask.sum()), int(events[mask].sum()))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--out",
        type=Path,
        default=Path("metadata/tcga_auxiliary/tcga_coadread_cox_os.tsv"),
        help="TSV de salida.",
    )
    p.add_argument(
        "--cache",
        type=Path,
        default=Path("metadata/tcga_auxiliary/.xena_cache"),
        help="Cache de matrices descargadas.",
    )
    p.add_argument("--max-genes", type=int, default=0, help="Limita a N genes para pruebas (0 = todos).")
    p.add_argument(
        "--only-genes-tsv",
        type=Path,
        default=None,
        help="TSV con columna gene_symbol; restringir Cox a esos genes (recomendado: matriz proteómica curada).",
    )
    p.add_argument(
        "--min-events",
        type=int,
        default=10,
        help="Mínimo de eventos (muertes) para considerar el gen.",
    )
    args = p.parse_args()

    cache = args.cache
    cache.mkdir(parents=True, exist_ok=True)

    expr_gz = cache / "HiSeqV2.gz"
    expr_tsv = cache / "HiSeqV2.tsv"
    surv_tsv = cache / "COADREAD_survival.txt"

    _download_to(EXPR_URL, expr_gz)
    _maybe_gunzip(expr_gz, expr_tsv)
    _download_to(SURV_URL, surv_tsv)

    print(f"[xena-fetch] leyendo supervivencia {surv_tsv.name} …", file=sys.stderr)
    surv = pd.read_csv(surv_tsv, sep="\t", low_memory=False)
    cols = {c.lower(): c for c in surv.columns}
    sample_c = cols.get("sample") or cols.get("sampleid") or list(surv.columns)[0]
    os_event_c = cols.get("os") or cols.get("os_status")
    os_time_c = cols.get("os.time") or cols.get("os_time") or cols.get("os_days")
    if not os_event_c or not os_time_c:
        print(f"[xena-fetch] no encuentro columnas OS; columnas: {surv.columns.tolist()}", file=sys.stderr)
        return 1
    surv = surv[[sample_c, os_event_c, os_time_c]].rename(
        columns={sample_c: "sample_id", os_event_c: "OS", os_time_c: "OS_time"}
    )
    surv["sample_id"] = surv["sample_id"].astype(str).str.strip()
    surv["OS"] = pd.to_numeric(surv["OS"], errors="coerce")
    surv["OS_time"] = pd.to_numeric(surv["OS_time"], errors="coerce")
    surv = surv.dropna(subset=["OS", "OS_time"])
    print(f"[xena-fetch] supervivencia: {len(surv)} muestras con OS válido.", file=sys.stderr)

    print(f"[xena-fetch] leyendo expresión {expr_tsv.name} …", file=sys.stderr)
    expr = pd.read_csv(expr_tsv, sep="\t", low_memory=False)
    gene_col = expr.columns[0]
    expr = expr.rename(columns={gene_col: "gene_symbol"})
    expr["gene_symbol"] = expr["gene_symbol"].astype(str).str.strip().str.upper()
    # Alinear muestras a las que tienen OS
    keep = ["gene_symbol"] + [c for c in expr.columns[1:] if c in set(surv["sample_id"])]
    expr = expr[keep].copy()
    print(f"[xena-fetch] matriz alineada: {expr.shape[1] - 1} muestras × {len(expr)} genes.", file=sys.stderr)

    only_genes = _load_only_genes(args.only_genes_tsv)
    if only_genes:
        before = len(expr)
        expr = expr[expr["gene_symbol"].isin(only_genes)].copy()
        print(
            f"[xena-fetch] restringido a {len(expr)} genes (de {before}) por --only-genes-tsv.",
            file=sys.stderr,
        )

    if args.max_genes > 0:
        expr = expr.head(args.max_genes).copy()
        print(f"[xena-fetch] (modo prueba) limitando a {len(expr)} genes.", file=sys.stderr)

    # Reordenar matriz a las muestras de survival
    sample_cols = [c for c in expr.columns[1:] if c in set(surv["sample_id"])]
    surv2 = surv.set_index("sample_id").loc[sample_cols].reset_index()
    durations = surv2["OS_time"].to_numpy(dtype=float)
    events = surv2["OS"].to_numpy(dtype=float)

    if events.sum() < args.min_events:
        print(
            f"[xena-fetch] eventos insuficientes ({int(events.sum())}<{args.min_events}); abortando.",
            file=sys.stderr,
        )
        return 1

    print(
        f"[xena-fetch] ajustando Cox por gen sobre {len(sample_cols)} muestras y {int(events.sum())} eventos …",
        file=sys.stderr,
    )

    expr_mat = expr[sample_cols].to_numpy(dtype=float)
    results: list[dict[str, float | str | int]] = []
    total = len(expr)
    log_every = max(1, total // 20)
    for i, gene in enumerate(expr["gene_symbol"].values):
        hr, log_hr, lo, hi, pv, n, ev = _cox_per_gene(expr_mat[i], durations, events)
        results.append(
            {
                "gene_symbol": gene,
                "cox_hr": hr,
                "cox_log_hr": log_hr,
                "cox_ci_low": lo,
                "cox_ci_high": hi,
                "cox_p": pv,
                "cox_n": n,
                "cox_events": ev,
            }
        )
        if (i + 1) % log_every == 0:
            print(f"[xena-fetch]   procesados {i + 1}/{total}", file=sys.stderr)

    out = pd.DataFrame(results).dropna(subset=["cox_hr"])
    out = out.sort_values("cox_p")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    n_sig = int((out["cox_p"] < 0.05).sum())
    print(
        f"[xena-fetch] escrito: {args.out.resolve()} ({len(out)} genes; p<0.05: {n_sig})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
