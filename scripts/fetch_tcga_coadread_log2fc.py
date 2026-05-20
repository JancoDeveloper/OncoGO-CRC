#!/usr/bin/env python3
"""Descarga TCGA COAD+READ (RNA-seq) desde UCSC Xena y calcula log2FC tumor vs normal.

Fuentes:
    - Expresión: ``https://tcga.xenahubs.net/download/TCGA.COADREAD.sampleMap/HiSeqV2.gz``
      (HiSeq V2 log2(RSEM+1), ~20 530 genes × 434 muestras, ~150 MB descomprimido)
    - Phenotype: ``https://tcga.xenahubs.net/download/TCGA.COADREAD.sampleMap/COADREAD_clinicalMatrix``
      (incluye columna ``sample_type`` con valores ``Primary Tumor`` / ``Solid Tissue Normal``)

Salida (TSV, columnas):
    ``gene_symbol  log2fc_rnaseq  n_tumor  n_normal  mean_tumor  mean_normal  t_stat  p_value``

El TSV es directamente consumible por ``external_evidence.tcga_rnaseq_log2fc_tsv``
(`config.yaml`). Una sola corrida (descarga + cómputo + escritura) toma ~3-5 min y ~150 MB de disco.

Uso típico::

    python scripts/fetch_tcga_coadread_log2fc.py \\
        --out metadata/tcga_auxiliary/tcga_coadread_log2fc.tsv \\
        --cache metadata/tcga_auxiliary/.xena_cache
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
from scipy import stats


USER_AGENT = "OncoGO-CRC/xena-fetch (https://github.com; academic research)"
EXPR_URL = "https://tcga.xenahubs.net/download/TCGA.COADREAD.sampleMap/HiSeqV2.gz"
PHENO_URL = "https://tcga.xenahubs.net/download/TCGA.COADREAD.sampleMap/COADREAD_clinicalMatrix"


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


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--out",
        type=Path,
        default=Path("metadata/tcga_auxiliary/tcga_coadread_log2fc.tsv"),
        help="TSV de salida (gene_symbol, log2fc_rnaseq, …).",
    )
    p.add_argument(
        "--cache",
        type=Path,
        default=Path("metadata/tcga_auxiliary/.xena_cache"),
        help="Directorio para cache de matrices descargadas.",
    )
    p.add_argument(
        "--min-samples-per-group",
        type=int,
        default=5,
        help="Mínimo de muestras tumor y normal por gen (default 5).",
    )
    p.add_argument(
        "--max-genes",
        type=int,
        default=0,
        help="Para pruebas: limita a N primeros genes (0 = todos).",
    )
    args = p.parse_args()

    cache = args.cache
    cache.mkdir(parents=True, exist_ok=True)

    expr_gz = cache / "HiSeqV2.gz"
    expr_tsv = cache / "HiSeqV2.tsv"
    pheno_tsv = cache / "COADREAD_clinicalMatrix.tsv"

    _download_to(EXPR_URL, expr_gz)
    _maybe_gunzip(expr_gz, expr_tsv)
    _download_to(PHENO_URL, pheno_tsv)

    print(f"[xena-fetch] leyendo phenotype {pheno_tsv.name} …", file=sys.stderr)
    ph = pd.read_csv(pheno_tsv, sep="\t", low_memory=False, dtype=str)
    sid_col = "sampleID" if "sampleID" in ph.columns else ph.columns[0]
    type_col = next((c for c in ("sample_type", "_sample_type") if c in ph.columns), None)
    if not type_col:
        print(f"[xena-fetch] no encuentro columna sample_type; columnas: {ph.columns.tolist()[:20]}", file=sys.stderr)
        return 1
    ph = ph[[sid_col, type_col]].rename(columns={sid_col: "sample_id", type_col: "sample_type"})
    ph["sample_id"] = ph["sample_id"].astype(str).str.strip()
    ph["sample_type"] = ph["sample_type"].astype(str).str.strip()

    tumor_set = set(ph.loc[ph["sample_type"].isin(["Primary Tumor", "Metastatic"]), "sample_id"])
    normal_set = set(ph.loc[ph["sample_type"].isin(["Solid Tissue Normal", "Normal Tissue"]), "sample_id"])
    print(
        f"[xena-fetch] phenotype: {len(tumor_set)} tumor + {len(normal_set)} normal (de {len(ph)} totales).",
        file=sys.stderr,
    )
    if not tumor_set or not normal_set:
        print("[xena-fetch] sin tumores o normales suficientes; abortando.", file=sys.stderr)
        return 1

    print(f"[xena-fetch] leyendo expresión {expr_tsv.name} (puede tardar)…", file=sys.stderr)
    # La matriz Xena es ~20 530 × 434; con bajo RAM, leer en chunks por filas.
    # Aquí usamos read_csv directo (cabe en ~120 MB de RAM).
    expr = pd.read_csv(expr_tsv, sep="\t", low_memory=False)
    gene_col = expr.columns[0]
    sample_cols = [c for c in expr.columns[1:] if c in tumor_set or c in normal_set]
    expr = expr[[gene_col] + sample_cols].copy()
    expr = expr.rename(columns={gene_col: "gene_symbol"})
    expr["gene_symbol"] = expr["gene_symbol"].astype(str).str.strip().str.upper()
    print(f"[xena-fetch] matriz reducida a {expr.shape[1] - 1} muestras con phenotype conocido.", file=sys.stderr)

    if args.max_genes > 0:
        expr = expr.head(args.max_genes).copy()
        print(f"[xena-fetch] (modo prueba) limitando a {len(expr)} genes.", file=sys.stderr)

    tumor_cols = [c for c in expr.columns[1:] if c in tumor_set]
    normal_cols = [c for c in expr.columns[1:] if c in normal_set]
    if len(tumor_cols) < args.min_samples_per_group or len(normal_cols) < args.min_samples_per_group:
        print(
            f"[xena-fetch] muestras insuficientes: tumor={len(tumor_cols)}, normal={len(normal_cols)}; "
            f"min={args.min_samples_per_group}.",
            file=sys.stderr,
        )
        return 1

    tum_mat = expr[tumor_cols].astype(float).values
    nor_mat = expr[normal_cols].astype(float).values

    print("[xena-fetch] calculando log2FC y t-test por gen …", file=sys.stderr)
    mean_t = np.nanmean(tum_mat, axis=1)
    mean_n = np.nanmean(nor_mat, axis=1)
    n_t = (~np.isnan(tum_mat)).sum(axis=1)
    n_n = (~np.isnan(nor_mat)).sum(axis=1)
    # Datos Xena HiSeqV2 ya están en log2(RSEM+1) → la diferencia de medias YA es log2FC.
    log2fc = mean_t - mean_n
    # Welch t-test ignorando NaN (vectorizado por fila).
    t_stat, p_value = stats.ttest_ind(
        tum_mat, nor_mat, axis=1, equal_var=False, nan_policy="omit"
    )

    out = pd.DataFrame(
        {
            "gene_symbol": expr["gene_symbol"].values,
            "log2fc_rnaseq": log2fc,
            "n_tumor": n_t,
            "n_normal": n_n,
            "mean_tumor": mean_t,
            "mean_normal": mean_n,
            "t_stat": t_stat,
            "p_value": p_value,
        }
    )
    out = out.dropna(subset=["log2fc_rnaseq"])
    # Drop duplicate gene symbols keeping the one with largest |t_stat|
    out["_abs_t"] = out["t_stat"].abs().fillna(0.0)
    out = out.sort_values("_abs_t", ascending=False).drop_duplicates("gene_symbol").drop(columns="_abs_t")
    out = out.sort_values("gene_symbol")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    n_up = int((out["log2fc_rnaseq"] > 1).sum())
    n_dn = int((out["log2fc_rnaseq"] < -1).sum())
    print(
        f"[xena-fetch] escrito: {args.out.resolve()} "
        f"({len(out)} genes; UP |log2FC|>1: {n_up}; DOWN |log2FC|<-1: {n_dn})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
