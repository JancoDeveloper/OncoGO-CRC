#!/usr/bin/env python3
"""Reduce **DepMap** (CRISPR gene effect, Chronos) a esencialidad en líneas de CCR → snapshot Capa C.

Los archivos de DepMap son grandes y sus URLs cambian por release, así que este script parte de
**archivos ya descargados** desde https://depmap.org/portal/download/ (release reciente):

- ``CRISPRGeneEffect.csv`` — filas = ModelID, columnas = ``SÍMBOLO (entrez)``, valores = Chronos.
- ``Model.csv``            — metadatos de líneas (incluye ``OncotreeLineage`` / ``OncotreePrimaryDisease``).

Filtra líneas de CCR (``OncotreeLineage == Bowel`` o enfermedad con «colorectal») y calcula por gen:
la media de Chronos y la fracción de líneas esenciales (Chronos < umbral).

Salida (``metadata/evidence/depmap_crc_gene_effect.tsv``):
    ``gene_symbol   chronos_mean   essential_fraction   n_crc_lines``

Uso::

    python scripts/fetch_evidence_depmap.py \\
        --gene-effect-csv ~/Downloads/CRISPRGeneEffect.csv \\
        --model-csv ~/Downloads/Model.csv \\
        --out metadata/evidence/depmap_crc_gene_effect.tsv
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _evidence_common as ec  # noqa: E402

_SYMBOL_RE = re.compile(r"^([A-Za-z0-9\-\.]+)")


def _crc_model_ids(model_csv: Path) -> set[str]:
    m = pd.read_csv(model_csv, low_memory=False)
    id_col = next((c for c in ("ModelID", "DepMap_ID", "depmap_id") if c in m.columns), m.columns[0])
    mask = pd.Series(False, index=m.index)
    for col in ("OncotreeLineage", "OncotreePrimaryDisease", "primary_disease", "lineage"):
        if col in m.columns:
            s = m[col].astype(str).str.lower()
            mask = mask | s.str.contains("bowel") | s.str.contains("colorectal") | s.str.contains("colon")
    ids = set(m.loc[mask, id_col].astype(str))
    return ids


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gene-effect-csv", type=Path, required=True, help="CRISPRGeneEffect.csv de DepMap.")
    p.add_argument("--model-csv", type=Path, required=True, help="Model.csv de DepMap.")
    p.add_argument("--out", type=Path, default=Path("metadata/evidence/depmap_crc_gene_effect.tsv"))
    p.add_argument("--essential-threshold", type=float, default=-0.5, help="Chronos < umbral ⇒ esencial.")
    p.add_argument("--min-lines", type=int, default=3, help="Mín. líneas CCR para reportar un gen.")
    args = p.parse_args()

    if not args.gene_effect_csv.exists() or not args.model_csv.exists():
        ec.log("[depmap] faltan archivos. Descárgalos de https://depmap.org/portal/download/")
        return 1

    crc_ids = _crc_model_ids(args.model_csv)
    ec.log(f"[depmap] líneas CCR detectadas en Model.csv: {len(crc_ids)}")
    if not crc_ids:
        ec.log("[depmap] no se detectaron líneas CCR (revisa columnas OncotreeLineage/Disease).")
        return 1

    ec.log("[depmap] leyendo CRISPRGeneEffect.csv (grande) …")
    ge = pd.read_csv(args.gene_effect_csv, low_memory=False)
    id_col = ge.columns[0]
    ge[id_col] = ge[id_col].astype(str)
    sub = ge[ge[id_col].isin(crc_ids)]
    n_lines = len(sub)
    ec.log(f"[depmap] líneas CCR con datos CRISPR: {n_lines}")
    if n_lines < args.min_lines:
        ec.log(f"[depmap] muy pocas líneas ({n_lines} < {args.min_lines}); abortando.")
        return 1

    gene_cols = [c for c in ge.columns if c != id_col]
    mat = sub[gene_cols].apply(pd.to_numeric, errors="coerce")
    chronos_mean = mat.mean(axis=0)
    essential_fraction = (mat < args.essential_threshold).sum(axis=0) / n_lines

    symbols = [(_SYMBOL_RE.match(c).group(1).strip().upper() if _SYMBOL_RE.match(c) else c.upper()) for c in gene_cols]
    out = pd.DataFrame(
        {
            "gene_symbol": symbols,
            "chronos_mean": np.round(chronos_mean.values, 5),
            "essential_fraction": np.round(essential_fraction.values, 4),
            "n_crc_lines": n_lines,
        }
    )
    out = out.dropna(subset=["chronos_mean"]).sort_values("chronos_mean").drop_duplicates("gene_symbol")
    ec.write_snapshot(out, args.out, source="DepMap CRISPRGeneEffect (Chronos), CCR lines")
    n_ess = int((out["chronos_mean"] < args.essential_threshold).sum())
    ec.log(f"[depmap] genes: {len(out)}; esenciales medios (Chronos<{args.essential_threshold}): {n_ess}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
