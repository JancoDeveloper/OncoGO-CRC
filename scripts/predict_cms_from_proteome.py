#!/usr/bin/env python3
"""Asignación aproximada de CMS (Consensus Molecular Subtypes) desde matriz proteómica CPTAC.

⚠ Importante (limitación honesta):
    Esta NO es la clasificación oficial CMScaller (Guinney et al. 2015), que requiere
    expresión transcriptómica completa y centroides entrenados. Esta es una **aproximación
    de proteoma** que usa firmas canónicas curadas de literatura. Para tesis, declarar
    explícitamente como "CMS inferido por firma simplificada proteómica" y NO como
    "clasificación CMScaller".

    El estándar de oro para reemplazar este TSV es:
      Vasaikar et al. 2019, Cell 177(4):1035-1049 — Supplementary Table S1
      DOI: 10.1016/j.cell.2019.03.030
      Contiene los CMS oficiales para los 110 participantes de CPTAC2 prospective colon (PDC000116).

Firmas (curadas a partir de Guinney et al. 2015, Vasaikar et al. 2019, y revisiones CMS):
    - CMS1 (MSI / Immune):     CIITA, IDO1, CXCL9, CXCL10, CD8A, HLA-DRA, GBP1, GZMB, PRF1, IRF1
    - CMS2 (Canonical WNT/MYC): MYC, AXIN2, CDX2, REG4, LGR5, ASCL2, EREG, TFF3, MUC2
    - CMS3 (Metabolic):        IGF2, NRG1, IGFBP3, MGAM, ALDOB, GCK, SLC2A2, PCK1
    - CMS4 (Mesenchymal/EMT):  ZEB1, ZEB2, TGFB1, COL1A1, COL3A1, FAP, THBS1, SPARC, VIM, ACTA2

Salida (TSV): ``participant_id  cms_subtype  cms1_score  cms2_score  cms3_score  cms4_score  cms_score_max  confidence``

Uso::

    python scripts/predict_cms_from_proteome.py \\
        --matrix data/processed/cptac_colorectal_log2norm.tsv \\
        --metadata data/processed/cptac_colorectal_metadata.tsv \\
        --out metadata/clinical_layers/cptac_cms_predicted.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SIGNATURES: dict[str, list[str]] = {
    "CMS1": [
        "CIITA", "IDO1", "CXCL9", "CXCL10", "CD8A",
        "HLA-DRA", "HLA-DRB1", "GBP1", "GZMB", "PRF1", "IRF1",
        "CXCL11", "STAT1", "PSMB8", "PSMB9", "TAP1", "TAP2",
    ],
    "CMS2": [
        "MYC", "AXIN2", "CDX2", "REG4", "LGR5", "ASCL2",
        "EREG", "TFF3", "MUC2", "VAV3", "EPHB2", "PROM1",
    ],
    "CMS3": [
        "IGF2", "NRG1", "IGFBP3", "MGAM", "ALDOB", "GCK",
        "SLC2A2", "PCK1", "MUC1", "FUT3", "FUT2",
    ],
    "CMS4": [
        "ZEB1", "ZEB2", "TGFB1", "COL1A1", "COL3A1", "FAP",
        "THBS1", "SPARC", "VIM", "ACTA2", "POSTN", "INHBA",
        "SNAI2", "FN1", "TWIST1", "SERPINH1",
    ],
}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--matrix", type=Path, required=True, help="TSV log2norm (protein_id + columnas muestra).")
    p.add_argument("--metadata", type=Path, required=True, help="TSV de metadata (sample_id, participant_id, condition).")
    p.add_argument(
        "--out",
        type=Path,
        default=Path("metadata/clinical_layers/cptac_cms_predicted.tsv"),
        help="TSV de salida con asignación CMS por participante.",
    )
    p.add_argument(
        "--min-signature-coverage",
        type=int,
        default=4,
        help="Mínimo de genes detectados por firma para que sea considerada (default 4).",
    )
    p.add_argument(
        "--confidence-margin",
        type=float,
        default=0.3,
        help="Margen mínimo entre z-score top y segundo para alta confianza (default 0.3).",
    )
    args = p.parse_args()

    if not args.matrix.is_file() or not args.metadata.is_file():
        print(f"[predict_cms] requiere {args.matrix} y {args.metadata}", file=sys.stderr)
        return 1

    print(f"[predict_cms] leyendo {args.matrix.name} …", file=sys.stderr)
    mx = pd.read_csv(args.matrix, sep="\t", low_memory=False)
    # Identificar columna gene_symbol; si no existe, usar protein_id (no ideal pero servirá si harmonize ya lo resolvió).
    gene_col = "gene_symbol" if "gene_symbol" in mx.columns else "protein_id"
    mx[gene_col] = mx[gene_col].astype(str).str.strip().str.upper()
    sample_cols = [c for c in mx.columns if c not in {"protein_id", "gene_symbol", "uniprot_accession"}]
    if not sample_cols:
        print("[predict_cms] no encuentro columnas de muestra.", file=sys.stderr)
        return 1

    md = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    # Restringir a muestras tumor con participant_id
    if "condition" in md.columns:
        md_tum = md[md["condition"].astype(str).str.lower() == "tumor"].copy()
    else:
        md_tum = md.copy()
    if "participant_id" not in md_tum.columns or "sample_id" not in md_tum.columns:
        print("[predict_cms] metadata requiere sample_id y participant_id.", file=sys.stderr)
        return 1

    # Mantener primera muestra tumor por participante
    md_tum = md_tum.drop_duplicates("participant_id", keep="first")
    tumor_samples = [s for s in md_tum["sample_id"].astype(str) if s in sample_cols]
    if len(tumor_samples) < 5:
        print(f"[predict_cms] solo {len(tumor_samples)} muestras tumor cruzadas; insuficiente.", file=sys.stderr)
        return 1
    print(f"[predict_cms] {len(tumor_samples)} muestras tumor con participant_id.", file=sys.stderr)

    # Construir subset matriz [genes_de_firmas × muestras_tumor], z-score por gen sobre muestras.
    all_sig_genes = sorted({g for sig in SIGNATURES.values() for g in sig})
    sub = mx[mx[gene_col].isin(all_sig_genes)][[gene_col] + tumor_samples].copy()
    sub = sub.drop_duplicates(gene_col, keep="first").set_index(gene_col)
    print(f"[predict_cms] genes de firma detectados en matriz: {len(sub)} de {len(all_sig_genes)}.", file=sys.stderr)

    # Diagnóstico de cobertura por firma
    for sig_name, sig_genes in SIGNATURES.items():
        found = [g for g in sig_genes if g in sub.index]
        print(f"  {sig_name}: {len(found)}/{len(sig_genes)} genes — {found[:5]}…", file=sys.stderr)

    # z-score por gen (across tumor samples)
    arr = sub.to_numpy(dtype=float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, ddof=1, keepdims=True)
    sd = np.where(sd < 1e-12, np.nan, sd)
    z = (arr - mu) / sd
    z_df = pd.DataFrame(z, index=sub.index, columns=sub.columns)

    # Score por firma = media de z-scores de genes presentes
    scores: dict[str, pd.Series] = {}
    coverage: dict[str, int] = {}
    for sig_name, sig_genes in SIGNATURES.items():
        present = [g for g in sig_genes if g in z_df.index]
        coverage[sig_name] = len(present)
        if len(present) < args.min_signature_coverage:
            print(
                f"[predict_cms] firma {sig_name} con cobertura {len(present)} (< {args.min_signature_coverage}); marcaré NA.",
                file=sys.stderr,
            )
            scores[sig_name] = pd.Series(np.nan, index=tumor_samples)
        else:
            scores[sig_name] = z_df.loc[present].mean(axis=0)

    score_df = pd.DataFrame(scores).reindex(tumor_samples)
    # Asignar CMS por argmax
    sample_to_part = dict(zip(md_tum["sample_id"].astype(str), md_tum["participant_id"].astype(str)))

    out_rows: list[dict] = []
    for s in tumor_samples:
        row = score_df.loc[s]
        if row.isna().all():
            cms = ""
            top = float("nan")
            conf = "low"
        else:
            cms = str(row.idxmax())
            top = float(row.max())
            sorted_vals = row.sort_values(ascending=False).dropna().to_numpy()
            margin = float(sorted_vals[0] - sorted_vals[1]) if len(sorted_vals) >= 2 else 0.0
            conf = "high" if margin >= args.confidence_margin and top >= 0.0 else "low"
        out_rows.append(
            {
                "sample_id": s,
                "participant_id": sample_to_part.get(s, ""),
                "cms_subtype": cms,
                "cms1_score": float(row.get("CMS1", np.nan)),
                "cms2_score": float(row.get("CMS2", np.nan)),
                "cms3_score": float(row.get("CMS3", np.nan)),
                "cms4_score": float(row.get("CMS4", np.nan)),
                "cms_score_max": top,
                "confidence": conf,
                "source": "proteome_signature_simplified_pipeline",
            }
        )
    out = pd.DataFrame(out_rows)
    out = out.drop_duplicates("participant_id", keep="first")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    dist = out["cms_subtype"].value_counts(dropna=False).to_dict()
    n_high = int((out["confidence"] == "high").sum())
    print(f"[predict_cms] escrito: {args.out.resolve()} ({len(out)} participantes)", file=sys.stderr)
    print(f"[predict_cms] distribución CMS: {dist}", file=sys.stderr)
    print(f"[predict_cms] alta confianza: {n_high}/{len(out)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
