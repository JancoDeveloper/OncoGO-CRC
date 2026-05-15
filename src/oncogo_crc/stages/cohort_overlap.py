"""Intersección CPTAC ∩ PRIDE en tablas diferenciales (firma proteómica núcleo reproducible).

CPTAC suele usar símbolos de gen como `protein_id`; PRIDE/MaxQuant suele usar accesiones UniProt
(o grupos `P123;Q456`). El cruce literal por `protein_id` suele dar 0 filas: aquí se construye una
clave común **UniProt canónico** (primera accesión, sin isoforma `-2`) usando `protein_mapping_master.tsv`.
"""

from __future__ import annotations

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.ids import effective_protein_mapping
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.reproducibility_overlap import (
    cptac_pride_overlap_tables,
    load_gene_uniprot_from_mapping,
)


def run(settings: Settings) -> None:
    rcfg = settings.config.get("reproducibility", {}) or {}
    block = rcfg.get("cptac_pride_overlap", {}) or {}
    if not bool(block.get("enabled", True)):
        print("[cohort_overlap] deshabilitado en config (reproducibility.cptac_pride_overlap.enabled=false).")
        return
    res = settings.resolve(settings.config["paths"]["results_dir"])
    diff_dir = res / "differential"
    out_dir = res / "reproducibility"
    out_dir.mkdir(parents=True, exist_ok=True)

    ckey = str(block.get("cptac_comparison", "cptac_tumor_vs_normal"))
    pkey = str(block.get("pride_comparison", "pride_validation_tumor_vs_normal"))
    alpha = float(block.get("padj_alpha", settings.config.get("differential", {}).get("alpha", 0.05)))

    cptac_path = diff_dir / f"{ckey}_diff.tsv"
    pride_path = diff_dir / f"{pkey}_diff.tsv"
    if not cptac_path.exists() or not pride_path.exists():
        print(
            f"[cohort_overlap] omitido (falta {cptac_path.name} o {pride_path.name}); "
            "ejecuta antes `differential` para ambas cohortes."
        )
        return

    c = pd.read_csv(cptac_path, sep="\t")
    p = pd.read_csv(pride_path, sep="\t")
    if "protein_id" not in c.columns or "protein_id" not in p.columns:
        print("[cohort_overlap] ambos TSV deben tener columna protein_id.")
        return

    map_path = settings.protein_id_map_path()
    map_df = effective_protein_mapping(settings)
    gene_strict = load_gene_uniprot_from_mapping(map_df, exclude_ambiguous=True)
    gene_full = load_gene_uniprot_from_mapping(map_df, exclude_ambiguous=False)
    if not gene_full:
        print(
            "[cohort_overlap] aviso: mapa vacío o sin columnas gene_symbol/uniprot_accession en "
            f"{settings.relative_to_project(map_path)}; amplía el mapa para cruzar CPTAC (genes) con PRIDE (UniProt)."
        )
    else:
        print(
            f"[cohort_overlap] genes con UniProt (mapa efectivo): estricto 1:1={len(gene_strict)}, lenient={len(gene_full)}"
        )

    merged, summary = cptac_pride_overlap_tables(c, p, gene_to_uni=gene_strict, alpha=alpha)
    summary = summary.copy()
    summary["scope"] = "strict"

    out_merged = out_dir / "cptac_pride_differential_overlap.tsv"
    out_strict = out_dir / "cptac_pride_strict_overlap.tsv"
    sum_path = out_dir / "cptac_pride_differential_overlap_summary.tsv"
    write_tsv(merged, out_merged)
    write_tsv(merged, out_strict)

    summaries = [summary]
    if gene_full != gene_strict:
        merged_l, summary_l = cptac_pride_overlap_tables(c, p, gene_to_uni=gene_full, alpha=alpha)
        summary_l = summary_l.copy()
        summary_l["scope"] = "lenient"
        summaries.append(summary_l)
        write_tsv(merged_l, out_dir / "cptac_pride_differential_overlap_lenient.tsv")

    write_tsv(pd.concat(summaries, ignore_index=True), sum_path)

    sum_s = summaries[0]
    n_bridge = int(sum_s.loc[sum_s["metric"] == "n_intersection_uniprot_bridge", "value"].iloc[0])
    n_core = int(sum_s.loc[sum_s["metric"] == f"n_both_significant_padj_lt_{alpha}", "value"].iloc[0])
    n_raw = int(sum_s.loc[sum_s["metric"] == "n_intersection_raw_protein_id", "value"].iloc[0])
    print(
        f"[cohort_overlap] UniProt-bridge (estricto): {n_bridge} filas comunes ({n_core} doble significativas padj<{alpha}) "
        f"| cruce ingenuo protein_id={n_raw}. -> {out_dir.name}/"
    )
