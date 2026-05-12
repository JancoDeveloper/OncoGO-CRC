"""Enriquecimiento KEGG vía Enrichr (gseapy), opcional — complementa GO."""

from __future__ import annotations

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv, write_tsv


def run(settings: Settings) -> None:
    cfg = settings.config.get("kegg_enrichment", {}) or {}
    if not bool(cfg.get("enabled", False)):
        print("[kegg] deshabilitado (kegg_enrichment.enabled=false en config.yaml).")
        return
    try:
        import gseapy as gp  # type: ignore[import-untyped]
    except ImportError as e:
        raise ImportError(
            "KEGG requiere gseapy. Instala: pip install -e \".[go-remote]\""
        ) from e

    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")
    alpha = float(settings.config.get("differential", {}).get("alpha", 0.05))
    comp = str(cfg.get("comparison", "cptac_tumor_vs_normal"))
    diff_path = settings.resolve(settings.config["paths"]["results_dir"]) / "differential" / f"{comp}_diff.tsv"
    if not diff_path.exists():
        print(f"[kegg] sin diferencial {diff_path.name}; omitiendo.")
        return
    diff = read_tsv(diff_path)
    if "padj" not in diff.columns or "protein_id" not in diff.columns:
        print("[kegg] diferencial sin columnas esperadas; omitiendo.")
        return
    sig = diff.loc[diff["padj"].astype(float) < alpha, "protein_id"].astype(str).str.upper().unique().tolist()
    if len(sig) < int(cfg.get("min_genes", 10)):
        print(f"[kegg] pocos genes significativos ({len(sig)}); omitiendo.")
        return
    lib = str(cfg.get("enrichr_library", "KEGG_2021_Human"))
    out_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "kegg"
    out_dir.mkdir(parents=True, exist_ok=True)
    enr = gp.enrichr(gene_list=sig, gene_sets=lib, organism="human", outdir=None, no_plot=True)
    res = pd.DataFrame(enr.results)
    out_tsv = out_dir / f"{primary}_kegg_enrichr.tsv"
    write_tsv(res, out_tsv)
    print(f"[kegg] Enrichr {lib} -> {out_tsv}")
