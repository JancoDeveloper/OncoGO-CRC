"""Construcción del bundle JSON para ProteinsLanguage / contexto externo."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from urllib.parse import quote

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv
from oncogo_crc.utils.uniprot_lookup import enrich_accessions


def build_protein_context_bundle(
    settings: Settings,
    *,
    comparison: str,
    top: int = 250,
    uniprot_top_n: int = 0,
    uniprot_cache_path: Path | None = None,
) -> dict[str, Any]:
    res = settings.resolve(settings.config["paths"]["results_dir"])
    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")

    mapping = read_tsv(settings.protein_id_map_path())
    map_cols = [c for c in ["protein_id", "gene_symbol", "uniprot_accession"] if c in mapping.columns]

    rows: list[dict[str, Any]] = []
    diff_path = res / "differential" / f"{comparison}_diff.tsv"
    if diff_path.exists():
        diff = read_tsv(diff_path).sort_values("padj", na_position="last").head(top)
        merged = (
            diff.merge(mapping[map_cols], on="protein_id", how="left")
            if "protein_id" in mapping.columns
            else diff
        )
        for _, r in merged.iterrows():
            gn = str(r.get("gene_symbol", "") or "").strip()
            ot_url = (
                f"https://platform.opentargets.org/search?q={quote(gn)}&page=1&entity=target" if gn else None
            )
            rows.append(
                {
                    "source": "differential",
                    "comparison": comparison,
                    "protein_id": str(r.get("protein_id", "")),
                    "gene_symbol": gn if pd.notna(r.get("gene_symbol")) else "",
                    "uniprot_accession": str(r.get("uniprot_accession", ""))
                    if pd.notna(r.get("uniprot_accession"))
                    else "",
                    "log2fc": float(r["log2fc"]) if pd.notna(r.get("log2fc")) else None,
                    "padj": float(r["padj"]) if pd.notna(r.get("padj")) else None,
                    "t_stat": float(r["t_stat"]) if "t_stat" in r and pd.notna(r.get("t_stat")) else None,
                    "druggability_open_targets_query_url": ot_url,
                }
            )

    ml_path = res / "ml" / f"{primary}_feature_ranking.tsv"
    if ml_path.exists():
        ml = read_tsv(ml_path).head(80)
        ml = ml.merge(mapping[map_cols], on="protein_id", how="left") if "protein_id" in ml.columns else ml
        for _, r in ml.iterrows():
            rows.append(
                {
                    "source": "ml_ranking",
                    "comparison": None,
                    "protein_id": str(r.get("protein_id", "")),
                    "gene_symbol": str(r.get("gene_symbol", "")) if pd.notna(r.get("gene_symbol")) else "",
                    "uniprot_accession": str(r.get("uniprot_accession", ""))
                    if pd.notna(r.get("uniprot_accession"))
                    else "",
                    "combined_score": float(r["combined_score"])
                    if "combined_score" in r and pd.notna(r.get("combined_score"))
                    else None,
                }
            )

    uniprot_by_acc: dict[str, Any] = {}
    if uniprot_top_n > 0:
        accs: list[str] = []
        for row in rows:
            u = str(row.get("uniprot_accession") or "").strip()
            if not u or ";" in u:
                continue
            base = u.split("-")[0]
            if base and base not in accs:
                accs.append(base)
            if len(accs) >= uniprot_top_n:
                break
        cache = uniprot_cache_path or settings.resolve("data/interim/uniprot_lookup_cache.json")
        uniprot_by_acc = enrich_accessions(accs, cache_path=cache)

    gsea_path = res / "gsea" / f"{comparison}_gsea_go.tsv"
    gsea_leading: list[str] = []
    if gsea_path.exists():
        try:
            gdf = read_tsv(gsea_path)
            fdr_c = "FDR q-val" if "FDR q-val" in gdf.columns else None
            if fdr_c and "Term" in gdf.columns and "Lead_genes" in gdf.columns:
                sub = gdf[gdf[fdr_c].astype(float) < 0.25].head(5)
                for _, gr in sub.iterrows():
                    gsea_leading.append(f"{gr['Term']}: {str(gr['Lead_genes'])[:120]}")
        except Exception:  # noqa: BLE001
            pass

    overlap_summary = res / "reproducibility" / "cptac_pride_differential_overlap_summary.tsv"
    repro_metrics: dict[str, float | int] = {}
    if overlap_summary.exists():
        odf = read_tsv(overlap_summary)
        if not odf.empty and "metric" in odf.columns and "value" in odf.columns:
            for _, rr in odf.iterrows():
                try:
                    repro_metrics[str(rr["metric"])] = float(rr["value"])
                except (TypeError, ValueError):
                    repro_metrics[str(rr["metric"])] = rr["value"]

    bundle: dict[str, Any] = {
        "project": settings.config.get("project", {}).get("name", "OncoGO-CRC"),
        "primary_dataset": primary,
        "comparison": comparison,
        "n_entries": len(rows),
        "note": (
            "Capa 1 (OncoGO-CRC): estadísticos + ranking ML. "
            "Capa 2 (ProteinsLanguage): GO profundo, STRING, AlphaFold, grafos — usar gene_symbol / UniProt. "
            "URLs: STRING https://string-db.org/ ; UniProt https://www.uniprot.org/ ; "
            "AlphaFold https://alphafold.ebi.ac.uk/ ."
        ),
        "external_links": {
            "string_search_gene": "https://string-db.org/cgi/input?input_identifier={gene}&species=9606",
            "uniprot_lookup": "https://www.uniprot.org/uniprotkb/{accession}",
            "alphafold_search": "https://alphafold.ebi.ac.uk/search/{gene}",
        },
        "druggability_hypothesis": {
            "disclaimer_es": (
                "Solo hipótesis de diana o reposicionamiento in silico; no afirmar respuesta clínica. "
                "Consultar bases curatoriales (Open Targets, ChEMBL, DrugBank) con criterio experto."
            ),
            "open_targets_search_template": "https://platform.opentargets.org/search?q={gene_symbol}&page=1&entity=target",
            "chembl_gene_search_template": "https://www.ebi.ac.uk/chembl/gene_search_results/{gene_symbol}",
        },
        "gsea_leading_edges_preview": gsea_leading,
        "cptac_pride_reproducibility": repro_metrics or None,
        "uniprot_summaries": uniprot_by_acc or None,
        "proteins": rows,
    }
    return bundle


def write_protein_context_bundle(
    settings: Settings,
    out_path: Path,
    *,
    comparison: str,
    top: int = 250,
    uniprot_top_n: int = 0,
) -> Path:
    exp_cfg = settings.config.get("export_bundle", {}) or {}
    top_n = int(exp_cfg.get("uniprot_top_n", uniprot_top_n))
    b = build_protein_context_bundle(
        settings,
        comparison=comparison,
        top=int(exp_cfg.get("diff_top", top)),
        uniprot_top_n=top_n,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(b, indent=2, ensure_ascii=False), encoding="utf-8")
    return out_path
