"""Evidencia in silico: matriz de literatura / módulos, concordancia con RNA-seq TCGA y tabla Cox opcional."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv, write_tsv


def _resolve_repo_or_workspace(settings: Settings, rel: str) -> Path:
    p_ws = settings.resolve(rel)
    if p_ws.exists():
        return p_ws
    p_repo = settings.root / rel
    if p_repo.exists():
        return p_repo
    return p_ws


def _load_gene_flags_from_layers(layers: dict[str, Any]) -> dict[str, set[str]]:
    """Devuelve nombre_columna_limpio -> conjunto de genes (mayúsculas)."""
    out: dict[str, set[str]] = {}

    et = layers.get("evidence_tiers") or {}
    for tier_key, block in et.items():
        if not isinstance(block, dict):
            continue
        genes = block.get("gene_symbols") or []
        col = f"hit_tier_{tier_key}"
        out[col] = {str(g).strip().upper() for g in genes if str(g).strip()}

    cms = layers.get("cms_pathway_hints") or {}
    for name, block in cms.items():
        if not isinstance(block, dict):
            continue
        genes = block.get("gene_symbols") or []
        safe = "".join(c if c.isalnum() else "_" for c in str(name))
        out[f"hit_cms_{safe}"] = {str(g).strip().upper() for g in genes if str(g).strip()}

    mods = layers.get("functional_modules") or {}
    for name, block in mods.items():
        if not isinstance(block, dict):
            continue
        genes = block.get("gene_symbols") or []
        safe = "".join(c if c.isalnum() else "_" for c in str(name))
        out[f"hit_module_{safe}"] = {str(g).strip().upper() for g in genes if str(g).strip()}

    hb = layers.get("hallmarks_liquid_biopsy_context") or {}
    for name, block in hb.items():
        if not isinstance(block, dict):
            continue
        genes = block.get("gene_symbols") or []
        safe = "".join(c if c.isalnum() else "_" for c in str(name))
        out[f"hit_hallmark_{safe}"] = {str(g).strip().upper() for g in genes if str(g).strip()}

    return out


def _purge_stale_external_evidence_outputs(out_dir: Path) -> None:
    """Borra TSVs de TCGA/Cox previos: si la corrida actual no los regenera, no deben mostrarse."""
    for name in ("tcga_rnaseq_concordance.tsv", "survival_gene_cox_merged.tsv"):
        p = out_dir / name
        if p.exists():
            try:
                p.unlink()
            except OSError:
                pass


def run(settings: Settings) -> None:
    cfg = settings.config.get("external_evidence") or {}
    if not bool(cfg.get("enabled", False)):
        print("[external_evidence] deshabilitado (external_evidence.enabled=false).")
        return
    comp = str(cfg.get("primary_comparison", "cptac_tumor_vs_normal"))
    res = settings.resolve(settings.config["paths"]["results_dir"])
    out_dir = res / "external_evidence"
    out_dir.mkdir(parents=True, exist_ok=True)
    _purge_stale_external_evidence_outputs(out_dir)

    diff_path = res / "differential" / f"{comp}_diff.tsv"
    if not diff_path.exists():
        print(f"[external_evidence] sin {diff_path.name}; omitiendo.")
        return
    diff = read_tsv(diff_path)
    if "protein_id" not in diff.columns:
        print("[external_evidence] diff sin protein_id.")
        return

    ctx_path = settings.resolve(
        str(settings.config.get("translation", {}).get("context_layers_path", ""))
        or "metadata/clinical_layers/crc_context_layers.yaml"
    )
    if not ctx_path.exists():
        ctx_path = settings.root / "metadata/clinical_layers/crc_context_layers.yaml"
    layers: dict[str, Any] = {}
    if ctx_path.exists():
        layers = yaml.safe_load(ctx_path.read_text(encoding="utf-8")) or {}
    gene_flags = _load_gene_flags_from_layers(layers)

    diff = diff.copy()
    diff["_gene"] = diff["protein_id"].astype(str).str.strip().str.upper()
    for col, gset in gene_flags.items():
        diff[col] = diff["_gene"].isin(gset)

    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")
    trans_path = res / "clinical_translation" / f"{primary}_translation_ranked.tsv"
    if bool(cfg.get("merge_translation", True)) and trans_path.exists():
        tr = read_tsv(trans_path)
        key = "feature_id" if "feature_id" in tr.columns else None
        if key:
            tr["_gene"] = tr[key].astype(str).str.strip().str.upper()
            keep = [
                c
                for c in tr.columns
                if c
                in (
                    "translational_priority_score",
                    "evidence_tier",
                    "evidence_score",
                    "combined_score",
                )
            ]
            if keep:
                m = tr[["_gene"] + keep].drop_duplicates("_gene")
                diff = diff.merge(m, on="_gene", how="left")

    diff = diff.drop(columns=["_gene"], errors="ignore")
    write_tsv(diff, out_dir / "evidence_matrix.tsv")
    print(f"[external_evidence] matriz de evidencia -> {settings.relative_to_project(out_dir / 'evidence_matrix.tsv')}")

    summary: dict[str, Any] = {
        "primary_comparison": comp,
        "n_rows_evidence_matrix": int(len(diff)),
        "n_annotation_columns": len(gene_flags),
    }

    # --- TCGA / RNA-seq direccional (tabla opcional) ---
    tcga_rel = str(cfg.get("tcga_rnaseq_log2fc_tsv", "") or "").strip()
    if tcga_rel:
        tcga_p = _resolve_repo_or_workspace(settings, tcga_rel)
        if tcga_p.exists():
            tcga = read_tsv(tcga_p)
            gcol = "gene_symbol" if "gene_symbol" in tcga.columns else tcga.columns[0]
            fcol = None
            for cand in ("log2fc_rnaseq", "log2fc_tcga", "log2fc", "rna_log2fc"):
                if cand in tcga.columns:
                    fcol = cand
                    break
            if fcol:
                tcga["_g"] = tcga[gcol].astype(str).str.strip().str.upper()
                d2 = read_tsv(diff_path)
                d2["_g"] = d2["protein_id"].astype(str).str.strip().str.upper()
                m = d2.merge(
                    tcga[["_g", fcol]].rename(columns={fcol: "log2fc_rnaseq_ref"}),
                    on="_g",
                    how="inner",
                )
                sub = m.dropna(subset=["log2fc", "log2fc_rnaseq_ref"])
                if len(sub) >= 3:
                    r = float(sub["log2fc"].astype(float).corr(sub["log2fc_rnaseq_ref"].astype(float), method="spearman"))
                    same = float(
                        (
                            np.sign(sub["log2fc"].astype(float))
                            == np.sign(sub["log2fc_rnaseq_ref"].astype(float))
                        ).mean()
                    )
                    m["direction_concordant"] = (
                        np.sign(m["log2fc"].astype(float)) == np.sign(m["log2fc_rnaseq_ref"].astype(float))
                    )
                    write_tsv(m.drop(columns=["_g"], errors="ignore"), out_dir / "tcga_rnaseq_concordance.tsv")
                    summary["tcga_rnaseq"] = {
                        "n_matched_genes": int(len(m)),
                        "spearman_r_log2fc": r,
                        "fraction_same_sign": same,
                        "reference_file": str(tcga_p),
                    }
                else:
                    summary["tcga_rnaseq"] = {"status": "insufficient_overlap_or_na", "n_matched": int(len(m))}
            else:
                summary["tcga_rnaseq"] = {"status": "missing_log2fc_column"}
        else:
            summary["tcga_rnaseq"] = {"status": "file_not_found", "path": str(tcga_p)}
    else:
        summary["tcga_rnaseq"] = {"status": "skipped_no_path"}

    # --- Supervivencia (Cox por gen precomputado; p. ej. export GEPIA / cBioPortal) ---
    surv_rel = str(cfg.get("survival_gene_cox_tsv", "") or "").strip()
    if surv_rel:
        sp = _resolve_repo_or_workspace(settings, surv_rel)
        if sp.exists():
            sv = read_tsv(sp)
            gcol = "gene_symbol" if "gene_symbol" in sv.columns else sv.columns[0]
            hr_c = "cox_hr" if "cox_hr" in sv.columns else None
            p_c = "cox_p" if "cox_p" in sv.columns else ("cox_p_value" if "cox_p_value" in sv.columns else None)
            if hr_c and p_c:
                sv["_g"] = sv[gcol].astype(str).str.strip().str.upper()
                d2 = read_tsv(diff_path)
                d2["_g"] = d2["protein_id"].astype(str).str.strip().str.upper()
                ms = d2.merge(sv[["_g", hr_c, p_c]], on="_g", how="inner")
                ms = ms.rename(columns={hr_c: "cox_hr_external", p_c: "cox_p_external"})
                write_tsv(ms.drop(columns=["_g"], errors="ignore"), out_dir / "survival_gene_cox_merged.tsv")
                summary["survival_cox"] = {
                    "n_matched": int(len(ms)),
                    "reference_file": str(sp),
                }
            else:
                summary["survival_cox"] = {"status": "missing_cox_hr_or_p_columns"}
        else:
            summary["survival_cox"] = {"status": "file_not_found", "path": str(sp)}
    else:
        summary["survival_cox"] = {"status": "skipped_no_path"}

    (out_dir / "external_evidence_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    flat_rows: list[dict[str, Any]] = []
    for k, v in summary.items():
        if isinstance(v, dict):
            for k2, v2 in v.items():
                flat_rows.append(
                    {
                        "section": k,
                        "metric": k2,
                        "value": json.dumps(v2) if isinstance(v2, (dict, list)) else v2,
                    }
                )
        else:
            flat_rows.append({"section": "_top", "metric": k, "value": v})
    write_tsv(pd.DataFrame(flat_rows), out_dir / "external_evidence_summary.tsv")
    print(f"[external_evidence] resumen -> {settings.relative_to_project(out_dir / 'external_evidence_summary.json')}")
