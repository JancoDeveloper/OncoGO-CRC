from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


def _load_context_layers(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.is_file():
        return None
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    return data if isinstance(data, dict) else None


def _annotate_context_layers(base: pd.DataFrame, layers: dict[str, Any] | None) -> pd.DataFrame:
    base = base.copy()
    base["literature_evidence_tier"] = ""
    base["cms_pathway_hints"] = ""
    base["mapk_liquid_biomarker_context"] = False
    if not layers:
        return base
    gene_col = "gene_lookup"
    if gene_col not in base.columns:
        base[gene_col] = base["feature_id"].astype(str).str.upper()
    tiers_cfg = layers.get("evidence_tiers") or {}
    tier_priority = ["tier1_validated_clinical", "tier2_emerging", "tier3_pathway_anchor"]
    gene_label: dict[str, str] = {}
    for tkey in tier_priority:
        spec = tiers_cfg.get(tkey) or {}
        lbl = str(spec.get("label", tkey))
        for g in spec.get("gene_symbols", []) or []:
            gu = str(g).strip().upper()
            if gu and gu not in gene_label:
                gene_label[gu] = lbl
    base["literature_evidence_tier"] = base[gene_col].map(lambda x: gene_label.get(str(x).strip().upper(), ""))

    cms_cfg = layers.get("cms_pathway_hints") or {}
    cms_tags: list[str] = []
    for _, row in base.iterrows():
        gu = str(row[gene_col]).strip().upper()
        hits = [
            short
            for short, spec in cms_cfg.items()
            if gu in {str(x).strip().upper() for x in (spec or {}).get("gene_symbols", []) or []}
        ]
        cms_tags.append(";".join(hits))
    base["cms_pathway_hints"] = cms_tags

    liq = layers.get("hallmarks_liquid_biopsy_context") or {}
    liq_genes = {str(x).strip().upper() for x in (liq.get("MAPK_EGFR_axis") or {}).get("gene_symbols", []) or []}
    base["mapk_liquid_biomarker_context"] = base[gene_col].astype(str).str.upper().isin(liq_genes)
    return base


def _annotate_functional_modules(base: pd.DataFrame, layers: dict[str, Any] | None) -> pd.DataFrame:
    base = base.copy()
    base["functional_module_hints"] = ""
    if not layers:
        return base
    fm = layers.get("functional_modules") or {}
    gene_col = "gene_lookup"
    tags: list[str] = []
    for _, row in base.iterrows():
        gu = str(row[gene_col]).strip().upper()
        hits = [
            short
            for short, spec in fm.items()
            if gu in {str(x).strip().upper() for x in (spec or {}).get("gene_symbols", []) or []}
        ]
        tags.append(";".join(hits))
    base["functional_module_hints"] = tags
    return base


def build_translation_table(
    diff: pd.DataFrame,
    hubs: pd.DataFrame | None,
    ml_scores: pd.DataFrame | None,
    *,
    tiers: list[str],
    context_layers_path: Path | None = None,
) -> pd.DataFrame:
    """Integra evidencia diferencial, topología y ML en un ranking traslacional."""
    base = diff.rename(columns={"protein_id": "feature_id"}).copy()
    if "gene_symbol" in diff.columns:
        base["gene_lookup"] = diff["gene_symbol"].astype(str).str.upper()
    else:
        base["gene_lookup"] = base["feature_id"].astype(str).str.upper()
    base["neglog10_padj"] = -np.log10(np.clip(base["padj"].astype(float), 1e-300, None))
    base["evidence_score"] = (
        np.abs(base["log2fc"]).rank(pct=True) * 0.35 + base["neglog10_padj"].rank(pct=True) * 0.35
    )
    if hubs is not None and not hubs.empty:
        h = hubs.rename(columns={"protein_id": "feature_id"})
        base = base.merge(h[["feature_id", "degree", "is_hub"]], on="feature_id", how="left")
        base["degree"] = base["degree"].fillna(0)
        base["is_hub"] = base["is_hub"].fillna(False)
        base["evidence_score"] += base["degree"].rank(pct=True) * 0.15
    else:
        base["degree"] = 0
        base["is_hub"] = False
    if ml_scores is not None and not ml_scores.empty:
        m = ml_scores.reset_index()
        first = m.columns[0]
        m = m.rename(columns={first: "feature_id"})
        cols = ["feature_id"] + [c for c in m.columns if c == "combined_score"]
        if "combined_score" in m.columns:
            base = base.merge(m[cols], on="feature_id", how="left")
            base["combined_score"] = base["combined_score"].fillna(0)
            base["evidence_score"] += base["combined_score"].rank(pct=True) * 0.15
    if "combined_score" not in base.columns:
        base["combined_score"] = 0.0

    def tier(row: pd.Series) -> str:
        s = row["evidence_score"]
        padj = row.get("padj", 1.0)
        if len(tiers) >= 4:
            q75, q50, q25 = np.quantile(base["evidence_score"], [0.75, 0.5, 0.25])
            if s >= q75 and pd.notna(padj) and padj < 0.05:
                return tiers[0]
            if s >= q50:
                return tiers[1] if len(tiers) > 1 else tiers[0]
            if s >= q25:
                return tiers[2] if len(tiers) > 2 else tiers[-1]
            return tiers[3] if len(tiers) > 3 else tiers[-1]
        q1, q2 = np.quantile(base["evidence_score"], [0.66, 0.33])
        if s >= q1 and pd.notna(padj) and padj < 0.05:
            return tiers[0] if tiers else "validated"
        if s >= q2:
            return tiers[1] if len(tiers) > 1 else "hypothesis_generating"
        return tiers[2] if len(tiers) > 2 else "future_validation"

    base["evidence_tier"] = base.apply(tier, axis=1)
    layers = _load_context_layers(context_layers_path)
    base = _annotate_context_layers(base, layers)
    base = _annotate_functional_modules(base, layers)
    lit = (base["literature_evidence_tier"].astype(str).str.len() > 0).astype(float)
    cms = (base["cms_pathway_hints"].astype(str).str.len() > 0).astype(float)
    fm = (base["functional_module_hints"].astype(str).str.len() > 0).astype(float)
    mapk = base["mapk_liquid_biomarker_context"].astype(float)
    z_fc = np.abs(base["log2fc"]).rank(pct=True)
    z_p = base["neglog10_padj"].rank(pct=True)
    z_deg = base["degree"].rank(pct=True)
    z_ml = base["combined_score"].rank(pct=True)
    base["translational_priority_score"] = (
        0.22 * z_fc + 0.22 * z_p + 0.12 * z_deg + 0.14 * z_ml + 0.15 * lit + 0.07 * cms + 0.05 * fm + 0.03 * mapk
    ).clip(0, 1)
    return base.sort_values(["translational_priority_score", "evidence_score"], ascending=[False, False])
