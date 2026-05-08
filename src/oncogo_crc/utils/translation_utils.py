from __future__ import annotations

import numpy as np
import pandas as pd


def build_translation_table(
    diff: pd.DataFrame,
    hubs: pd.DataFrame | None,
    ml_scores: pd.DataFrame | None,
    *,
    tiers: list[str],
) -> pd.DataFrame:
    """Integra evidencia diferencial, topología y ML en un ranking traslacional."""
    base = diff.rename(columns={"protein_id": "feature_id"}).copy()
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
        q1, q2 = np.quantile(base["evidence_score"], [0.66, 0.33])
        padj = row.get("padj", 1.0)
        if s >= q1 and pd.notna(padj) and padj < 0.05:
            return tiers[0] if tiers else "validated"
        if s >= q2:
            return tiers[1] if len(tiers) > 1 else "hypothesis_generating"
        return tiers[2] if len(tiers) > 2 else "future_validation"

    base["evidence_tier"] = base.apply(tier, axis=1)
    return base.sort_values("evidence_score", ascending=False)
