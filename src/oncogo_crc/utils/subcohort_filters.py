from __future__ import annotations

from typing import Any

import pandas as pd


def apply_subcohort_filters(
    matrix_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    *,
    group_col: str,
    case_label: str,
    control_label: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    report: dict[str, Any] = {
        "applied": False,
        "group_col": group_col,
        "case": case_label,
        "control": control_label,
        "n_case": 0,
        "n_control": 0,
        "warnings": [],
    }
    if group_col not in meta_df.columns:
        report["warnings"].append(f"group_col '{group_col}' no existe")
        return matrix_df, meta_df, report

    keep_meta = meta_df[meta_df[group_col].isin([case_label, control_label])].copy()
    report["n_case"] = int((keep_meta[group_col] == case_label).sum())
    report["n_control"] = int((keep_meta[group_col] == control_label).sum())
    report["applied"] = True

    keep_samples = set(keep_meta["sample_id"].astype(str))
    keep_cols = [matrix_df.columns[0]] + [c for c in matrix_df.columns[1:] if str(c) in keep_samples]
    keep_matrix = matrix_df.loc[:, keep_cols].copy()
    return keep_matrix, keep_meta, report


def passes_group_minimums(meta_df: pd.DataFrame, group_col: str, min_n: int) -> dict[str, Any]:
    if group_col not in meta_df.columns:
        return {"group_col": group_col, "min_n": min_n, "counts": {}, "passes": False}
    counts = meta_df[group_col].value_counts(dropna=False).to_dict()
    return {"group_col": group_col, "min_n": min_n, "counts": counts, "passes": all(v >= min_n for v in counts.values())}
