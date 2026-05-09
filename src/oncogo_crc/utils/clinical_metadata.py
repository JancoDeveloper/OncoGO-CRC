from __future__ import annotations

from typing import Any

import pandas as pd

SEX_MAP = {
    "male": "male",
    "m": "male",
    "female": "female",
    "f": "female",
    "unknown": "unknown",
    "not reported": "unknown",
    "na": "unknown",
    "none": "unknown",
}


def normalize_sex(value: Any) -> str:
    if value is None or pd.isna(value):
        return "unknown"
    return SEX_MAP.get(str(value).strip().lower(), "unknown")


def normalize_age(value: Any) -> int | None:
    if value is None or pd.isna(value):
        return None
    try:
        return int(float(value))
    except Exception:  # noqa: BLE001
        return None


def bin_age(age: Any, bins_cfg: list[dict[str, Any]] | None = None) -> str:
    a = normalize_age(age)
    if a is None:
        return "unknown"
    if bins_cfg:
        for item in bins_cfg:
            if a >= int(item.get("min", -999)) and a <= int(item.get("max", 999)):
                return str(item.get("name", "unknown"))
    return "older" if a >= 50 else "young"


def normalize_stage(value: Any) -> str:
    if value is None or pd.isna(value):
        return "unknown"
    return str(value).strip().lower()


def bin_stage(stage: Any, stage_mapping: dict[str, list[str]] | None) -> str:
    s = normalize_stage(stage)
    if s == "unknown":
        return "unknown"
    for label, values in (stage_mapping or {}).items():
        vals = {str(v).strip().lower() for v in values}
        if s in vals:
            return str(label)
    return "unknown"


def normalize_msi(value: Any) -> str:
    if value is None or pd.isna(value):
        return "unknown"
    s = str(value).strip().upper()
    if s in {"MSI-H", "MSIH", "HIGH"}:
        return "MSI-H"
    if s in {"MSS", "STABLE"}:
        return "MSS"
    if s in {"MSI-L", "MSIL", "LOW"}:
        return "MSI-L"
    return "unknown"


def harmonize_clinical_metadata(df: pd.DataFrame, global_cfg: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    bins_cfg = global_cfg.get("clinical_filters", {}).get("age_bins", {}).get("bins", [])
    stage_map = global_cfg.get("clinical_filters", {}).get("stage_mapping", {})

    if "sex" not in out.columns and "gender" in out.columns:
        out["sex"] = out["gender"]
    if "sex" not in out.columns and "gender_raw" in out.columns:
        out["sex"] = out["gender_raw"]

    out["sex"] = out.get("sex", pd.Series(index=out.index, dtype="object")).apply(normalize_sex)
    out["age"] = out.get("age", pd.Series(index=out.index, dtype="object")).apply(normalize_age)
    out["age_bin"] = out["age"].apply(lambda x: bin_age(x, bins_cfg))
    out["stage"] = out.get("stage", pd.Series(index=out.index, dtype="object")).apply(normalize_stage)
    out["stage_bin"] = out["stage"].apply(lambda x: bin_stage(x, stage_map))
    out["msi_status"] = out.get("msi_status", pd.Series(index=out.index, dtype="object")).apply(normalize_msi)
    if "metastatic_status" not in out.columns:
        out["metastatic_status"] = "unknown"
    if "treatment_response" not in out.columns:
        out["treatment_response"] = "unknown"
    if "condition" not in out.columns:
        out["condition"] = "unknown"

    ordered = [
        "sample_id",
        "condition",
        "sex",
        "age",
        "age_bin",
        "stage",
        "stage_bin",
        "msi_status",
        "metastatic_status",
        "treatment_response",
    ]
    existing = [c for c in ordered if c in out.columns]
    rest = [c for c in out.columns if c not in existing]
    return out[existing + rest]
