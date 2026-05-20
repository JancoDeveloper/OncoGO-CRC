from __future__ import annotations

import json
import urllib.parse
from pathlib import Path
from typing import Any

import pandas as pd


def validate_matrix(matrix_path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "exists": matrix_path.exists(),
        "first_column_ok": False,
        "duplicate_protein_id": False,
        "duplicate_sample_id": False,
        "n_rows": 0,
        "n_samples": 0,
        "errors": [],
    }
    if not matrix_path.exists():
        result["errors"].append("matrix_missing")
        return result
    df = pd.read_csv(matrix_path, sep="\t", low_memory=False)
    result["n_rows"] = int(df.shape[0])
    result["n_samples"] = int(max(df.shape[1] - 1, 0))
    if len(df.columns) > 0 and df.columns[0] == "protein_id":
        result["first_column_ok"] = True
    else:
        result["errors"].append("first_column_not_protein_id")
    if "protein_id" in df.columns and bool(df["protein_id"].duplicated().any()):
        result["duplicate_protein_id"] = True
        result["errors"].append("duplicate_protein_id")
    sample_cols = list(df.columns[1:])
    if len(sample_cols) != len(set(sample_cols)):
        result["duplicate_sample_id"] = True
        result["errors"].append("duplicate_sample_id_columns")
    return result


def validate_metadata(metadata_path: Path, recommended_vars: list[str]) -> dict[str, Any]:
    result: dict[str, Any] = {
        "exists": metadata_path.exists(),
        "sample_id_ok": False,
        "condition_ok": False,
        "n_rows": 0,
        "missing_recommended_vars": [],
        "coverage": {},
        "errors": [],
    }
    if not metadata_path.exists():
        result["errors"].append("metadata_missing")
        return result
    df = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    result["n_rows"] = int(df.shape[0])
    if "sample_id" in df.columns:
        result["sample_id_ok"] = True
    else:
        result["errors"].append("metadata_missing_sample_id")
    if "condition" in df.columns:
        result["condition_ok"] = True
    else:
        result["errors"].append("metadata_missing_condition")
    sentinel_missing = {"", "unknown", "nan", "none", "na", "n/a"}
    for var in recommended_vars:
        if var not in df.columns:
            result["missing_recommended_vars"].append(var)
        else:
            col = df[var]
            informative = col.notna() & ~col.astype(str).str.strip().str.lower().isin(sentinel_missing)
            result["coverage"][var] = float(int(informative.sum()) / max(len(df), 1))
    return result


def validate_sample_matching(matrix_path: Path, metadata_path: Path) -> dict[str, Any]:
    out: dict[str, Any] = {"exact_match": False, "matrix_only_samples": [], "metadata_only_samples": []}
    if not (matrix_path.exists() and metadata_path.exists()):
        return out
    mdf = pd.read_csv(matrix_path, sep="\t", nrows=0, low_memory=False)
    cdf = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    matrix_samples = set(map(str, mdf.columns[1:]))
    meta_samples = set(map(str, cdf["sample_id"].astype(str).tolist())) if "sample_id" in cdf.columns else set()
    out["matrix_only_samples"] = sorted(list(matrix_samples - meta_samples))
    out["metadata_only_samples"] = sorted(list(meta_samples - matrix_samples))
    out["exact_match"] = (len(out["matrix_only_samples"]) == 0) and (len(out["metadata_only_samples"]) == 0)
    return out


def write_validation_reports(report: dict[str, Any], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "input_validation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    rows: list[dict[str, Any]] = []
    for dataset, payload in report.get("datasets", {}).items():
        row: dict[str, Any] = {
            "dataset": dataset,
            "matrix_exists": payload["matrix"].get("exists"),
            "metadata_exists": payload["metadata"].get("exists"),
            "exact_match": payload["matching"].get("exact_match"),
            "n_samples_matrix": payload["matrix"].get("n_samples"),
            "n_rows_matrix": payload["matrix"].get("n_rows"),
            "n_rows_metadata": payload["metadata"].get("n_rows"),
        }
        for k, v in payload["metadata"].get("coverage", {}).items():
            row[f"coverage_{k}"] = v
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_dir / "input_validation_report.tsv", sep="\t", index=False)


def validate_dataset_intent(
    *,
    dataset_name: str,
    matrix_path: Path,
    metadata_path: Path,
    comparisons_cfg: dict[str, Any],
) -> dict[str, Any]:
    """
    Verifica si lo descargado sirve para la consulta analítica esperada.
    """
    out: dict[str, Any] = {
        "dataset": dataset_name,
        "comparisons_checked": [],
        "all_core_comparisons_supported": True,
        "reasons": [],
    }
    if not matrix_path.exists() or not metadata_path.exists():
        out["all_core_comparisons_supported"] = False
        out["reasons"].append("matrix_or_metadata_missing")
        return out

    mdf = pd.read_csv(matrix_path, sep="\t", nrows=0, low_memory=False)
    cdf = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    matrix_samples = set(map(str, mdf.columns[1:]))

    comparisons = (comparisons_cfg or {}).get("comparisons", {})
    for comp_name, comp in comparisons.items():
        if comp.get("dataset") != dataset_name:
            continue
        group_col = comp.get("condition_column") or comp.get("group_col")
        case = comp.get("group_a") or comp.get("case")
        control = comp.get("group_b") or comp.get("control")
        optional = bool(comp.get("optional", False))
        entry = {
            "comparison": comp_name,
            "optional": optional,
            "group_col": group_col,
            "case": case,
            "control": control,
            "supported": True,
            "reason": "ok",
        }
        if group_col not in cdf.columns:
            entry["supported"] = False
            entry["reason"] = f"group_col_missing:{group_col}"
        else:
            vals = set(map(str, cdf[group_col].dropna().astype(str)))
            if str(case) not in vals or str(control) not in vals:
                entry["supported"] = False
                entry["reason"] = "case_or_control_not_present"
        if "sample_id" in cdf.columns:
            meta_samples = set(map(str, cdf["sample_id"].astype(str)))
            if len(matrix_samples & meta_samples) == 0:
                entry["supported"] = False
                entry["reason"] = "no_sample_overlap"

        if (not optional) and (not entry["supported"]):
            out["all_core_comparisons_supported"] = False
        out["comparisons_checked"].append(entry)

    if out["all_core_comparisons_supported"]:
        out["reasons"].append("core_comparisons_supported")
    else:
        out["reasons"].append("some_core_comparisons_not_supported")
    return out


def validate_remote_sources(manifest_json_path: Path) -> dict[str, Any]:
    out: dict[str, Any] = {"exists": manifest_json_path.exists(), "datasets": {}, "has_file_scheme": False, "has_invalid_remote_url": False}
    if not manifest_json_path.exists():
        return out
    payload = json.loads(manifest_json_path.read_text(encoding="utf-8"))
    datasets = payload.get("datasets", payload if isinstance(payload, dict) else {})
    if isinstance(datasets, list):
        datasets = {row.get("dataset", f"dataset_{i}"): row for i, row in enumerate(datasets)}

    def _is_useful_remote(u: str) -> bool:
        ul = (u or "").lower()
        if ul == "":
            return False
        if ul.startswith("file://"):
            return False
        if not ul.startswith(("http://", "https://", "ftp://")):
            return False
        bad = [".css", ".js", ".svg", ".png", ".jpg", ".jpeg", ".gif", "fonts.gstatic", "/cgi/", "view=datasets"]
        if any(b in ul for b in bad):
            return False
        parsed = urllib.parse.urlparse(ul)
        path = parsed.path or ""
        if path in {"", "/"}:
            return False
        good_ext = [".tsv", ".txt", ".csv", ".xlsx", ".xls", ".zip", ".gz"]
        if any(path.endswith(ext) for ext in good_ext):
            return True
        return "/download/" in path or "/files/" in path

    for ds, row in datasets.items():
        matrix_url = str((row or {}).get("matrix_url") or "")
        metadata_url = str((row or {}).get("metadata_url") or "")
        matrix_status = str((row or {}).get("matrix_status") or "")
        metadata_status = str((row or {}).get("metadata_status") or "")
        has_file = matrix_url.startswith("file://") or metadata_url.startswith("file://")
        if has_file:
            out["has_file_scheme"] = True
        matrix_useful = _is_useful_remote(matrix_url) if matrix_status == "ok" else True
        metadata_useful = _is_useful_remote(metadata_url) if metadata_status == "ok" else True
        invalid = not (matrix_useful and metadata_useful)
        if invalid:
            out["has_invalid_remote_url"] = True
        out["datasets"][ds] = {
            "matrix_url": matrix_url,
            "metadata_url": metadata_url,
            "matrix_remote_http": matrix_url.startswith(("http://", "https://", "ftp://")) or matrix_url == "",
            "metadata_remote_http": metadata_url.startswith(("http://", "https://", "ftp://")) or metadata_url == "",
            "has_file_scheme": has_file,
            "matrix_status": matrix_status,
            "metadata_status": metadata_status,
            "matrix_url_useful": matrix_useful,
            "metadata_url_useful": metadata_useful,
            "invalid_remote_url_for_ok_status": invalid,
        }
    return out
