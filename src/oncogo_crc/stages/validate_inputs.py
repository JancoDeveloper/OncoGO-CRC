from __future__ import annotations

from oncogo_crc.settings import Settings
from oncogo_crc.utils.input_validation import (
    validate_dataset_intent,
    validate_remote_sources,
    validate_matrix,
    validate_metadata,
    validate_sample_matching,
    write_validation_reports,
)


def run(settings: Settings) -> None:
    val_cfg = settings.config.get("validation", {})
    recommended_vars = list(val_cfg.get("recommended_clinical_vars", []))
    report: dict[str, dict] = {"datasets": {}}
    for dataset_name, dcfg in settings.datasets.get("datasets", {}).items():
        if dcfg.get("enabled") is False:
            continue
        matrix_path = settings.resolve(dcfg["matrix_path"])
        metadata_path = settings.resolve(dcfg["sample_metadata"])
        report["datasets"][dataset_name] = {
            "matrix": validate_matrix(matrix_path),
            "metadata": validate_metadata(metadata_path, recommended_vars),
            "matching": validate_sample_matching(matrix_path, metadata_path),
            "intent_validation": validate_dataset_intent(
                dataset_name=dataset_name,
                matrix_path=matrix_path,
                metadata_path=metadata_path,
                comparisons_cfg=settings.comparisons,
            ),
        }
    out_dir = settings.resolve(settings.config.get("download_real", {}).get("audit_dir", "results/validation"))
    src_report = validate_remote_sources(settings.resolve("data/manifests/_summary.json"))
    report["remote_sources"] = src_report
    write_validation_reports(report, out_dir)
    fail_on_error = bool(settings.config.get("download_real", {}).get("fail_on_validation_error", False))
    if fail_on_error:
        bad = []
        for dataset, payload in report["datasets"].items():
            m = payload["matrix"]
            md = payload["metadata"]
            mt = payload["matching"]
            iv = payload.get("intent_validation", {})
            if m.get("errors") or md.get("errors") or not mt.get("exact_match", False):
                bad.append(dataset)
            if not iv.get("all_core_comparisons_supported", True):
                bad.append(f"{dataset}:core_comparison_not_supported")
        if report.get("remote_sources", {}).get("has_file_scheme", False):
            bad.append("remote_sources_file_scheme")
        if report.get("remote_sources", {}).get("has_invalid_remote_url", False):
            bad.append("remote_sources_invalid_url")
        if bad:
            raise RuntimeError(f"validate_inputs failed for datasets: {', '.join(bad)}")
    print(f"[validate_inputs] report -> {settings.relative_to_project(out_dir)}")
