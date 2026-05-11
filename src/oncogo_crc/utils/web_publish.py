from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from oncogo_crc.settings import Settings


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return {}


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def build_web_pages(settings: Settings) -> Path:
    global_cfg = settings.config
    web_cfg = global_cfg.get("web_publish", {})
    web_dir = settings.resolve(str(web_cfg.get("web_reports_dir", "reports/web")))
    docs_web_dir = settings.root / "docs" / "web"
    runs_dir = settings.resolve(str(web_cfg.get("runs_dir", "reports/runs")))
    latest = web_cfg.get("latest_run_alias", "latest")
    latest_dir = runs_dir / latest
    latest_dir.mkdir(parents=True, exist_ok=True)

    summary = _read_json(settings.resolve("data/manifests/_summary.json"))
    if not summary:
        summary = _read_json(settings.resolve("data/raw/download_real_manifest.json"))
        if isinstance(summary, list):
            summary = {"datasets": {row.get("dataset", f"dataset_{i}"): row for i, row in enumerate(summary)}}
    validation = _read_json(settings.resolve("results/validation/input_validation_report.json"))
    stability_json = _read_json(settings.resolve("results/stability/stability_summary.json"))

    run_manifest = {
        "datasets": summary.get("datasets", {}),
        "validation": {"datasets": validation.get("datasets", {})},
        "stability": stability_json,
    }
    (latest_dir / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2), encoding="utf-8")

    index_md = """# OncoGO-CRC

## Estado general

Este sitio resume el avance reproducible del pipeline:

- estado de datasets
- cobertura clínica
- resultados base
- resultados estratificados
- estabilidad entre estratos
- limitaciones metodológicas
"""
    datasets_md = "# Datasets\n\n"
    for ds, payload in summary.get("datasets", {}).items():
        datasets_md += f"## {ds}\n"
        datasets_md += f"- status: {payload.get('status', payload.get('matrix_status', 'unknown'))}\n"
        if "provider" in payload:
            datasets_md += f"- provider: {payload.get('provider')}\n"
        if "n_files" in payload:
            datasets_md += f"- n_files: {payload.get('n_files')}\n"
        if payload.get("matrix_url"):
            datasets_md += f"- matrix_url: {payload.get('matrix_url')}\n"
        if payload.get("metadata_url"):
            datasets_md += f"- metadata_url: {payload.get('metadata_url')}\n"
        datasets_md += "\n"

    coverage_md = "# Cobertura clínica\n\n"
    for ds, payload in validation.get("datasets", {}).items():
        coverage_md += f"## {ds}\n"
        cov = payload.get("metadata", {}).get("coverage", {})
        if not cov:
            coverage_md += "- Sin cobertura clínica consolidada todavía.\n\n"
            continue
        for k, v in cov.items():
            coverage_md += f"- {k}: {v:.2%}\n"
        coverage_md += "\n"

    stability_md = "# Estabilidad\n\n"
    if isinstance(stability_json, list) and stability_json:
        stability_md += "| dataset | core | estratificado | n_core | n_strat | intersección | jaccard |\n"
        stability_md += "|---|---|---:|---:|---:|---:|---:|\n"
        for r in stability_json:
            stability_md += (
                f"| {r.get('dataset','')} | {r.get('core_comparison','')} | {r.get('stratified_comparison','')} "
                f"| {r.get('n_sig_core',0)} | {r.get('n_sig_stratified',0)} | {r.get('n_intersection',0)} | {r.get('jaccard',0):.3f} |\n"
            )
    else:
        stability_md += "Sin resultados de estabilidad todavía.\n"

    pages = {
        "index.md": index_md,
        "datasets.md": datasets_md,
        "clinical_coverage.md": coverage_md,
        "results_base.md": "# Resultados base\n\nEn construcción automática desde outputs de pipeline.\n",
        "results_stratified.md": "# Resultados estratificados\n\nEn construcción automática desde comparaciones stratified.\n",
        "stability.md": stability_md,
        "limitations.md": "# Limitaciones\n\nEn construcción automática.\n",
    }
    for name, content in pages.items():
        _write(web_dir / name, content)
        _write(docs_web_dir / name, content)
    return web_dir
