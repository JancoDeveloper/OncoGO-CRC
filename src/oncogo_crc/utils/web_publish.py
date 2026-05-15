from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

from oncogo_crc.settings import Settings
from oncogo_crc.utils.run_manifest import load_doc as load_pipeline_run_manifest
from oncogo_crc.utils.thesis_web_digest import (
    build_external_evidence_markdown,
    build_limitations_markdown,
    build_results_base_markdown,
    build_results_inventory_markdown,
    build_thesis_synthesis_markdown,
    build_web_index_markdown,
)


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
    docs_web_dir = settings.resolve("docs/web")
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

    pipeline_run = load_pipeline_run_manifest(settings)
    run_manifest_doc = {
        "datasets": summary.get("datasets", {}),
        "validation": {"datasets": validation.get("datasets", {})},
        "stability": stability_json,
        "pipeline_run": pipeline_run,
    }
    (latest_dir / "run_manifest.json").write_text(json.dumps(run_manifest_doc, indent=2), encoding="utf-8")

    index_md = build_web_index_markdown(settings)
    datasets_md = "# Datasets\n\n"
    for ds, payload in summary.get("datasets", {}).items():
        datasets_md += f"## {ds}\n"
        raw_status = payload.get("status", payload.get("matrix_status", "unknown"))
        display_status = raw_status
        if ds in settings.datasets.get("datasets", {}) and settings.processed_matrix_path(ds, "log2norm").exists():
            if str(raw_status).lower() in ("failed", "error", "download_failed", "partial"):
                display_status = f"{raw_status} (override: matriz log2norm local presente)"
        datasets_md += f"- status: {display_status}\n"
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
        "thesis_synthesis.md": build_thesis_synthesis_markdown(settings),
        "results_base.md": build_results_base_markdown(settings),
        "results_inventory.md": build_results_inventory_markdown(settings),
        "results_stratified.md": "# Resultados estratificados\n\n"
        "Comparaciones con `optional: true` en `config/comparisons.yaml`: "
        "si hay estratos y muestras suficientes, sus salidas aparecerán bajo `results/` con el prefijo de comparación.\n",
        "stability.md": stability_md,
        "limitations.md": build_limitations_markdown(settings),
        "external_evidence.md": build_external_evidence_markdown(settings),
    }
    for name, content in pages.items():
        _write(web_dir / name, content)
        _write(docs_web_dir / name, content)
    # MkDocs lee `docs/` bajo la raíz del repo; si el workspace es externo, espejamos docs/web.
    if settings.data_workspace is not None:
        repo_docs_web = settings.root / "docs" / "web"
        if docs_web_dir.resolve() != repo_docs_web.resolve():
            repo_docs_web.mkdir(parents=True, exist_ok=True)
            for name in pages:
                shutil.copy2(docs_web_dir / name, repo_docs_web / name)
    return web_dir
