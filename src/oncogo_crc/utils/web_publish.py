from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv
from oncogo_crc.utils.run_manifest import load_doc as load_pipeline_run_manifest
from oncogo_crc.utils.thesis_web_digest import (
    build_external_evidence_markdown,
    build_limitations_markdown,
    build_results_base_markdown,
    build_results_inventory_markdown,
    build_thesis_synthesis_markdown,
    build_web_index_markdown,
)
from oncogo_crc.utils.web_interpretation import (
    render_glossary_markdown,
    render_interpretation_markdown,
)
from oncogo_crc.utils.clinical_demographics import (
    render_clinical_demographics_markdown,
    write_demographics_artifacts,
)
from oncogo_crc.utils.web_visualizations import (
    plot_clinical_coverage,
    plot_clinical_demographics_bars,
    plot_cohort_overlap_scatter,
    plot_gsea_top_nes,
    plot_nested_cv_auc,
    plot_pca_by_condition,
    plot_top_features_heatmap,
    plot_top_go_ora,
    plot_transfer_auc_bars,
    plot_volcano_with_labels,
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

    try:
        write_demographics_artifacts(settings)
    except OSError:
        pass

    coverage_md = (
        "# Cobertura clínica\n\n"
        "Fracción de muestras con valor **informativo** por variable (≠ vacío / unknown). "
        "Para **conteos y distribución** (estadio, hispano/latino, región geográfica), ver "
        "[Demografía y contexto clínico](clinical_demographics.md).\n\n"
    )
    for ds, payload in validation.get("datasets", {}).items():
        coverage_md += f"## {ds}\n"
        cov = payload.get("metadata", {}).get("coverage", {})
        if not cov:
            coverage_md += "- Sin cobertura clínica consolidada todavía.\n\n"
            continue
        for k, v in sorted(cov.items()):
            coverage_md += f"- `{k}`: {v:.2%}\n"
        coverage_md += "\n"

    demographics_md = render_clinical_demographics_markdown(settings)

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

    # === Generación de figuras (PNG) que las páginas Markdown embeben ===
    figures_subdir = "figures"
    figs_web = web_dir / figures_subdir
    figs_docs = docs_web_dir / figures_subdir
    generated_figs = _generate_web_figures(
        settings, figs_web, figs_docs, coverage_payload=validation.get("datasets", {})
    )

    # === Interpretación biológica por comparación ===
    interpretation_md = _build_interpretation_page(settings)

    glossary_md = render_glossary_markdown()
    figures_index_md = _build_figures_index_md(generated_figs, figures_subdir)

    pages = {
        "index.md": index_md,
        "datasets.md": datasets_md,
        "clinical_coverage.md": _augment_coverage_md(coverage_md, generated_figs, figures_subdir),
        "clinical_demographics.md": _augment_demographics_md(
            demographics_md, generated_figs, figures_subdir
        ),
        "thesis_synthesis.md": build_thesis_synthesis_markdown(settings),
        "results_base.md": _augment_results_base_md(
            build_results_base_markdown(settings), generated_figs, figures_subdir
        ),
        "results_inventory.md": build_results_inventory_markdown(settings),
        "results_stratified.md": "# Resultados estratificados\n\n"
        "Comparaciones con `optional: true` en `config/comparisons.yaml`: "
        "si hay estratos y muestras suficientes, sus salidas aparecerán bajo `results/` con el prefijo de comparación.\n",
        "stability.md": stability_md,
        "limitations.md": build_limitations_markdown(settings),
        "external_evidence.md": build_external_evidence_markdown(settings),
        "glosario.md": glossary_md,
        "interpretacion_biologica.md": interpretation_md,
        "figuras.md": figures_index_md,
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
            # Espejamos también las figuras.
            repo_figs = repo_docs_web / figures_subdir
            repo_figs.mkdir(parents=True, exist_ok=True)
            if figs_docs.exists():
                for fig in figs_docs.glob("*.png"):
                    shutil.copy2(fig, repo_figs / fig.name)
    return web_dir


# ---------------------------------------------------------------------------
# Helpers — figuras e interpretación biológica
# ---------------------------------------------------------------------------
def _generate_web_figures(
    settings: Settings,
    figs_web: Path,
    figs_docs: Path,
    *,
    coverage_payload: dict[str, Any],
) -> dict[str, Path]:
    """Genera todos los PNG e indexa por nombre lógico → ruta dentro de ``figs_docs``.

    Las claves usan el patrón ``<comp_or_ds>__<kind>`` para que las páginas Markdown puedan
    embeberlas con `![](figures/<nombre>.png)` sin colisiones.
    """
    figs_web.mkdir(parents=True, exist_ok=True)
    figs_docs.mkdir(parents=True, exist_ok=True)
    res = settings.resolve(settings.config["paths"]["results_dir"])
    proc = settings.resolve(settings.config["paths"]["processed_dir"])
    out: dict[str, Path] = {}

    def _write_both(key: str, name: str, generator) -> None:
        web_path = figs_web / name
        result = generator(web_path)
        if result is not None:
            docs_path = figs_docs / name
            try:
                shutil.copy2(web_path, docs_path)
            except OSError:
                return
            out[key] = docs_path

    # Cobertura clínica (consolidada todos los datasets).
    cov_dict: dict[str, dict[str, float]] = {}
    for ds, payload in coverage_payload.items():
        cov = payload.get("metadata", {}).get("coverage", {})
        if cov:
            cov_dict[ds] = {k: float(v) for k, v in cov.items()}
    if cov_dict:
        _write_both("clinical_coverage", "clinical_coverage.png",
                    lambda p: plot_clinical_coverage(cov_dict, p))

    demo_json = settings.results_dir() / "clinical" / "demographics_summary.json"
    if demo_json.is_file():
        summary_demo = json.loads(demo_json.read_text(encoding="utf-8"))
        _write_both(
            "clinical_demographics_cptac",
            "clinical_demographics_cptac.png",
            lambda p: plot_clinical_demographics_bars(summary_demo, p, dataset="cptac_colorectal"),
        )

    # Por comparación principal: volcano enhanced + top GO + GSEA NES.
    comps = settings.comparisons.get("comparisons", {})
    alpha = float(settings.config.get("differential", {}).get("alpha", 0.05))
    gsea_fdr = float(settings.config.get("gsea", {}).get("fdr_threshold", 0.25))
    for comp_key, spec in comps.items():
        diff = res / "differential" / f"{comp_key}_diff.tsv"
        go = res / "go" / f"{comp_key}_go.tsv"
        gsea = res / "gsea" / f"{comp_key}_gsea_go.tsv"
        label = str(spec.get("label", comp_key))
        if diff.exists():
            _write_both(f"{comp_key}__volcano", f"{comp_key}_volcano_labeled.png",
                        lambda p, d=diff, lab=label: plot_volcano_with_labels(d, p, title=lab, alpha=alpha))
        if go.exists():
            _write_both(f"{comp_key}__ora", f"{comp_key}_ora_top_go.png",
                        lambda p, g=go, lab=label: plot_top_go_ora(g, p, title=f"ORA — {lab}"))
        if gsea.exists():
            _write_both(f"{comp_key}__gsea_nes", f"{comp_key}_gsea_top_nes.png",
                        lambda p, g=gsea, lab=label, f=gsea_fdr: plot_gsea_top_nes(
                            g, p, fdr_threshold=f, title=f"GSEA — {lab}"
                        ))

    # Cohort overlap (CPTAC ∩ PRIDE) scatter.
    overlap_tsv = res / "reproducibility" / "cptac_pride_differential_overlap.tsv"
    if overlap_tsv.exists():
        _write_both("cohort_overlap", "cohort_overlap_scatter.png",
                    lambda p: plot_cohort_overlap_scatter(overlap_tsv, p))

    # PCA por dataset (log2norm).
    for ds in settings.datasets.get("datasets", {}):
        pca = res / "qc" / f"{ds}_pca_log2norm.tsv"
        meta = settings.resolve(settings.datasets["datasets"][ds].get("sample_metadata", ""))
        if pca.exists() and meta.exists():
            _write_both(f"{ds}__pca", f"{ds}_pca.png",
                        lambda p, q=pca, m=meta, lab=ds: plot_pca_by_condition(q, m, p, title=f"PCA — {lab}"))

    # Heatmap top features para la cohorte primaria.
    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")
    primary_log2 = proc / f"{primary}_log2norm.tsv"
    primary_meta = settings.resolve(settings.datasets["datasets"].get(primary, {}).get("sample_metadata", ""))
    primary_diff = res / "differential" / f"cptac_tumor_vs_normal_diff.tsv"
    if primary_log2.exists() and primary_diff.exists():
        _write_both(f"{primary}__heatmap_top", f"{primary}_top_features_heatmap.png",
                    lambda p, m=primary_log2, md=primary_meta, d=primary_diff:
                        plot_top_features_heatmap(m, md, d, p, top_n=40,
                                                  title=f"Top 40 proteínas — {primary} tumor vs normal"))

    # ML: nested CV boxplot + transfer AUC bars.
    per_fold = res / "ml" / f"{primary}_nested_cv_per_fold.tsv"
    if per_fold.exists():
        _write_both("ml__nested_cv", f"{primary}_nested_cv_auc.png",
                    lambda p, q=per_fold: plot_nested_cv_auc(q, p, title="Nested CV — ROC AUC por fold externo"))
    # Auto-detectar JSON de transfer (puede haber varios).
    transfer_dir = res / "ml"
    if transfer_dir.exists():
        for tj in sorted(transfer_dir.glob("*_train_*_test_transfer.json")):
            key = f"ml__transfer_{tj.stem}"
            _write_both(key, f"{tj.stem}_auc_bars.png",
                        lambda p, t=tj: plot_transfer_auc_bars(t, p))

    return out


def _build_interpretation_page(settings: Settings) -> str:
    """Página `interpretacion_biologica.md` con módulos detectados por comparación."""
    res = settings.resolve(settings.config["paths"]["results_dir"])
    gsea_fdr = float(settings.config.get("gsea", {}).get("fdr_threshold", 0.25))
    comps = settings.comparisons.get("comparisons", {})
    parts = [
        "# Interpretación biológica de los módulos detectados\n\n",
        "Esta página **traduce los términos GO/GSEA significativos a contexto biológico de cáncer "
        "colorrectal**, agrupando los _leading edges_ por familia funcional y emparejándolos con un "
        "catálogo curado a partir de literatura proteogenómica reciente (Cancer Cell 2023, Nature 2014, "
        "PLOS Comp Biol 2024, García-Foncillas et al. 2026).\n\n",
        "_Los textos siguientes son **hipótesis generadas automáticamente**, no diagnóstico. La "
        "discusión final de la tesis debe contrastar cada módulo con la literatura específica._\n\n",
        "Para definiciones de NES, leading edge, ORA vs GSEA, etc. ver el [Glosario](glosario.md).\n\n",
        "---\n\n",
    ]
    for comp_key, spec in comps.items():
        if spec.get("optional") and not (res / "gsea" / f"{comp_key}_gsea_go.tsv").exists():
            continue
        label = str(spec.get("label", comp_key))
        gsea_path = res / "gsea" / f"{comp_key}_gsea_go.tsv"
        df = None
        if gsea_path.exists():
            try:
                df = read_tsv(gsea_path)
            except Exception:  # noqa: BLE001
                df = None
        parts.append(render_interpretation_markdown(label, comp_key, df, fdr_threshold=gsea_fdr))
    return "".join(parts)


def _build_figures_index_md(figs: dict[str, Path], subdir: str) -> str:
    if not figs:
        return (
            "# Figuras\n\nLas figuras se generan en `publish_results`. Si esta página está vacía, "
            "revisa que las etapas correspondientes (`differential`, `go`, `gsea`, `cohort_overlap`, `ml`) "
            "hayan corrido y que matplotlib esté disponible.\n"
        )
    parts = [
        "# Figuras (generadas automáticamente)\n\n",
        "Cada figura se regenera en cada `publish_results`. Las versiones PNG residen bajo "
        f"`reports/web/{subdir}/` y se espejan en `docs/web/{subdir}/` para MkDocs.\n\n",
        "Para definiciones técnicas (volcano, NES, AUC, padj…) ver [Glosario](glosario.md).\n\n",
    ]
    groups: dict[str, list[tuple[str, Path]]] = {
        "Cobertura clínica": [],
        "Volcanos (proteínas diferenciales)": [],
        "ORA — top términos GO": [],
        "GSEA — top NES (up/down)": [],
        "Reproducibilidad multi-cohorte": [],
        "PCA por dataset": [],
        "Heatmaps top features": [],
        "ML — nested CV y transferencia": [],
    }
    for key, path in figs.items():
        if "clinical_coverage" in key:
            groups["Cobertura clínica"].append((key, path))
        elif "__volcano" in key:
            groups["Volcanos (proteínas diferenciales)"].append((key, path))
        elif "__ora" in key:
            groups["ORA — top términos GO"].append((key, path))
        elif "__gsea_nes" in key:
            groups["GSEA — top NES (up/down)"].append((key, path))
        elif "cohort_overlap" in key:
            groups["Reproducibilidad multi-cohorte"].append((key, path))
        elif "__pca" in key:
            groups["PCA por dataset"].append((key, path))
        elif "__heatmap_top" in key:
            groups["Heatmaps top features"].append((key, path))
        elif key.startswith("ml__"):
            groups["ML — nested CV y transferencia"].append((key, path))
    for section, items in groups.items():
        if not items:
            continue
        parts.append(f"## {section}\n\n")
        for key, path in sorted(items):
            parts.append(f"![{key}]({subdir}/{path.name})\n\n")
    return "".join(parts)


def _augment_results_base_md(md: str, figs: dict[str, Path], subdir: str) -> str:
    """Inserta enlaces a figuras relevantes al inicio de `results_base.md`."""
    extras = ["\n\n## Figuras de un vistazo\n\n"]
    overview_keys = [
        ("cohort_overlap", "Reproducibilidad CPTAC ∩ PRIDE (concordancia direccional de log2FC)"),
        ("clinical_coverage", "Cobertura clínica por variable y dataset"),
        ("ml__nested_cv", "Nested CV — distribución de AUC por fold externo"),
    ]
    any_added = False
    for k, caption in overview_keys:
        if k in figs:
            extras.append(f"**{caption}**\n\n![{k}]({subdir}/{figs[k].name})\n\n")
            any_added = True
    if not any_added:
        return md
    extras.append(
        "Ver también [**Figuras**](figuras.md) (volcanos, ORA top GO, GSEA NES, PCA, heatmaps por "
        "comparación) e [**Interpretación biológica**](interpretacion_biologica.md) para los módulos "
        "funcionales detectados.\n\n---\n\n"
    )
    return md + "".join(extras)


def _augment_demographics_md(md: str, figs: dict[str, Path], subdir: str) -> str:
    if "clinical_demographics_cptac" not in figs:
        return md
    return md + (
        "\n\n## Distribución CPTAC (muestras tumor)\n\n"
        f"![clinical_demographics_cptac]({subdir}/{figs['clinical_demographics_cptac'].name})\n\n"
        "_Conteos detallados en `results/clinical/demographics_counts.tsv`._\n"
    )


def _augment_coverage_md(md: str, figs: dict[str, Path], subdir: str) -> str:
    if "clinical_coverage" not in figs:
        return md
    return md + (
        "\n\n## Vista gráfica\n\n"
        f"![clinical_coverage]({subdir}/{figs['clinical_coverage'].name})\n\n"
        "_Cada barra reporta la fracción de muestras con valor informativo (no NaN ni `unknown`). "
        "Variables con cobertura cero indican que la metadata aún no se ha curado para ese dataset._\n"
    )
