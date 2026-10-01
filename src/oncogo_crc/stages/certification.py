"""Etapa de certificación: informe de brecha y subíndices R calculables.

Lee el cruce CPTAC ∩ PRIDE si existe. Si no existe, igual escribe el contrato
y la lista de lo que el pipeline actual todavía no cumple. No modifica el
diferencial ni la priorización.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone

import pandas as pd

from oncogo_crc.certification.contract import (
    H1_PATHWAY_CONCORDANCE_MIN,
    H1_PROTEIN_CONCORDANCE_MAX,
    H4_MEASURABLE_FRACTION_REDUNDANT,
    TREAT_LFC,
    build_certification_report,
)
from oncogo_crc.settings import Settings
from oncogo_crc.utils import run_manifest
from oncogo_crc.utils.io import read_tsv, write_tsv


def _fmt(value: float | None, digits: int = 3) -> str:
    if value is None:
        return "no calculado"
    return f"{value:.{digits}f}"


def _report_markdown(report_dict: dict, overlap_path: str | None) -> str:
    lines = [
        "# Informe de certificación (generado por la etapa `certification`)",
        "",
        f"Fecha UTC: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        "",
        "Este informe no es un capítulo de resultados. Dice qué se puede afirmar hoy",
        "con los artefactos del pipeline y qué queda fuera hasta ejecutar el diseño 3.0.",
        "",
        f"Cruce usado: `{overlap_path or 'ninguno'}`",
        "",
        "## Índice R",
        "",
        "| Subíndice | Valor | Estado |",
        "|---|---|---|",
    ]
    for sub in report_dict["subindices"]:
        val = "—" if sub["value"] is None else sub["value"]
        lines.append(f"| {sub['code']} | {val} | {sub['status']} |")
    lines += [
        "",
        "## Cláusulas evaluables",
        "",
        f"- Concordancia de signo (descubrimiento CPTAC, q < alpha): {_fmt(report_dict['sign_concordance'])}",
        f"- IC 95 % exacto: {_fmt(report_dict['sign_ci_low'])} – {_fmt(report_dict['sign_ci_high'])}",
        f"- Proteínas evaluables: {report_dict['n_sign_evaluable']} "
        f"(mismo signo: {report_dict['n_same_sign']})",
        f"- Cláusula de proteína de H1 (umbral < {H1_PROTEIN_CONCORDANCE_MAX}): "
        f"**{report_dict['h1_protein_clause']}**",
        f"- Cláusula de vía de H1 (umbral > {H1_PATHWAY_CONCORDANCE_MIN}): no evaluable en esta etapa.",
        f"- Proxy de medibilidad H4 (se falsifica si ≥ {H4_MEASURABLE_FRACTION_REDUNDANT}): "
        f"**{report_dict['h4_proxy']}** (fracción {_fmt(report_dict['measurable_fraction'])}).",
        f"- Umbral treat declarado y aún no aplicado al diferencial: log2FC = {TREAT_LFC}.",
        "",
        "## Lo que este pipeline todavía no cumple",
        "",
    ]
    for gap in report_dict["gaps"]:
        lines.append(f"- {gap}")
    lines.append("")
    return "\n".join(lines)


def run(settings: Settings) -> None:
    cfg = settings.config.get("certification", {}) or {}
    alpha = float(cfg.get("alpha", settings.config.get("differential", {}).get("alpha", 0.05)))
    res = settings.resolve(settings.config["paths"]["results_dir"])
    overlap_rel = str(
        cfg.get(
            "overlap_tsv",
            "results/reproducibility/cptac_pride_differential_overlap.tsv",
        )
    )
    overlap_path = settings.resolve(overlap_rel)
    overlap = read_tsv(overlap_path) if overlap_path.exists() else None
    report = build_certification_report(overlap, alpha=alpha)
    payload = report.to_dict()
    payload["alpha"] = alpha
    payload["overlap_tsv"] = str(overlap_path) if overlap is not None else None
    payload["thresholds"] = {
        "h1_protein_concordance_max": H1_PROTEIN_CONCORDANCE_MAX,
        "h1_pathway_concordance_min": H1_PATHWAY_CONCORDANCE_MIN,
        "h4_measurable_fraction_redundant": H4_MEASURABLE_FRACTION_REDUNDANT,
        "treat_lfc": TREAT_LFC,
    }

    out_dir = res / "certification"
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "certificado_resumen.json"
    md_path = out_dir / "informe_brecha.md"
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    md_path.write_text(
        _report_markdown(payload, str(overlap_path) if overlap is not None else None),
        encoding="utf-8",
    )

    rows = []
    for sub in report.subindices:
        rows.append(
            {
                "code": sub.code,
                "name": sub.name,
                "value": sub.value if sub.value is not None else "",
                "status": sub.status,
                "definition": sub.definition,
                "note": sub.note,
            }
        )
    write_tsv(pd.DataFrame(rows), out_dir / "indice_R.tsv")

    run_manifest.record_stage_meta(
        settings,
        "certification",
        executed=True,
        overlap_found=overlap is not None,
        h1_protein_clause=report.h1_protein_clause,
        h4_proxy=report.h4_proxy,
        outputs=["certificado_resumen.json", "informe_brecha.md", "indice_R.tsv"],
    )
    print(
        f"[certification] H1 proteína={report.h1_protein_clause}; "
        f"H4 proxy={report.h4_proxy} -> {settings.relative_to_project(md_path)}"
    )
