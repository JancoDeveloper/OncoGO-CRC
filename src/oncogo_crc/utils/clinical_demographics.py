"""Resúmenes demográficos y clínicos para metadata de cohortes (tesis / sitio web)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from oncogo_crc.settings import Settings

# Variables categóricas a tabular (conteo por nivel, solo muestras tumor si existe condition).
_DEMO_VARS = [
    "stage",
    "stage_bin",
    "stage_group",
    "msi_status",
    "age_bin",
    "sex",
    "ethnicity",
    "race",
    "hispanic_latino_bin",
    "country",
    "geographic_region_bin",
    "cms_subtype",
    "cms_subtype_bin",
    "tumor_site",
]


def _informative_mask(series: pd.Series) -> pd.Series:
    sentinel = {"", "unknown", "nan", "none", "na", "n/a"}
    return series.notna() & ~series.astype(str).str.strip().str.lower().isin(sentinel)


def summarize_dataset_metadata(
    df: pd.DataFrame,
    *,
    dataset: str,
    tumor_only: bool = True,
) -> dict[str, Any]:
    """Conteos por variable; opcionalmente solo muestras ``condition=tumor``."""
    out: dict[str, Any] = {"dataset": dataset, "n_samples_total": int(len(df)), "variables": {}}
    if df.empty:
        return out

    work = df.copy()
    if tumor_only and "condition" in work.columns:
        work = work[work["condition"].astype(str).str.lower() == "tumor"]
    out["n_samples_tumor"] = int(len(work))

    # Participantes únicos (evita duplicar TMT shared/unshared)
    if "participant_id" in work.columns:
        pids = work["participant_id"].astype(str).str.strip()
        pids = pids[pids.notna() & (pids != "") & (pids.str.lower() != "nan")]
        out["n_participants_unique"] = int(pids.nunique())
    else:
        out["n_participants_unique"] = out["n_samples_tumor"]

    for var in _DEMO_VARS:
        if var not in work.columns:
            continue
        col = work[var]
        mask = _informative_mask(col)
        n_inf = int(mask.sum())
        counts: dict[str, int] = {}
        if n_inf > 0:
            vc = col[mask].astype(str).str.strip().value_counts()
            counts = {str(k): int(v) for k, v in vc.items()}
        out["variables"][var] = {
            "n_informative": n_inf,
            "coverage_fraction": round(n_inf / max(len(work), 1), 4),
            "counts": counts,
        }
    return out


def build_all_demographics_summaries(settings: Settings) -> dict[str, Any]:
    """Lee metadata TSV de cada dataset declarado y devuelve resumen agregado."""
    summaries: dict[str, Any] = {"datasets": {}, "limitations_note": _limitations_note()}
    for ds_name, ds_cfg in (settings.datasets.get("datasets") or {}).items():
        meta_rel = ds_cfg.get("sample_metadata")
        if not meta_rel:
            continue
        meta_path = settings.resolve(str(meta_rel))
        if not meta_path.is_file():
            summaries["datasets"][ds_name] = {"dataset": ds_name, "status": "metadata_missing"}
            continue
        df = pd.read_csv(meta_path, sep="\t", low_memory=False)
        summaries["datasets"][ds_name] = summarize_dataset_metadata(df, dataset=ds_name)
    return summaries


def write_demographics_artifacts(settings: Settings, out_dir: Path | None = None) -> Path:
    """Escribe ``demographics_summary.json`` y ``demographics_counts.tsv`` bajo results/clinical/."""
    out_dir = out_dir or (settings.results_dir() / "clinical")
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = build_all_demographics_summaries(settings)

    json_path = out_dir / "demographics_summary.json"
    json_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    rows: list[dict[str, Any]] = []
    for ds_name, block in summary.get("datasets", {}).items():
        if not isinstance(block, dict) or "variables" not in block:
            continue
        for var, payload in block.get("variables", {}).items():
            for level, cnt in (payload.get("counts") or {}).items():
                rows.append(
                    {
                        "dataset": ds_name,
                        "variable": var,
                        "level": level,
                        "n_samples": cnt,
                        "coverage_fraction": payload.get("coverage_fraction"),
                    }
                )
    tsv_path = out_dir / "demographics_counts.tsv"
    pd.DataFrame(rows).to_csv(tsv_path, sep="\t", index=False)
    return json_path


def _limitations_note() -> str:
    return (
        "CPTAC colon (PDC000116) recluta principalmente en centros de EE.UU. La variable "
        "'hispanic_latino_bin' proviene de la categoría NIH 'Hispanic or Latino' (autorreportada), "
        "NO de nacionalidad colombiana ni de residencia en Latinoamérica. 'geographic_region_bin' "
        "usa país/sitio cuando el Excel clínico lo incluye; si falta, queda 'unknown'. "
        "No sustituye un registro clínico local (p. ej. cohorte CTIC Colombia)."
    )


def render_clinical_demographics_markdown(settings: Settings) -> str:
    """Página MkDocs con tablas de conteo y advertencias metodológicas."""
    summary_path = settings.results_dir() / "clinical" / "demographics_summary.json"
    if summary_path.is_file():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    else:
        write_demographics_artifacts(settings)
        summary = json.loads(summary_path.read_text(encoding="utf-8"))

    lines = [
        "# Demografía y contexto clínico (estadio, etnia, geografía)\n",
        _limitations_note(),
        "\n\n",
        "Esta página resume **qué información clínica está realmente disponible** en cada cohorte "
        "pública integrada. Los conteos se calculan sobre muestras **tumor** (cuando `condition` existe) "
        "y, cuando aplica, se reportan participantes únicos (`participant_id`).\n",
        "\n## Cómo interpretar los estratos\n",
        "\n| Variable | Significado operativo | Uso en el pipeline |\n",
        "|----------|----------------------|-------------------|\n",
        "| `stage` / `stage_bin` | Estadio patológico crudo y binario temprano vs avanzado | Comparación `cptac_stage_advanced_vs_early` |\n",
        "| `stage_group` | Grupo AJCC I–IV extraído del texto de estadio | Tablas y futuras comparaciones por estadio |\n",
        "| `hispanic_latino_bin` | Categoría NIH: hispano/latino vs no (autorreportada, contexto US) | Comparación opcional `cptac_hispanic_latino_vs_non` |\n",
        "| `geographic_region_bin` | País/sitio → `latin_america`, `united_states`, `other`, `unknown` | Exploratorio; cobertura suele ser baja en CPTAC |\n",
        "| `ethnicity` / `race` | Texto normalizado del Excel PDC | Transparencia; no implica origen colombiano |\n",
        "\n",
    ]

    for ds_name, block in summary.get("datasets", {}).items():
        if not isinstance(block, dict):
            continue
        lines.append(f"\n## {ds_name}\n\n")
        if block.get("status") == "metadata_missing":
            lines.append("- Metadata no encontrada en esta corrida.\n")
            continue
        lines.append(
            f"- Muestras tumor en metadata: **{block.get('n_samples_tumor', 0)}** "
            f"(participantes únicos: **{block.get('n_participants_unique', 0)}**)\n"
        )
        vars_block = block.get("variables") or {}
        if not vars_block:
            lines.append("- Sin variables demográficas tabuladas.\n")
            continue
        lines.append("\n| Variable | Cobertura | Distribución (conteo) |\n")
        lines.append("|---|---:|---|\n")
        for var, payload in sorted(vars_block.items()):
            cov = float(payload.get("coverage_fraction", 0.0))
            counts = payload.get("counts") or {}
            if counts:
                dist = ", ".join(f"`{k}`: {v}" for k, v in sorted(counts.items(), key=lambda x: -x[1])[:12])
                if len(counts) > 12:
                    dist += f" … (+{len(counts) - 12} niveles)"
            else:
                dist = "_sin valores informativos_"
            lines.append(f"| `{var}` | {cov:.1%} | {dist} |\n")

    lines.append(
        "\n## Comparaciones opcionales habilitadas por estas variables\n\n"
        "- `cptac_stage_advanced_vs_early` — requiere `stage_bin` con suficientes `early` y `advanced`.\n"
        "- `cptac_hispanic_latino_vs_non` — requiere `hispanic_latino_bin` (derivado de `ethnicity`).\n"
        "- `cptac_eocrc_vs_lorc`, `cptac_msi_vs_mss`, `cptac_cms4_vs_rest` — ver [Cobertura clínica](clinical_coverage.md).\n"
        "\nSi una comparación queda omitida, el log de `differential` indica la columna o el N mínimo.\n"
        "\n## Relevancia para Colombia / Latinoamérica (tesis)\n\n"
        "La revisión de García-Foncillas et al. enfatiza heterogeneidad del CCR y tendencias en población joven "
        "en la región. **Los datos públicos CPTAC no sustituyen** una cohorte colombiana curada: aquí la capa "
        "`hispanic_latino_bin` solo permite explorar, de forma prudente, si participantes autorreportados como "
        "hispanos/latinos en un estudio US muestran perfiles proteómicos distinguibles — como **hipótesis "
        "generadora**, no como generalización a Colombia.\n"
        "\nPara validación regional futura: integrar metadata local (edad, estadio TNM, MSI, origen) en el mismo "
        "esquema de columnas (`metadata/sample_metadata/`) y re-ejecutar `harmonize_public_metadata` + `differential`.\n"
    )
    return "".join(lines)
