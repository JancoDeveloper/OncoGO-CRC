"""Resúmenes cuantitativos listos para tesis / MkDocs a partir de TSV en results/."""

from __future__ import annotations

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv


def _df_to_md_table(df: pd.DataFrame, *, max_rows: int = 15) -> str:
    if df is None or df.empty:
        return "_Sin filas._\n"
    sub = df.head(max_rows).copy()
    sub = sub.astype(object).where(pd.notna(sub), "")
    header = "| " + " | ".join(map(str, sub.columns)) + " |"
    sep = "| " + " | ".join(["---"] * len(sub.columns)) + " |"
    body_rows = ["| " + " | ".join(map(str, row)) + " |" for row in sub.values.tolist()]
    return "\n".join([header, sep, *body_rows]) + "\n"


def _n_significant_proteins(df: pd.DataFrame, alpha: float) -> int:
    if df is None or df.empty:
        return 0
    if "significant" in df.columns:
        return int(df["significant"].fillna(False).sum())
    if "padj" in df.columns:
        return int((df["padj"].astype(float) < alpha).sum())
    return 0


def _top_diff(df: pd.DataFrame, alpha: float, n: int) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    out = df.copy()
    if "padj" not in out.columns:
        return out.head(n)
    out = out.sort_values("padj", na_position="last")
    sig = out[out["padj"].astype(float) < alpha] if "padj" in out.columns else out
    if len(sig) >= 3:
        return sig.head(n)
    return out.head(n)


def _top_go(df: pd.DataFrame, n: int) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    out = df.copy()
    if "padj" in out.columns:
        out = out.sort_values("padj", na_position="last")
        if "significant" in out.columns:
            sig = out[out["significant"] == True]  # noqa: E712
            if not sig.empty:
                return sig.head(n)
    if "p_value" in out.columns:
        return out.sort_values("p_value").head(n)
    return out.head(n)


def build_quantitative_bullets(settings: Settings) -> str:
    """Párrafos cortos para informe ejecutivo (Markdown plano)."""
    alpha = float(settings.config.get("differential", {}).get("alpha", 0.05))
    go_alpha = float(settings.config.get("go", {}).get("fdr_threshold", 0.05))
    res = settings.resolve(settings.config["paths"]["results_dir"])
    comps = settings.comparisons.get("comparisons", {})
    lines: list[str] = []
    for comp_key, spec in comps.items():
        if spec.get("optional"):
            continue
        label = spec.get("label", comp_key)
        diff_path = res / "differential" / f"{comp_key}_diff.tsv"
        go_path = res / "go" / f"{comp_key}_go.tsv"
        if not diff_path.exists():
            lines.append(f"- **{label}** (`{comp_key}`): sin archivo diferencial.\n")
            continue
        diff = read_tsv(diff_path)
        n_sig = _n_significant_proteins(diff, alpha)
        n_tot = len(diff)
        lines.append(
            f"- **{label}** (`{comp_key}`): **{n_sig}** proteínas con FDR por debajo de **{alpha:g}** "
            f"(de **{n_tot}** filas en el contraste).\n"
        )
        if go_path.exists():
            go = read_tsv(go_path)
            if "significant" in go.columns:
                n_go = int(go["significant"].fillna(False).sum())
            elif "padj" in go.columns:
                n_go = int((go["padj"].astype(float) < go_alpha).sum())
            else:
                n_go = 0
            lines.append(f"  - GO (ORA): **{n_go}** términos enriquecidos (FDR por debajo de **{go_alpha:g}**).\n")
    return "".join(lines) if lines else "_Sin comparaciones núcleo ejecutadas._\n"


def build_results_base_markdown(settings: Settings) -> str:
    """Página web: tablas de proteínas y GO por comparación."""
    alpha = float(settings.config.get("differential", {}).get("alpha", 0.05))
    go_alpha = float(settings.config.get("go", {}).get("fdr_threshold", 0.05))
    res = settings.resolve(settings.config["paths"]["results_dir"])
    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")
    comps = settings.comparisons.get("comparisons", {})

    parts: list[str] = [
        "# Resultados base (generado automáticamente)\n\n",
        "Resumen **cuantitativo** del último `oncogo-crc run`: proteínas ordenadas por `padj` y términos GO por enriquecimiento. ",
        "La narrativa biológica y las conclusiones para capítulos de tesis van en ",
        "[**Síntesis para tesis**](thesis_synthesis.md) y en el prompt ",
        "[`docs/PROMPT_SINTESIS_TESIS.md`](../PROMPT_SINTESIS_TESIS.md) para usar con un LLM.\n\n",
        f"Umbral diferencial (FDR): **{alpha:g}**. Umbral GO (ORA): **{go_alpha:g}**.\n\n",
    ]

    for comp_key, spec in comps.items():
        if spec.get("optional"):
            continue
        label = spec.get("label", comp_key)
        diff_path = res / "differential" / f"{comp_key}_diff.tsv"
        go_path = res / "go" / f"{comp_key}_go.tsv"
        go_red = res / "go" / f"{comp_key}_go_reduced.tsv"
        parts.append(f"## {label}\n\n`{comp_key}`\n\n")
        if not diff_path.exists():
            parts.append("_No hay TSV diferencial._\n\n")
            continue
        diff = read_tsv(diff_path)
        n_sig = _n_significant_proteins(diff, alpha)
        parts.append(
            f"**Proteínas significativas** (criterio FDR **{alpha:g}**): **{n_sig}** / **{len(diff)}**\n\n"
        )
        parts.append("### Top filas diferenciales (ordenadas por `padj` ascendente)\n\n")
        topd = _top_diff(diff, alpha, 15)
        cols = [c for c in ["protein_id", "log2fc", "t_stat", "p_value", "padj", "significant"] if c in topd.columns]
        parts.append(_df_to_md_table(topd[cols] if cols else topd, max_rows=15))

        if go_path.exists():
            go = read_tsv(go_path)
            if "significant" in go.columns:
                n_go = int(go["significant"].fillna(False).sum())
            elif "padj" in go.columns:
                n_go = int((go["padj"].astype(float) < go_alpha).sum())
            else:
                n_go = len(go)
            parts.append(f"\n### GO — términos más enriquecidos (FDR **{go_alpha:g}**; **{n_go}** términos pasan el umbral)\n\n")
            top_go = _top_go(go, 15)
            gcols = [c for c in ["go_id", "overlap", "term_size", "odds_ratio", "p_value", "padj"] if c in top_go.columns]
            parts.append(_df_to_md_table(top_go[gcols] if gcols else top_go, max_rows=15))
            if go_red.exists():
                parts.append(
                    f"\n_Reducción semántica:_ `{settings.relative_to_project(go_red)}`\n\n"
                )
        else:
            parts.append("\n_Sin tabla GO (falta diferencial o universo)._\n\n")

    parts.append("## ML — priorización (cohorte primaria)\n\n")
    ml_rank = res / "ml" / f"{primary}_feature_ranking.tsv"
    ml_met = res / "ml" / f"{primary}_model_metrics.tsv"
    if ml_rank.exists():
        rank = read_tsv(ml_rank)
        rcols = [c for c in ["protein_id", "combined_score", "rf_importance", "l1_abs_coef"] if c in rank.columns]
        parts.append(_df_to_md_table(rank[rcols] if rcols else rank, max_rows=15))
    else:
        parts.append("_No se generó ranking ML (revisar muestras tumor/normal en metadata primaria)._\n\n")
    if ml_met.exists():
        parts.append("\n### Métricas test (hold-out)\n\n")
        parts.append(_df_to_md_table(read_tsv(ml_met), max_rows=10))

    parts.append("\n## Traslación clínica (extracto)\n\n")
    clin = res / "clinical_translation" / f"{primary}_translation_ranked.tsv"
    if clin.exists():
        cdf = read_tsv(clin)
        parts.append(_df_to_md_table(cdf.head(15), max_rows=15))
    else:
        parts.append("_No disponible._\n")

    parts.append("\n## Red — hubs (extracto)\n\n")
    hubs = res / "networks" / f"{primary}_hub_proteins.tsv"
    if hubs.exists():
        parts.append(_df_to_md_table(read_tsv(hubs).head(15), max_rows=15))
    else:
        parts.append("_No disponible._\n")

    return "".join(parts)


def build_thesis_synthesis_markdown(settings: Settings) -> str:
    """Guía explícita de qué dice el pipeline y qué falta redactar."""
    bullets = build_quantitative_bullets(settings)
    return (
        "# Síntesis para tesis\n\n"
        "## Qué hace este pipeline (en una frase)\n\n"
        "Integra matrices proteómicas públicas (CPTAC, validación PRIDE, plasma), "
        "normaliza e imputa, contrasta **tumor vs normal** (y otras comparaciones núcleo), "
        "enriquece genes significativos en **Gene Ontology (ORA offline)**, resume redes y "
        "prioriza proteínas con modelos interpretables.\n\n"
        "## Qué está **demostrado** aquí (evidencia cuantitativa)\n\n"
        f"{bullets}\n"
        "---\n\n"
        "## Contexto clínico‑epidemiológico (Colombia) y subtipos CMS\n\n"
        "La revisión **García-Foncillas *et al.*, 2026** (*Rev. colomb. hematol. oncol.*, DOI [10.51643/22562915.853](https://doi.org/10.51643/22562915.853)) "
        "subraya el aumento de CCR en población joven y el papel de clasificaciones moleculares (**CMS**, MSI, MAPK/EGFR, biopsia líquida). "
        "OncoGO-CRC incorpora:\n\n"
        "- **Edad / EOCRC:** si el Excel clínico PDC trae edad, la metadata CPTAC incluye `age`, `age_bin` (`young`/`older` según `clinical_filters.age_bins`) "
        "y la comparación opcional `cptac_eocrc_vs_lorc`.\n"
        "- **MSI vs MSS:** comparación opcional `cptac_msi_vs_mss` cuando `msi_status` se rellena desde clínica.\n"
        "- **CMS4 vs resto:** columna `cms_subtype_bin` si aportas un TSV suplementario (`harmonize_public_metadata.cptac.cms_subtype_supplement_tsv`) "
        "y la comparación opcional `cptac_cms4_vs_rest`.\n"
        "- **Traslación:** `metadata/clinical_layers/crc_context_layers.yaml` etiqueta genes priorizados (biomarcadores, pistas CMS, eje MAPK/EGFR para discusión de plasma).\n\n"
        "Interpretación sugerida: hubs en matriz extracelular / EMT pueden discutirse como **coherente con fenotipo CMS4** (TGF‑β, angiogénesis) "
        "si el contraste y la literatura lo sustentan — citar siempre el subtipo molecular **medido o inferido con criterio explícito** en tesis.\n\n"
        "---\n\n"
        "## Qué **no** sustituye este documento\n\n"
        "- No interpreta mecanismos ni discute contraste con la bibliografía clínica: eso es texto de tesis.\n"
        "- Los avisos de SciPy («catastrophic cancellation») afectan filas con varianza casi nula; "
        "no invalidan todo el contraste.\n\n"
        "## Cómo cerrar el capítulo de resultados\n\n"
        "1. Lee [Resultados base](results_base.md) y valida que las tablas coinciden con tus preguntas de investigación.\n"
        "2. Abre los TSV completos en `results/` (misma jerarquía bajo tu `data_workspace` si lo usas).\n"
        "3. Usa el prompt en [`docs/PROMPT_SINTESIS_TESIS.md`](../PROMPT_SINTESIS_TESIS.md) "
        "con un LLM pegando extractos de tablas y el contexto de la cohorte.\n\n"
        "## Archivos clave para citar en la memoria\n\n"
        "- `results/differential/*_diff.tsv` — estadísticos por proteína.\n"
        "- `results/go/*_go.tsv` — enriquecimiento por término.\n"
        "- `results/ml/<cohorte>_feature_ranking.tsv` — priorización multivariante.\n"
        "- `results/validation/input_validation_report.json` — trazabilidad de calidad.\n"
    )


def build_limitations_markdown(settings: Settings) -> str:
    """Limitaciones explícitas enlazadas a config."""
    dc = settings.config.get("differential", {})
    go = settings.config.get("go", {})
    ml = settings.config.get("ml", {})
    return (
        "# Limitaciones (auto)\n\n"
        "- **Modelo estadístico:** t-test por proteína + FDR Benjamini–Hochberg; asume razonablemente "
        "muestreo independiente entre grupos; no sustituye modelos mixtos si hay batches fuertes no corregidos.\n"
        "- **GO:** ORA con universo = proteínas detectadas en la matriz; términos muy grandes o redundantes "
        f"pueden dominar; FDR término a {go.get('fdr_threshold', 0.05)}.\n"
        "- **ML:** hold-out simple (`test_size={ml.get('test_size', 0.25)}`); métricas en conjunto de prueba, "
        "no validación externa salvo que añadas cohortes.\n"
        "- **Datos públicos:** versiones de matrices y metadata dependen de PDC/PRIDE en la fecha de descarga; "
        "documenta `data/versions/` y manifiestos.\n"
        f"- **α diferencial:** {dc.get('alpha', 0.05)} — ajustar solo con criterio metodológico defendible en tesis.\n"
        "- **Tejido vs plasma:** el contraste tumor/normal en tejido suele separar mejor que plasma circulante; "
        "AUC más altos en CPTAC frente a capas líquidas débiles es **esperable** y defendible frente a revisiones que citan AUC modestos en cohortes sanguíneas.\n"
        "- **CMS / KEGG:** CMS suele provenir de transcriptómica; aquí solo se usa si aportas suplemento explícito. "
        "KEGG (`kegg_enrichment`, gseapy) es opcional y depende de conectividad Enrichr.\n"
    )


def build_results_inventory_markdown(settings: Settings) -> str:
    """Lista compacta de artefactos bajo results/."""
    res = settings.resolve(settings.config["paths"]["results_dir"])
    lines = ["# Inventario de archivos en `results/`\n\n"]
    if not res.is_dir():
        lines.append("_No existe el directorio results._\n")
        return "".join(lines)
    for sub in sorted([p for p in res.iterdir() if p.is_dir()]):
        files = sorted(sub.glob("*.tsv")) + sorted(sub.glob("*.json")) + sorted(sub.glob("*.png"))
        if not files:
            continue
        lines.append(f"## `{sub.name}/`\n\n")
        for f in files[:40]:
            rel = settings.relative_to_project(f)
            try:
                kb = max(1, f.stat().st_size // 1024)
            except OSError:
                kb = "?"
            lines.append(f"- `{rel}` (~{kb} KiB)\n")
        if len(files) > 40:
            lines.append(f"- _… y {len(files) - 40} archivos más_\n")
        lines.append("\n")
    return "".join(lines)


def build_web_index_markdown(settings: Settings) -> str:
    """Portada del bloque Proyecto Web con enlaces accionables."""
    return (
        "# Proyecto web — OncoGO-CRC\n\n"
        "Navegación sugerida para **revisión de tesis**:\n\n"
        "1. [**Síntesis para tesis**](thesis_synthesis.md) — qué demuestra el pipeline y qué redactar aparte.\n"
        "2. [**Resultados base**](results_base.md) — tablas de **proteínas**, **GO** y extractos ML/red/traslación.\n"
        "3. [**Inventario**](results_inventory.md) — lista de archivos generados.\n"
        "4. [Datasets](datasets.md) · [Cobertura clínica](clinical_coverage.md) · [Estabilidad](stability.md) · "
        "[Limitaciones](limitations.md)\n\n"
        "Prompt para redacción interpretativa con LLM: "
        "[`docs/PROMPT_SINTESIS_TESIS.md`](../PROMPT_SINTESIS_TESIS.md).\n"
    )
