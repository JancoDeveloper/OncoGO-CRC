"""Resúmenes cuantitativos listos para tesis / MkDocs a partir de TSV en results/."""

from __future__ import annotations

import json

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv
from oncogo_crc.utils.run_manifest import gsea_comparison_trustworthy, ml_transfer_trustworthy


def _workspace_doc_note(settings: Settings) -> str:
    """Párrafo Markdown opcional: deja explícito el volumen externo para informes generados."""
    if settings.data_workspace is None:
        return ""
    w = settings.data_workspace.resolve()
    return (
        "> **Workspace de datos:** `data/`, `results/`, `reports/`, `metadata/` y `docs/web/` "
        f"apuntan a `{w}` (el clon conserva `config/`, `src/` y el resto de `docs/`). "
        "Trazabilidad de corrida: `results/run_manifest.json`.\n\n"
    )


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


def _gsea_fdr_column(df: pd.DataFrame) -> str | None:
    for c in ("FDR q-val", "FDR", "fdr"):
        if c in df.columns:
            return c
    return None


def _n_sig_gsea(df: pd.DataFrame, fdr: float) -> int:
    col = _gsea_fdr_column(df)
    if not col:
        return 0
    return int((df[col].astype(float) < fdr).sum())


def _top_gsea(df: pd.DataFrame, fdr: float, n: int) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    col = _gsea_fdr_column(df)
    if not col:
        return df.head(n)
    sub = df.loc[df[col].astype(float) < fdr].copy()
    if sub.empty:
        sub = df.copy()
    if "NES" in sub.columns:
        sub = sub.iloc[sub["NES"].abs().argsort()[::-1]]
    return sub.head(n)


def _gsea_ranking_metrics_list(gsea_cfg: dict) -> list[str]:
    raw = gsea_cfg.get("ranking_metrics")
    if isinstance(raw, list) and raw:
        return [str(x).strip() for x in raw if str(x).strip()]
    return [str(gsea_cfg.get("ranking_metric", "t_stat")).strip()]


def build_quantitative_bullets(settings: Settings) -> str:
    """Párrafos cortos para informe ejecutivo (Markdown plano)."""
    alpha = float(settings.config.get("differential", {}).get("alpha", 0.05))
    go_alpha = float(settings.config.get("go", {}).get("fdr_threshold", 0.05))
    gsea_cfg = settings.config.get("gsea", {}) or {}
    gsea_fdr = float(gsea_cfg.get("fdr_threshold", 0.25))
    rk_metrics = _gsea_ranking_metrics_list(gsea_cfg)
    res = settings.resolve(settings.config["paths"]["results_dir"])
    comps = settings.comparisons.get("comparisons", {})
    lines: list[str] = []

    def _lines_for_comparison(comp_key: str, label: str) -> None:
        diff_path = res / "differential" / f"{comp_key}_diff.tsv"
        go_path = res / "go" / f"{comp_key}_go.tsv"
        gsea_path = res / "gsea" / f"{comp_key}_gsea_go.tsv"
        if not diff_path.exists():
            lines.append(f"- **{label}** (`{comp_key}`): sin archivo diferencial.\n")
            return
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
            lines.append(
                f"  - GO (ORA, universo = proteínas detectadas en la matriz): **{n_go}** términos "
                f"(FDR término < **{go_alpha:g}**).\n"
            )
        else:
            lines.append("  - GO (ORA): sin tabla.\n")
        if gsea_comparison_trustworthy(settings, comp_key) and gsea_path.exists():
            gdf = read_tsv(gsea_path)
            n_ge = _n_sig_gsea(gdf, gsea_fdr)
            lines.append(
                f"  - GSEA (GO offline, ranking continuo; métrica principal `{rk_metrics[0]}`): "
                f"**{n_ge}** términos con FDR < **{gsea_fdr:g}**.\n"
            )
            for m in rk_metrics[1:]:
                safe = m.replace(" ", "_")
                alt = res / "gsea" / f"{comp_key}_gsea_go_{safe}.tsv"
                if alt.exists() and gsea_comparison_trustworthy(settings, comp_key):
                    gdf2 = read_tsv(alt)
                    n2 = _n_sig_gsea(gdf2, gsea_fdr)
                    lines.append(
                        f"  - GSEA (métrica `{m}`): **{n2}** términos con FDR < **{gsea_fdr:g}** "
                        f"(`results/gsea/{alt.name}`).\n"
                    )
        elif gsea_path.exists() and not gsea_comparison_trustworthy(settings, comp_key):
            lines.append(
                "  - GSEA: hay archivos en `results/gsea/` pero **no son de la corrida actual** "
                "(o esta comparación se omitió). No se muestran cifras; ver `results/run_manifest.json`.\n"
            )
        else:
            lines.append(
                "  - GSEA: sin tabla (etapa `gsea` o `pip install -e \".[go-remote]\"` si falta gseapy).\n"
            )

    for comp_key, spec in comps.items():
        if spec.get("optional"):
            continue
        _lines_for_comparison(comp_key, str(spec.get("label", comp_key)))

    opt_written = False
    for comp_key, spec in comps.items():
        if not spec.get("optional"):
            continue
        diff_path = res / "differential" / f"{comp_key}_diff.tsv"
        if not diff_path.exists():
            continue
        if not opt_written:
            lines.append("\n**Comparaciones opcionales** (solo si `differential.auto_run_optional_comparisons` y datos):\n\n")
            opt_written = True
        _lines_for_comparison(comp_key, str(spec.get("label", comp_key)))

    rep = res / "reproducibility" / "cptac_pride_differential_overlap_summary.tsv"
    if rep.exists():
        odf = read_tsv(rep)
        if not odf.empty and "metric" in odf.columns and "value" in odf.columns:
            kv = dict(zip(odf["metric"].astype(str), odf["value"]))
            n_core = None
            for k, v in kv.items():
                if str(k).startswith("n_both_significant"):
                    n_core = v
                    break
            n_bridge = kv.get("n_intersection_uniprot_bridge", kv.get("n_intersection_ids"))
            n_raw = kv.get("n_intersection_raw_protein_id", "")
            pr = kv.get("log2fc_pearson_r_cptac_pride")
            ss = kv.get("log2fc_same_sign_fraction")
            corr_txt = ""
            if pr is not None and ss is not None:
                try:
                    pvf = float(pr)
                    ssf = float(ss)
                    if pvf == pvf and ssf == ssf:
                        corr_txt = (
                            f" Correlación Pearson log2FC CPTAC–PRIDE (UniProt común): **{pvf:.3f}**; "
                            f"mismo signo: **{ssf:.0%}**."
                        )
                except (TypeError, ValueError):
                    pass
            lines.append(
                f"\n**Firma núcleo CPTAC ∩ PRIDE** (clave común **UniProt** vía mapa `gene_symbol`↔`uniprot_accession`; "
                f"filas con UniProt en ambas cohortes: **{n_bridge}**; doble significativo `padj`: **{n_core}**). "
                f"Cruce ingenuo por el mismo texto en `protein_id` suele ser **{n_raw}** (CPTAC=genes, PRIDE=UniProt).{corr_txt} "
                f"Tablas: `results/reproducibility/cptac_pride_differential_overlap.tsv` y `..._summary.tsv`.\n"
            )

    return "".join(lines) if lines else "_Sin comparaciones núcleo ejecutadas._\n"


def build_results_base_markdown(settings: Settings) -> str:
    """Página web: tablas de proteínas, ORA, GSEA y firma multi-cohorte."""
    alpha = float(settings.config.get("differential", {}).get("alpha", 0.05))
    go_alpha = float(settings.config.get("go", {}).get("fdr_threshold", 0.05))
    dc = settings.config.get("differential", {})
    go_cfg = settings.config.get("go", {})
    gsea_cfg = settings.config.get("gsea", {}) or {}
    gsea_fdr = float(gsea_cfg.get("fdr_threshold", 0.25))
    min_fc = float(dc.get("min_abs_log2fc", 0) or 0)
    ora_fc = float(go_cfg.get("ora_min_abs_log2fc", 0) or 0)
    rk_metrics = _gsea_ranking_metrics_list(gsea_cfg)
    rk_display = ", ".join(rk_metrics)
    res = settings.resolve(settings.config["paths"]["results_dir"])
    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")
    comps = settings.comparisons.get("comparisons", {})

    parts: list[str] = [
        "# Resultados base (generado automáticamente)\n\n",
        _workspace_doc_note(settings),
        "Resumen **cuantitativo** del último `oncogo-crc run`: proteínas por `padj`, **ORA** (lista binaria), "
        "**GSEA pre-ranked** (ranking continuo, mejor tolerancia cuando casi todo el proteoma cambia), "
        "firma **CPTAC ∩ PRIDE** y extractos ML/red/traslación. La narrativa interpretativa va en "
        "[**Síntesis para tesis**](thesis_synthesis.md).\n\n",
        "**Umbrales en esta corrida:** FDR diferencial **{:.3g}**; columna `significant` incluye |log2fc| ≥ **{:.3g}** "
        "si `differential.min_abs_log2fc` > 0; ORA puede exigir además |log2fc| ≥ **{:.3g}** (`go.ora_min_abs_log2fc`); "
        "GSEA usa métrica(s) de ranking **{}** y FDR término **{:.3g}** (exploratorio). "
        "Cada métrica adicional genera `*_gsea_go_<métrica>.tsv`; el primero se refleja también en `*_gsea_go.tsv`. "
        "**Universo ORA/GSEA:** genes mapeados desde la matriz normalizada ∩ `metadata/gene_ontology/go_gene_long.tsv`.\n\n".format(
            alpha, min_fc, ora_fc, rk_display, gsea_fdr
        ),
    ]

    def append_one_comparison(comp_key: str, label: str, *, heading_level: str = "##") -> None:
        diff_path = res / "differential" / f"{comp_key}_diff.tsv"
        go_path = res / "go" / f"{comp_key}_go.tsv"
        go_red = res / "go" / f"{comp_key}_go_reduced.tsv"
        gsea_path = res / "gsea" / f"{comp_key}_gsea_go.tsv"
        parts.append(f"{heading_level} {label}\n\n`{comp_key}`\n\n")
        if not diff_path.exists():
            parts.append("_No hay TSV diferencial._\n\n")
            return
        diff = read_tsv(diff_path)
        n_sig = _n_significant_proteins(diff, alpha)
        parts.append(
            f"**Proteínas significativas** (FDR **{alpha:g}**; ver `significant` y filtros en config): "
            f"**{n_sig}** / **{len(diff)}**\n\n"
        )
        parts.append("### Top filas diferenciales (por `padj`)\n\n")
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
            parts.append(
                f"\n### GO — ORA offline (FDR término **{go_alpha:g}**; **{n_go}** términos pasan)\n\n"
            )
            top_go = _top_go(go, 15)
            gcols = [c for c in ["go_id", "overlap", "term_size", "odds_ratio", "p_value", "padj"] if c in top_go.columns]
            parts.append(_df_to_md_table(top_go[gcols] if gcols else top_go, max_rows=15))
            if go_red.exists():
                parts.append(
                    f"\n_Reducción semántica:_ `{settings.relative_to_project(go_red)}`\n\n"
                )
        else:
            parts.append("\n_Sin tabla GO._\n\n")

        if gsea_comparison_trustworthy(settings, comp_key) and gsea_path.exists():
            gdf = read_tsv(gsea_path)
            n_ge = _n_sig_gsea(gdf, gsea_fdr)
            parts.append(
                f"\n### GSEA — GO offline pre-ranked (FDR **{gsea_fdr:g}**; **{n_ge}** términos; "
                f"ranking principal `{rk_metrics[0]}`)\n\n"
            )
            tcols = [c for c in ["Term", "NES", "NOM p-val", "FDR q-val", "Tag %", "Lead_genes"] if c in gdf.columns]
            topg = _top_gsea(gdf, gsea_fdr, 12)
            parts.append(_df_to_md_table(topg[tcols] if tcols else topg, max_rows=12))
            parts.append(
                f"\n_Fichero completo:_ `{settings.relative_to_project(gsea_path)}`\n\n"
            )
            for m in rk_metrics[1:]:
                safe = m.replace(" ", "_")
                alt = res / "gsea" / f"{comp_key}_gsea_go_{safe}.tsv"
                if not alt.exists() or not gsea_comparison_trustworthy(settings, comp_key):
                    continue
                gdf2 = read_tsv(alt)
                n2 = _n_sig_gsea(gdf2, gsea_fdr)
                parts.append(
                    f"\n#### GSEA — métrica alternativa `{m}` (FDR **{gsea_fdr:g}**; **{n2}** términos)\n\n"
                )
                top2 = _top_gsea(gdf2, gsea_fdr, 12)
                parts.append(_df_to_md_table(top2[tcols] if tcols else top2, max_rows=12))
                parts.append(
                    f"\n_Fichero completo:_ `{settings.relative_to_project(alt)}`\n\n"
                )
        elif gsea_path.exists() and not gsea_comparison_trustworthy(settings, comp_key):
            parts.append(
                "\n### GSEA\n\n"
                "_Tablas presentes en `results/gsea/` pero **no verificadas para esta corrida** "
                "(GSEA omitido o comparación sin salida). No se incrustan tablas antiguas. "
                "Ver `results/run_manifest.json`._\n\n"
            )
        else:
            parts.append(
                "\n_Sin GSEA (etapa `gsea`, requiere `gseapy` en `pip install -e \".[go-remote]\"`)._\n\n"
            )

    for comp_key, spec in comps.items():
        if spec.get("optional"):
            continue
        append_one_comparison(comp_key, str(spec.get("label", comp_key)))

    parts.append("## Comparaciones opcionales ejecutadas\n\n")
    opt_any = False
    for comp_key, spec in comps.items():
        if not spec.get("optional"):
            continue
        if not (res / "differential" / f"{comp_key}_diff.tsv").exists():
            continue
        opt_any = True
        append_one_comparison(comp_key, str(spec.get("label", comp_key)), heading_level="###")
    if not opt_any:
        parts.append("_Ninguna comparación opcional generó `*_diff.tsv` en esta corrida._\n\n")

    parts.append("## Firma núcleo CPTAC ∩ PRIDE (reproducibilidad multi-cohorte)\n\n")
    rep_sum = res / "reproducibility" / "cptac_pride_differential_overlap_summary.tsv"
    rep_full = res / "reproducibility" / "cptac_pride_differential_overlap.tsv"
    if rep_sum.exists():
        parts.append(
            "Cruce **por UniProt canónico** (primera accesión del grupo MaxQuant; isoformas `-n` colapsadas) "
            "usando el mapa maestro `metadata/protein_id_maps/protein_mapping_master.tsv` para traducir símbolos "
            "CPTAC → UniProt. Sin mapa amplio, `n_cptac_with_join_uni` será bajo y la firma núcleo se subestima. "
            "Métricas `n_intersection_raw_protein_id` vs `n_intersection_uniprot_bridge` diagnostican el problema "
            "histórico «0 intersección» por IDs incompatibles.\n\n"
        )
        parts.append(_df_to_md_table(read_tsv(rep_sum), max_rows=20))
        parts.append(
            f"\n_Detalle ancho:_ `{settings.relative_to_project(rep_full)}` (columna `_join_uni`, etapa `cohort_overlap`).\n\n"
        )
    else:
        parts.append(
            "_Aún no generado._ Ejecuta `cohort_overlap` tras tener diferenciales CPTAC y PRIDE.\n\n"
        )

    parts.append("## ML — priorización (cohorte primaria)\n\n")
    ml_rank = res / "ml" / f"{primary}_feature_ranking.tsv"
    ml_met = res / "ml" / f"{primary}_model_metrics.tsv"
    if ml_rank.exists():
        rank = read_tsv(ml_rank)
        rcols = [
            c
            for c in [
                "protein_id",
                "combined_score",
                "rf_importance",
                "l1_abs_coef",
                "elasticnet_abs_coef",
            ]
            if c in rank.columns
        ]
        parts.append(_df_to_md_table(rank[rcols] if rcols else rank, max_rows=15))
    else:
        parts.append("_No se generó ranking ML (revisar muestras tumor/normal en metadata primaria)._\n\n")
    if ml_met.exists():
        parts.append("\n### Métricas test (hold-out)\n\n")
        parts.append(_df_to_md_table(read_tsv(ml_met), max_rows=16))

    nest = res / "ml" / f"{primary}_nested_cv_per_fold.tsv"
    if nest.exists():
        parts.append(
            "\n### Nested cross-validation (sesgo reducido frente a un solo hold-out)\n\n"
            "Pliegues **externos** estiman AUC ROC; en cada entrenamiento externo, **GridSearchCV** con pliegues "
            "**internos** elige hiperparámetros. Con `ml.nested_cv.only_rf: true` solo aparece **random_forest_nested** "
            "(más rápido); con `only_rf: false` se incluyen también L1 y Elastic Net. Los pliegues se **recortan** "
            "automáticamente si la clase minoritaria es pequeña (`split_adjustment` en el JSON de resumen). "
            "Ver `ml.nested_cv` en `config.yaml`.\n\n"
        )
        parts.append(_df_to_md_table(read_tsv(nest), max_rows=25))
        parts.append(
            f"\n_Resumen numérico:_ `{settings.relative_to_project(res / 'ml' / f'{primary}_nested_cv_summary.tsv')}`.\n\n"
        )

    xfer = list(res.glob("ml/*_train_*_test_transfer.json"))
    if ml_transfer_trustworthy(settings) and xfer:
        parts.append("\n### Transferencia entre cohortes (entrenar en CPTAC, evaluar en PRIDE)\n\n")
        parts.append(
            "Genes alineados por símbolo post-`harmonize` (filas colapsadas a `gene_symbol`). "
            "Las métricas en PRIDE pueden ser **más modestas** que el hold-out dentro de CPTAC: es un resultado "
            "válido sobre **dominio distinto** y batch/procedimiento distinto. "
            "Si `ml.external_transfer.calibration_method` es `sigmoid` o `isotonic`, se listan también bloques "
            "`*_calibrated` (calibración con CV **solo sobre CPTAC** antes de proyectar a PRIDE).\n\n"
        )
        for jf in sorted(xfer)[:3]:
            payload = json.loads(jf.read_text(encoding="utf-8"))
            parts.append(f"**`{jf.name}`** — estado: `{payload.get('status', '')}`; genes: **{payload.get('n_genes_aligned', '—')}**; "
            f"muestras train/test: **{payload.get('n_train_samples', '—')}** / **{payload.get('n_test_samples', '—')}**.\n\n")
            cal = payload.get("calibration") or {}
            if cal.get("applied"):
                parts.append(
                    f"- _Calibración_: **{cal.get('method')}**, `calibration_cv`={cal.get('calibration_cv')} (ajuste en CPTAC).\n"
                )
            for model in ("sparse_linear_l1", "elastic_net_sgd", "random_forest"):
                block = payload.get(model)
                if isinstance(block, dict):
                    auc_v = block.get("roc_auc")
                    acc_v = block.get("accuracy")
                    try:
                        auc_s = f"{float(auc_v):.4f}"
                    except (TypeError, ValueError):
                        auc_s = "—"
                    try:
                        acc_s = f"{float(acc_v):.4f}"
                    except (TypeError, ValueError):
                        acc_s = "—"
                    parts.append(
                        f"- _{model}_ (sin calibrar): ROC-AUC **{auc_s}**, accuracy **{acc_s}**\n"
                    )
            for model in ("sparse_linear_l1_calibrated", "elastic_net_sgd_calibrated", "random_forest_calibrated"):
                block = payload.get(model)
                if isinstance(block, dict):
                    auc_v = block.get("roc_auc")
                    acc_v = block.get("accuracy")
                    try:
                        auc_s = f"{float(auc_v):.4f}"
                    except (TypeError, ValueError):
                        auc_s = "—"
                    try:
                        acc_s = f"{float(acc_v):.4f}"
                    except (TypeError, ValueError):
                        acc_s = "—"
                    parts.append(
                        f"- _{model}_: ROC-AUC **{auc_s}**, accuracy **{acc_s}**\n"
                    )
            parts.append("\n")

    elif xfer and not ml_transfer_trustworthy(settings):
        parts.append("\n### Transferencia entre cohortes (entrenar en CPTAC, evaluar en PRIDE)\n\n")
        parts.append(
            "_Hay `*_transfer.json` en `results/ml/` pero **no corresponden a la corrida actual** "
            "(transferencia omitida o etapa `ml` sin salida verificada). Ver `results/run_manifest.json` "
            "(`stages.ml_transfer`)._\n\n"
        )

    boot = res / "ml" / f"{primary}_feature_stability_bootstrap.tsv"
    if boot.exists():
        parts.append("\n### Estabilidad bootstrap de features (L1)\n\n")
        parts.append(_df_to_md_table(read_tsv(boot), max_rows=20))
        parts.append(
            "\n_Frecuencia en la que cada gen entra en el top‑k por |coef| en réplicas bootstrap._ "
            "Config: `ml.feature_stability_bootstrap_n`.\n\n"
        )

    parts.append("\n## Traslación clínica (extracto)\n\n")
    clin = res / "clinical_translation" / f"{primary}_translation_ranked.tsv"
    if clin.exists():
        cdf = read_tsv(clin)
        pref = [
            "feature_id",
            "translational_priority_score",
            "evidence_tier",
            "functional_module_hints",
            "cms_pathway_hints",
            "log2fc",
            "padj",
            "combined_score",
            "evidence_score",
        ]
        tcols = [c for c in pref if c in cdf.columns]
        parts.append(_df_to_md_table(cdf[tcols].head(15) if tcols else cdf.head(15), max_rows=15))
    else:
        parts.append("_No disponible._\n")

    parts.append("\n## Red — hubs (extracto)\n\n")
    hubs = res / "networks" / f"{primary}_hub_proteins.tsv"
    if hubs.exists():
        parts.append(_df_to_md_table(read_tsv(hubs).head(15), max_rows=15))
    else:
        parts.append("_No disponible._\n")

    parts.append("\n## Export contextual (ProteinsLanguage / STRING / AlphaFold)\n\n")
    bundle = res / "export" / "protein_context_bundle.json"
    if bundle.exists():
        parts.append(
            f"Bundle JSON generado por la etapa `export_bundle`: `{settings.relative_to_project(bundle)}`. "
            "Incluye `external_links` (STRING, UniProt, AlphaFold), vista previa de *leading edges* GSEA si existen, "
            "métricas de reproducibilidad CPTAC∩PRIDE y, si activas `export_bundle.uniprot_top_n`, resúmenes UniProt cacheados.\n\n"
        )
    else:
        parts.append(
            "_Sin bundle._ Ejecuta `export_bundle` (o `python scripts/export_protein_context_bundle.py`).\n\n"
        )

    return "".join(parts)


def build_thesis_synthesis_markdown(settings: Settings) -> str:
    """Guía explícita de qué dice el pipeline y qué falta redactar."""
    bullets = build_quantitative_bullets(settings)
    return (
        "# Síntesis para tesis\n\n"
        f"{_workspace_doc_note(settings)}"
        "## Qué hace este pipeline (en una frase)\n\n"
        "Integra matrices proteómicas públicas (CPTAC, validación PRIDE, plasma), normaliza e imputa, "
        "contrasta grupos (tumor vs normal y comparaciones opcionales MSI/EOCRC si hay datos), "
        "interpreta con **ORA + GSEA pre-ranked** sobre Gene Ontology offline (mismo `go_gene_long.tsv`), "
        "cuantifica **firma núcleo CPTAC ∩ PRIDE**, resume redes, prioriza con ML explicable y exporta un "
        "**bundle** listo para capas contextuales (STRING, UniProt, AlphaFold, ProteinsLanguage).\n\n"
        "## Arquitectura conceptual (3 niveles)\n\n"
        "1. **Nivel 1 — OncoGO-CRC:** pipeline estadístico reproducible (diferencial, ORA, GSEA, KEGG opcional, ML, redes, estabilidad).\n"
        "2. **Nivel 2 — ProteinsLanguage / recursos externos:** GO profundo, estructura, embeddings, grafos, visualización "
        "(consume `results/export/protein_context_bundle.json` + tablas en `results/`).\n"
        "3. **Nivel 3 — Traslación clínica:** CMS / MSI / ECM / plasma / EOCRC vía metadata y `crc_context_layers.yaml` "
        "(interpretación siempre acotada a lo medido o inferido con criterio explícito).\n\n"
        "## Enfoque metodológico: ORA, GSEA y señal distribuida\n\n"
        "En proteómica, **listas muy grandes o heterogéneas** pueden hacer que la **ORA (sobrerrepresentación binaria)** "
        "sea **poco informativa** (pocos o ningún término GO con FDR estricto), **sin que ello invalide** el contraste diferencial. "
        "Por ello el proyecto **prioriza también una estrategia basada en ranking continuo (GSEA pre-ranked)** y en "
        "**reducción semántica** de términos GO (`*_go_reduced.tsv`), para capturar **cambios distribuidos a nivel de proceso**. "
        "Se reportan **ORA y GSEA juntos**: la ORA describe conjuntos discretos bajo reglas explícitas; el GSEA explora "
        "si la ranking completa enriquece vías cuando la ORA colapsa. En la memoria se puede enmarcar un **0 términos ORA** "
        "como **decisión metodológica defendible** (lista no específica / universo próximo al foreground), no como “fallo”.\n\n"
        "## Preguntas de investigación: módulos funcionales, no solo “proteína top”\n\n"
        "Además de identificar proteínas punteras, el marco sugiere preguntar **qué procesos o módulos son consistentes** "
        "entre cohortes o con la literatura (p. ej. EMT/ECM, metabolismo lipídico/glicolítico, MAPK/EGFR, inmunidad/MSI-like, "
        "reparación del ADN, angiogénesis/hipoxia). Las etiquetas `functional_module_hints` y `cms_pathway_hints` en "
        "`results/clinical_translation/*.tsv` orientan esa lectura; la síntesis clínica debe citar **García-Foncillas *et al.*, 2026** "
        "y tu revisión narrativa donde corresponda.\n\n"
        "## Validación externa y transferencia entre cohortes\n\n"
        "Diseño explícito: **descubrimiento en CPTAC**, **validación de dirección en PRIDE**, **apoyo traslacional en plasma**. "
        "El resumen `results/reproducibility/cptac_pride_differential_overlap_summary.tsv` incluye **correlación de `log2fc`** "
        "(Pearson/Spearman) y **fracción de mismo signo** entre cohortes unidas por **UniProt canónico**. "
        "Con `ml.external_transfer.enabled`, el pipeline escribe **`results/ml/*_train_*_test_transfer.json`**: AUC/accuracy en PRIDE "
        "tras alinear **genes** (post-`harmonize`) y entrenar en CPTAC. "
        "Un **rendimiento ML alto dentro de CPTAC** pero **transferencia modesta a PRIDE** es defendible científicamente: "
        "subraya batch/tecnología y la necesidad de **priorizar señales biológicas robustas** frente a clasificadores sobreajustados.\n\n"
        "## ML: ranking estable y explicable (no “el mejor AUC”)\n\n"
        "Los modelos **L1**, **Elastic Net (SGD)** y **RF** con métricas de hold-out sirven para **ordenar genes con evidencia multivariante**, no para "
        "vender un clasificador cerrado. Con `ml.nested_cv.enabled`, se generan **`results/ml/<cohorte>_nested_cv_per_fold.tsv`** "
        "(AUC por pliegue externo; `GridSearchCV` interno para `alpha` o hiperparámetros de RF). Opcional: **`ml.feature_stability_bootstrap_n`** "
        "para frecuencias bootstrap del top‑k en L1. La narrativa honesta favorece **nested CV + transferencia** frente a un único split optimista.\n\n"
        "## Translational Priority Score (TPS)\n\n"
        "La tabla de traslación incluye **`translational_priority_score`**: combinación explícita de magnitud/FDR, red, ML, "
        "contexto CMS/MAPK, módulos funcionales y capas literarias (`crc_context_layers.yaml`). Sirve como **marco cuantitativo** "
        "para priorizar hipótesis, no como decisión clínica.\n\n"
        "## Druggability e hipótesis de reposicionamiento (in silico)\n\n"
        "Sin experimentos de laboratorio, el **bundle** JSON incluye plantillas de consulta (p. ej. Open Targets / ChEMBL) "
        "para **hipótesis de diana o reposicionamiento**. La redacción debe evitar afirmaciones tipo “este fármaco sirve en CRC”; "
        "el formulación correcta es **“priorización de diana / hipótesis generadora”**, alineada con revisiones de farmacoproteómica.\n\n"
        "## Validaciones in silico de bajo costo (opcionales)\n\n"
        "La etapa **`external_evidence`** (tras `clinical`) escribe **`results/external_evidence/evidence_matrix.tsv`** "
        "(columnas `hit_*` desde `crc_context_layers.yaml` + fusión con traslación). "
        "Si configuras `external_evidence.tcga_rnaseq_log2fc_tsv` o `survival_gene_cox_tsv` (plantillas en `metadata/tcga_auxiliary/`), "
        "se generan **`tcga_rnaseq_concordance.tsv`** y **`survival_gene_cox_merged.tsv`** con Spearman/mismo signo o Cox fusionado. "
        "Rutas adicionales típicas: soporte en **TCGA/GEO**, **supervivencia** pública, **meta-análisis direccional** (≥2 cohortes).\n\n"
        "## Qué está **demostrado** aquí (evidencia cuantitativa)\n\n"
        f"{bullets}\n"
        "---\n\n"
        "## Contexto clínico‑epidemiológico (Colombia) y subtipos CMS\n\n"
        "La revisión **García-Foncillas *et al.*, 2026** (*Rev. colomb. hematol. oncol.*, DOI [10.51643/22562915.853](https://doi.org/10.51643/22562915.853)) "
        "subraya el aumento de CCR en población joven y el papel de clasificaciones moleculares (**CMS**, MSI, MAPK/EGFR, biopsia líquida). "
        "OncoGO-CRC incorpora:\n\n"
        "- **Edad / EOCRC:** si el Excel clínico PDC trae edad, la metadata CPTAC incluye `age`, `age_bin` (`young`/`older` según `clinical_filters.age_bins`) "
        "y la comparación opcional `cptac_eocrc_vs_lorc` (se ejecuta sola si `differential.auto_run_optional_comparisons` y hay estratos con N mínimo).\n"
        "- **MSI vs MSS:** comparación opcional `cptac_msi_vs_mss`; con `auto_run_optional_comparisons` corre si `msi_status` tiene suficientes MSI y MSS.\n"
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
        "- `results/differential/*_diff.tsv` — estadísticos por proteína (todas las comparaciones corridas).\n"
        "- `results/gsea/*_gsea_go.tsv` y `results/gsea/*_gsea_go_<métrica>.tsv` — GSEA pre-ranked (p. ej. `t_stat` y `signed_neglog10_p`).\n"
        "- `results/go/*_go.tsv` — ORA por término (lista significativos + universo detectado).\n"
        "- `results/reproducibility/cptac_pride_differential_overlap*.tsv` — firma núcleo multi-cohorte.\n"
        "- `results/export/protein_context_bundle.json` — puente a ProteinsLanguage / STRING / UniProt / AlphaFold.\n"
        "- `results/ml/<cohorte>_feature_ranking.tsv` — priorización multivariante.\n"
        "- `results/ml/<cohorte>_nested_cv_*.tsv` — nested CV (AUC por pliegue externo; hiperparámetros por grid interno).\n"
        "- `results/ml/*_train_*_test_transfer.json` — evaluación **CPTAC → PRIDE** (genes alineados).\n"
        "- `results/ml/<cohorte>_feature_stability_bootstrap.tsv` — estabilidad de ranking L1 (si `ml.feature_stability_bootstrap_n>0`).\n"
        "- `results/external_evidence/evidence_matrix.tsv` — matriz de evidencia bibliográfica / módulos.\n"
        "- `results/external_evidence/tcga_rnaseq_concordance.tsv` — concordancia direccional proteómica vs RNA-seq (opcional).\n"
        "- `results/validation/input_validation_report.json` — trazabilidad de calidad.\n"
    )


def build_limitations_markdown(settings: Settings) -> str:
    """Limitaciones explícitas enlazadas a config."""
    dc = settings.config.get("differential", {})
    go = settings.config.get("go", {})
    gsea = settings.config.get("gsea", {})
    ml = settings.config.get("ml", {})
    return (
        "# Limitaciones (auto)\n\n"
        "- **Modelo estadístico:** t-test por proteína + FDR Benjamini–Hochberg; asume razonablemente "
        "muestreo independiente entre grupos; **no** modela emparejamiento tumor-normal por sujeto "
        "(mejora futura: t pareado o modelos mixtos si hay `participant_id` confiable).\n"
        "- **GO ORA:** sensible a listas enormes de significativos; usa `go.ora_min_abs_log2fc` y "
        f"`ora_foreground_max_genes`; universo = proteínas detectadas; FDR término a {go.get('fdr_threshold', 0.05)}.\n"
        "- **GSEA:** depende de la calidad del ranking y del archivo `go_gene_long.tsv`; FDR término habitualmente "
        f"más laxo ({gsea.get('fdr_threshold', 0.25)} en config); permutaciones finitas ({gsea.get('permutations', 200)}).\n"
        f"- **ML:** hold-out simple (`test_size={ml.get('test_size', 0.25)}`) **además** de `results/ml/*_nested_cv_*.tsv` "
        "si `ml.nested_cv.enabled=true` (estimación menos optimista que un único split). "
        "**Transferencia CPTAC→PRIDE** (`ml.external_transfer.enabled`) reporta AUC en cohorte externa con genes alineados; "
        "opcional **`calibration_method`** (`sigmoid`/`isotonic`) recalibra probabilidades con CV **solo en CPTAC** antes de PRIDE. "
        "AUC modesta es interpretable (cambio de dominio). Bootstrap de features: `ml.feature_stability_bootstrap_n`.\n"
        "- **Datos públicos:** versiones de matrices y metadata dependen de PDC/PRIDE en la fecha de descarga; "
        "documenta `data/versions/` y manifiestos.\n"
        f"- **α diferencial:** {dc.get('alpha', 0.05)}; efecto mínimo opcional `min_abs_log2fc` en la escala **log2** de la matriz **log2norm**.\n"
        "- **Tejido vs plasma:** el contraste tumor/normal en tejido suele separar mejor que plasma circulante; "
        "AUC más altos en CPTAC frente a capas líquidas débiles es **esperable** y defendible frente a revisiones que citan AUC modestos en cohortes sanguíneas.\n"
        "- **CMS / KEGG:** CMS suele provenir de transcriptómica; aquí solo se usa si aportas suplemento explícito. "
        "KEGG (`kegg_enrichment`, gseapy) es opcional y depende de conectividad Enrichr.\n"
        "- **Druggability / reposicionamiento:** solo enlaces in silico en el bundle; no equivalen a evidencia de respuesta clínica.\n"
        "- **UniProt REST:** el enriquecimiento en el bundle requiere red y respeta límites de frecuencia; resultados cacheados en `data/interim/`.\n"
    )


def build_results_inventory_markdown(settings: Settings) -> str:
    """Lista compacta de artefactos bajo results/."""
    res = settings.resolve(settings.config["paths"]["results_dir"])
    lines = ["# Inventario de archivos en `results/`\n\n", _workspace_doc_note(settings)]
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


def _external_evidence_status(summary: dict) -> dict[str, bool]:
    """Devuelve flags: ¿debe renderizarse tcga / cox como resultado real?"""
    def _ok(block: dict | None) -> bool:
        if not isinstance(block, dict):
            return False
        st = str(block.get("status", "")).lower()
        if st in ("skipped_no_path", "file_not_found", "missing_log2fc_column",
                  "missing_cox_hr_or_p_columns", "insufficient_overlap_or_na"):
            return False
        return bool(block.get("n_matched_genes") or block.get("n_matched"))

    return {
        "tcga_ok": _ok(summary.get("tcga_rnaseq")),
        "cox_ok": _ok(summary.get("survival_cox")),
    }


def build_external_evidence_markdown(settings: Settings) -> str:
    """Documenta matriz de evidencia, concordancia RNA-seq y Cox opcional.

    Solo incrusta extractos cuando ``external_evidence_summary.json`` reporta una unión real;
    si la etapa estuvo ``skipped_no_path`` / sin archivo / sin solapamiento, declara explícitamente
    que no hay evidencia externa cargada (evita publicar plantillas o TSVs antiguos).
    """
    res = settings.resolve(settings.config["paths"]["results_dir"])
    ex = res / "external_evidence"
    parts: list[str] = [
        "# Evidencia externa e in silico\n\n",
        "Esta etapa (`external_evidence`) **no sustituye** experimentos de laboratorio: integra el diferencial núcleo con "
        "**capas YAML** (`crc_context_layers.yaml`), tablas **opcionales** de referencia transcriptómica (TCGA u otras) "
        "y **Cox por gen** si aportas un TSV (p. ej. export GEPIA / cBioPortal).\n\n",
        "## Qué se genera\n\n",
        "1. **`evidence_matrix.tsv`**: cada fila = feature del contraste `external_evidence.primary_comparison`; columnas "
        "`hit_*` = pertenencia a listas de la revisión (tiers, CMS, módulos funcionales, eje MAPK/plasma). "
        "Si existe traslación clínica, se fusionan columnas de scoring (`translational_priority_score`, etc.).\n",
        "2. **`tcga_rnaseq_concordance.tsv`** (si config apunta a un TSV con `gene_symbol` y `log2fc_rnaseq` o alias): "
        "unión por gen con el `log2fc` proteómico; columnas de concordancia de signo.\n",
        "3. **`survival_gene_cox_merged.tsv`** (si config apunta a `gene_symbol`, `cox_hr`, `cox_p`): unión con el diferencial.\n",
        "4. **`external_evidence_summary.json`** + **`external_evidence_summary.tsv`**: métricas agregadas (Spearman signo, n matched, …).\n\n",
    ]
    summ_path = ex / "external_evidence_summary.json"
    summary: dict = {}
    if summ_path.exists():
        try:
            summary = json.loads(summ_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            summary = {}
        parts.append("## Resumen de la última corrida\n\n")
        txt = summ_path.read_text(encoding="utf-8")
        if len(txt) > 12000:
            parts.append(f"```json\n{txt[:12000]}\n…\n```\n\n")
        else:
            parts.append(f"```json\n{txt}\n```\n\n")
    else:
        parts.append(
            "_Sin salidas todavía._ Activa `external_evidence.enabled: true` en `config.yaml` y ejecuta el pipeline "
            "después de `differential` y `clinical`.\n\n"
        )

    flags = _external_evidence_status(summary)
    em = ex / "evidence_matrix.tsv"
    if em.exists():
        parts.append("## Extracto — matriz de evidencia (capas YAML internas)\n\n")
        parts.append(_df_to_md_table(read_tsv(em).head(12), max_rows=12))

    tc = ex / "tcga_rnaseq_concordance.tsv"
    if flags["tcga_ok"] and tc.exists():
        parts.append("\n## Extracto — concordancia proteómica vs RNA-seq de referencia\n\n")
        parts.append(_df_to_md_table(read_tsv(tc).head(12), max_rows=12))
    else:
        parts.append(
            "\n## Concordancia RNA-seq externa (TCGA u otro)\n\n"
            "_No cargada en esta corrida._ Para activarla, fija `external_evidence.tcga_rnaseq_log2fc_tsv` "
            "(formato: columnas `gene_symbol`, `log2fc_rnaseq` o alias). "
            "Mientras tanto, **no se publican extractos**: cualquier `tcga_rnaseq_concordance.tsv` previo "
            "se borra al inicio de la etapa para no arrastrar plantillas.\n\n"
        )

    sv = ex / "survival_gene_cox_merged.tsv"
    if flags["cox_ok"] and sv.exists():
        parts.append("\n## Extracto — Cox externo fusionado\n\n")
        parts.append(_df_to_md_table(read_tsv(sv).head(12), max_rows=12))
    else:
        parts.append(
            "\n## Supervivencia externa (Cox por gen)\n\n"
            "_No cargada en esta corrida._ Fija `external_evidence.survival_gene_cox_tsv` con columnas "
            "`gene_symbol`, `cox_hr`, `cox_p` (export GEPIA / cBioPortal). Sin TSV real, no se renderiza tabla.\n\n"
        )
    parts.append(
        "\n---\n\n**Cómo preparar las tablas de referencia:** coloca TSV bajo `metadata/tcga_auxiliary/` "
        "(o ruta absoluta) y fija `external_evidence.tcga_rnaseq_log2fc_tsv` / `survival_gene_cox_tsv`. "
        "Los archivos `*_template.tsv` en el repo son **plantillas de columnas**; sustituye valores por datos "
        "curados de tu fuente (GDC, Xena, GEPIA, etc.).\n"
    )
    return "".join(parts)


def build_web_index_markdown(settings: Settings) -> str:
    """Portada del bloque Proyecto Web con enlaces accionables."""
    return (
        "# Proyecto web — OncoGO-CRC\n\n"
        f"{_workspace_doc_note(settings)}"
        "Sitio generado automáticamente por `oncogo-crc run --stage publish_results` a partir de los "
        "TSV bajo `results/`. Las páginas se reescriben en cada corrida; los hashes y `run_id` para "
        "trazabilidad están en `results/run_manifest.json`.\n\n"
        "## Lectura recomendada (orden para tesis)\n\n"
        "1. [**Síntesis para tesis**](thesis_synthesis.md) — qué demuestra el pipeline y qué redactar aparte (narrativa).\n"
        "2. [**Interpretación biológica**](interpretacion_biologica.md) — módulos funcionales detectados "
        "por comparación (OXPHOS, mitorribosoma, proteasoma, ECM/CMS4, fase aguda…) con contexto CRC.\n"
        "3. [**Figuras**](figuras.md) — volcanos etiquetados, ORA/GSEA top, scatter CPTAC↔PRIDE, "
        "PCA, heatmaps, AUC nested CV y transferencia.\n"
        "4. [**Resultados base**](results_base.md) — tablas cuantitativas: proteínas, ORA, GSEA, "
        "firma CPTAC∩PRIDE, ML, red, traslación, export bundle.\n"
        "5. [**Glosario**](glosario.md) — definiciones operativas: padj/FDR, NES, leading edge, "
        "nested CV, transferencia, calibración, TPS, etc.\n\n"
        "## Referencia secundaria\n\n"
        "- [Datasets](datasets.md) · [Cobertura clínica](clinical_coverage.md) · "
        "[Demografía (estadio, etnia)](clinical_demographics.md) · "
        "[Evidencia externa](external_evidence.md) · [Estabilidad](stability.md) · "
        "[Inventario](results_inventory.md) · [Limitaciones](limitations.md)\n"
        "- Documento de limitaciones honestas y plan de validación: "
        "[`docs/RESULTADOS_Y_LIMITACIONES_TESIS.md`](../RESULTADOS_Y_LIMITACIONES_TESIS.md)\n"
        "- Prompt para redacción interpretativa con LLM: "
        "[`docs/PROMPT_SINTESIS_TESIS.md`](../PROMPT_SINTESIS_TESIS.md)\n"
    )
