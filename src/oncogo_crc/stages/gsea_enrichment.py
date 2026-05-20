"""GSEA pre-ranked (GO offline) con gseapy — ranking continuo sin umbral duro de significación."""

from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils import run_manifest
from oncogo_crc.utils.go_utils import (
    _build_term_genes,
    assert_go_annotation_not_synthetic,
    demo_like_go_fraction,
    expression_protein_to_gene_map,
    gene_symbol_set_from_strings,
    load_go_long,
)
from oncogo_crc.utils.io import read_tsv, write_tsv


def _gene_sets_dict_for_prerank(
    term2g: dict[str, set[str]],
    universe: set[str],
    *,
    min_size: int,
    max_size: int,
) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for term, genes in term2g.items():
        g2 = genes & universe
        n = len(g2)
        if n < min_size or n > max_size:
            continue
        out[str(term)] = sorted(g2)
    return out


def _n_terms_with_universe_overlap(term2g: dict[str, set[str]], universe: set[str], *, min_genes: int) -> int:
    return sum(1 for _t, gs in term2g.items() if len(gs & universe) >= min_genes)


def _rank_dataframe(
    diff: pd.DataFrame,
    pid_to_gene: dict[str, str],
    *,
    ranking_metric: str,
) -> pd.DataFrame:
    rows: list[tuple[str, float]] = []
    for _, r in diff.iterrows():
        gene = pid_to_gene.get(str(r["protein_id"]), str(r["protein_id"]))
        gene = str(gene).strip()
        if not gene or gene.lower() == "nan":
            continue
        if ranking_metric == "signed_neglog10_p":
            p = r.get("p_value")
            fc = r.get("log2fc")
            if p is None or fc is None or pd.isna(p) or pd.isna(fc):
                continue
            pv = float(max(float(p), 1e-300))
            score = float(np.sign(float(fc)) * (-np.log10(pv)))
        else:
            v = r.get("t_stat")
            if v is None or pd.isna(v):
                continue
            score = float(v)
        rows.append((gene, score))
    if not rows:
        return pd.DataFrame(columns=["gene", "score"])
    df = pd.DataFrame(rows, columns=["gene", "score"])
    df["_ab"] = df["score"].abs()
    df = df.sort_values("_ab", ascending=False).drop_duplicates("gene", keep="first").drop(
        columns="_ab"
    )
    return df.sort_values("score", ascending=False).reset_index(drop=True)


def _purge_gsea_output_dir(out_dir: Path) -> None:
    if not out_dir.exists():
        return
    for p in out_dir.iterdir():
        if p.is_file():
            p.unlink()


def run(settings: Settings) -> None:
    cfg = settings.config.get("gsea", {}) or {}
    res_root = settings.resolve(settings.config["paths"]["results_dir"])
    out_dir = res_root / "gsea"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not bool(cfg.get("enabled", False)):
        print("[gsea] deshabilitado (gsea.enabled=false o ausente en config.yaml).")
        run_manifest.record_stage_meta(
            settings,
            "gsea",
            executed=False,
            reason="disabled_in_config",
            comparisons={},
        )
        return
    try:
        import gseapy as gp  # type: ignore[import-untyped]
    except ImportError as e:
        print(
            "[gsea] omitido: falta gseapy. Instala con: pip install -e \".[go-remote]\" "
            f"({e!r})."
        )
        _purge_gsea_output_dir(out_dir)
        run_manifest.record_stage_meta(
            settings,
            "gsea",
            executed=False,
            reason="missing_gseapy",
            comparisons={},
        )
        return

    go_cfg = settings.config.get("go", {})
    go_path = settings.go_gene_long_path()
    go_long = load_go_long(go_path)
    assert_go_annotation_not_synthetic(go_long, demo=settings.demo)
    demo_frac = demo_like_go_fraction(go_long)
    namespaces = list(cfg.get("ontologies", go_cfg.get("ontologies", ["BP", "MF", "CC"])))
    term2g = _build_term_genes(go_long, namespaces)

    if bool(cfg.get("clear_stale_outputs_before_run", True)):
        _purge_gsea_output_dir(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

    diff_dir = res_root / "differential"

    min_size = int(cfg.get("min_size", 15))
    max_size = int(cfg.get("max_size", 500))
    permutations = int(cfg.get("permutations", 200))
    fdr_thr = float(cfg.get("fdr_threshold", 0.25))
    metrics = cfg.get("ranking_metrics")
    if isinstance(metrics, list) and metrics:
        ranking_metrics = [str(m).strip() for m in metrics if str(m).strip()]
    else:
        ranking_metrics = [str(cfg.get("ranking_metric", "t_stat")).strip()]

    seed = int(cfg.get("seed", 42))

    comps = settings.comparisons.get("comparisons", {})
    comp_state: dict[str, dict[str, object]] = {}
    any_prerank_ok = False

    for comp_key, spec in comps.items():
        diff_path = diff_dir / f"{comp_key}_diff.tsv"
        if not diff_path.exists():
            comp_state[comp_key] = {"status": "skipped", "reason": "missing_differential"}
            continue
        diff = read_tsv(diff_path)
        if "protein_id" not in diff.columns:
            comp_state[comp_key] = {"status": "skipped", "reason": "missing_protein_id"}
            continue
        ds_key = spec["dataset"]
        mat_path = settings.processed_matrix_path(ds_key, "log2norm")
        if not mat_path.exists():
            comp_state[comp_key] = {"status": "skipped", "reason": "missing_log2norm_matrix"}
            continue
        expr = read_tsv(mat_path)
        pid_to_gene = expression_protein_to_gene_map(expr)
        universe = gene_symbol_set_from_strings(
            {pid_to_gene.get(str(r["protein_id"]), str(r["protein_id"])) for _, r in expr.iterrows()}
        )
        all_term_genes: set[str] = set()
        for genes in term2g.values():
            all_term_genes |= genes
        genes_annotated_in_universe = universe & all_term_genes

        n_terms_before = _n_terms_with_universe_overlap(term2g, universe, min_genes=3)
        gsets = _gene_sets_dict_for_prerank(term2g, universe, min_size=min_size, max_size=max_size)
        n_terms_after = len(gsets)
        if len(gsets) < 3:
            print(f"[gsea] {comp_key}: pocos conjuntos GO tras filtrar tamaño ({len(gsets)}); omitiendo.")
            comp_state[comp_key] = {
                "status": "skipped_no_sets",
                "n_terms_before_filter": n_terms_before,
                "n_terms_after_filter": n_terms_after,
                "n_genes_universe": len(universe),
                "n_genes_in_go_universe": len(genes_annotated_in_universe),
            }
            continue

        comp_ok_for_comp = False
        for idx, ranking_metric in enumerate(ranking_metrics):
            rnk = _rank_dataframe(diff, pid_to_gene, ranking_metric=ranking_metric)
            if len(rnk) < min_size:
                print(
                    f"[gsea] {comp_key} metric={ranking_metric}: ranking insuficiente ({len(rnk)} genes); omitiendo."
                )
                continue
            safe_tag = ranking_metric.replace(" ", "_")
            rank_path = out_dir / f"{comp_key}_gsea_ranking_{safe_tag}.tsv"
            write_tsv(rnk.rename(columns={"gene": "gene_symbol", "score": "rank_metric"}), rank_path)
            demo_ids = bool(rnk["gene"].astype(str).str.match(r"^GENE\d{3,}$", na=False).any())
            sample_genes = rnk["gene"].astype(str).head(15).tolist()
            ranked_in_go = sorted(set(rnk["gene"].astype(str).str.upper()) & all_term_genes)
            ranked_missing_go = sorted(set(rnk["gene"].astype(str).str.upper()) - all_term_genes)
            audit = {
                "comparison_id": comp_key,
                "ranking_metric": ranking_metric,
                "diff_path": str(settings.relative_to_project(diff_path)),
                "sha256_diff_input": sha256(diff_path.read_bytes()).hexdigest(),
                "sha256_ranking_file": sha256(rank_path.read_bytes()).hexdigest(),
                "matrix_log2norm": str(settings.relative_to_project(mat_path)),
                "go_annotation_path": str(settings.relative_to_project(go_path)),
                "n_ranked_genes": int(len(rnk)),
                "n_ranked_genes_total": int(len(rnk)),
                "n_ranked_genes_with_go": int(len(ranked_in_go)),
                "n_genes_universe_from_matrix": len(universe),
                "n_genes_universe_intersect_go_annotation": len(genes_annotated_in_universe),
                "n_genes_in_go": len(genes_annotated_in_universe),
                "n_terms_pre_filter": n_terms_before,
                "n_terms_before_filter": n_terms_before,
                "n_terms_post_filter": n_terms_after,
                "n_terms_after_filter": n_terms_after,
                "n_gene_sets_passing_size_filter": len(gsets),
                "min_size": min_size,
                "max_size": max_size,
                "demo_like_go_row_fraction": round(demo_frac, 4),
                "contains_demo_gene_ids": demo_ids,
                "sample_ranked_genes": sample_genes,
                "sample_ranked_genes_with_go": ranked_in_go[:20],
                "top_missing_genes_from_go": ranked_missing_go[:20],
                "n_ranked_genes_missing_from_go": int(len(ranked_missing_go)),
            }
            (out_dir / f"{comp_key}_gsea_input_{safe_tag}.json").write_text(
                json.dumps(audit, indent=2), encoding="utf-8"
            )
            pre = gp.prerank(
                rnk=rnk,
                gene_sets=gsets,
                permutation_num=permutations,
                min_size=min_size,
                max_size=max_size,
                seed=seed,
                verbose=False,
                outdir=None,
                no_plot=True,
            )
            res2 = pre.res2d.copy()
            fdr_col = "FDR q-val"
            if fdr_col not in res2.columns:
                for cand in ("FDR q-val", "FDR", "fdr", "padj"):
                    if cand in res2.columns:
                        fdr_col = cand
                        break
                else:
                    fdr_col = res2.columns[-1]
            res2["significant"] = res2[fdr_col].astype(float) < fdr_thr
            out_named = out_dir / f"{comp_key}_gsea_go_{safe_tag}.tsv"
            write_tsv(res2, out_named)
            if idx == 0:
                write_tsv(res2, out_dir / f"{comp_key}_gsea_go.tsv")
            n_sig = int(res2["significant"].sum()) if "significant" in res2.columns else 0
            print(
                f"[gsea] {comp_key} metric={ranking_metric} sets={len(gsets)} "
                f"NES/FDR<{fdr_thr}: {n_sig} términos -> {out_named.name}"
            )
            comp_ok_for_comp = True
            any_prerank_ok = True

        comp_state[comp_key] = {
            "status": "ok" if comp_ok_for_comp else "skipped_no_ranking",
            "n_terms_before_filter": n_terms_before,
            "n_terms_after_filter": n_terms_after,
            "n_genes_universe": len(universe),
            "n_genes_in_go_universe": len(genes_annotated_in_universe),
        }
        if comp_ok_for_comp:
            main = out_dir / f"{comp_key}_gsea_go.tsv"
            comp_state[comp_key]["outputs"] = [str(settings.relative_to_project(main))]
            if main.exists():
                comp_state[comp_key]["sha256_main_output"] = sha256(main.read_bytes()).hexdigest()

    run_manifest.record_stage_meta(
        settings,
        "gsea",
        executed=bool(any_prerank_ok),
        comparisons=comp_state,
    )
