from __future__ import annotations

import json
from hashlib import sha256

from oncogo_crc.settings import Settings
from oncogo_crc.utils.go_utils import (
    _build_term_genes,
    assert_go_annotation_not_synthetic,
    demo_like_go_fraction,
    expression_protein_to_gene_map,
    gene_symbol_set_from_strings,
    load_go_long,
    ora_fisher,
    semantic_reduce_go,
)
from oncogo_crc.utils.io import read_tsv, write_tsv
from oncogo_crc.utils import run_manifest


def run(settings: Settings) -> None:
    """GO como eje central: ORA offline + reducción semántica por solapamiento de genes."""
    go_cfg = settings.config.get("go", {})
    go_path = settings.go_gene_long_path()
    go_long = load_go_long(go_path)
    assert_go_annotation_not_synthetic(go_long, demo=settings.demo)
    demo_frac = demo_like_go_fraction(go_long)
    if settings.demo and demo_frac > 0.08:
        print(
            f"[go] aviso: ~{demo_frac:.0%} de anotaciones parecen sintéticas (GENE###). "
            f"Sustituye {settings.relative_to_project(go_path)} por GOA/UniProt+ real para resultados publicables."
        )
    namespaces = go_cfg.get("ontologies", ["BP", "MF", "CC"])
    term2g = _build_term_genes(go_long, namespaces)

    proc = settings.resolve(settings.config["paths"]["processed_dir"])
    diff_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "differential"
    out_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "go"

    comps = settings.comparisons.get("comparisons", {})
    n_written = 0
    for comp_key, spec in comps.items():
        diff_path = diff_dir / f"{comp_key}_diff.tsv"
        if not diff_path.exists():
            continue
        diff = read_tsv(diff_path)
        fdr_go = float(go_cfg.get("fdr_threshold", 0.05))
        padj = diff["padj"].astype(float)
        sig_mask = padj.notna() & (padj < fdr_go)
        if "significant" in diff.columns:
            sig_mask = sig_mask | diff["significant"].fillna(False)
        ora_fc_min = float(go_cfg.get("ora_min_abs_log2fc", 0.0) or 0.0)
        if ora_fc_min > 0 and "log2fc" in diff.columns:
            sig_mask = sig_mask & diff["log2fc"].notna() & (diff["log2fc"].abs() >= ora_fc_min)

        ds_key = spec["dataset"]
        mat_path = settings.processed_matrix_path(ds_key, "log2norm")
        if not mat_path.exists():
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

        sub = diff.loc[sig_mask].copy()
        sub = sub.sort_values("padj", na_position="last")
        max_fg = int(go_cfg.get("ora_foreground_max_genes", 800) or 0)
        n_sig_raw = len(sub)
        if max_fg > 0 and len(sub) > max_fg:
            sub = sub.head(max_fg)
            print(
                f"[go] {comp_key}: ORA usa top {max_fg} genes por padj "
                f"(había {n_sig_raw} filas significativas; evita dilución ORA)."
            )
        sig_genes = gene_symbol_set_from_strings(
            {pid_to_gene.get(str(r["protein_id"]), str(r["protein_id"])) for _, r in sub.iterrows()}
        )
        sig_genes = sig_genes & universe

        res = ora_fisher(
            sig_genes,
            universe,
            term2g,
            fdr_threshold=float(go_cfg.get("fdr_threshold", 0.05)),
            min_genes_in_term=int(go_cfg.get("min_genes_per_term", 5)),
        )
        write_tsv(res, out_dir / f"{comp_key}_go.tsv")
        if go_cfg.get("semantic_reduction", True) and not res.empty:
            red = semantic_reduce_go(res, term2g)
            write_tsv(red, out_dir / f"{comp_key}_go_reduced.tsv")

        min_t = max(3, int(go_cfg.get("min_genes_per_term", 5)))
        n_terms_ok = sum(1 for _t, gs in term2g.items() if len(gs & universe) >= min_t)
        audit = {
            "comparison_id": comp_key,
            "diff_path": str(settings.relative_to_project(diff_path)),
            "sha256_diff_input": sha256(diff_path.read_bytes()).hexdigest(),
            "matrix_log2norm": str(settings.relative_to_project(mat_path)),
            "go_annotation_path": str(settings.relative_to_project(go_path)),
            "n_genes_universe_from_matrix": len(universe),
            "n_genes_universe_intersect_go_annotation": len(genes_annotated_in_universe),
            "n_genes_foreground_ora": len(sig_genes),
            "n_go_terms_passing_size_filter": int(n_terms_ok),
            "demo_like_go_row_fraction": round(demo_frac, 4),
        }
        (out_dir / f"{comp_key}_go_input_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
        print(f"[go] {comp_key} -> results/go/{comp_key}_go.tsv")
        n_written += 1

    run_manifest.record_stage_meta(
        settings,
        "go",
        executed=n_written > 0,
        comparisons_written=int(n_written),
    )
