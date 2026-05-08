from __future__ import annotations


from oncogo_crc.settings import Settings
from oncogo_crc.utils.go_utils import load_go_long, ora_fisher, semantic_reduce_go, _build_term_genes
from oncogo_crc.utils.io import read_tsv, write_tsv


def run(settings: Settings) -> None:
    """GO como eje central: ORA offline + reducción semántica por solapamiento de genes."""
    go_cfg = settings.config.get("go", {})
    go_long = load_go_long(settings.go_gene_long_path())
    namespaces = go_cfg.get("ontologies", ["BP", "MF", "CC"])
    term2g = _build_term_genes(go_long, namespaces)

    proc = settings.resolve(settings.config["paths"]["processed_dir"])
    diff_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "differential"
    out_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "go"

    comps = settings.comparisons.get("comparisons", {})
    for comp_key, spec in comps.items():
        if spec.get("optional"):
            continue
        diff_path = diff_dir / f"{comp_key}_diff.tsv"
        if not diff_path.exists():
            continue
        diff = read_tsv(diff_path)
        sig_genes = set(diff.loc[diff["padj"] < float(go_cfg.get("fdr_threshold", 0.05)), "protein_id"].astype(str))

        ds_key = spec["dataset"]
        mat_path = proc / f"{ds_key}_normalized.tsv"
        if not mat_path.exists():
            continue
        expr = read_tsv(mat_path)
        universe = set(expr["protein_id"].astype(str))

        res = ora_fisher(sig_genes, universe, term2g, fdr_threshold=float(go_cfg.get("fdr_threshold", 0.05)))
        write_tsv(res, out_dir / f"{comp_key}_go.tsv")
        if go_cfg.get("semantic_reduction", True) and not res.empty:
            red = semantic_reduce_go(res, term2g)
            write_tsv(red, out_dir / f"{comp_key}_go_reduced.tsv")
        print(f"[go] {comp_key} -> results/go/{comp_key}_go.tsv")
