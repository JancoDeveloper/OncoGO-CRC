from __future__ import annotations


import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.go_utils import load_go_long, ora_fisher, semantic_reduce_go, _build_term_genes
from oncogo_crc.utils.io import read_tsv, write_tsv


def _pid_to_gene_map(expr: pd.DataFrame) -> dict[str, str]:
    """Mapea protein_id de la matriz normalizada a gene_symbol para coincidir con go_gene_long."""
    m: dict[str, str] = {}
    for _, row in expr.iterrows():
        pid = str(row["protein_id"]).strip()
        gs = row["gene_symbol"] if "gene_symbol" in row.index else None
        if gs is None or (isinstance(gs, float) and pd.isna(gs)) or str(gs).strip() == "":
            m[pid] = pid
        else:
            m[pid] = str(gs).strip()
    return m


def _gene_set(raw: set[str]) -> set[str]:
    out: set[str] = set()
    for x in raw:
        s = str(x).strip()
        if s and s.lower() != "nan":
            out.add(s)
    return out


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
        fdr_go = float(go_cfg.get("fdr_threshold", 0.05))
        padj = diff["padj"].astype(float)
        sig_mask = padj.notna() & (padj < fdr_go)
        if "significant" in diff.columns:
            sig_mask = sig_mask | diff["significant"].fillna(False)

        ds_key = spec["dataset"]
        mat_path = proc / f"{ds_key}_normalized.tsv"
        if not mat_path.exists():
            continue
        expr = read_tsv(mat_path)
        pid_to_gene = _pid_to_gene_map(expr)
        universe = _gene_set(
            {pid_to_gene.get(str(r["protein_id"]), str(r["protein_id"])) for _, r in expr.iterrows()}
        )

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
        sig_genes = _gene_set(
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
        print(f"[go] {comp_key} -> results/go/{comp_key}_go.tsv")
