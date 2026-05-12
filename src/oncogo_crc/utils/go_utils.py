from __future__ import annotations

from collections import defaultdict
from typing import Iterable

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def load_go_long(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _build_term_genes(go_long: pd.DataFrame, namespaces: Iterable[str] | None) -> dict[str, set[str]]:
    dfg = go_long.copy()
    if namespaces is not None:
        dfg = dfg[dfg["namespace"].isin(list(namespaces))]
    term_to_genes: dict[str, set[str]] = defaultdict(set)
    for _, r in dfg.iterrows():
        term_to_genes[str(r["go_id"])].add(str(r["gene_symbol"]))
    return dict(term_to_genes)


def ora_fisher(
    significant_genes: set[str],
    universe_genes: set[str],
    term_to_genes: dict[str, set[str]],
    *,
    fdr_threshold: float = 0.05,
    min_genes_in_term: int = 5,
) -> pd.DataFrame:
    rows = []
    sig = significant_genes & universe_genes
    unsig = universe_genes - sig
    min_t = max(3, int(min_genes_in_term))
    for term, genes in term_to_genes.items():
        in_term = genes & universe_genes
        if len(in_term) < min_t:
            continue
        a = len(sig & genes)
        b = len(sig - genes)
        c = len(unsig & genes)
        d = len(unsig - genes)
        if a + c == 0 or b + d == 0:
            continue
        oddsratio, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append((term, len(in_term), a, oddsratio, p))
    if not rows:
        return pd.DataFrame(columns=["go_id", "term_size", "overlap", "odds_ratio", "p_value", "padj"])
    res = pd.DataFrame(rows, columns=["go_id", "term_size", "overlap", "odds_ratio", "p_value"])
    _, q, _, _ = multipletests(res["p_value"], method="fdr_bh")
    res["padj"] = q
    res["significant"] = res["padj"] < fdr_threshold
    return res.sort_values("p_value")


def semantic_reduce_go(
    res: pd.DataFrame,
    term_to_genes: dict[str, set[str]],
    *,
    jaccard_min: float = 0.45,
) -> pd.DataFrame:
    terms = res["go_id"].astype(str).tolist()

    def jaccard(a: set[str], b: set[str]) -> float:
        u = len(a | b)
        if u == 0:
            return 0.0
        return len(a & b) / u

    assigned: set[str] = set()
    clusters: list[list[str]] = []
    porder = res.set_index("go_id")["p_value"].to_dict()
    for t in sorted(terms, key=lambda x: float(porder.get(x, 1.0))):
        if t in assigned:
            continue
        cl = [t]
        genes_t = term_to_genes.get(t, set())
        for u in terms:
            if u == t or u in assigned:
                continue
            if jaccard(genes_t, term_to_genes.get(u, set())) >= jaccard_min:
                cl.append(u)
        for x in cl:
            assigned.add(x)
        clusters.append(cl)
    reduced_rows = []
    for cl in clusters:
        sub = res[res["go_id"].isin(cl)].sort_values("p_value")
        rep = sub.iloc[0]
        reduced_rows.append(
            {
                "representative_go_id": rep["go_id"],
                "cluster_size": len(cl),
                "member_go_ids": ";".join(cl),
                "p_value": rep["p_value"],
                "padj": rep["padj"],
                "overlap": rep["overlap"],
                "term_size": rep["term_size"],
            }
        )
    return pd.DataFrame(reduced_rows)
