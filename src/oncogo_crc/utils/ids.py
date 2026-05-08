from __future__ import annotations

import pandas as pd


def harmonize_protein_table(
    expr: pd.DataFrame,
    mapping: pd.DataFrame,
    *,
    protein_col: str = "protein_id",
) -> pd.DataFrame:
    """Une la matriz con el mapa maestro y colapsa filas duplicadas por símbolo de gen (mediana por muestra)."""
    merge_cols = [protein_col] + [c for c in ("gene_symbol", "uniprot_accession") if c in mapping.columns]
    merged = expr.merge(mapping[[c for c in merge_cols if c in mapping.columns]], on=protein_col, how="left")
    if "gene_symbol" not in merged.columns:
        merged["gene_symbol"] = merged[protein_col].astype(str)
    merged["gene_symbol"] = merged["gene_symbol"].fillna(merged[protein_col].astype(str))
    meta = {protein_col, "gene_symbol", "uniprot_accession"}
    sample_cols = [c for c in merged.columns if c not in meta]

    def _collapse(g: pd.DataFrame) -> pd.Series:
        sym = str(g.name)
        med = g[sample_cols].median(axis=0)
        row = pd.Series({protein_col: sym, "gene_symbol": sym})
        if "uniprot_accession" in g.columns:
            row["uniprot_accession"] = g["uniprot_accession"].iloc[0]
        for c in sample_cols:
            row[c] = med[c]
        return row

    gb = merged.groupby("gene_symbol", sort=False)
    try:
        out = gb.apply(_collapse, include_groups=False)
    except TypeError:
        out = gb.apply(_collapse)
    return out.reset_index(drop=True)
