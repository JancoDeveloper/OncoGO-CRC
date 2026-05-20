"""Diagnósticos de cobertura del mapa maestro de IDs (por dataset, post-harmonize)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from oncogo_crc.utils.io import read_tsv
from oncogo_crc.utils.reproducibility_overlap import (
    ambiguous_gene_symbols,
    load_gene_uniprot_from_mapping,
    uniprot_base,
)


def mapping_file_stats(mapping: pd.DataFrame) -> dict[str, int | float]:
    """Resumen del TSV de mapa (no depende de una matriz concreta)."""
    out: dict[str, int | float] = {
        "n_map_rows": int(len(mapping)),
        "n_unique_protein_id": 0,
        "n_rows_with_gene_symbol": 0,
        "n_rows_with_uniprot": 0,
        "n_ambiguous_genes": 0,
    }
    if mapping is None or mapping.empty:
        return out
    if "protein_id" in mapping.columns:
        out["n_unique_protein_id"] = int(mapping["protein_id"].nunique())
    if "gene_symbol" in mapping.columns:
        gs = mapping["gene_symbol"].astype(str).str.strip()
        out["n_rows_with_gene_symbol"] = int(((gs != "") & (gs.str.upper() != "NAN")).sum())
    if "uniprot_accession" in mapping.columns:
        ub = mapping["uniprot_accession"].map(lambda x: uniprot_base(str(x)))
        out["n_rows_with_uniprot"] = int(ub.notna().sum())
    out["n_ambiguous_genes"] = len(ambiguous_gene_symbols(mapping))
    return out


def harmonized_matrix_stats(harmonized: pd.DataFrame) -> dict[str, int | float]:
    """Cobertura de columnas de ID en una matriz ya harmonizada (interim).

    Métricas:
    - ``n_proteins_harmonized``: filas totales tras harmonize.
    - ``n_nonempty_gene_symbol_harmonized``: filas con ``gene_symbol`` no vacío/NaN.
    - ``n_resolved_gene_symbol_harmonized``: filas con ``gene_symbol`` **distinto** del ``protein_id``
      (es decir, donde el merge con el mapa aportó información nueva). Mide cuán útil fue el mapa.
    - ``n_uniprot_like_gene_symbol_harmonized``: filas donde el ``gene_symbol`` parece un accession UniProt
      (heurística de longitud/forma) → indica filas que no se resolvieron a HGNC (fallback a UniProt).
    - ``n_nonempty_uniprot_harmonized``: filas con columna ``uniprot_accession`` con base resolvable.
    """
    from oncogo_crc.utils.ids import _is_uniprot_accession

    out: dict[str, int | float] = {"n_proteins_harmonized": int(len(harmonized))}
    if "gene_symbol" in harmonized.columns:
        gs = harmonized["gene_symbol"].astype(str).str.strip()
        nonempty = (gs != "") & (gs.str.upper() != "NAN")
        out["n_nonempty_gene_symbol_harmonized"] = int(nonempty.sum())
        if "protein_id" in harmonized.columns:
            pid = harmonized["protein_id"].astype(str).str.strip()
            out["n_resolved_gene_symbol_harmonized"] = int((nonempty & (gs != pid)).sum())
        else:
            out["n_resolved_gene_symbol_harmonized"] = int(nonempty.sum())
        out["n_uniprot_like_gene_symbol_harmonized"] = int(
            (nonempty & gs.map(_is_uniprot_accession)).sum()
        )
    else:
        out["n_nonempty_gene_symbol_harmonized"] = 0
        out["n_resolved_gene_symbol_harmonized"] = 0
        out["n_uniprot_like_gene_symbol_harmonized"] = 0
    if "uniprot_accession" in harmonized.columns:
        ub = harmonized["uniprot_accession"].map(lambda x: uniprot_base(str(x)))
        out["n_nonempty_uniprot_harmonized"] = int(ub.notna().sum())
    else:
        out["n_nonempty_uniprot_harmonized"] = 0
    return out


def build_mapping_diagnostics_row(
    *,
    dataset_key: str,
    harmonized_path: Path,
    mapping: pd.DataFrame,
) -> dict[str, int | float | str]:
    row: dict[str, int | float | str] = {"dataset": dataset_key}
    row.update(mapping_file_stats(mapping))
    row["n_gene_uniprot_map_strict_1to1"] = len(load_gene_uniprot_from_mapping(mapping, exclude_ambiguous=True))
    row["n_gene_uniprot_map_lenient"] = len(load_gene_uniprot_from_mapping(mapping, exclude_ambiguous=False))
    if harmonized_path.exists():
        row.update(harmonized_matrix_stats(read_tsv(harmonized_path)))
    else:
        row["n_proteins_harmonized"] = 0
        row["n_nonempty_gene_symbol_harmonized"] = 0
        row["n_nonempty_uniprot_harmonized"] = 0
    return row


def write_mapping_diagnostics_rows(rows: list[dict[str, int | float | str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
