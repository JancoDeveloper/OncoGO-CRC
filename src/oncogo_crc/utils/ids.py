from __future__ import annotations

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv


def effective_protein_mapping(settings: Settings) -> pd.DataFrame:
    """Mapa maestro + TSV suplementarios en ``mapping.supplementary_map_paths`` (prioridad al maestro)."""
    mapping = read_tsv(settings.protein_id_map_path())
    map_cfg = settings.config.get("mapping", {}) or {}
    extras: list[pd.DataFrame] = []
    for rel in map_cfg.get("supplementary_map_paths", []) or []:
        p = settings.resolve(str(rel).strip())
        if p.exists():
            extras.append(read_tsv(p))
    if extras:
        return merge_protein_mapping_tables(mapping, *extras)
    return mapping


def _dedupe_mapping_rows(mapping: pd.DataFrame, *, protein_col: str) -> pd.DataFrame:
    """Una fila por ``protein_col`` (prioriza filas con UniProt) para evitar multiplicación en el merge."""
    if mapping is None or mapping.empty or protein_col not in mapping.columns:
        return mapping
    m = mapping.copy()
    if "uniprot_accession" in m.columns:
        m["_u_pri"] = m["uniprot_accession"].notna() & (m["uniprot_accession"].astype(str).str.strip() != "")
    else:
        m["_u_pri"] = False
    m = m.sort_values("_u_pri", ascending=False).drop_duplicates(subset=[protein_col], keep="first")
    return m.drop(columns=["_u_pri"], errors="ignore")


def merge_protein_mapping_tables(
    primary: pd.DataFrame,
    *extras: pd.DataFrame,
    protein_col: str = "protein_id",
) -> pd.DataFrame:
    """Concatena mapas; el ``primary`` gana en conflicto de ``protein_id`` (primer bloque conservado)."""
    frames: list[pd.DataFrame] = []
    if primary is not None and not primary.empty:
        frames.append(primary.copy())
    for ex in extras:
        if ex is not None and not ex.empty:
            frames.append(ex.copy())
    if not frames:
        return pd.DataFrame()
    all_cols: set[str] = set()
    for f in frames:
        all_cols.update(f.columns)
    aligned: list[pd.DataFrame] = []
    for i, f in enumerate(frames):
        fc = f.copy()
        fc["_map_rank"] = i
        if "uniprot_accession" in fc.columns:
            fc["_u_pri"] = fc["uniprot_accession"].notna() & (fc["uniprot_accession"].astype(str).str.strip() != "")
        else:
            fc["_u_pri"] = False
        for c in sorted(all_cols):
            if c not in fc.columns:
                fc[c] = pd.NA
        aligned.append(fc[list(sorted(all_cols)) + ["_map_rank", "_u_pri"]])
    out = pd.concat(aligned, ignore_index=True)
    out = out.sort_values([protein_col, "_map_rank", "_u_pri"], ascending=[True, True, False])
    out = out.drop_duplicates(subset=[protein_col], keep="first")
    return out.drop(columns=["_map_rank", "_u_pri"], errors="ignore")


def harmonize_protein_table(
    expr: pd.DataFrame,
    mapping: pd.DataFrame,
    *,
    protein_col: str = "protein_id",
) -> pd.DataFrame:
    """Une la matriz con el mapa maestro y colapsa filas duplicadas por símbolo de gen (mediana por muestra)."""
    mapping = _dedupe_mapping_rows(mapping, protein_col=protein_col)
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
