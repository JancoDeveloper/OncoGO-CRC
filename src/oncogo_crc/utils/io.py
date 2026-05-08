from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def write_tsv(df: pd.DataFrame, path: Path, *, mkdir: bool = True) -> None:
    if mkdir:
        path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def read_expression_matrix(path: Path, protein_column: str = "protein_id") -> pd.DataFrame:
    df = read_tsv(path)
    if protein_column not in df.columns:
        raise ValueError(f"Falta columna '{protein_column}' en {path}")
    return df
