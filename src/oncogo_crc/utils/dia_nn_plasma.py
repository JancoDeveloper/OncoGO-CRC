"""Post-procesado de informes DIA-NN (formato largo) a matriz ancho proteína × muestra."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def pivot_dia_nn_report_to_wide(matrix_path: Path) -> tuple[pd.DataFrame, list[str]]:
    """
    Convierte DIA-NN_report.tsv (precursor/protein long) en tabla ancha:
    protein_id + una columna por Run (mediana de Genes.MaxLFQ o PG.MaxLFQ).
    """
    peek = pd.read_csv(matrix_path, sep="\t", nrows=1)
    if "Run" not in peek.columns:
        raise ValueError("DIA-NN: falta columna Run")
    value_col = "Genes.MaxLFQ" if "Genes.MaxLFQ" in peek.columns else "PG.MaxLFQ"
    idx_col = "Protein.Group" if "Protein.Group" in peek.columns else "Protein.Ids"
    if idx_col not in peek.columns or value_col not in peek.columns:
        raise ValueError(f"DIA-NN: faltan columnas esperadas ({idx_col}, {value_col})")
    usecols = [idx_col, "Run", value_col]
    df = pd.read_csv(matrix_path, sep="\t", usecols=usecols, low_memory=False)
    sub = df.dropna(subset=[idx_col, "Run"])
    sub[value_col] = pd.to_numeric(sub[value_col], errors="coerce")
    sub = sub.dropna(subset=[value_col])
    wide = sub.pivot_table(index=idx_col, columns="Run", values=value_col, aggfunc="median")
    wide = wide.reset_index().rename(columns={idx_col: "protein_id"})
    sample_cols = [c for c in wide.columns if c != "protein_id"]
    return wide, sample_cols


def build_plasma_metadata_for_pxd055218_runs(run_ids: list[str]) -> pd.DataFrame:
    """PXD055218 trae tres grupos en el nombre de la corrida.

    CP_* es cáncer colorrectal, CN_* es colonoscopia negativa (otras enfermedades
    intestinales) y Normal_* es control sano. Los tres se conservan. El contraste
    binario crc frente a sano sigue filtrando por esas dos etiquetas; tirar CN_*
    aquí borraba el único grupo que permite preguntar especificidad.
    """
    rows: list[dict[str, str]] = []
    for r in sorted(set(run_ids)):
        if not r or r == "Run":
            continue
        stem = str(r).replace("\\", "/").split("/")[-1]
        if stem.startswith("CP_"):
            cond = "crc"
        elif stem.startswith("CN_"):
            cond = "colonoscopy_negative"
        elif stem.startswith("Normal_"):
            cond = "healthy"
        else:
            continue
        rows.append({"sample_id": r, "condition": cond, "batch": "PXD055218"})
    return pd.DataFrame(rows)


def postprocess_dia_nn_plasma_download(matrix_path: Path, metadata_path: Path) -> None:
    wide, runs = pivot_dia_nn_report_to_wide(matrix_path)
    wide.to_csv(matrix_path, sep="\t", index=False)
    md = build_plasma_metadata_for_pxd055218_runs(runs)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    md.to_csv(metadata_path, sep="\t", index=False)
