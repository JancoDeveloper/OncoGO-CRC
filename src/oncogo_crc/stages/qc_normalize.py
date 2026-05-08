from __future__ import annotations

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv, write_tsv
from oncogo_crc.utils.qc import filter_missingness, knn_impute, log2_transform, zscore_by_protein


def run(settings: Settings) -> None:
    """QC: missingness, log2, imputación KNN, normalización z-score por proteína."""
    pp = settings.config.get("preprocessing", {})
    interim = settings.resolve(settings.config["paths"]["interim_dir"])
    proc = settings.resolve(settings.config["paths"]["processed_dir"])
    meta_cols = {"protein_id", "gene_symbol", "uniprot_accession"}

    for key in settings.datasets["datasets"]:
        path = interim / f"{key}_harmonized.tsv"
        if not path.exists():
            continue
        df = read_tsv(path)
        meta = [c for c in df.columns if c in meta_cols]
        sample_cols = [c for c in df.columns if c not in meta_cols]
        df = filter_missingness(
            df,
            meta,
            protein_missing_max=float(pp.get("missing_threshold_protein", 0.3)),
            sample_missing_max=float(pp.get("missing_threshold_sample", 0.25)),
        )
        sample_cols = [c for c in df.columns if c not in meta_cols]
        mat = df[sample_cols].astype(float)
        if pp.get("log2_transform", True):
            mat = log2_transform(mat)
        mat = knn_impute(mat, n_neighbors=int(pp.get("knn_neighbors", 5)))
        mat = zscore_by_protein(mat)
        out = pd.concat([df[meta].reset_index(drop=True), mat.reset_index(drop=True)], axis=1)
        out_path = proc / f"{key}_normalized.tsv"
        write_tsv(out, out_path)
        print(f"[qc] {key} -> {out_path.relative_to(settings.root)}")
