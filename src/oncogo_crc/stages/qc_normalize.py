from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv, write_tsv
from oncogo_crc.utils.qc import (
    filter_missingness,
    knn_impute,
    log2_transform_masked,
    zscore_by_protein,
)
from oncogo_crc.utils.quant_columns import columns_for_family, is_summary_protein


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _pca_sample_coords(mat: pd.DataFrame, *, n_comp: int = 2) -> pd.DataFrame | None:
    """PCA muestras × componentes (solo diagnóstico; NaN→0 solo para este cálculo)."""
    try:
        from sklearn.decomposition import PCA
    except ImportError:
        return None
    X = mat.astype(float).T.fillna(0.0).values
    if X.shape[0] < 3 or X.shape[1] < 2:
        return None
    k = min(int(n_comp), X.shape[0], X.shape[1])
    if k < 1:
        return None
    pca = PCA(n_components=k)
    Z = pca.fit_transform(X)
    cols = [f"PC{i + 1}" for i in range(Z.shape[1])]
    return pd.DataFrame(Z, index=mat.columns, columns=cols)


def run(settings: Settings) -> None:
    """QC: missingness → log2 (solo x>0) + KNN → producto log2norm; z-score por proteína → zscore (ML/red)."""
    pp = settings.config.get("preprocessing", {})
    interim = settings.resolve(settings.config["paths"]["interim_dir"])
    proc = settings.resolve(settings.config["paths"]["processed_dir"])
    proc.mkdir(parents=True, exist_ok=True)
    qc_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    meta_cols = {"protein_id", "gene_symbol", "uniprot_accession"}
    legacy_alias = bool(pp.get("write_legacy_normalized_alias", True))
    write_pca = bool(pp.get("write_qc_pca", True))

    for key in settings.datasets["datasets"]:
        path = interim / f"{key}_harmonized.tsv"
        if not path.exists():
            continue
        df = read_tsv(path)
        if "protein_id" in df.columns:
            df = df.loc[~df["protein_id"].map(is_summary_protein)].copy()
        meta = [c for c in df.columns if c in meta_cols]
        sc0 = [c for c in df.columns if c not in meta_cols]
        ds_cfg = settings.datasets["datasets"].get(key, {})
        family = ds_cfg.get("quant_column_family")
        picked = columns_for_family(sc0, str(family) if family else None)
        if family and picked and len(picked) < len(sc0):
            df = df[meta + picked].copy()
            sc0 = picked
        abundance_scale = str(ds_cfg.get("abundance_scale") or "")
        prot_frac_before = df[sc0].isna().mean(axis=1)
        samp_frac_before = df[sc0].isna().mean(axis=0)

        df0 = filter_missingness(
            df,
            meta,
            protein_missing_max=float(pp.get("missing_threshold_protein", 0.3)),
            sample_missing_max=float(pp.get("missing_threshold_sample", 0.25)),
        )
        sample_cols = [c for c in df0.columns if c not in meta_cols]
        mat_raw = df0[sample_cols].astype(float)

        prot_frac_after = df0[sample_cols].isna().mean(axis=1)
        samp_frac_after = df0[sample_cols].isna().mean(axis=0)

        prot_rows = pd.DataFrame(
            {
                "protein_id": df["protein_id"].astype(str),
                "missing_frac_before_filter": prot_frac_before.values,
            }
        )
        keep_set = set(df0["protein_id"].astype(str))
        prot_rows["kept_after_filter"] = prot_rows["protein_id"].isin(keep_set)
        prot_rows = prot_rows.merge(
            pd.DataFrame(
                {
                    "protein_id": df0["protein_id"].astype(str),
                    "missing_frac_after_filter": prot_frac_after.values,
                }
            ),
            on="protein_id",
            how="left",
        )
        write_tsv(prot_rows, qc_dir / f"{key}_protein_missingness.tsv")

        samp_df = pd.DataFrame(
            {
                "sample_id": list(sc0),
                "missing_frac_before_filter": [float(samp_frac_before.get(c, np.nan)) for c in sc0],
            }
        )
        samp_df["kept_after_filter"] = samp_df["sample_id"].isin(sample_cols)
        samp_df["missing_frac_after_filter"] = samp_df["sample_id"].map(
            lambda sid: float(samp_frac_after[sid]) if sid in samp_frac_after.index else float("nan")
        )
        write_tsv(samp_df, qc_dir / f"{key}_sample_missingness.tsv")

        n_nonpositive = 0
        # log2_ratio: la matriz pública ya es un log2. Transformar otra vez no es un fold change.
        if pp.get("log2_transform", True) and abundance_scale != "log2_ratio":
            mat_log, n_nonpositive = log2_transform_masked(mat_raw)
        else:
            mat_log = mat_raw

        # Matriz observada: las ausencias siguen siendo ausencias. El diseño 3.0
        # exige usarla para inferencia. El diferencial actual sigue leyendo log2norm
        # (imputada) hasta que ese cambio se preregistre.
        out_obs = pd.concat([df0[meta].reset_index(drop=True), mat_log.reset_index(drop=True)], axis=1)
        p_obs = proc / f"{key}_log2_observed.tsv"
        write_tsv(out_obs, p_obs)

        mat_log_imp = knn_impute(mat_log, n_neighbors=int(pp.get("knn_neighbors", 5)))
        out_log = pd.concat([df0[meta].reset_index(drop=True), mat_log_imp.reset_index(drop=True)], axis=1)
        p_log = proc / f"{key}_log2norm.tsv"
        write_tsv(out_log, p_log)

        mat_z = zscore_by_protein(mat_log_imp)
        out_z = pd.concat([df0[meta].reset_index(drop=True), mat_z.reset_index(drop=True)], axis=1)
        p_z = proc / f"{key}_zscore.tsv"
        write_tsv(out_z, p_z)

        if legacy_alias:
            write_tsv(out_z, proc / f"{key}_normalized.tsv")

        summary = {
            "dataset": key,
            "n_proteins_input": int(len(df)),
            "n_proteins_after_missingness_filter": int(len(df0)),
            "n_samples_input": len(sc0),
            "n_samples_after_missingness_filter": len(sample_cols),
            "n_nonpositive_values_masked_for_log2": int(n_nonpositive),
            "outputs": {
                "log2_observed": settings.relative_to_project(p_obs),
                "log2norm": settings.relative_to_project(p_log),
                "zscore": settings.relative_to_project(p_z),
                **(
                    {
                        "legacy_normalized_equals_zscore": settings.relative_to_project(
                            proc / f"{key}_normalized.tsv"
                        )
                    }
                    if legacy_alias
                    else {}
                ),
            },
            "sha256_log2norm": _sha256_file(p_log),
            "sha256_zscore": _sha256_file(p_z),
        }
        (qc_dir / f"{key}_qc_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

        if write_pca:
            pc_log = _pca_sample_coords(mat_log_imp)
            if pc_log is not None:
                pc_log_out = qc_dir / f"{key}_pca_log2norm.tsv"
                write_tsv(pc_log.reset_index().rename(columns={"index": "sample_id"}), pc_log_out)
            pc_z = _pca_sample_coords(mat_z)
            if pc_z is not None:
                pc_z_out = qc_dir / f"{key}_pca_zscore.tsv"
                write_tsv(pc_z.reset_index().rename(columns={"index": "sample_id"}), pc_z_out)

        print(
            f"[qc] {key} -> {settings.relative_to_project(p_obs)} (log2 observado, sin imputar), "
            f"{settings.relative_to_project(p_log)} (log2norm imputado), "
            f"{settings.relative_to_project(p_z)} (zscore)"
            + (f", legacy {key}_normalized.tsv" if legacy_alias else "")
        )
