from __future__ import annotations

import json

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv, write_tsv
from oncogo_crc.utils.ml_utils import train_prioritization_models


def run(settings: Settings) -> None:
    """Modelos explicables (L1-logit + RandomForest) para priorizar proteínas."""
    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")
    proc = settings.resolve(settings.config["paths"]["processed_dir"])
    mat_path = proc / f"{primary}_normalized.tsv"
    meta_path = settings.sample_metadata_path(primary)
    if not mat_path.exists() or not meta_path.exists():
        print("[ml] faltan matriz o metadata; omitiendo.")
        return
    df = read_tsv(mat_path)
    sm = read_tsv(meta_path)
    meta = {"protein_id", "gene_symbol", "uniprot_accession"}
    sample_cols = [c for c in df.columns if c not in meta]
    mat = df.set_index("protein_id")[sample_cols].astype(float)
    labels = sm.set_index("sample_id")["condition"]
    samples = mat.columns.intersection(labels.index)
    lab = labels.loc[samples]
    mask = lab.isin(["tumor", "normal"])
    samples = lab.index[mask]
    if len(samples) < 6:
        print("[ml] pocas muestras etiquetadas tumor/normal; omitiendo.")
        return
    mat_sub = mat[samples]
    lab = lab.loc[samples]
    X = mat_sub.T
    y = (lab == "tumor").astype(int).to_numpy()
    if len(np.unique(y)) < 2:
        print("[ml] se requieren clases tumor y normal en metadata; omitiendo.")
        return
    ml_cfg = settings.config.get("ml", {})
    metrics, scores = train_prioritization_models(
        X,
        y,
        test_size=float(ml_cfg.get("test_size", 0.25)),
        seed=int(settings.config.get("project", {}).get("seed", 42)),
        rf_n_estimators=int(ml_cfg.get("random_forest_n_estimators", 200)),
    )
    out_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "ml"
    ranking = scores.reset_index()
    ranking = ranking.rename(columns={ranking.columns[0]: "protein_id"})
    write_tsv(ranking, out_dir / f"{primary}_feature_ranking.tsv")
    (out_dir / f"{primary}_model_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    metrics_flat = [{"model": k, **v} for k, v in metrics.items()]
    write_tsv(pd.DataFrame(metrics_flat), out_dir / f"{primary}_model_metrics.tsv")
    print(f"[ml] métricas y ranking -> results/ml/{primary}_*")
