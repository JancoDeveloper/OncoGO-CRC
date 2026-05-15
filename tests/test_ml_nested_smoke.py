"""Smoke tests for nested CV / transfer helpers (no full pipeline)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from oncogo_crc.utils.ml_utils import (
    nested_cross_validate_models,
    normalized_matrix_to_samples_x_genes,
    recommend_nested_cv_splits,
    transfer_train_test_eval,
)
from oncogo_crc.utils.reproducibility_overlap import cptac_pride_overlap_tables


def test_transfer_train_test_eval_ok():
    rng = np.random.default_rng(0)
    Xa = pd.DataFrame(rng.standard_normal((40, 10)), columns=[f"g{i}" for i in range(10)])
    Xb = pd.DataFrame(rng.standard_normal((30, 10)), columns=[f"g{i}" for i in range(10)])
    ya = np.array([0] * 20 + [1] * 20)
    yb = np.array([0] * 15 + [1] * 15)
    out = transfer_train_test_eval(Xa, ya, Xb, yb, seed=2, rf_n_estimators=50)
    assert out["status"] == "ok"
    assert "sparse_linear_l1" in out
    assert out["calibration"]["applied"] is False


def test_transfer_with_calibration():
    rng = np.random.default_rng(4)
    Xa = pd.DataFrame(rng.standard_normal((44, 8)), columns=[f"g{i}" for i in range(8)])
    Xb = pd.DataFrame(rng.standard_normal((26, 8)), columns=[f"g{i}" for i in range(8)])
    ya = np.array([0] * 22 + [1] * 22)
    yb = np.array([0] * 13 + [1] * 13)
    out = transfer_train_test_eval(
        Xa, ya, Xb, yb, seed=5, rf_n_estimators=40, calibration_method="sigmoid", calibration_cv=3
    )
    assert out["status"] == "ok"
    assert out["calibration"]["applied"] is True
    assert "sparse_linear_l1_calibrated" in out


def test_nested_cv_runs_small_grid():
    rng = np.random.default_rng(1)
    X = rng.standard_normal((36, 12))
    y = np.array([0] * 18 + [1] * 18)
    per_fold, summary = nested_cross_validate_models(
        X,
        y,
        outer_splits=3,
        inner_splits=2,
        seed=3,
        rf_n_estimators=40,
        linear_alpha_grid=[0.02, 0.08],
        rf_max_depth_grid=[None, 6],
        rf_min_samples_leaf_grid=[1],
        run_l1=True,
        run_elastic_net=True,
        run_rf=True,
    )
    assert len(per_fold) == 9
    assert "sparse_linear_l1_nested" in summary


def test_nested_cv_only_rf():
    rng = np.random.default_rng(6)
    X = rng.standard_normal((30, 8))
    y = np.array([0] * 15 + [1] * 15)
    per_fold, summary = nested_cross_validate_models(
        X,
        y,
        outer_splits=3,
        inner_splits=2,
        seed=7,
        rf_n_estimators=35,
        linear_alpha_grid=[0.02],
        rf_max_depth_grid=[4],
        rf_min_samples_leaf_grid=[1],
        run_l1=False,
        run_elastic_net=False,
        run_rf=True,
    )
    assert len(per_fold) == 3
    assert set(per_fold["model"].unique()) == {"random_forest_nested"}


def test_recommend_nested_splits_caps_by_class():
    y = np.array([0] * 5 + [1] * 8)
    outer, inner, _note = recommend_nested_cv_splits(y, outer_requested=5, inner_requested=4)
    assert outer <= 5
    assert inner <= 4
    assert outer >= 2 and inner >= 2


def test_cptac_pride_overlap_merge_skips_nan_join_key():
    """Sin clave UniProt válida, el merge no debe emparejar todos los NaN entre sí (producto cartesiano)."""
    c = pd.DataFrame(
        {
            "protein_id": ["TP53", "KRAS", "EGFR"],
            "padj": [0.01, 0.02, 0.5],
            "log2fc": [1.0, -0.5, 0.1],
        }
    )
    p = pd.DataFrame(
        {
            "protein_id": ["O001", "O002"],
            "padj": [0.03, 0.04],
            "log2fc": [0.8, -0.3],
        }
    )
    merged, summary = cptac_pride_overlap_tables(c, p, gene_to_uni={}, alpha=0.05)
    assert len(merged) == 0
    n_bridge = int(summary.loc[summary["metric"] == "n_intersection_uniprot_bridge", "value"].iloc[0])
    assert n_bridge == 0


def test_cptac_pride_overlap_uniprot_bridge_one_row():
    gene_to_uni = {"TP53": "P04637"}
    c = pd.DataFrame(
        {"protein_id": ["TP53"], "padj": [0.01], "log2fc": [1.2]},
    )
    p = pd.DataFrame(
        {"protein_id": ["P04637"], "padj": [0.02], "log2fc": [1.0]},
    )
    merged, summary = cptac_pride_overlap_tables(c, p, gene_to_uni=gene_to_uni, alpha=0.05)
    assert len(merged) == 1
    n_bridge = int(summary.loc[summary["metric"] == "n_intersection_uniprot_bridge", "value"].iloc[0])
    assert n_bridge == 1


def test_normalized_matrix_resolves_uniprot_via_map():
    df = pd.DataFrame(
        {
            "protein_id": ["P04637", "P04637-2"],
            "gene_symbol": [pd.NA, pd.NA],
            "s1": [1.0, 3.0],
            "s2": [2.0, 4.0],
        }
    )
    utg = {"P04637": "TP53"}
    X = normalized_matrix_to_samples_x_genes(df, uniprot_to_gene=utg)
    assert list(X.columns) == ["TP53"]
    assert np.allclose(X.loc["s1", "TP53"], 2.0)
    assert np.allclose(X.loc["s2", "TP53"], 3.0)
