"""t moderada, treat, y que un log-ratio no se vuelva a transformar."""

from __future__ import annotations

import numpy as np

from oncogo_crc.certification.contract import h1_pathway_clause, h1_status
from oncogo_crc.utils.moderated import (
    paired_covariate_stats,
    paired_observed_stats,
    treat_pvalues,
)
from oncogo_crc.utils.quant_columns import columns_for_family, is_summary_protein, quant_column_kind


def test_treat_rejects_only_effects_past_the_threshold():
    p_null = treat_pvalues(np.array([0.0]), np.array([0.2]), np.array([40.0]), 0.585)
    p_small = treat_pvalues(np.array([0.2]), np.array([0.05]), np.array([40.0]), 0.585)
    p_large = treat_pvalues(np.array([2.0]), np.array([0.2]), np.array([40.0]), 0.585)
    assert p_null[0] > 0.5
    assert p_small[0] > 0.05
    assert p_large[0] < 1e-6


def test_paired_difference_keeps_a_negative_log_ratio():
    diffs = np.array([[-0.5, -0.4, -0.6]])
    coef, _s2, _df, _stdev, n = paired_observed_stats(diffs, min_pairs=3)
    assert n[0] == 3
    assert coef[0] == np.mean([-0.5, -0.4, -0.6])
    assert coef[0] < 0


def test_missing_pair_is_not_imputed():
    diffs = np.array([[1.0, np.nan, 1.0, 1.0]])
    coef, _s2, _df, _stdev, n = paired_observed_stats(diffs, min_pairs=3)
    assert n[0] == 3
    assert coef[0] == 1.0


def test_unshared_columns_are_not_mixed_with_log_ratio():
    cols = [
        "s1 Unshared Log Ratio",
        "s1 Log Ratio",
        "s2 Unshared Log Ratio",
    ]
    assert quant_column_kind("s1 Log Ratio") == "log_ratio"
    assert columns_for_family(cols, "unshared_log_ratio") == [
        "s1 Unshared Log Ratio",
        "s2 Unshared Log Ratio",
    ]
    assert is_summary_protein("Mean")
    assert not is_summary_protein("A1BG")


def test_effect_fully_explained_by_covariate_collapses_at_zero():
    # La diferencia es exactamente 2*x: sin diferencia de composición no queda efecto.
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    diffs = (2.0 * x).reshape(1, -1)
    coef, _s2, _df, _stdev, n, slope, _t = paired_covariate_stats(diffs, x, min_pairs=3, at=0.0)
    assert n[0] == 5
    assert abs(slope[0] - 2.0) < 1e-9
    assert abs(coef[0]) < 1e-9
    # Evaluar en la media devolvería la media de las diferencias y no ajustaría nada.
    coef_mean, *_ = paired_covariate_stats(diffs, x, min_pairs=3, at=float(x.mean()))
    assert abs(coef_mean[0] - float(np.mean(2.0 * x))) < 1e-9


def test_effect_independent_of_covariate_is_preserved():
    x = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    diffs = np.array([[1.5, 1.5, 1.5, 1.5, 1.5]])
    coef, _s2, _df, _stdev, _n, slope, _t = paired_covariate_stats(diffs, x, min_pairs=3, at=0.0)
    assert abs(slope[0]) < 1e-9
    assert abs(coef[0] - 1.5) < 1e-9


def test_h1_needs_both_clauses():
    assert h1_pathway_clause(0.80) == "rechaza_clausula_via"
    assert h1_pathway_clause(0.81) == "sostiene_clausula_via"
    assert h1_status("sostiene_clausula_proteina", "sostiene_clausula_via") == "sostiene_h1"
    assert h1_status("sostiene_clausula_proteina", "rechaza_clausula_via") == "rechaza_h1"
    assert h1_status("sostiene_clausula_proteina", "no_evaluable") == "parcial_solo_proteina"
