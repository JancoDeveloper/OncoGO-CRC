"""QC log2 enmascarado y diferencial pareado (Fase 1)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from oncogo_crc.utils.qc import log2_transform_masked
from oncogo_crc.utils.stats import differential_ttest


def test_log2_masks_nonpositive():
    mat = pd.DataFrame([[4.0, 0.0, -1.0], [8.0, 2.0, 1.0]])
    out, n_bad = log2_transform_masked(mat)
    assert n_bad == 2
    assert np.isnan(out.iloc[0, 1]) and np.isnan(out.iloc[0, 2])
    assert out.iloc[0, 0] == 2.0


def test_paired_differential_mean_of_paired_diffs():
    mat = pd.DataFrame(
        [[1.0, 2.0, 3.0, 4.0], [10.0, 12.0, 9.0, 11.0]],
        index=["A", "B"],
        columns=["s1", "s2", "s3", "s4"],
    )
    labels = pd.Series({"s1": "tumor", "s2": "normal", "s3": "tumor", "s4": "normal"})
    pairs = [("s1", "s2"), ("s3", "s4")]
    res = differential_ttest(
        mat,
        labels,
        group_a="tumor",
        group_b="normal",
        paired_sample_order=pairs,
    )
    assert res.loc[res["protein_id"] == "A", "log2fc"].iloc[0] == -1.0
    assert res["stat_method"].iloc[0] == "paired_ttest_mean_diff"
