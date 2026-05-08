from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def volcano_plot(
    diff: pd.DataFrame,
    out_path: Path,
    *,
    fc_col: str = "log2fc",
    p_col: str = "padj",
    label_col: str = "protein_id",
    top_n: int = 12,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = diff.copy()
    df["neglog10q"] = -np.log10(np.clip(df[p_col].astype(float), 1e-300, None))
    plt.figure(figsize=(8, 6))
    sig = df["padj"] < 0.05
    sns.scatterplot(data=df, x=fc_col, y="neglog10q", hue=sig, palette={True: "tab:red", False: "gray"}, s=26)
    top = df.nsmallest(top_n, p_col)
    for _, r in top.iterrows():
        plt.text(r[fc_col], r["neglog10q"], str(r[label_col]), fontsize=7, ha="left")
    plt.axhline(-np.log10(0.05), color="k", ls="--", lw=0.8)
    plt.title("Volcano (t-test + FDR)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=160)
    plt.close()
