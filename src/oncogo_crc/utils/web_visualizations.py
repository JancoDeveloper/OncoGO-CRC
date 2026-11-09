"""Gráficos para el sitio MkDocs (matplotlib headless).

Todas las funciones aceptan rutas concretas y devuelven la ruta donde quedó el PNG (o ``None``
si el input no tenía suficiente señal). Diseñadas para que ``web_publish`` las invoque y embeba
los PNG en los Markdown generados con `![alt](figures/<archivo>.png)`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# matplotlib en modo no interactivo (sin DISPLAY); seguro en CI y headless.
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _safe_savefig(fig: "plt.Figure", path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return path


def _short_label(s: str, n: int = 36) -> str:
    s = str(s)
    return s if len(s) <= n else (s[: n - 1] + "…")


# ---------------------------------------------------------------------------
# 1) Top GO (ORA) por -log10(padj)
# ---------------------------------------------------------------------------
def plot_top_go_ora(go_tsv: Path, out_png: Path, *, top_n: int = 15, title: str | None = None) -> Path | None:
    if not go_tsv.exists():
        return None
    df = pd.read_csv(go_tsv, sep="\t")
    if df.empty or "padj" not in df.columns:
        return None
    df = df.copy()
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df = df.dropna(subset=["padj"])
    if df.empty:
        return None
    df["neglog10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df = df.sort_values("neglog10_padj", ascending=False).head(top_n)
    if df.empty:
        return None

    fig, ax = plt.subplots(figsize=(8, max(3, 0.32 * len(df) + 1)))
    y = np.arange(len(df))
    ax.barh(y, df["neglog10_padj"].values, color="#3a6ea5")
    label_col = "go_id"
    if "term_size" in df.columns and "overlap" in df.columns:
        labels = [
            f"{r['go_id']}  (n={int(r['overlap'])}/{int(r['term_size'])})"
            for _, r in df.iterrows()
        ]
    else:
        labels = df[label_col].astype(str).tolist()
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(r"$-\log_{10}(\mathrm{padj})$")
    if title:
        ax.set_title(title, fontsize=10)
    ax.grid(axis="x", linestyle=":", alpha=0.6)
    return _safe_savefig(fig, out_png)


# ---------------------------------------------------------------------------
# 2) GSEA NES — top up y down
# ---------------------------------------------------------------------------
def plot_gsea_top_nes(
    gsea_tsv: Path,
    out_png: Path,
    *,
    top_n: int = 10,
    fdr_threshold: float = 0.25,
    title: str | None = None,
) -> Path | None:
    if not gsea_tsv.exists():
        return None
    df = pd.read_csv(gsea_tsv, sep="\t")
    if df.empty:
        return None
    fdr_col = next((c for c in ("FDR q-val", "FDR", "fdr") if c in df.columns), None)
    nes_col = next((c for c in ("NES", "nes") if c in df.columns), None)
    term_col = next((c for c in ("Term", "term", "go_id") if c in df.columns), None)
    if not (fdr_col and nes_col and term_col):
        return None
    df = df.copy()
    df["_nes"] = pd.to_numeric(df[nes_col], errors="coerce")
    df["_fdr"] = pd.to_numeric(df[fdr_col], errors="coerce")
    df = df.dropna(subset=["_nes", "_fdr"])
    df = df[df["_fdr"] < fdr_threshold]
    if df.empty:
        return None
    df_up = df[df["_nes"] > 0].nlargest(top_n, "_nes")
    df_dn = df[df["_nes"] < 0].nsmallest(top_n, "_nes")
    if df_up.empty and df_dn.empty:
        return None
    sub = pd.concat([df_up, df_dn]).sort_values("_nes")
    colors = ["#b3001b" if x > 0 else "#1a4e8a" for x in sub["_nes"].values]
    fig, ax = plt.subplots(figsize=(8, max(3, 0.32 * len(sub) + 1)))
    y = np.arange(len(sub))
    ax.barh(y, sub["_nes"].values, color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels([_short_label(t, 44) for t in sub[term_col].astype(str)], fontsize=8)
    ax.axvline(0, color="black", linewidth=0.7)
    ax.set_xlabel("NES (Normalized Enrichment Score)")
    if title:
        ax.set_title(title, fontsize=10)
    ax.grid(axis="x", linestyle=":", alpha=0.6)
    return _safe_savefig(fig, out_png)


# ---------------------------------------------------------------------------
# 3) Overlap CPTAC vs PRIDE — scatter log2fc
# ---------------------------------------------------------------------------
def plot_cohort_overlap_scatter(merged_tsv: Path, out_png: Path) -> Path | None:
    if not merged_tsv.exists():
        return None
    df = pd.read_csv(merged_tsv, sep="\t")
    if df.empty or "log2fc_cptac" not in df.columns or "log2fc_pride" not in df.columns:
        return None
    sub = df[["log2fc_cptac", "log2fc_pride"]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(sub) < 5:
        return None
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(sub["log2fc_cptac"], sub["log2fc_pride"], s=8, alpha=0.45, color="#3a6ea5")
    lim = float(np.nanmax(np.abs(sub.values))) * 1.05 if sub.size else 4.0
    if not np.isfinite(lim) or lim <= 0:
        lim = 4.0
    ax.plot([-lim, lim], [-lim, lim], "--", color="grey", linewidth=0.8)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2FC — CPTAC tumor vs normal")
    ax.set_ylabel("log2FC — PRIDE tumor vs normal")
    pr = float(sub["log2fc_cptac"].corr(sub["log2fc_pride"], method="pearson"))
    sr = float(sub["log2fc_cptac"].corr(sub["log2fc_pride"], method="spearman"))
    same_sign = float((np.sign(sub["log2fc_cptac"]) == np.sign(sub["log2fc_pride"])).mean())
    ax.set_title(
        f"Reproducibilidad multi-cohorte (n={len(sub)})\n"
        f"Pearson={pr:.2f} · Spearman={sr:.2f} · mismo signo={same_sign:.2f}",
        fontsize=9,
    )
    ax.grid(linestyle=":", alpha=0.5)
    return _safe_savefig(fig, out_png)


# ---------------------------------------------------------------------------
# 4) PCA por muestra coloreado por condition
# ---------------------------------------------------------------------------
def plot_pca_by_condition(pca_tsv: Path, meta_tsv: Path, out_png: Path, *, title: str | None = None) -> Path | None:
    if not pca_tsv.exists() or not meta_tsv.exists():
        return None
    pca = pd.read_csv(pca_tsv, sep="\t")
    meta = pd.read_csv(meta_tsv, sep="\t")
    if "sample_id" not in pca.columns or "sample_id" not in meta.columns:
        return None
    if not {"PC1", "PC2"}.issubset(pca.columns):
        return None
    m = pca.merge(meta[["sample_id", "condition"]], on="sample_id", how="left")
    m["condition"] = m["condition"].fillna("unknown").astype(str)
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    palette = {"tumor": "#b3001b", "normal": "#1a4e8a", "crc": "#b3001b", "healthy": "#1a4e8a"}
    for cond, sub in m.groupby("condition"):
        color = palette.get(cond.lower(), "#888888")
        ax.scatter(sub["PC1"], sub["PC2"], s=22, alpha=0.7, color=color, label=f"{cond} (n={len(sub)})", edgecolors="white", linewidths=0.5)
    ax.axhline(0, color="black", linewidth=0.3)
    ax.axvline(0, color="black", linewidth=0.3)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    if title:
        ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(linestyle=":", alpha=0.4)
    return _safe_savefig(fig, out_png)


# ---------------------------------------------------------------------------
# 5) Nested CV AUC por fold
# ---------------------------------------------------------------------------
def plot_nested_cv_auc(per_fold_tsv: Path, out_png: Path, *, title: str | None = None) -> Path | None:
    if not per_fold_tsv.exists():
        return None
    df = pd.read_csv(per_fold_tsv, sep="\t")
    if df.empty or "roc_auc" not in df.columns:
        return None
    models = df["model"].unique().tolist() if "model" in df.columns else ["model"]
    fig, ax = plt.subplots(figsize=(max(4.0, 0.9 * len(models) + 2), 4.0))
    box_data = []
    for m in models:
        if "model" in df.columns:
            box_data.append(df[df["model"] == m]["roc_auc"].astype(float).dropna().values)
        else:
            box_data.append(df["roc_auc"].astype(float).dropna().values)
    bp = ax.boxplot(box_data, labels=[_short_label(m, 20) for m in models], showmeans=True, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("#9ec5e8")
        patch.set_alpha(0.75)
    for i, data in enumerate(box_data, start=1):
        ax.scatter(np.full(len(data), i) + (np.random.RandomState(0).uniform(-0.07, 0.07, len(data))),
                   data, color="#1a4e8a", s=20, alpha=0.85, zorder=3)
    ax.set_ylabel("ROC AUC (nested outer fold)")
    ax.set_ylim(0.0, 1.05)
    ax.axhline(0.5, color="grey", linestyle="--", linewidth=0.7)
    if title:
        ax.set_title(title, fontsize=10)
    ax.grid(axis="y", linestyle=":", alpha=0.5)
    return _safe_savefig(fig, out_png)


# ---------------------------------------------------------------------------
# 6) Transfer comparativo CPTAC→PRIDE (uncalibrated vs calibrated)
# ---------------------------------------------------------------------------
def plot_transfer_auc_bars(transfer_json_path: Path, out_png: Path) -> Path | None:
    import json as _json
    if not transfer_json_path.exists():
        return None
    try:
        payload = _json.loads(transfer_json_path.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return None
    rows: list[tuple[str, float]] = []
    for k, v in payload.items():
        if isinstance(v, dict) and "roc_auc" in v:
            try:
                rows.append((str(k), float(v["roc_auc"])))
            except (ValueError, TypeError):
                continue
    if not rows:
        return None
    rows.sort(key=lambda x: x[1], reverse=True)
    labels = [r[0] for r in rows]
    aucs = [r[1] for r in rows]
    colors = ["#b3001b" if "calibrated" in name else "#3a6ea5" for name in labels]
    fig, ax = plt.subplots(figsize=(max(5.0, 0.55 * len(labels) + 2), 3.6))
    y = np.arange(len(labels))
    ax.barh(y, aucs, color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels([_short_label(etiqueta, 32) for etiqueta in labels], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.axvline(0.5, color="grey", linestyle="--", linewidth=0.7)
    ax.set_xlabel("ROC AUC en PRIDE (train: CPTAC)")
    ax.set_title("Transferencia CPTAC → PRIDE", fontsize=10)
    for i, v in enumerate(aucs):
        ax.text(min(0.98, v + 0.01), i, f"{v:.3f}", va="center", fontsize=8)
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    return _safe_savefig(fig, out_png)


# ---------------------------------------------------------------------------
# 7) Volcano enhanced (etiqueta top-k por -log10 padj)
# ---------------------------------------------------------------------------
def plot_volcano_with_labels(
    diff_tsv: Path,
    out_png: Path,
    *,
    top_n_labels: int = 12,
    alpha: float = 0.05,
    title: str | None = None,
) -> Path | None:
    if not diff_tsv.exists():
        return None
    df = pd.read_csv(diff_tsv, sep="\t")
    needed = {"protein_id", "log2fc", "padj"}
    if df.empty or not needed.issubset(df.columns):
        return None
    df = df.copy()
    df["log2fc"] = pd.to_numeric(df["log2fc"], errors="coerce")
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df = df.dropna(subset=["log2fc", "padj"])
    df["nlp"] = -np.log10(df["padj"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sig = df["padj"] < alpha
    ax.scatter(df.loc[~sig, "log2fc"], df.loc[~sig, "nlp"], s=6, alpha=0.35, color="#cccccc")
    ax.scatter(
        df.loc[sig & (df["log2fc"] > 0), "log2fc"],
        df.loc[sig & (df["log2fc"] > 0), "nlp"],
        s=8, alpha=0.55, color="#b3001b", label=f"UP padj<{alpha}",
    )
    ax.scatter(
        df.loc[sig & (df["log2fc"] < 0), "log2fc"],
        df.loc[sig & (df["log2fc"] < 0), "nlp"],
        s=8, alpha=0.55, color="#1a4e8a", label=f"DOWN padj<{alpha}",
    )
    ax.axhline(-np.log10(alpha), color="grey", linestyle="--", linewidth=0.7)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel(r"$\log_2$ FC (tumor vs normal)")
    ax.set_ylabel(r"$-\log_{10}$ padj")
    if title:
        ax.set_title(title, fontsize=10)
    # Labels: top-k by nlp, balanced between up and down.
    df_sig = df[sig].copy()
    top_up = df_sig[df_sig["log2fc"] > 0].nlargest(top_n_labels // 2, "nlp")
    top_dn = df_sig[df_sig["log2fc"] < 0].nlargest(top_n_labels - len(top_up), "nlp")
    for _, r in pd.concat([top_up, top_dn]).iterrows():
        ax.annotate(
            str(r["protein_id"]),
            (r["log2fc"], r["nlp"]),
            fontsize=7, color="black",
            xytext=(3, 3), textcoords="offset points",
        )
    ax.legend(fontsize=8, loc="upper left", frameon=False)
    ax.grid(linestyle=":", alpha=0.4)
    return _safe_savefig(fig, out_png)


# ---------------------------------------------------------------------------
# 8) Heatmap top diff features (samples × top-k genes)
# ---------------------------------------------------------------------------
def plot_top_features_heatmap(
    matrix_tsv: Path,
    meta_tsv: Path,
    diff_tsv: Path,
    out_png: Path,
    *,
    top_n: int = 40,
    title: str | None = None,
) -> Path | None:
    if not matrix_tsv.exists() or not diff_tsv.exists():
        return None
    diff = pd.read_csv(diff_tsv, sep="\t")
    if "padj" not in diff.columns or "protein_id" not in diff.columns:
        return None
    diff["padj"] = pd.to_numeric(diff["padj"], errors="coerce")
    diff = diff.dropna(subset=["padj"]).sort_values("padj").head(top_n)
    if diff.empty:
        return None
    top_ids = diff["protein_id"].astype(str).tolist()
    mat = pd.read_csv(matrix_tsv, sep="\t")
    meta_cols = {"protein_id", "gene_symbol", "uniprot_accession"}
    if "protein_id" not in mat.columns:
        return None
    sub = mat[mat["protein_id"].astype(str).isin(top_ids)].copy()
    if sub.empty:
        return None
    sample_cols = [c for c in sub.columns if c not in meta_cols]
    if not sample_cols:
        return None
    # Orden de muestras por condition si hay metadata.
    sample_order = sample_cols
    col_labels = sample_cols
    if meta_tsv.exists():
        meta = pd.read_csv(meta_tsv, sep="\t")
        if "sample_id" in meta.columns and "condition" in meta.columns:
            meta = meta[["sample_id", "condition"]].copy()
            meta["sample_id"] = meta["sample_id"].astype(str)
            meta = meta[meta["sample_id"].isin(sample_cols)]
            meta = meta.sort_values(["condition", "sample_id"])
            sample_order = meta["sample_id"].tolist()
            col_labels = [f"{sid} ({cond[:1]})" for sid, cond in zip(meta["sample_id"], meta["condition"])]
    sub = sub.set_index("protein_id")[sample_order]
    # Z-score por fila para visualización (no afecta TSV original).
    z = sub.subtract(sub.mean(axis=1), axis=0).divide(sub.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    fig, ax = plt.subplots(figsize=(max(6.0, 0.06 * len(sample_order) + 4), max(4.0, 0.18 * len(top_ids) + 2)))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
    ax.set_yticks(np.arange(len(top_ids)))
    ax.set_yticklabels(top_ids, fontsize=7)
    if len(sample_order) <= 60:
        ax.set_xticks(np.arange(len(sample_order)))
        ax.set_xticklabels(col_labels, rotation=90, fontsize=6)
    else:
        ax.set_xticks([])
    fig.colorbar(im, ax=ax, shrink=0.6, label="z-score (por gen)")
    if title:
        ax.set_title(title, fontsize=10)
    return _safe_savefig(fig, out_png)


# ---------------------------------------------------------------------------
# 9) Cobertura clínica como barplot horizontal por dataset
# ---------------------------------------------------------------------------
def plot_clinical_demographics_bars(
    summary: dict[str, Any],
    out_png: Path,
    *,
    dataset: str = "cptac_colorectal",
    max_vars: int = 6,
) -> Path | None:
    """Barplot horizontal de conteos categóricos (estadio, etnia, geografía) para un dataset."""
    block = (summary.get("datasets") or {}).get(dataset)
    if not isinstance(block, dict):
        return None
    vars_block = block.get("variables") or {}
    priority = [
        "stage_group",
        "stage_bin",
        "hispanic_latino_bin",
        "geographic_region_bin",
        "msi_status",
        "age_bin",
    ]
    chosen: list[tuple[str, dict[str, int]]] = []
    for v in priority:
        if v in vars_block and vars_block[v].get("counts"):
            chosen.append((v, vars_block[v]["counts"]))
        if len(chosen) >= max_vars:
            break
    if not chosen:
        return None

    n = len(chosen)
    fig, axes = plt.subplots(1, n, figsize=(max(4.0, 2.2 * n), 3.8))
    if n == 1:
        axes = [axes]
    palette = ["#3a6ea5", "#b3001b", "#2a9d8f", "#e9c46a", "#9c89b8", "#6d6875"]
    for ax, (var, counts) in zip(axes, chosen):
        items = sorted(counts.items(), key=lambda x: -x[1])
        labels = [str(k) for k, _ in items]
        vals = [float(v) for _, v in items]
        y = np.arange(len(labels))
        ax.barh(y, vals, color=palette[: len(labels)], height=0.65)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("n (muestras tumor)", fontsize=8)
        ax.set_title(var, fontsize=9)
        ax.grid(axis="x", linestyle=":", alpha=0.4)
    fig.suptitle(f"Demografía clínica — {dataset} (muestras tumor)", fontsize=10, y=1.02)
    fig.tight_layout()
    return _safe_savefig(fig, out_png)


def plot_clinical_coverage(coverage: dict[str, dict[str, float]], out_png: Path) -> Path | None:
    if not coverage:
        return None
    rows: list[tuple[str, str, float]] = []
    for ds, varmap in coverage.items():
        for var, frac in varmap.items():
            try:
                rows.append((ds, str(var), float(frac)))
            except (ValueError, TypeError):
                continue
    if not rows:
        return None
    df = pd.DataFrame(rows, columns=["dataset", "var", "coverage"])
    pivot = df.pivot(index="var", columns="dataset", values="coverage").fillna(0.0)
    fig, ax = plt.subplots(figsize=(max(5.0, 1.0 * pivot.shape[1] + 2), max(3.0, 0.32 * pivot.shape[0] + 1)))
    y = np.arange(pivot.shape[0])
    width = 0.9 / max(1, pivot.shape[1])
    palette = ["#3a6ea5", "#b3001b", "#2a9d8f", "#e9c46a", "#9c89b8"]
    for i, ds in enumerate(pivot.columns):
        ax.barh(y + (i - (pivot.shape[1] - 1) / 2) * width, pivot[ds].values, height=width,
                color=palette[i % len(palette)], label=str(ds))
    ax.set_yticks(y)
    ax.set_yticklabels(pivot.index.tolist(), fontsize=8)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Fracción de muestras con valor informativo (≠ NaN / unknown)")
    ax.set_title("Cobertura clínica por variable y dataset", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    return _safe_savefig(fig, out_png)
