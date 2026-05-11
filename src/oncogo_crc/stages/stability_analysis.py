from __future__ import annotations

from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings


def _sig_set(path: Path) -> set[str]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep="\t")
    if "protein_id" not in df.columns or "padj" not in df.columns:
        return set()
    sdf = df[df["padj"].notna() & (df["padj"] < 0.05)]
    return set(sdf["protein_id"].astype(str))


def run(settings: Settings) -> None:
    comps = settings.comparisons.get("comparisons", {})
    diff_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "differential"
    out_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "stability"
    out_dir.mkdir(parents=True, exist_ok=True)

    core_by_dataset: dict[str, list[str]] = {}
    stratified: list[str] = []
    for name, spec in comps.items():
        if spec.get("optional", False):
            stratified.append(name)
        else:
            core_by_dataset.setdefault(spec["dataset"], []).append(name)

    rows: list[dict[str, object]] = []
    for strat in stratified:
        s_spec = comps[strat]
        ds = s_spec["dataset"]
        s_set = _sig_set(diff_dir / f"{strat}_diff.tsv")
        for core in core_by_dataset.get(ds, []):
            c_set = _sig_set(diff_dir / f"{core}_diff.tsv")
            inter = len(s_set & c_set)
            union = len(s_set | c_set)
            jaccard = (inter / union) if union else 0.0
            rows.append(
                {
                    "dataset": ds,
                    "core_comparison": core,
                    "stratified_comparison": strat,
                    "n_sig_core": len(c_set),
                    "n_sig_stratified": len(s_set),
                    "n_intersection": inter,
                    "jaccard": jaccard,
                }
            )

    df = pd.DataFrame(rows)
    tsv = out_dir / "stability_summary.tsv"
    js = out_dir / "stability_summary.json"
    if not df.empty:
        df.to_csv(tsv, sep="\t", index=False)
        js.write_text(df.to_json(orient="records", indent=2), encoding="utf-8")
    else:
        pd.DataFrame(
            columns=[
                "dataset",
                "core_comparison",
                "stratified_comparison",
                "n_sig_core",
                "n_sig_stratified",
                "n_intersection",
                "jaccard",
            ]
        ).to_csv(tsv, sep="\t", index=False)
        js.write_text("[]", encoding="utf-8")
    print(f"[stability] -> {settings.relative_to_project(tsv)}")
