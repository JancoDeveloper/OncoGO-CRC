from __future__ import annotations


from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv, write_tsv
from oncogo_crc.utils.plots import volcano_plot
from oncogo_crc.utils.stats import differential_ttest


def run(settings: Settings) -> None:
    """Modelado estadístico simple (t-test + FDR) por cada comparación declarada."""
    proc = settings.resolve(settings.config["paths"]["processed_dir"])
    res_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "differential"
    fig_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "figures"
    dc = settings.config.get("differential", {})
    meta_cols = {"protein_id", "gene_symbol", "uniprot_accession"}

    comps = settings.comparisons.get("comparisons", {})
    for comp_key, spec in comps.items():
        if spec.get("optional"):
            continue
        ds_key = spec["dataset"]
        mat_path = proc / f"{ds_key}_normalized.tsv"
        if not mat_path.exists():
            continue
        df = read_tsv(mat_path)
        meta_path = settings.sample_metadata_path(ds_key)
        if not meta_path.exists():
            print(f"[differential] sin metadata: {meta_path}")
            continue
        sm = read_tsv(meta_path)
        if spec["condition_column"] not in sm.columns:
            print(f"[differential] columna {spec['condition_column']} ausente en {meta_path}")
            continue
        sample_cols = [c for c in df.columns if c not in meta_cols]
        labels = sm.set_index("sample_id")[spec["condition_column"]]
        mat = df.set_index("protein_id")[sample_cols]
        lab = labels.loc[labels.isin([spec["group_a"], spec["group_b"]])].dropna()
        mat = mat[[c for c in mat.columns if c in lab.index]]
        diff = differential_ttest(
            mat,
            lab,
            group_a=spec["group_a"],
            group_b=spec["group_b"],
            alpha=float(dc.get("alpha", 0.05)),
            fdr_method=dc.get("fdr_method", "fdr_bh"),
        )
        out_tsv = res_dir / f"{comp_key}_diff.tsv"
        write_tsv(diff, out_tsv)
        volcano_plot(diff, fig_dir / f"{comp_key}_volcano.png")
        print(f"[differential] {comp_key} -> {settings.relative_to_project(out_tsv)}")
