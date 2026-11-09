from __future__ import annotations

import json
from hashlib import sha256

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv, write_tsv
from oncogo_crc.utils.plots import volcano_plot
from oncogo_crc.utils.stats import build_paired_sample_order, differential_ttest


def _optional_comparison_viable(spec: dict, sm: pd.DataFrame, default_min_n: int) -> bool:
    """True si la columna existe y hay suficientes muestras en cada grupo."""
    col = spec["condition_column"]
    if col not in sm.columns:
        return False
    min_n = int(spec.get("min_samples_per_group", default_min_n))
    s = sm[col].astype(str)
    ga = (s == str(spec["group_a"])).sum()
    gb = (s == str(spec["group_b"])).sum()
    return int(ga) >= min_n and int(gb) >= min_n


def run(settings: Settings) -> None:
    """Modelado estadístico simple (t-test + FDR) por cada comparación declarada (matriz log2norm)."""
    res_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "differential"
    fig_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "figures"
    dc = settings.config.get("differential", {})
    meta_cols = {"protein_id", "gene_symbol", "uniprot_accession"}
    min_fc = float(dc.get("min_abs_log2fc", 0.0) or 0.0)
    auto_opt = bool(dc.get("auto_run_optional_comparisons", False))
    default_paired_col = str(dc.get("paired_column", "participant_id") or "").strip()

    comps = settings.comparisons.get("comparisons", {})
    for comp_key, spec in comps.items():
        if spec.get("optional") and not auto_opt:
            continue
        ds_key = spec["dataset"]
        mat_path = settings.processed_matrix_path(ds_key, "log2norm")
        if not mat_path.exists():
            continue
        df = read_tsv(mat_path)
        meta_path = settings.sample_metadata_path(ds_key)
        if not meta_path.exists():
            print(f"[differential] sin metadata: {meta_path}")
            continue
        sm = read_tsv(meta_path)
        if spec["condition_column"] not in sm.columns:
            if spec.get("optional"):
                print(
                    f"[differential] {comp_key}: omitido (columna {spec['condition_column']} ausente)."
                )
            else:
                print(f"[differential] columna {spec['condition_column']} ausente en {meta_path}")
            continue
        if spec.get("optional") and auto_opt and not _optional_comparison_viable(
            spec, sm, int(dc.get("optional_comparison_min_n_per_group", 3))
        ):
            print(
                f"[differential] {comp_key}: omitido (pocas muestras en "
                f"{spec['group_a']}/{spec['group_b']} o estratos ausentes)."
            )
            continue
        sample_cols = [c for c in df.columns if c not in meta_cols]
        labels = sm.set_index("sample_id")[spec["condition_column"]]
        mat = df.set_index("protein_id")[sample_cols]
        lab = labels.loc[labels.isin([spec["group_a"], spec["group_b"]])].dropna()
        mat = mat[[c for c in mat.columns if c in lab.index]]

        paired_order: list[tuple[str, str]] | None = None
        pcol = str(spec.get("paired_by") or "").strip() or default_paired_col
        use_paired = bool(spec.get("use_paired", True))
        if use_paired and pcol:
            paired_order = build_paired_sample_order(
                sm,
                condition_column=spec["condition_column"],
                group_a=spec["group_a"],
                group_b=spec["group_b"],
                participant_col=pcol,
            )
            if paired_order and len(paired_order) >= 2:
                print(
                    f"[differential] {comp_key}: diseño pareado ({len(paired_order)} pares, columna {pcol})."
                )
            else:
                paired_order = None

        diff = differential_ttest(
            mat,
            lab,
            group_a=spec["group_a"],
            group_b=spec["group_b"],
            alpha=float(dc.get("alpha", 0.05)),
            fdr_method=dc.get("fdr_method", "fdr_bh"),
            min_abs_log2fc=min_fc,
            paired_sample_order=paired_order,
        )
        diff["comparison_id"] = comp_key
        diff["matrix_input"] = str(settings.relative_to_project(mat_path))
        out_tsv = res_dir / f"{comp_key}_diff.tsv"
        write_tsv(diff, out_tsv)
        raw = out_tsv.read_bytes()
        audit = {
            "comparison_id": comp_key,
            "dataset": ds_key,
            "matrix_kind": "log2norm",
            "matrix_path": str(settings.relative_to_project(mat_path)),
            "sha256_diff_output": sha256(raw).hexdigest(),
            "n_rows_diff": int(len(diff)),
            "n_samples_in_contrast": int(len(lab)),
            "stat_method": str(diff["stat_method"].iloc[0]) if len(diff) else None,
        }
        (res_dir / f"{comp_key}_diff_input_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
        volcano_plot(diff, fig_dir / f"{comp_key}_volcano.png")
        print(f"[differential] {comp_key} -> {settings.relative_to_project(out_tsv)}")
