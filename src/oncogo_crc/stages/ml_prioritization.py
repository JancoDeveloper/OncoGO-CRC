from __future__ import annotations

import json

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.ids import effective_protein_mapping
from oncogo_crc.utils.io import read_tsv, write_tsv
from oncogo_crc.utils import run_manifest
from oncogo_crc.utils.reproducibility_overlap import uniprot_to_gene_dict
from oncogo_crc.utils.ml_utils import (
    align_two_cohort_matrices,
    bootstrap_linear_feature_stability,
    nested_cross_validate_models,
    nested_cv_model_flags,
    normalized_matrix_to_samples_x_genes,
    recommend_nested_cv_splits,
    train_prioritization_models,
    transfer_train_test_eval,
)


def run(settings: Settings) -> None:
    """Modelos explicables (L1-logit + Elastic Net SGD + RF) y validación rigurosa (nested CV, transferencia)."""
    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")
    mat_path = settings.processed_matrix_path(primary, "zscore")
    meta_path = settings.sample_metadata_path(primary)
    if not mat_path.exists() or not meta_path.exists():
        print("[ml] faltan matriz o metadata; omitiendo.")
        if bool((settings.config.get("ml") or {}).get("external_transfer", {}).get("enabled", False)):
            run_manifest.record_stage_meta(
                settings,
                "ml_transfer",
                executed=False,
                reason="ml_stage_skipped_missing_primary_matrix",
            )
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
        if bool((settings.config.get("ml") or {}).get("external_transfer", {}).get("enabled", False)):
            run_manifest.record_stage_meta(
                settings,
                "ml_transfer",
                executed=False,
                reason="ml_stage_skipped_few_labeled_samples",
            )
        return
    mat_sub = mat[samples]
    lab = lab.loc[samples]
    X = mat_sub.T
    y = (lab == "tumor").astype(int).to_numpy()
    if len(np.unique(y)) < 2:
        print("[ml] se requieren clases tumor y normal en metadata; omitiendo.")
        if bool((settings.config.get("ml") or {}).get("external_transfer", {}).get("enabled", False)):
            run_manifest.record_stage_meta(
                settings,
                "ml_transfer",
                executed=False,
                reason="ml_stage_skipped_single_class",
            )
        return
    ml_cfg = settings.config.get("ml", {}) or {}
    seed = int(settings.config.get("project", {}).get("seed", 42))
    use_shap = bool(ml_cfg.get("use_shap", False))
    cv_folds = int(ml_cfg.get("cv_folds", 0) or 0)
    shap_bg = int(ml_cfg.get("shap_background_max", 256) or 256)
    models = [str(m).lower() for m in (ml_cfg.get("models") or ["lasso", "random_forest"])]
    include_elastic_net = "elastic_net" in models or "elasticnet" in models
    linear_alpha = float(ml_cfg.get("linear_alpha", 0.02) or 0.02)

    def _train(with_shap: bool) -> tuple[dict, pd.DataFrame]:
        return train_prioritization_models(
            X,
            y,
            test_size=float(ml_cfg.get("test_size", 0.25)),
            seed=seed,
            rf_n_estimators=int(ml_cfg.get("random_forest_n_estimators", 200)),
            use_shap=with_shap,
            shap_background_max=shap_bg,
            cv_folds=cv_folds,
            include_elastic_net=include_elastic_net,
            linear_alpha=linear_alpha,
        )

    try:
        metrics, scores = _train(use_shap)
    except ImportError as e:
        if use_shap and "shap" in str(e).lower():
            print(f"[ml] SHAP no disponible ({e}); ejecutando sin SHAP. Instala: pip install -e \".[explain]\"")
            metrics, scores = _train(False)
        else:
            raise
    out_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "ml"
    out_dir.mkdir(parents=True, exist_ok=True)
    ranking = scores.reset_index()
    ranking = ranking.rename(columns={ranking.columns[0]: "protein_id"})
    write_tsv(ranking, out_dir / f"{primary}_feature_ranking.tsv")

    nest_cfg = ml_cfg.get("nested_cv") or {}
    if bool(nest_cfg.get("enabled", False)):
        outer_req = int(nest_cfg.get("outer_splits", 5))
        inner_req = int(nest_cfg.get("inner_splits", 3))
        outer, inner, split_note = recommend_nested_cv_splits(y, outer_requested=outer_req, inner_requested=inner_req)
        run_l1, run_en, run_rf = nested_cv_model_flags(nest_cfg)
        alphas = nest_cfg.get("linear_alpha_grid") or [0.005, 0.02, 0.08]
        if not isinstance(alphas, list):
            alphas = [0.005, 0.02, 0.08]
        depths_raw = nest_cfg.get("rf_max_depth_grid")
        if isinstance(depths_raw, list):
            rf_depths: list[int | None] = []
            for d in depths_raw:
                if d is None or (isinstance(d, str) and str(d).lower() in ("null", "none", "")):
                    rf_depths.append(None)
                else:
                    rf_depths.append(int(d))
        else:
            rf_depths = [None, 12, 20]
        leaves = nest_cfg.get("rf_min_samples_leaf_grid") or [1, 2]
        if not isinstance(leaves, list):
            leaves = [1, 2]
        leaves = [int(x) for x in leaves]
        rf_nested = int(ml_cfg.get("random_forest_n_estimators", 200))
        if nest_cfg.get("rf_n_estimators") is not None:
            rf_nested = int(nest_cfg["rf_n_estimators"])
        if run_rf and not run_l1 and not run_en and len(y) < 50:
            rf_nested = min(rf_nested, 120)
        per_fold, nest_summary = nested_cross_validate_models(
            X.values,
            y,
            outer_splits=outer,
            inner_splits=inner,
            seed=seed,
            rf_n_estimators=rf_nested,
            linear_alpha_grid=[float(a) for a in alphas],
            rf_max_depth_grid=rf_depths,
            rf_min_samples_leaf_grid=leaves,
            run_l1=run_l1,
            run_elastic_net=run_en,
            run_rf=run_rf,
            split_note=split_note,
        )
        write_tsv(per_fold, out_dir / f"{primary}_nested_cv_per_fold.tsv")
        summary_rows = [{"model": k, **v} for k, v in nest_summary.items() if isinstance(v, dict)]
        write_tsv(pd.DataFrame(summary_rows), out_dir / f"{primary}_nested_cv_summary.tsv")
        (out_dir / f"{primary}_nested_cv_summary.json").write_text(
            json.dumps(nest_summary, indent=2), encoding="utf-8"
        )
        metrics["nested_cv"] = nest_summary
        which = []
        if run_l1:
            which.append("L1")
        if run_en:
            which.append("EN")
        if run_rf:
            which.append("RF")
        print(
            f"[ml] nested CV ({outer}×outer, {inner}×inner; modelos: {','.join(which)}) -> {primary}_nested_cv_*.tsv"
        )

    xfer = ml_cfg.get("external_transfer") or {}
    xfer_enabled = bool(xfer.get("enabled", False))
    if not xfer_enabled:
        run_manifest.record_stage_meta(
            settings,
            "ml_transfer",
            executed=False,
            reason="disabled_in_config",
        )
    if xfer_enabled:
        train_ds = str(xfer.get("train_dataset", primary))
        test_ds = str(xfer.get("test_dataset", "pride_crc_meta"))
        m_tr = settings.processed_matrix_path(train_ds, "zscore")
        m_te = settings.processed_matrix_path(test_ds, "zscore")
        meta_tr = settings.sample_metadata_path(train_ds)
        meta_te = settings.sample_metadata_path(test_ds)
        for stale in out_dir.glob("*_train_*_test_transfer.json"):
            stale.unlink(missing_ok=True)
        for stale in out_dir.glob("*_train_*_test_sample_labels.tsv"):
            stale.unlink(missing_ok=True)
        if m_tr.exists() and m_te.exists() and meta_tr.exists() and meta_te.exists():
            dtr = read_tsv(m_tr)
            dte = read_tsv(m_te)
            sm_tr = read_tsv(meta_tr)
            sm_te = read_tsv(meta_te)
            map_df = effective_protein_mapping(settings)
            utg = uniprot_to_gene_dict(map_df)
            Xtr = normalized_matrix_to_samples_x_genes(dtr, uniprot_to_gene=utg or None)
            Xte = normalized_matrix_to_samples_x_genes(dte, uniprot_to_gene=utg or None)
            n_tr_feat = int(Xtr.shape[1])
            n_te_feat = int(Xte.shape[1])
            print(
                f"[ml] transfer genes (pre-intersección): train_cols={n_tr_feat}, test_cols={n_te_feat}"
            )
            Xtr_a, Xte_a, genes = align_two_cohort_matrices(Xtr, Xte)
            print(
                f"[ml] transfer genes (post-intersección canónica): common={len(genes)}"
            )
            lt = sm_tr.set_index("sample_id")["condition"]
            lte = sm_te.set_index("sample_id")["condition"]
            st = Xtr_a.index.intersection(lt.index)
            ste = Xte_a.index.intersection(lte.index)
            lt = lt.loc[st]
            lte = lte.loc[ste]
            mkt = lt.isin(["tumor", "normal"])
            mke = lte.isin(["tumor", "normal"])
            Xtr_a = Xtr_a.loc[st[mkt]]
            Xte_a = Xte_a.loc[ste[mke]]
            yt = (lt.loc[st[mkt]] == "tumor").astype(int).to_numpy()
            yte = (lte.loc[ste[mke]] == "tumor").astype(int).to_numpy()
            if len(genes) >= 5 and len(np.unique(yt)) == 2 and len(np.unique(yte)) == 2:
                xfer_metrics = transfer_train_test_eval(
                    Xtr_a,
                    yt,
                    Xte_a,
                    yte,
                    seed=seed,
                    rf_n_estimators=int(ml_cfg.get("random_forest_n_estimators", 200)),
                    linear_alpha=linear_alpha,
                    calibration_method=str(xfer.get("calibration_method") or "").strip() or None,
                    calibration_cv=int(xfer.get("calibration_cv", 3) or 3),
                )
                xfer_metrics["train_dataset"] = train_ds
                xfer_metrics["test_dataset"] = test_ds
                (out_dir / f"{train_ds}_train_{test_ds}_test_transfer.json").write_text(
                    json.dumps(xfer_metrics, indent=2), encoding="utf-8"
                )
                pred_rows = []
                for i, sid in enumerate(Xte_a.index):
                    pred_rows.append(
                        {
                            "sample_id": sid,
                            "y_true": int(yte[i]),
                            "cohort": test_ds,
                        }
                    )
                pred_df = pd.DataFrame(pred_rows)
                write_tsv(pred_df, out_dir / f"{train_ds}_train_{test_ds}_test_sample_labels.tsv")
                print(
                    f"[ml] transfer {train_ds}→{test_ds} ({len(genes)} genes) -> "
                    f"{train_ds}_train_{test_ds}_test_transfer.json"
                )
                run_manifest.record_stage_meta(
                    settings,
                    "ml_transfer",
                    executed=True,
                    n_genes_aligned=len(genes),
                    n_genes_train_pre_align=n_tr_feat,
                    n_genes_test_pre_align=n_te_feat,
                    train_dataset=train_ds,
                    test_dataset=test_ds,
                )
            else:
                print(
                    f"[ml] transfer omitido: genes comunes={len(genes)}, "
                    f"clases train={np.unique(yt)}, test={np.unique(yte)}."
                )
                run_manifest.record_stage_meta(
                    settings,
                    "ml_transfer",
                    executed=False,
                    reason="insufficient_genes_or_labels",
                    n_genes_aligned=len(genes),
                    n_genes_train_pre_align=n_tr_feat,
                    n_genes_test_pre_align=n_te_feat,
                    train_dataset=train_ds,
                    test_dataset=test_ds,
                )
        else:
            print(f"[ml] transfer: faltan matrices o metadata para {train_ds} / {test_ds}.")
            run_manifest.record_stage_meta(
                settings,
                "ml_transfer",
                executed=False,
                reason="missing_matrix_or_metadata",
                train_dataset=train_ds,
                test_dataset=test_ds,
            )

    boot_n = int(ml_cfg.get("feature_stability_bootstrap_n", 0) or 0)
    if boot_n >= 10:
        stab = bootstrap_linear_feature_stability(
            X, y, n_boot=boot_n, seed=seed, alpha=linear_alpha, top_k=int(ml_cfg.get("bootstrap_top_k", 30))
        )
        write_tsv(stab, out_dir / f"{primary}_feature_stability_bootstrap.tsv")
        print(f"[ml] bootstrap estabilidad L1 (n={boot_n}) -> {primary}_feature_stability_bootstrap.tsv")

    (out_dir / f"{primary}_model_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    metrics_flat: list[dict] = []
    for k, v in metrics.items():
        if k == "nested_cv" and isinstance(v, dict):
            for kk, vv in v.items():
                if isinstance(vv, dict) and "roc_auc_mean" in vv:
                    metrics_flat.append({"model": kk, **vv})
            continue
        if not isinstance(v, dict):
            metrics_flat.append({"model": k, "value": v})
            continue
        if "roc_auc_mean" in v or "accuracy_test" in v or "roc_auc_test" in v:
            metrics_flat.append({"model": k, **v})
    write_tsv(pd.DataFrame(metrics_flat), out_dir / f"{primary}_model_metrics.tsv")
    extra = []
    if "shap_mean_abs" in ranking.columns:
        extra.append("SHAP")
    if any("nested" in str(m.get("model", "")) for m in metrics_flat):
        extra.append("nested-CV")
    if boot_n >= 10:
        extra.append("bootstrap-L1")
    suf = f" ({', '.join(extra)})" if extra else ""
    print(f"[ml] métricas y ranking{suf} -> results/ml/{primary}_*")
