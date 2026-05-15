from __future__ import annotations

import json
from typing import Any

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import (
    GridSearchCV,
    StratifiedKFold,
    cross_val_score,
    train_test_split,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def normalized_matrix_to_samples_x_genes(
    df: pd.DataFrame,
    *,
    uniprot_to_gene: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Matriz normalizada (filas = features, columnas = muestras) → muestras × genes.

    Si ``uniprot_to_gene`` no está vacío, intenta resolver filas sin ``gene_symbol``
    útil (p. ej. PRIDE con UniProt en ``protein_id``) a símbolo HGNC para alinear
    con CPTAC. Varias filas que colapsan al mismo gen se agregan con mediana.
    """
    from oncogo_crc.utils.reproducibility_overlap import uniprot_base

    meta = {"protein_id", "gene_symbol", "uniprot_accession"}
    sample_cols = [c for c in df.columns if c not in meta]
    utg = uniprot_to_gene or {}

    if utg:

        def _row_gene(r: pd.Series) -> str:
            if "gene_symbol" in df.columns:
                g = r.get("gene_symbol")
                if g is not None and pd.notna(g):
                    s = str(g).strip().upper()
                    if s and s != "NAN":
                        return s
            pid_raw = str(r.get("protein_id", ""))
            pid = pid_raw.strip().upper()
            if pid in utg:
                return utg[pid]
            ub = uniprot_base(pid_raw)
            if ub and ub in utg:
                return utg[ub]
            if "uniprot_accession" in df.columns:
                ua = r.get("uniprot_accession")
                if ua is not None and pd.notna(ua):
                    u2 = uniprot_base(str(ua))
                    if u2 and u2 in utg:
                        return utg[u2]
            return pid

        genes = df.apply(_row_gene, axis=1)
    elif "gene_symbol" in df.columns:
        genes = (
            df["gene_symbol"]
            .fillna(df["protein_id"])
            .astype(str)
            .str.strip()
            .str.upper()
        )
    else:
        genes = df["protein_id"].astype(str).str.strip().str.upper()

    agg = (
        df.assign(_gene=genes)
        .groupby("_gene", sort=False)[sample_cols]
        .median()
        .reset_index()
    )
    X = agg.set_index("_gene").astype(float).T
    X.index.name = "sample_id"
    return X


def recommend_nested_cv_splits(
    y: np.ndarray, *, outer_requested: int, inner_requested: int
) -> tuple[int, int, str]:
    """Ajusta pliegues a límites de StratifiedKFold (n_splits ≤ mínimo de conteos por clase)."""
    y = np.asarray(y).astype(int).ravel()
    _, counts = np.unique(y, return_counts=True)
    mmin = int(counts.min()) if len(counts) else 0
    if mmin < 2:
        return 2, 2, "clase minoritaria < 2; nested CV no fiable"
    outer = min(max(2, int(outer_requested)), mmin)
    inner = min(max(2, int(inner_requested)), mmin)
    inner = min(inner, max(2, mmin - 1))
    notes: list[str] = []
    if outer != int(outer_requested):
        notes.append(f"outer_splits {outer_requested}→{outer} (mín. por clase={mmin})")
    if inner != int(inner_requested):
        notes.append(f"inner_splits {inner_requested}→{inner}")
    return outer, inner, "; ".join(notes) if notes else "sin recorte"


def nested_cv_model_flags(nest_cfg: dict[str, Any]) -> tuple[bool, bool, bool]:
    """Devuelve (run_l1, run_elastic_net, run_rf) para nested CV."""
    raw = nest_cfg.get("models")
    if isinstance(raw, list) and len(raw) > 0:
        sl = {str(m).strip().lower() for m in raw}
        run_l1 = bool(
            sl
            & {
                "l1",
                "lasso",
                "sparse_linear",
                "sparse_linear_l1",
                "sparse_linear_l1_nested",
                "sgd_l1",
            }
        )
        run_en = bool(
            sl
            & {
                "elastic",
                "elastic_net",
                "elasticnet",
                "elastic_net_sgd",
                "elastic_net_sgd_nested",
                "en",
            }
        )
        run_rf = bool(
            sl & {"rf", "random_forest", "randomforest", "random_forest_nested"}
        )
        if not (run_l1 or run_en or run_rf):
            return True, True, True
        return run_l1, run_en, run_rf
    if bool(nest_cfg.get("only_rf")):
        return False, False, True
    return True, True, True


def align_two_cohort_matrices(
    X_train: pd.DataFrame, X_test: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    X_train = X_train.copy()
    X_test = X_test.copy()
    X_train.columns = [str(c).strip().upper() for c in X_train.columns]
    X_test.columns = [str(c).strip().upper() for c in X_test.columns]
    common = sorted(set(X_train.columns) & set(X_test.columns))
    if not common:
        return (
            X_train.iloc[:, :0].copy(),
            X_test.iloc[:, :0].copy(),
            [],
        )
    return X_train[common].copy(), X_test[common].copy(), common


def _sparse_linear_pipe(seed: int, alpha: float) -> Pipeline:
    return Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "model",
                SGDClassifier(
                    loss="log_loss",
                    penalty="l1",
                    alpha=float(alpha),
                    max_iter=5000,
                    random_state=seed,
                    class_weight="balanced",
                ),
            ),
        ]
    )


def _elastic_net_pipe(seed: int, alpha: float, l1_ratio: float) -> Pipeline:
    return Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "model",
                SGDClassifier(
                    loss="log_loss",
                    penalty="elasticnet",
                    alpha=float(alpha),
                    l1_ratio=float(l1_ratio),
                    max_iter=5000,
                    random_state=seed,
                    class_weight="balanced",
                ),
            ),
        ]
    )


def nested_cross_validate_models(
    X: np.ndarray,
    y: np.ndarray,
    *,
    outer_splits: int,
    inner_splits: int,
    seed: int,
    rf_n_estimators: int,
    linear_alpha_grid: list[float],
    rf_max_depth_grid: list[int | None],
    rf_min_samples_leaf_grid: list[int],
    run_l1: bool = True,
    run_elastic_net: bool = True,
    run_rf: bool = True,
    split_note: str = "",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Nested stratified CV: inner GridSearchCV, outer hold-out AUC por pliegue."""
    outer_splits = max(2, int(outer_splits))
    inner_splits = max(2, int(inner_splits))
    outer = StratifiedKFold(n_splits=outer_splits, shuffle=True, random_state=seed)
    inner = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=seed)
    rows: list[dict[str, Any]] = []
    alphas = [float(a) for a in linear_alpha_grid] or [0.02]
    alpha0 = float(alphas[len(alphas) // 2])

    for fold_id, (tr_o, te_o) in enumerate(outer.split(X, y), start=1):
        X_tr, X_te = X[tr_o], X[te_o]
        y_tr, y_te = y[tr_o], y[te_o]

        if run_l1:
            base_l1 = _sparse_linear_pipe(seed, alpha0)
            grid_l1 = GridSearchCV(
                clone(base_l1),
                {"model__alpha": list(alphas)},
                cv=inner,
                scoring="roc_auc",
                n_jobs=-1,
                refit=True,
            )
            grid_l1.fit(X_tr, y_tr)
            p1 = grid_l1.predict_proba(X_te)[:, 1]
            rows.append(
                {
                    "outer_fold": fold_id,
                    "model": "sparse_linear_l1_nested",
                    "roc_auc": float(roc_auc_score(y_te, p1)),
                    "best_params_json": json.dumps(grid_l1.best_params_),
                }
            )

        if run_elastic_net:
            base_en = _elastic_net_pipe(seed, alpha0, 0.15)
            grid_en = GridSearchCV(
                clone(base_en),
                {"model__alpha": list(alphas)},
                cv=inner,
                scoring="roc_auc",
                n_jobs=-1,
                refit=True,
            )
            grid_en.fit(X_tr, y_tr)
            p_en = grid_en.predict_proba(X_te)[:, 1]
            rows.append(
                {
                    "outer_fold": fold_id,
                    "model": "elastic_net_sgd_nested",
                    "roc_auc": float(roc_auc_score(y_te, p_en)),
                    "best_params_json": json.dumps(grid_en.best_params_),
                }
            )

        if run_rf:
            base_rf = RandomForestClassifier(
                n_estimators=int(rf_n_estimators),
                random_state=seed,
                class_weight="balanced_subsample",
                n_jobs=-1,
            )
            depth_list: list[Any] = list(rf_max_depth_grid)
            grid_rf = GridSearchCV(
                clone(base_rf),
                {
                    "max_depth": depth_list,
                    "min_samples_leaf": [int(x) for x in rf_min_samples_leaf_grid],
                },
                cv=inner,
                scoring="roc_auc",
                n_jobs=-1,
                refit=True,
            )
            grid_rf.fit(X_tr, y_tr)
            prf = grid_rf.predict_proba(X_te)[:, 1]
            rows.append(
                {
                    "outer_fold": fold_id,
                    "model": "random_forest_nested",
                    "roc_auc": float(roc_auc_score(y_te, prf)),
                    "best_params_json": json.dumps(grid_rf.best_params_),
                }
            )

    per_fold = pd.DataFrame(rows)
    summary: dict[str, Any] = {}
    for mname in per_fold["model"].unique():
        sub = per_fold.loc[per_fold["model"] == mname, "roc_auc"].astype(float)
        summary[mname] = {
            "roc_auc_mean": float(sub.mean()),
            "roc_auc_std": float(sub.std(ddof=1)) if len(sub) > 1 else 0.0,
            "n_outer_folds": int(len(sub)),
        }
    summary["n_features"] = int(X.shape[1])
    summary["feature_names_note"] = "Columnas = genes (post-armonización gen-centric)."
    if split_note:
        summary["split_adjustment"] = split_note
    return per_fold, summary


def _effective_calibration_cv(y: np.ndarray, requested: int) -> int:
    y = np.asarray(y).astype(int).ravel()
    counts = np.bincount(y)
    mmin = int(counts.min()) if len(counts) else 1
    n = len(y)
    return max(2, min(int(requested), mmin, max(2, n // 2)))


def transfer_train_test_eval(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    X_test: pd.DataFrame,
    y_test: np.ndarray,
    *,
    seed: int,
    rf_n_estimators: int,
    linear_alpha: float = 0.02,
    calibration_method: str | None = None,
    calibration_cv: int = 3,
) -> dict[str, Any]:
    """Entrena en cohorte descubrimiento (p. ej. CPTAC), evalúa en cohorte externa (p. ej. PRIDE).

    Escala con StandardScaler en los pipelines lineales; RF usa scaler+RF. Opcionalmente aplica
    **calibración probabilística** (`sigmoid` / Platt o `isotonic`) con CV solo sobre CPTAC
    antes de proyectar probabilidades en PRIDE.
    """
    genes = list(X_train.columns)
    if not genes or X_test.shape[1] != len(genes):
        return {
            "status": "error",
            "detail": "Sin genes comunes alineados o dimensiones inconsistentes.",
        }
    Xt = X_train.values.astype(float)
    Xv = X_test.values.astype(float)
    pipe_l1 = _sparse_linear_pipe(seed, linear_alpha)
    pipe_l1.fit(Xt, y_train)
    pl1 = pipe_l1.predict_proba(Xv)[:, 1]

    en = _elastic_net_pipe(seed, linear_alpha, 0.15)
    en.fit(Xt, y_train)
    pen = en.predict_proba(Xv)[:, 1]

    rf = RandomForestClassifier(
        n_estimators=int(rf_n_estimators),
        random_state=seed,
        class_weight="balanced_subsample",
        n_jobs=-1,
    )
    scaler = StandardScaler()
    Xt_s = scaler.fit_transform(Xt)
    Xv_s = scaler.transform(Xv)
    rf.fit(Xt_s, y_train)
    prf = rf.predict_proba(Xv_s)[:, 1]

    out: dict[str, Any] = {
        "status": "ok",
        "n_genes_aligned": len(genes),
        "n_train_samples": int(len(y_train)),
        "n_test_samples": int(len(y_test)),
        "sparse_linear_l1": {
            "roc_auc": float(roc_auc_score(y_test, pl1)),
            "accuracy": float(accuracy_score(y_test, (pl1 >= 0.5).astype(int))),
        },
        "elastic_net_sgd": {
            "roc_auc": float(roc_auc_score(y_test, pen)),
            "accuracy": float(accuracy_score(y_test, (pen >= 0.5).astype(int))),
        },
        "random_forest": {
            "roc_auc": float(roc_auc_score(y_test, prf)),
            "accuracy": float(accuracy_score(y_test, (prf >= 0.5).astype(int))),
        },
        "note": (
            "El escalado L1/EN usa solo CPTAC; RF usa StandardScaler en CPTAC y la misma transformación en PRIDE. "
            "Distribuciones de cohorte distintas pueden reducir AUC aun con genes alineados (transferencia honesta)."
        ),
    }

    cal_raw = (calibration_method or "").strip().lower()
    if cal_raw in ("none", "", "false", "off", "0"):
        out["calibration"] = {"applied": False}
        return out

    cal_m = cal_raw if cal_raw in ("sigmoid", "isotonic") else "sigmoid"
    cv_eff = _effective_calibration_cv(y_train, calibration_cv)
    if len(y_train) < 2 * cv_eff or cv_eff < 2:
        out["calibration"] = {
            "applied": False,
            "reason": "muestras CPTAC insuficientes para CalibratedClassifierCV con el cv solicitado",
        }
        return out

    cal_l1 = CalibratedClassifierCV(
        clone(_sparse_linear_pipe(seed, linear_alpha)), method=cal_m, cv=cv_eff
    )
    cal_l1.fit(Xt, y_train)
    pl1c = cal_l1.predict_proba(Xv)[:, 1]
    out["sparse_linear_l1_calibrated"] = {
        "method": cal_m,
        "calibration_cv": cv_eff,
        "roc_auc": float(roc_auc_score(y_test, pl1c)),
        "accuracy": float(accuracy_score(y_test, (pl1c >= 0.5).astype(int))),
    }

    cal_en = CalibratedClassifierCV(
        clone(_elastic_net_pipe(seed, linear_alpha, 0.15)), method=cal_m, cv=cv_eff
    )
    cal_en.fit(Xt, y_train)
    penc = cal_en.predict_proba(Xv)[:, 1]
    out["elastic_net_sgd_calibrated"] = {
        "method": cal_m,
        "calibration_cv": cv_eff,
        "roc_auc": float(roc_auc_score(y_test, penc)),
        "accuracy": float(accuracy_score(y_test, (penc >= 0.5).astype(int))),
    }

    pipe_rf = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "model",
                RandomForestClassifier(
                    n_estimators=int(rf_n_estimators),
                    random_state=seed,
                    class_weight="balanced_subsample",
                    n_jobs=-1,
                ),
            ),
        ]
    )
    cal_rf = CalibratedClassifierCV(clone(pipe_rf), method=cal_m, cv=cv_eff)
    cal_rf.fit(Xt, y_train)
    prfc = cal_rf.predict_proba(Xv)[:, 1]
    out["random_forest_calibrated"] = {
        "method": cal_m,
        "calibration_cv": cv_eff,
        "roc_auc": float(roc_auc_score(y_test, prfc)),
        "accuracy": float(accuracy_score(y_test, (prfc >= 0.5).astype(int))),
    }
    out["calibration"] = {"applied": True, "method": cal_m, "calibration_cv": cv_eff}
    return out


def bootstrap_linear_feature_stability(
    X: pd.DataFrame,
    y: np.ndarray,
    *,
    n_boot: int,
    seed: int,
    alpha: float = 0.02,
    top_k: int = 30,
) -> pd.DataFrame:
    """Fracción de réplicas bootstrap en las que |coef| del modelo L1 está entre los top_k."""
    rng = np.random.default_rng(seed)
    names = list(X.columns)
    n = len(y)
    counts = np.zeros(len(names), dtype=float)
    for _ in range(int(n_boot)):
        idx = rng.choice(n, size=n, replace=True)
        if len(np.unique(y[idx])) < 2:
            continue
        pipe = _sparse_linear_pipe(seed, alpha)
        pipe.fit(X.values[idx], y[idx])
        coef = np.abs(pipe.named_steps["model"].coef_.ravel())
        if len(coef) != len(names):
            continue
        thr = np.partition(coef, -min(top_k, len(coef)))[-min(top_k, len(coef))]
        counts += (coef >= thr).astype(float)
    freq = counts / max(1.0, float(n_boot))
    return pd.DataFrame({"feature_id": names, "bootstrap_topk_frequency": freq}).sort_values(
        "bootstrap_topk_frequency", ascending=False
    )


def train_prioritization_models(
    X: pd.DataFrame,
    y: np.ndarray,
    *,
    test_size: float = 0.25,
    seed: int = 42,
    rf_n_estimators: int = 200,
    use_shap: bool = False,
    shap_background_max: int = 256,
    cv_folds: int = 0,
    include_elastic_net: bool = True,
    linear_alpha: float = 0.02,
) -> tuple[dict[str, dict], pd.DataFrame]:
    """
    X: muestras × proteínas (filas = muestras).
    y: etiquetas binarias (p. ej. tumor=1, normal=0).
    """
    strat = y if len(np.unique(y)) > 1 else None
    X_train, X_test, y_train, y_test = train_test_split(
        X.values, y, test_size=test_size, random_state=seed, stratify=strat
    )
    feature_names = list(X.columns)
    metrics: dict[str, dict] = {}

    pipe_sparse = _sparse_linear_pipe(seed, linear_alpha)
    pipe_sparse.fit(X_train, y_train)
    pred_l = pipe_sparse.predict(X_test)
    proba_l = pipe_sparse.predict_proba(X_test)[:, 1]
    metrics["sparse_linear_sgd"] = {
        "accuracy_test": float(accuracy_score(y_test, pred_l)),
        "roc_auc_test": float(roc_auc_score(y_test, proba_l)),
    }

    if include_elastic_net:
        pipe_en = _elastic_net_pipe(seed, linear_alpha, 0.15)
        pipe_en.fit(X_train, y_train)
        pred_en = pipe_en.predict(X_test)
        proba_en = pipe_en.predict_proba(X_test)[:, 1]
        metrics["elastic_net_sgd"] = {
            "accuracy_test": float(accuracy_score(y_test, pred_en)),
            "roc_auc_test": float(roc_auc_score(y_test, proba_en)),
        }

    clf = RandomForestClassifier(
        n_estimators=rf_n_estimators, random_state=seed, class_weight="balanced_subsample"
    )
    if cv_folds and cv_folds >= 2:
        skf = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=seed)
        aucs = cross_val_score(
            RandomForestClassifier(
                n_estimators=rf_n_estimators,
                random_state=seed,
                class_weight="balanced_subsample",
            ),
            X.values,
            y,
            cv=skf,
            scoring="roc_auc",
        )
        metrics["random_forest_cv"] = {
            "roc_auc_mean": float(np.mean(aucs)),
            "roc_auc_std": float(np.std(aucs)),
            "n_folds": int(cv_folds),
        }
    clf.fit(X_train, y_train)
    pred_rf = clf.predict(X_test)
    proba = clf.predict_proba(X_test)[:, 1]
    metrics["random_forest"] = {
        "accuracy_test": float(accuracy_score(y_test, pred_rf)),
        "roc_auc_test": float(roc_auc_score(y_test, proba)),
    }

    coefs = pipe_sparse.named_steps["model"].coef_.ravel()
    imp = pd.Series(clf.feature_importances_, index=feature_names, name="rf_importance")
    coef = pd.Series(np.abs(coefs), index=feature_names, name="l1_abs_coef")
    out = pd.concat([imp, coef], axis=1).fillna(0.0)
    if include_elastic_net:
        en_coef = np.abs(pipe_en.named_steps["model"].coef_.ravel())
        out["elasticnet_abs_coef"] = en_coef
    if use_shap:
        try:
            import shap  # type: ignore[import-untyped]
        except ImportError as e:
            raise ImportError(
                "ml.use_shap=true requiere el paquete 'shap'. Instala: pip install -e \".[explain]\""
            ) from e
        rng = np.random.default_rng(seed)
        n_tr = len(X_train)
        k = min(int(shap_background_max), n_tr)
        idx = rng.choice(n_tr, size=k, replace=False) if n_tr > k else np.arange(n_tr)
        bg = X_train[idx]
        explainer = shap.TreeExplainer(clf)
        sv = explainer.shap_values(bg)
        if isinstance(sv, list):
            shv = sv[1]
        else:
            shv = sv
        if shv.ndim == 3:
            shv = shv[:, :, 1]
        mean_abs = np.abs(shv).mean(axis=0)
        out["shap_mean_abs"] = mean_abs
    if "shap_mean_abs" in out.columns:
        if include_elastic_net and "elasticnet_abs_coef" in out.columns:
            out["combined_score"] = (
                out["rf_importance"].rank(pct=True) * 0.30
                + out["l1_abs_coef"].rank(pct=True) * 0.15
                + out["elasticnet_abs_coef"].rank(pct=True) * 0.15
                + out["shap_mean_abs"].rank(pct=True) * 0.40
            )
        else:
            out["combined_score"] = (
                out["rf_importance"].rank(pct=True) * 0.40
                + out["l1_abs_coef"].rank(pct=True) * 0.30
                + out["shap_mean_abs"].rank(pct=True) * 0.30
            )
    elif include_elastic_net and "elasticnet_abs_coef" in out.columns:
        out["combined_score"] = (
            out["rf_importance"].rank(pct=True) * 0.45
            + out["l1_abs_coef"].rank(pct=True) * 0.30
            + out["elasticnet_abs_coef"].rank(pct=True) * 0.25
        )
    else:
        out["combined_score"] = (
            out["rf_importance"].rank(pct=True) * 0.55 + out["l1_abs_coef"].rank(pct=True) * 0.45
        )
    return metrics, out.sort_values("combined_score", ascending=False)
