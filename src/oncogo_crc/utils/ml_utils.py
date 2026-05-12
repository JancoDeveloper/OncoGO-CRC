from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


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

    pipe_sparse = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "model",
                SGDClassifier(
                    loss="log_loss",
                    penalty="l1",
                    alpha=0.02,
                    max_iter=4000,
                    random_state=seed,
                    class_weight="balanced",
                ),
            ),
        ]
    )
    pipe_sparse.fit(X_train, y_train)
    pred_l = pipe_sparse.predict(X_test)
    proba_l = pipe_sparse.predict_proba(X_test)[:, 1]
    metrics["sparse_linear_sgd"] = {
        "accuracy_test": float(accuracy_score(y_test, pred_l)),
        "roc_auc_test": float(roc_auc_score(y_test, proba_l)),
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
        out["combined_score"] = (
            out["rf_importance"].rank(pct=True) * 0.40
            + out["l1_abs_coef"].rank(pct=True) * 0.30
            + out["shap_mean_abs"].rank(pct=True) * 0.30
        )
    else:
        out["combined_score"] = out["rf_importance"].rank(pct=True) * 0.55 + out["l1_abs_coef"].rank(pct=True) * 0.45
    return metrics, out.sort_values("combined_score", ascending=False)
