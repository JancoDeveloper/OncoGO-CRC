"""Ejes que el certificado primario no cubre, con los mismos datos públicos.

Cuatro preguntas, todas con la dirección esperada declarada antes de mirar:

1. Secuencia adenoma-carcinoma. PXD034575 trae pólipo, normal y tumor del mismo
   paciente. El pólipo nunca se había contrastado.
2. Brecha traslacional tejido→plasma. La literatura reciente muestra que la
   sobreexpresión en tejido no implica detectabilidad en circulación.
3. Panel externo predeclarado en plasma, con la dirección publicada.
4. Selección de modelo por AUC externo con intervalo, no por el AUC interno.

No se modifica ninguna matriz y no se fusionan cohortes.
"""

from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from oncogo_crc.certification.contract import TREAT_LFC
from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.moderated import moderate_table, paired_observed_stats
from oncogo_crc.utils.quant_columns import is_summary_protein
from oncogo_crc.utils.reproducibility_overlap import (
    load_gene_uniprot_from_mapping,
    uniprot_base,
)
from oncogo_crc.utils.ids import effective_protein_mapping
from oncogo_crc.stages.certified_inference import MIN_COMPLETE, SEED, _numeric_block

N_BOOT = 2000

# Dirección publicada en plasma o suero de CCR frente a control sano.
# Solo entran proteínas cuya dirección el artículo declara de forma explícita.
PLASMA_PANEL: dict[str, tuple[int, str]] = {
    "C5": (+1, "C5 elevado en plasma de CCR, confirmado por ELISA en cohorte de validación"),
    "C8A": (+1, "cascada del complemento aumentada en plasma de CCR"),
    "C9": (+1, "C9 aumentado en plasma de CCR"),
    "LRG1": (+1, "LRG aumentado en plasma de CCR"),
    "SERPINA1": (+1, "A1AT aumentado en plasma de CCR"),
    "ORM1": (+1, "AGP1 aumentado en plasma de CCR"),
    "FN1": (-1, "fibronectina disminuida en plasma de CCR"),
    "CEACAM5": (+1, "CEA/CEACAM5 aumentado en plasma de CCR, AUC 0,68"),
    "FTL": (-1, "ferritina disminuida en plasma de CCR, AUC 0,66"),
    "FTH1": (-1, "ferritina disminuida en plasma de CCR, AUC 0,66"),
    "SPP1": (+1, "osteopontina aumentada en los tres tumores del panel"),
    "PECAM1": (+1, "PECAM1/CD31 aumentado en plasma de CCR, AUC 0,69"),
}

_PRIDE_COL = re.compile(r"^LFQ intensity (CRC\d+)([MNP])$", re.IGNORECASE)


def _exact(same: int, n: int) -> dict:
    if n == 0:
        return {"n": 0, "same": 0, "concordance": None, "ci_low": None, "ci_high": None, "p_vs_half": None}
    test = binomtest(same, n, p=0.5, alternative="two-sided")
    ci = test.proportion_ci(confidence_level=0.95, method="exact")
    return {
        "n": n,
        "same": same,
        "concordance": same / n,
        "ci_low": float(ci.low),
        "ci_high": float(ci.high),
        "p_vs_half": float(test.pvalue),
    }


def _pride_log2(path) -> tuple[pd.DataFrame, dict[str, dict[str, str]]]:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    raw = raw.loc[~raw["protein_id"].map(is_summary_protein)].drop_duplicates("protein_id", keep="first")
    mat = raw.set_index("protein_id")
    patients: dict[str, dict[str, str]] = {}
    for col in mat.columns:
        m = _PRIDE_COL.match(str(col))
        if not m:
            continue
        patients.setdefault(m.group(1), {})[m.group(2).upper()] = str(col)
    values = _numeric_block(mat, list(mat.columns))
    with np.errstate(divide="ignore", invalid="ignore"):
        logged = np.where(values > 0, np.log2(values), np.nan)
    return pd.DataFrame(logged, index=mat.index, columns=mat.columns), patients


def _paired_arm(
    log2mat: pd.DataFrame,
    patients: dict[str, dict[str, str]],
    *,
    target: str,
    label: str,
    lfc: float,
    alpha: float,
) -> tuple[pd.DataFrame, int]:
    """Contraste pareado dentro del paciente: target frente a su propio tejido normal."""
    pairs = [(v[target], v["N"]) for v in patients.values() if target in v and "N" in v]
    diffs = log2mat[[a for a, _ in pairs]].to_numpy() - log2mat[[b for _, b in pairs]].to_numpy()
    coef, s2, df, stdev, n_obs = paired_observed_stats(diffs, min_pairs=MIN_COMPLETE)
    table, _d0, _s02 = moderate_table(
        protein_ids=log2mat.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method=f"moderated_t_treat_paired_lfq_{label}",
    )
    table["comparison_id"] = f"pride_{label}_vs_normal"
    return table, len(pairs)


def _adenoma_sequence(tumor: pd.DataFrame, polyp: pd.DataFrame, *, alpha: float) -> dict:
    """De las proteínas con efecto certificado en tumor, ¿el pólipo ya va en esa dirección?"""
    t = tumor.loc[tumor["significant"], ["protein_id", "log2fc"]].rename(columns={"log2fc": "lfc_tumor"})
    p = polyp[["protein_id", "log2fc"]].rename(columns={"log2fc": "lfc_polyp"})
    m = t.merge(p, on="protein_id", how="inner").dropna(subset=["lfc_tumor", "lfc_polyp"])
    m = m.loc[(m["lfc_tumor"] != 0) & (m["lfc_polyp"] != 0)]
    same = int((np.sign(m["lfc_tumor"]) == np.sign(m["lfc_polyp"])).sum())
    out = {"direction": _exact(same, int(len(m)))}
    agree = m.loc[np.sign(m["lfc_tumor"]) == np.sign(m["lfc_polyp"])]
    inter = int((agree["lfc_polyp"].abs() < agree["lfc_tumor"].abs()).sum())
    out["intermediate_given_same_direction"] = _exact(inter, int(len(agree)))
    out["median_abs_lfc_tumor"] = float(m["lfc_tumor"].abs().median()) if len(m) else None
    out["median_abs_lfc_polyp"] = float(m["lfc_polyp"].abs().median()) if len(m) else None
    return out


def _tissue_to_plasma(
    tissue: pd.DataFrame,
    plasma: pd.DataFrame,
    gene_to_uni: dict[str, str],
    *,
    alpha: float,
) -> tuple[dict, pd.DataFrame]:
    """¿La dirección en tejido se repite en suero? La expectativa declarada es que no."""
    t = tissue.loc[tissue["significant"], ["protein_id", "log2fc", "padj"]].copy()
    t["uniprot"] = t["protein_id"].astype(str).str.upper().map(gene_to_uni)
    t = t.dropna(subset=["uniprot"])
    p = plasma[["protein_id", "log2fc", "padj", "significant"]].copy()
    p["uniprot"] = p["protein_id"].map(lambda v: uniprot_base(str(v)))
    p = p.dropna(subset=["uniprot"]).drop_duplicates("uniprot")
    m = t.merge(
        p.rename(columns={"log2fc": "lfc_plasma", "padj": "padj_plasma", "significant": "sig_plasma"}),
        on="uniprot",
        how="inner",
        suffixes=("_tissue", "_plasma"),
    )
    m = m.loc[m["lfc_plasma"].notna() & (m["log2fc"] != 0) & (m["lfc_plasma"] != 0)]
    same = int((np.sign(m["log2fc"]) == np.sign(m["lfc_plasma"])).sum())
    stats = _exact(same, int(len(m)))
    stats["n_tissue_significant_with_uniprot"] = int(len(t))
    stats["n_measured_in_plasma"] = int(len(m))
    stats["fraction_tissue_hits_seen_in_plasma"] = (len(m) / len(t)) if len(t) else None
    stats["n_also_significant_in_plasma"] = int(m["sig_plasma"].fillna(False).sum())
    return stats, m


def _external_panel(
    plasma: pd.DataFrame,
    tissue: pd.DataFrame,
    gene_to_uni: dict[str, str],
) -> tuple[dict, pd.DataFrame]:
    p = plasma.copy()
    p["uniprot"] = p["protein_id"].map(lambda v: uniprot_base(str(v)))
    p = p.dropna(subset=["uniprot"]).drop_duplicates("uniprot").set_index("uniprot")
    ti = tissue.drop_duplicates("protein_id").set_index("protein_id")
    rows = []
    for gene, (expected, source) in PLASMA_PANEL.items():
        uni = gene_to_uni.get(gene)
        row = {
            "gene_symbol": gene,
            "uniprot": uni or "",
            "expected_direction_plasma": expected,
            "published_claim": source,
        }
        if uni and uni in p.index:
            r = p.loc[uni]
            row["plasma_log2fc"] = float(r["log2fc"]) if pd.notna(r["log2fc"]) else None
            row["plasma_padj"] = float(r["padj"]) if pd.notna(r["padj"]) else None
            row["plasma_significant"] = bool(r["significant"])
            row["plasma_measured"] = True
            row["direction_matches"] = (
                None
                if row["plasma_log2fc"] in (None, 0.0)
                else bool(np.sign(row["plasma_log2fc"]) == expected)
            )
        else:
            row["plasma_measured"] = False
            row["direction_matches"] = None
        if gene in ti.index:
            row["tissue_log2fc"] = float(ti.loc[gene, "log2fc"]) if pd.notna(ti.loc[gene, "log2fc"]) else None
            row["tissue_significant"] = bool(ti.loc[gene, "significant"])
        rows.append(row)
    frame = pd.DataFrame(rows)
    ev = frame.loc[frame["direction_matches"].notna()]
    stats = _exact(int(ev["direction_matches"].sum()), int(len(ev)))
    stats["n_panel"] = int(len(frame))
    stats["n_measured_in_plasma"] = int(frame["plasma_measured"].sum())
    stats["n_significant_in_plasma"] = int(frame.get("plasma_significant", pd.Series(dtype=bool)).fillna(False).sum())
    return stats, frame


def _zscore_rows(frame: pd.DataFrame) -> pd.DataFrame:
    mu = frame.mean(axis=1)
    sd = frame.std(axis=1, ddof=1).replace(0.0, np.nan)
    return frame.sub(mu, axis=0).div(sd, axis=0)


def _model_selection(
    settings: Settings,
    cptac_raw: pd.DataFrame,
    cptac_meta: pd.DataFrame,
    pride_log2: pd.DataFrame,
    patients: dict[str, dict[str, str]],
    gene_to_uni: dict[str, str],
) -> dict:
    """Regla declarada antes de mirar: la métrica primaria es el AUC en la cohorte externa.

    Las dos matrices están en escalas distintas (log-ratio frente a log2 de
    intensidad), así que cada cohorte se estandariza por proteína dentro de sí
    misma. Eso convierte la transferencia en una prueba de patrón relativo, no
    de calibración absoluta.
    """
    labels = cptac_meta.set_index("sample_id")["condition"].astype(str)
    cols = [
        c
        for c in cptac_raw.columns
        if str(c).endswith(" Unshared Log Ratio") and labels.get(c) in ("tumor", "normal")
    ]
    tissue = cptac_raw.loc[~cptac_raw["protein_id"].map(is_summary_protein)]
    tissue = tissue.drop_duplicates("protein_id", keep="first").set_index("protein_id")
    tissue = pd.DataFrame(_numeric_block(tissue, cols), index=tissue.index, columns=cols)
    tissue["_uni"] = [gene_to_uni.get(str(g).upper()) for g in tissue.index]
    tissue = tissue.dropna(subset=["_uni"]).drop_duplicates("_uni").set_index("_uni")

    pride = pride_log2.copy()
    pride["_uni"] = [uniprot_base(str(i)) for i in pride.index]
    pride = pride.dropna(subset=["_uni"]).drop_duplicates("_uni").set_index("_uni")
    pcols = [v[k] for v in patients.values() for k in ("M", "N") if k in v]
    pride = pride[pcols]

    common = sorted(set(tissue.index) & set(pride.index))
    tr = _zscore_rows(tissue.loc[common]).fillna(0.0)
    te = _zscore_rows(pride.loc[common]).fillna(0.0)
    y_tr = np.array([1 if labels.get(c) == "tumor" else 0 for c in cols])
    y_te = np.array([1 if str(c).upper().endswith("M") else 0 for c in pcols])

    x_tr = tr.T.to_numpy()
    x_te = te.T.to_numpy()
    models = {
        "logistica_l1": LogisticRegression(
            solver="saga", l1_ratio=1.0, C=0.1, max_iter=5000, random_state=SEED
        ),
        "elastic_net": LogisticRegression(
            solver="saga", l1_ratio=0.5, C=0.1, max_iter=5000, random_state=SEED
        ),
        "random_forest": RandomForestClassifier(n_estimators=500, random_state=SEED, n_jobs=-1),
    }
    rng = np.random.default_rng(SEED)
    out: dict = {
        "n_features_common_uniprot": len(common),
        "n_train": int(len(y_tr)),
        "n_test": int(len(y_te)),
        "primary_metric": "AUC en PXD034575, cohorte no vista en el entrenamiento",
        "declared_rule": (
            "Se elige un modelo solo si su intervalo del 95 % por bootstrap no se solapa con el "
            "del siguiente. Si se solapan, no hay ganador y se reportan los tres."
        ),
        "standardization": "z por proteína dentro de cada cohorte; prueba de patrón, no de calibración",
        "models": {},
    }
    for name, model in models.items():
        model.fit(x_tr, y_tr)
        score = model.predict_proba(x_te)[:, 1]
        auc = float(roc_auc_score(y_te, score))
        boots = []
        for _ in range(N_BOOT):
            idx = rng.integers(0, len(y_te), len(y_te))
            if len(np.unique(y_te[idx])) < 2:
                continue
            boots.append(roc_auc_score(y_te[idx], score[idx]))
        lo, hi = np.percentile(boots, [2.5, 97.5])
        out["models"][name] = {
            "auc_external": auc,
            "ci_low": float(lo),
            "ci_high": float(hi),
            "n_bootstrap": len(boots),
        }
        print(f"[axes] modelo {name}: AUC externo={auc:.3f} IC95 {lo:.3f}-{hi:.3f}")
    ranked = sorted(out["models"].items(), key=lambda kv: kv[1]["auc_external"], reverse=True)
    best, second = ranked[0], ranked[1]
    disjoint = best[1]["ci_low"] > second[1]["ci_high"]
    out["verdict"] = (
        f"gana {best[0]}" if disjoint else "sin ganador: los intervalos se solapan"
    )
    out["ranking"] = [name for name, _ in ranked]
    return out


def run(settings: Settings) -> None:
    cfg = settings.config.get("certification", {}) or {}
    alpha = float(cfg.get("alpha", 0.05))
    lfc = float(cfg.get("treat_lfc", TREAT_LFC))
    out = settings.resolve("results") / "axes"
    out.mkdir(parents=True, exist_ok=True)
    cert = settings.resolve("results") / "certified"

    cptac = pd.read_csv(cert / "cptac_unshared_treat.tsv", sep="\t")
    mapping = effective_protein_mapping(settings)
    gene_to_uni = load_gene_uniprot_from_mapping(mapping, exclude_ambiguous=True)

    pride_ds = settings.datasets["datasets"]["pride_crc_meta"]
    log2mat, patients = _pride_log2(settings.resolve(pride_ds["matrix_path"]))
    n_full = sum(1 for v in patients.values() if {"M", "N", "P"} <= set(v))
    print(f"[axes] PXD034575: {len(patients)} pacientes, {n_full} con tumor, normal y pólipo.")

    tumor_arm, n_t = _paired_arm(log2mat, patients, target="M", label="tumor", lfc=lfc, alpha=alpha)
    polyp_arm, n_p = _paired_arm(log2mat, patients, target="P", label="polyp", lfc=lfc, alpha=alpha)
    write_tsv(tumor_arm, out / "pride_paired_tumor_vs_normal_treat.tsv")
    write_tsv(polyp_arm, out / "pride_paired_polyp_vs_normal_treat.tsv")
    print(
        f"[axes] pareado PRIDE: tumor n={n_t} sig={int(tumor_arm['significant'].sum())}; "
        f"pólipo n={n_p} sig={int(polyp_arm['significant'].sum())}"
    )
    sequence = _adenoma_sequence(tumor_arm, polyp_arm, alpha=alpha)
    print(f"[axes] adenoma-carcinoma: {sequence['direction']}")

    plasma = pd.read_csv(settings.resolve("results") / "secondary" / "plasma_crc_vs_healthy_treat.tsv", sep="\t")
    gap, gap_table = _tissue_to_plasma(cptac, plasma, gene_to_uni, alpha=alpha)
    write_tsv(gap_table, out / "tejido_a_plasma.tsv")
    print(f"[axes] brecha tejido→plasma: {gap}")

    panel_stats, panel_table = _external_panel(plasma, cptac, gene_to_uni)
    write_tsv(panel_table, out / "panel_externo_plasma.tsv")
    print(f"[axes] panel externo predeclarado: {panel_stats}")

    cptac_ds = settings.datasets["datasets"]["cptac_colorectal"]
    cptac_raw = pd.read_csv(settings.resolve(cptac_ds["matrix_path"]), sep="\t", low_memory=False)
    cptac_meta = pd.read_csv(settings.sample_metadata_path("cptac_colorectal"), sep="\t")
    models = _model_selection(settings, cptac_raw, cptac_meta, log2mat, patients, gene_to_uni)
    print(f"[axes] veredicto de modelo: {models['verdict']}")

    audit = {
        "alpha": alpha,
        "treat_lfc": lfc,
        "seed": SEED,
        "pride_patients": len(patients),
        "pride_patients_with_three_tissues": n_full,
        "pride_paired_tumor": {"n_pairs": n_t, "n_significant": int(tumor_arm["significant"].sum())},
        "pride_paired_polyp": {"n_pairs": n_p, "n_significant": int(polyp_arm["significant"].sum())},
        "adenoma_carcinoma_sequence": sequence,
        "tissue_to_plasma_gap": gap,
        "external_plasma_panel": panel_stats,
        "model_selection": models,
        "caveats": [
            "El brazo de pólipo es adenoma del mismo paciente, no una cohorte independiente de adenomas.",
            "El plasma es suero DIA de otro estudio (PXD055218). No hay el mismo paciente en tejido y en sangre, "
            "así que la brecha tejido→plasma es entre cohortes, no intraindividual.",
            "El panel externo se declaró a partir de artículos publicados antes de mirar esta matriz; "
            "la dirección esperada no se ajustó después.",
            "La transferencia de modelo estandariza cada cohorte por separado: mide patrón relativo, no calibración.",
            "32 muestras de prueba dan intervalos anchos. Un AUC alto aquí no selecciona un modelo.",
        ],
    }
    (out / "axes_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[axes] escrito en {out}")
