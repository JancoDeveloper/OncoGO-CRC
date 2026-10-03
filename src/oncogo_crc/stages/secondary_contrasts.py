"""Contrastes secundarios sobre las mismas matrices, sin reescribir el certificado.

Solo entra una variable si los dos grupos tienen al menos ``MIN_PER_GROUP``
observaciones. Las que no llegan se inventarían en el inventario y no se
contrastan. El tumor frente a normal primario no se vuelve a calcular.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from oncogo_crc.certification.contract import TREAT_LFC
from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.moderated import moderate_table, unpaired_observed_stats
from oncogo_crc.utils.quant_columns import is_summary_protein
from oncogo_crc.stages.certified_inference import _cptac_contrast, _numeric_block

MIN_PER_GROUP = 15


def _sign_agreement(reference: pd.DataFrame, other: pd.DataFrame, *, alpha: float) -> dict:
    a = reference[["protein_id", "log2fc", "padj"]].rename(columns={"log2fc": "lfc_ref", "padj": "padj_ref"})
    b = other[["protein_id", "log2fc"]].rename(columns={"log2fc": "lfc_other"})
    m = a.merge(b, on="protein_id", how="inner")
    sig = m.loc[m["padj_ref"].astype(float) < alpha]
    ok = sig["lfc_ref"].notna() & sig["lfc_other"].notna() & (sig["lfc_ref"] != 0) & (sig["lfc_other"] != 0)
    ev = sig.loc[ok]
    n = int(len(ev))
    if n == 0:
        return {"n": 0, "same": 0, "concordance": None, "ci_low": None, "ci_high": None}
    same = int((np.sign(ev["lfc_ref"]) == np.sign(ev["lfc_other"])).sum())
    test = binomtest(same, n, p=0.5, alternative="two-sided")
    ci = test.proportion_ci(confidence_level=0.95, method="exact")
    return {
        "n": n,
        "same": same,
        "concordance": same / n,
        "ci_low": float(ci.low),
        "ci_high": float(ci.high),
    }


def _stage_contrast(raw: pd.DataFrame, meta: pd.DataFrame, *, lfc: float, alpha: float) -> pd.DataFrame | None:
    sm = meta.loc[meta["sample_id"].astype(str).str.endswith(" Unshared Log Ratio")].copy()
    sm = sm.loc[sm["condition"].astype(str) == "tumor"]
    sm = sm.loc[sm["stage_bin"].astype(str).isin(["advanced", "early"])]
    sm = sm.drop_duplicates("participant_id", keep="first")
    advanced = sm.loc[sm["stage_bin"].astype(str) == "advanced", "sample_id"].astype(str).tolist()
    early = sm.loc[sm["stage_bin"].astype(str) == "early", "sample_id"].astype(str).tolist()
    if len(advanced) < MIN_PER_GROUP or len(early) < MIN_PER_GROUP:
        return None
    mat = raw.loc[~raw["protein_id"].map(is_summary_protein)].drop_duplicates("protein_id", keep="first")
    mat = mat.set_index("protein_id")
    coef, s2, df, stdev, n_obs = unpaired_observed_stats(
        _numeric_block(mat, advanced),
        _numeric_block(mat, early),
        min_n=3,
    )
    table, _d0, _s02 = moderate_table(
        protein_ids=mat.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method="moderated_t_treat_unpaired_stage_unshared",
    )
    table["n_advanced"] = len(advanced)
    table["n_early"] = len(early)
    return table


def _plasma_contrast(path, meta: pd.DataFrame, *, lfc: float, alpha: float) -> pd.DataFrame | None:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    labels = meta.set_index("sample_id")["condition"].astype(str)
    crc = [c for c in raw.columns if c != "protein_id" and labels.get(c) == "crc"]
    healthy = [c for c in raw.columns if c != "protein_id" and labels.get(c) == "healthy"]
    if len(crc) < MIN_PER_GROUP or len(healthy) < MIN_PER_GROUP:
        return None
    mat = raw.drop_duplicates("protein_id", keep="first").set_index("protein_id")
    a = _numeric_block(mat, crc)
    b = _numeric_block(mat, healthy)
    with np.errstate(divide="ignore", invalid="ignore"):
        group_a = np.where(a > 0, np.log2(a), np.nan)
        group_b = np.where(b > 0, np.log2(b), np.nan)
    coef, s2, df, stdev, n_obs = unpaired_observed_stats(group_a, group_b, min_n=3)
    table, _d0, _s02 = moderate_table(
        protein_ids=mat.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method="moderated_t_treat_unpaired_plasma_log2_observed",
    )
    table["n_crc"] = len(crc)
    table["n_healthy"] = len(healthy)
    return table


def _inventory(meta: pd.DataFrame) -> list[dict]:
    """Qué variables de la clínica pública alcanzan para un contraste. El umbral es 15 por lado."""
    u = meta.loc[meta["sample_id"].astype(str).str.endswith(" Unshared Log Ratio")].copy()
    tumor = u.loc[u["condition"].astype(str) == "tumor"]

    def n_pairs(sex: str) -> int:
        sub = u.loc[u["sex"].astype(str) == sex]
        n = 0
        for _, g in sub.groupby("participant_id"):
            if (g["condition"].astype(str) == "tumor").sum() == 1 and (g["condition"].astype(str) == "normal").sum() == 1:
                n += 1
        return n

    rows = [
        {"variable": "tumor_vs_normal_pareado", "aplica": True, "detalle": "Contraste primario ya certificado. 94 pares. No se repite aquí."},
        {"variable": "sexo_estratificado", "aplica": True, "detalle": f"Pares mujer={n_pairs('female')}, hombre={n_pairs('male')}. Réplica interna del mismo contraste."},
        {"variable": "estadio_avanzado_vs_temprano", "aplica": True, "detalle": f"Tumores advanced={(tumor['stage_bin'].astype(str)=='advanced').sum()}, early={(tumor['stage_bin'].astype(str)=='early').sum()}."},
        {"variable": "plasma_crc_vs_sano", "aplica": True, "detalle": "Otro proteoma (suero DIA). No entra en H1."},
        {"variable": "edad", "aplica": False, "detalle": "Edad vacía en las 198 filas unshared."},
        {"variable": "msi", "aplica": False, "detalle": "msi_status = unknown en todas las filas."},
        {"variable": "cms", "aplica": False, "detalle": "cms_subtype_bin vacío."},
        {"variable": "sitio_tumor", "aplica": False, "detalle": "Todos los especímenes con sitio informado son colon. No hay recto para contrastar."},
        {"variable": "metastasis", "aplica": False, "detalle": "metastatic_status = unknown."},
        {"variable": "respuesta_tratamiento", "aplica": False, "detalle": "treatment_response = unknown."},
        {"variable": "pais_centro", "aplica": False, "detalle": "country y enrollment_site vacíos."},
        {"variable": "lote", "aplica": False, "detalle": "Un solo lote, pdc_tmt10. No hay lote que ajustar."},
        {
            "variable": "pureza_tumoral",
            "aplica": True,
            "detalle": "El porcentaje del patólogo no está en la clínica pública. Se usa un sustituto estromal "
            "del propio proteoma (matriz extracelular, GO:0031012) como covariable de par en R4.",
        },
        {
            "variable": "polipo_adenoma",
            "aplica": True,
            "detalle": "PXD034575 trae 16 pólipos emparejados con el normal y el tumor del mismo paciente. "
            "Se contrasta en la etapa translational_axes.",
        },
        {
            "variable": "etnicidad_hispana",
            "aplica": False,
            "detalle": f"En tumor: hispanic={(tumor['hispanic_latino_bin'].astype(str)=='hispanic_latino').sum()} frente a non={(tumor['hispanic_latino_bin'].astype(str)=='non_hispanic_latino').sum()}. Por debajo de {MIN_PER_GROUP}.",
        },
        {
            "variable": "raza_census",
            "aplica": False,
            "detalle": "Es una categoría de censo, no ancestría genética. En tumor, asian=15 y black=7. No se contrasta.",
        },
    ]
    return rows


def run(settings: Settings) -> None:
    cfg = settings.config.get("certification", {}) or {}
    alpha = float(cfg.get("alpha", 0.05))
    lfc = float(cfg.get("treat_lfc", TREAT_LFC))
    out = settings.resolve("results") / "secondary"
    out.mkdir(parents=True, exist_ok=True)

    cptac_ds = settings.datasets["datasets"]["cptac_colorectal"]
    raw = pd.read_csv(settings.resolve(cptac_ds["matrix_path"]), sep="\t", low_memory=False)
    meta = pd.read_csv(settings.sample_metadata_path("cptac_colorectal"), sep="\t")
    inventory = _inventory(meta)
    write_tsv(pd.DataFrame(inventory), out / "variables_inventario.tsv")

    audit: dict = {"min_per_group": MIN_PER_GROUP, "treat_lfc": lfc, "alpha": alpha, "contrasts": {}}

    for sex, label in (("female", "mujer"), ("male", "hombre")):
        sm = meta.loc[meta["sex"].astype(str) == sex].copy()
        table, _diffs, pairs, d0 = _cptac_contrast(
            raw,
            sm,
            family="unshared_log_ratio",
            stat_method=f"moderated_t_treat_paired_unshared_{sex}",
            lfc=lfc,
            alpha=alpha,
        )
        write_tsv(table, out / f"cptac_unshared_treat_{sex}.tsv")
        n_pairs = len(pairs)
        audit["contrasts"][sex] = {
            "label": label,
            "n_pairs": n_pairs,
            "prior_df": d0,
            "n_significant": int(table["significant"].sum()),
            "n_tested": int(table["p_value"].notna().sum()),
            "n_up": int(((table["significant"]) & (table["log2fc"] > 0)).sum()),
            "n_down": int(((table["significant"]) & (table["log2fc"] < 0)).sum()),
        }
        print(f"[secondary] sexo {label}: pares={n_pairs} significativas={audit['contrasts'][sex]['n_significant']}")

    female = pd.read_csv(out / "cptac_unshared_treat_female.tsv", sep="\t")
    male = pd.read_csv(out / "cptac_unshared_treat_male.tsv", sep="\t")
    agreement = _sign_agreement(female, male, alpha=alpha)
    audit["sex_sign_agreement_female_reference"] = agreement
    print(f"[secondary] concordancia de signo mujer→hombre: {agreement}")

    stage = _stage_contrast(raw, meta, lfc=lfc, alpha=alpha)
    if stage is not None:
        write_tsv(stage, out / "cptac_stage_advanced_vs_early_treat.tsv")
        audit["contrasts"]["stage"] = {
            "n_advanced": int(stage["n_advanced"].iloc[0]),
            "n_early": int(stage["n_early"].iloc[0]),
            "n_significant": int(stage["significant"].sum()),
            "n_tested": int(stage["p_value"].notna().sum()),
        }
        print(f"[secondary] estadio: {audit['contrasts']['stage']}")

    plasma_ds = settings.datasets["datasets"]["plasma_crc_2024"]
    plasma_meta = pd.read_csv(settings.sample_metadata_path("plasma_crc_2024"), sep="\t")
    plasma = _plasma_contrast(settings.resolve(plasma_ds["matrix_path"]), plasma_meta, lfc=lfc, alpha=alpha)
    if plasma is not None:
        write_tsv(plasma, out / "plasma_crc_vs_healthy_treat.tsv")
        audit["contrasts"]["plasma"] = {
            "n_crc": int(plasma["n_crc"].iloc[0]),
            "n_healthy": int(plasma["n_healthy"].iloc[0]),
            "n_significant": int(plasma["significant"].sum()),
            "n_tested": int(plasma["p_value"].notna().sum()),
        }
        print(f"[secondary] plasma: {audit['contrasts']['plasma']}")

    (out / "secondary_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[secondary] escrito en {out}")
