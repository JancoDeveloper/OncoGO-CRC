"""Contraste certificable sobre las matrices públicas ya descargadas.

No reescribe los ficheros crudos, no imputa ausencias y no fusiona cohortes.
CPTAC entra como log-ratio (ya está en log2): no se vuelve a transformar.
La cuantificación primaria es ``Unshared Log Ratio`` (péptidos no compartidos).
``Log Ratio`` completo se guarda como sensibilidad y no decide las hipótesis.
PRIDE entra como log2 de intensidades LFQ positivas; el cero de MaxQuant es ausencia.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from oncogo_crc.certification.contract import TREAT_LFC
from oncogo_crc.settings import Settings
from oncogo_crc.utils.go_utils import _build_term_genes, load_go_long, ora_fisher
from oncogo_crc.utils.ids import effective_protein_mapping
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.moderated import (
    moderate_table,
    n_significant_after_signs,
    paired_covariate_stats,
    paired_observed_stats,
    posterior_scale,
    fit_f_prior,
    unpaired_observed_stats,
)
from oncogo_crc.utils.quant_columns import columns_for_family, is_summary_protein
from oncogo_crc.utils.reproducibility_overlap import (
    cptac_pride_overlap_tables,
    load_gene_uniprot_from_mapping,
    uniprot_base,
    uniprot_to_gene_dict,
)
from oncogo_crc.utils.stats import build_paired_sample_order

# Fijados con el diseño, antes de mirar este contraste. No se mueven después.
SEED = 20260914
N_PERM = 100
MIN_COMPLETE = 3
MIN_TERM_GENES = 10
MIN_DIRECTION_GENES = 5
UNIQUE_PEPTIDES_MIN = 1
ECM_GO_TERM = "GO:0031012"  # matriz extracelular; el eje estromal de la proteómica bulk
PROTEIN_GROUPS_URL = (
    "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/01/PXD034575/Cohort1_proteinGroups.txt"
)


def _numeric_block(frame: pd.DataFrame, columns: list[str]) -> np.ndarray:
    block = frame.loc[:, columns].apply(pd.to_numeric, errors="coerce")
    return block.to_numpy(dtype=float)


def _cptac_contrast(
    raw: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    family: str,
    stat_method: str,
    lfc: float,
    alpha: float,
) -> tuple[pd.DataFrame, np.ndarray, list[tuple[str, str]], float]:
    sample_cols = [c for c in raw.columns if c != "protein_id"]
    cols = columns_for_family(sample_cols, family)
    sm = meta.loc[meta["sample_id"].astype(str).isin(cols)].copy()
    pairs = build_paired_sample_order(
        sm,
        condition_column="condition",
        group_a="tumor",
        group_b="normal",
        participant_col="participant_id",
    )
    if not pairs:
        raise RuntimeError(f"CPTAC {family}: no hay pares tumor-normal con un espécimen por lado.")
    mat = raw.loc[~raw["protein_id"].map(is_summary_protein)].copy()
    mat = mat.drop_duplicates("protein_id", keep="first").set_index("protein_id")
    a_cols = [a for a, _b in pairs]
    b_cols = [b for _a, b in pairs]
    diffs = _numeric_block(mat, a_cols) - _numeric_block(mat, b_cols)
    coef, s2, df, stdev, n_obs = paired_observed_stats(diffs, min_pairs=MIN_COMPLETE)
    table, d0, _s02 = moderate_table(
        protein_ids=mat.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method=stat_method,
    )
    table["comparison_id"] = f"cptac_{family}"
    table["quant_family"] = family
    table["n_pairs_design"] = len(pairs)
    return table, diffs, pairs, d0


def stromal_score(raw: pd.DataFrame, columns: list[str], ecm_genes: set[str]) -> pd.Series:
    """Contenido estromal por espécimen: media de proteínas de matriz extracelular, z por proteína.

    Sustituto del porcentaje de pureza, que la clínica pública de PDC no trae.
    La proteómica bulk de CCR se separa por contenido estromal y mesenquimal
    (CAVIN1, PRELP), así que un efecto puede ser composición de tejido.
    """
    mat = raw.loc[~raw["protein_id"].map(is_summary_protein)].drop_duplicates("protein_id", keep="first")
    mat = mat.set_index("protein_id")
    present = [g for g in mat.index.astype(str) if g.upper() in ecm_genes]
    block = pd.DataFrame(_numeric_block(mat.loc[present], columns), index=present, columns=columns)
    mu = block.mean(axis=1)
    sd = block.std(axis=1, ddof=1).replace(0.0, np.nan)
    z = block.sub(mu, axis=0).div(sd, axis=0)
    return z.mean(axis=0, skipna=True)


def _pride_contrast(path: Path, meta: pd.DataFrame, *, lfc: float, alpha: float) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    raw = raw.loc[~raw["protein_id"].map(is_summary_protein)].drop_duplicates("protein_id", keep="first")
    labels = meta.set_index("sample_id")["condition"].astype(str)
    tumor = [c for c in raw.columns if c != "protein_id" and labels.get(c) == "tumor"]
    normal = [c for c in raw.columns if c != "protein_id" and labels.get(c) == "normal"]
    mat = raw.set_index("protein_id")
    # LFQ 0 es ausencia en MaxQuant. No se suma un offset.
    a = _numeric_block(mat, tumor)
    b = _numeric_block(mat, normal)
    with np.errstate(divide="ignore", invalid="ignore"):
        group_a = np.where(a > 0, np.log2(a), np.nan)
        group_b = np.where(b > 0, np.log2(b), np.nan)
    coef, s2, df, stdev, n_obs = unpaired_observed_stats(group_a, group_b, min_n=MIN_COMPLETE)
    table, _d0, _s02 = moderate_table(
        protein_ids=mat.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method="moderated_t_treat_unpaired_lfq_observed",
        )
    table["comparison_id"] = "pride_tumor_vs_normal"
    return table


def _unique_peptide_counts(settings: Settings) -> dict[str, int] | None:
    """Recupera Unique peptides del proteinGroups público, sin tocar la matriz LFQ ya reducida."""
    slim = settings.resolve("data/raw/pride/cohort1_unique_peptides.tsv")
    if slim.exists() and slim.stat().st_size > 0:
        tab = pd.read_csv(slim, sep="\t")
    else:
        dest = settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt")
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists() or dest.stat().st_size < 1_000_000:
            print("[certified] descargando proteinGroups de PXD034575 para Unique peptides.")
            try:
                urllib.request.urlretrieve(PROTEIN_GROUPS_URL, dest)
            except Exception as exc:
                print(f"[certified] sin péptidos únicos: {exc}")
                return None
        header = pd.read_csv(dest, sep="\t", nrows=0).columns.tolist()
        id_col = next((c for c in ("Protein IDs", "Majority protein IDs") if c in header), None)
        pep_col = next((c for c in header if c.lower() == "unique peptides"), None)
        if id_col is None or pep_col is None:
            print("[certified] proteinGroups no trae Protein IDs y Unique peptides.")
            return None
        tab = pd.read_csv(dest, sep="\t", usecols=[id_col, pep_col], low_memory=False)
        tab = tab.rename(columns={id_col: "protein_id", pep_col: "unique_peptides"})
        tab["uniprot"] = tab["protein_id"].map(lambda v: uniprot_base(str(v).split(";")[0]))
        tab = tab.dropna(subset=["uniprot"])
        tab["unique_peptides"] = pd.to_numeric(tab["unique_peptides"], errors="coerce")
        tab = tab.groupby("uniprot", as_index=False)["unique_peptides"].max()
        write_tsv(tab, slim)
    if "uniprot" not in tab.columns:
        return None
    out: dict[str, int] = {}
    for _, row in tab.iterrows():
        if pd.isna(row["unique_peptides"]):
            continue
        out[str(row["uniprot"])] = int(row["unique_peptides"])
    return out


def _peptide_fraction(
    cptac: pd.DataFrame,
    gene_to_uni: dict[str, str],
    unique_peptides: dict[str, int],
    *,
    alpha: float,
) -> tuple[float | None, int, int]:
    sig = cptac.loc[cptac["padj"].astype(float) < alpha, "protein_id"].astype(str)
    n = int(len(sig))
    if n == 0:
        return None, 0, 0
    ok = 0
    for gene in sig:
        uni = gene_to_uni.get(gene.strip().upper())
        if uni is None:
            continue
        n_pep = unique_peptides.get(uni)
        if n_pep is not None and n_pep >= UNIQUE_PEPTIDES_MIN:
            ok += 1
    return ok / n, ok, n


def _pathway_concordance(
    cptac: pd.DataFrame,
    pride: pd.DataFrame,
    term_to_genes: dict[str, set[str]],
    uni_to_gene: dict[str, str],
    *,
    alpha: float,
) -> tuple[float | None, int, int, float | None, float | None]:
    """Vías BP significativas en CPTAC: ¿el signo medio de sus genes coincide en PRIDE?

    El universo de cada cohorte son los genes con t moderada finita.
    Un término entra en el denominador solo si hay al menos 5 genes con
    log2FC en las dos cohortes. No se elige el término mirando PRIDE.
    """
    c = cptac.copy()
    c["gene"] = c["protein_id"].astype(str).str.upper()
    p = pride.copy()
    p["gene"] = p["protein_id"].map(lambda v: uni_to_gene.get(uniprot_base(str(v)) or "", ""))
    p = p.loc[p["gene"].astype(str).str.len() > 0]
    fc_c = c.dropna(subset=["log2fc"]).drop_duplicates("gene").set_index("gene")["log2fc"]
    fc_p = p.dropna(subset=["log2fc"]).drop_duplicates("gene").set_index("gene")["log2fc"]
    tested = set(c.loc[c["t_stat"].notna(), "gene"])
    sig = set(c.loc[c["padj"].astype(float) < alpha, "gene"])
    ora = ora_fisher(sig, tested, term_to_genes, min_genes_in_term=MIN_TERM_GENES)
    if ora.empty or "significant" not in ora.columns:
        return None, 0, 0, None, None
    same = 0
    n = 0
    n_sig_terms = int(ora["significant"].sum())
    for term in ora.loc[ora["significant"], "go_id"]:
        members = term_to_genes.get(str(term), set())
        fg = [g for g in members if g in sig and g in fc_c.index and float(fc_c[g]) != 0]
        if len(fg) < MIN_DIRECTION_GENES:
            continue
        fg_p = [g for g in fg if g in fc_p.index and float(fc_p[g]) != 0]
        if len(fg_p) < MIN_DIRECTION_GENES:
            continue
        dir_c = np.sign(float(np.median([fc_c[g] for g in fg])))
        dir_p = np.sign(float(np.median([fc_p[g] for g in fg_p])))
        if dir_c == 0 or dir_p == 0:
            continue
        n += 1
        same += int(dir_c == dir_p)
    if n == 0:
        return None, 0, n_sig_terms, None, None
    test = binomtest(same, n, p=0.5, alternative="two-sided")
    ci = test.proportion_ci(confidence_level=0.95, method="exact")
    return same / n, n, n_sig_terms, float(ci.low), float(ci.high)


def _permute_pair_signs(diffs: np.ndarray, *, lfc: float, alpha: float) -> tuple[float, int, float]:
    coef, s2, df, stdev, _n = paired_observed_stats(diffs, min_pairs=MIN_COMPLETE)
    d0, s02 = fit_f_prior(s2, df)
    s2_post, df_total = posterior_scale(s2, df, d0, s02)
    observed = n_significant_after_signs(
        diffs,
        np.ones(diffs.shape[1]),
        s2_post=s2_post,
        df_total=df_total,
        stdev_unscaled=stdev,
        lfc=lfc,
        alpha=alpha,
    )
    rng = np.random.default_rng(SEED)
    extreme = 0
    null_counts: list[int] = []
    for _ in range(N_PERM):
        signs = rng.choice(np.array([-1.0, 1.0]), size=diffs.shape[1])
        n_sig = n_significant_after_signs(
            diffs,
            signs,
            s2_post=s2_post,
            df_total=df_total,
            stdev_unscaled=stdev,
            lfc=lfc,
            alpha=alpha,
        )
        null_counts.append(n_sig)
        extreme += int(n_sig >= observed)
    empirical_p = (1 + extreme) / (N_PERM + 1)
    return empirical_p, observed, float(np.median(null_counts))


def purity_adjusted_fraction(
    raw: pd.DataFrame,
    pairs: list[tuple[str, str]],
    diffs: np.ndarray,
    protein_ids: np.ndarray,
    certified_sig: set[str],
    ecm_genes: set[str],
    *,
    lfc: float,
    alpha: float,
) -> tuple[pd.DataFrame, dict]:
    """R4: qué fracción del efecto certificado sobrevive al ajuste por contenido estromal.

    El efecto se evalúa en diferencia estromal cero, es decir, el efecto que
    quedaría si el tumor y su normal tuvieran la misma composición.
    """
    cols = [c for pair in pairs for c in pair]
    score = stromal_score(raw, cols, ecm_genes)
    delta = np.array([float(score.get(a, np.nan)) - float(score.get(b, np.nan)) for a, b in pairs])
    coef, s2, df, stdev, n_used, slope, slope_t = paired_covariate_stats(
        diffs, delta, min_pairs=MIN_COMPLETE, at=0.0
    )
    table, d0, _s02 = moderate_table(
        protein_ids=protein_ids,
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_used,
        lfc=lfc,
        alpha=alpha,
        stat_method="moderated_t_treat_paired_unshared_stroma_adjusted",
    )
    table["stroma_slope"] = slope
    table["stroma_slope_t"] = slope_t
    table["certified_unadjusted"] = table["protein_id"].astype(str).isin(certified_sig)
    cert_rows = table.loc[table["certified_unadjusted"]]
    kept = cert_rows.loc[cert_rows["significant"]]
    fraction = (len(kept) / len(certified_sig)) if certified_sig else None
    # |t| > 2 sobre 92 grados de libertad: el efecto del par se mueve con su contraste estromal.
    stroma_driven = int((cert_rows["stroma_slope_t"].abs() > 2.0).sum())
    stats = {
        "ecm_proteins_used": int(sum(1 for g in raw["protein_id"].astype(str) if g.upper() in ecm_genes)),
        "stroma_delta_mean": float(np.nanmean(delta)),
        "stroma_delta_sd": float(np.nanstd(delta, ddof=1)),
        "evaluated_at_stroma_delta": 0.0,
        "prior_df": d0,
        "n_certified": len(certified_sig),
        "n_certified_surviving": int(len(kept)),
        "fraction_surviving": fraction,
        "n_certified_stroma_dependent_slope": stroma_driven,
        "fraction_stroma_dependent_slope": (stroma_driven / len(certified_sig)) if certified_sig else None,
        "median_abs_log2fc_adjusted": float(cert_rows["log2fc"].abs().median()),
        "n_significant_adjusted_total": int(table["significant"].sum()),
    }
    return table, stats


def _gaps(*, peptides_ok: bool, pathway_ok: bool) -> list[str]:
    gaps = [
        "Las cohortes no se fusionaron y no hay meta-análisis REML. R3 no es calculable.",
        "R4 no usa el porcentaje de pureza del patólogo: la clínica pública de PDC no lo trae. "
        "Usa un sustituto estromal calculado del propio proteoma (proteínas de matriz extracelular, "
        f"{ECM_GO_TERM}) como covariable de par. Es un ajuste por composición de tejido, no una pureza medida.",
        "La matriz CPTAC pública es un log-ratio de proteína, no las intensidades de cada canal TMT. msqrob2 con el plex como efecto aleatorio no se puede ajustar sin inventar esos canales.",
        "PXD034575 es el control de independencia del diseño, no la réplica primaria PXD019504. Este cruce no sustituye ese diseño.",
        "H0 nivel A (checksum PXD019504 × PXD031705) no se repitió en esta etapa.",
        "El contraste anterior (segundo log2 sobre el log-ratio, Log Ratio y Unshared como muestras distintas, t sobre matriz imputada) no entra en este certificado.",
        f"R5 son {N_PERM} permutaciones del signo de los pares tumor-normal en CPTAC, semilla {SEED}. No permuta la réplica ni el enriquecimiento.",
    ]
    if peptides_ok:
        gaps.append(
            "La ventana de péptido usa Unshared Log Ratio en CPTAC y Unique peptides ≥ 1 en el proteinGroups de PRIDE. "
            "No demuestra que la misma secuencia se haya visto en las dos plataformas."
        )
    else:
        gaps.append(
            "No se recuperó la columna Unique peptides. La ventana de péptido proteotípico queda sin calcular."
        )
    if not pathway_ok:
        gaps.append(
            "La concordancia de vía no se pudo calcular: ningún término GO-BP significativo en CPTAC "
            "tuvo al menos 5 genes con log2FC en las dos cohortes."
        )
    return gaps


def run(settings: Settings) -> None:
    cfg = settings.config.get("certification", {}) or {}
    alpha = float(cfg.get("alpha", 0.05))
    lfc = float(cfg.get("treat_lfc", TREAT_LFC))
    out = settings.resolve("results") / "certified"
    out.mkdir(parents=True, exist_ok=True)

    cptac_ds = settings.datasets["datasets"]["cptac_colorectal"]
    pride_ds = settings.datasets["datasets"]["pride_crc_meta"]
    raw = pd.read_csv(settings.resolve(cptac_ds["matrix_path"]), sep="\t", low_memory=False)
    meta = pd.read_csv(settings.sample_metadata_path("cptac_colorectal"), sep="\t")
    print("[certified] CPTAC: Unshared Log Ratio, pares tumor-normal, sin segundo log2 y sin imputar.")
    cptac, diffs, pairs, d0 = _cptac_contrast(
        raw,
        meta,
        family="unshared_log_ratio",
        stat_method="moderated_t_treat_paired_unshared_logratio",
        lfc=lfc,
        alpha=alpha,
    )
    n_pairs = len(pairs)
    print(
        f"[certified] CPTAC unshared: pares={n_pairs} d0={d0:.2f} "
        f"significativas={int(cptac['significant'].sum())}/{len(cptac)}"
    )
    sensitivity, _diffs_all, pairs_all, _d0_all = _cptac_contrast(
        raw,
        meta,
        family="log_ratio",
        stat_method="moderated_t_treat_paired_all_peptides_sensitivity",
        lfc=lfc,
        alpha=alpha,
    )
    print(
        f"[certified] sensibilidad Log Ratio (no decide H1): "
        f"pares={len(pairs_all)} significativas={int(sensitivity['significant'].sum())}"
    )

    pride_meta = pd.read_csv(settings.sample_metadata_path("pride_crc_meta"), sep="\t")
    pride = _pride_contrast(settings.resolve(pride_ds["matrix_path"]), pride_meta, lfc=lfc, alpha=alpha)
    print(f"[certified] PRIDE LFQ observado: significativas={int(pride['significant'].sum())}/{len(pride)}")

    write_tsv(cptac, out / "cptac_unshared_treat.tsv")
    write_tsv(sensitivity, out / "cptac_all_peptides_treat_sensitivity.tsv")
    write_tsv(pride, out / "pride_lfq_treat.tsv")

    mapping = effective_protein_mapping(settings)
    gene_to_uni = load_gene_uniprot_from_mapping(mapping, exclude_ambiguous=True)
    uni_to_gene = uniprot_to_gene_dict(mapping)
    merged, summary = cptac_pride_overlap_tables(cptac, pride, gene_to_uni=gene_to_uni, alpha=alpha)
    write_tsv(merged, out / "cptac_pride_overlap.tsv")
    write_tsv(summary, out / "cptac_pride_overlap_summary.tsv")

    peptides = _unique_peptide_counts(settings)
    peptide_fraction = None
    peptide_ok_n = 0
    peptide_den = 0
    if peptides:
        peptide_fraction, peptide_ok_n, peptide_den = _peptide_fraction(
            cptac, gene_to_uni, peptides, alpha=alpha
        )
        print(
            f"[certified] péptido único en PRIDE entre significativas CPTAC: "
            f"{peptide_ok_n}/{peptide_den} = {peptide_fraction}"
        )

    pathway = (None, 0, 0, None, None)
    purity_stats: dict = {}
    go_path = settings.go_gene_long_path()
    if go_path.exists():
        go_long = load_go_long(go_path)
        term2g = {
            term: {g.upper() for g in genes}
            for term, genes in _build_term_genes(go_long, ["BP"]).items()
        }
        pathway = _pathway_concordance(cptac, pride, term2g, uni_to_gene, alpha=alpha)
        print(f"[certified] vía BP: concordancia={pathway[0]} n={pathway[1]} términos_sig={pathway[2]}")
        ecm_genes = {
            str(g).strip().upper()
            for g in go_long.loc[go_long["go_id"].astype(str) == ECM_GO_TERM, "gene_symbol"]
        }
        certified_sig = set(cptac.loc[cptac["significant"], "protein_id"].astype(str))
        adjusted, purity_stats = purity_adjusted_fraction(
            raw,
            pairs,
            diffs,
            cptac["protein_id"].to_numpy(),
            certified_sig,
            ecm_genes,
            lfc=lfc,
            alpha=alpha,
        )
        write_tsv(adjusted, out / "cptac_unshared_treat_stroma_adjusted.tsv")
        print(
            f"[certified] R4 estroma: sobreviven {purity_stats['n_certified_surviving']}/"
            f"{purity_stats['n_certified']} = {purity_stats['fraction_surviving']:.4f} "
            f"(proteínas ECM usadas={purity_stats['ecm_proteins_used']})"
        )
    else:
        print(f"[certified] sin anotación GO en {go_path}")

    empirical_p, observed_nsig, null_median = _permute_pair_signs(diffs, lfc=lfc, alpha=alpha)
    print(
        f"[certified] R5: significativas observadas={observed_nsig} "
        f"mediana nula={null_median:.0f} p empírico={empirical_p:.4f}"
    )

    conc, n_terms, n_sig_terms, ci_low, ci_high = pathway
    gaps = _gaps(peptides_ok=peptides is not None, pathway_ok=conc is not None)
    audit = {
        "method": "moderated_t_treat_paired_unshared_logratio",
        "alpha": alpha,
        "treat_lfc": lfc,
        "seed": SEED,
        "n_permutations": N_PERM,
        "min_complete_observations": MIN_COMPLETE,
        "cptac_pairs": n_pairs,
        "cptac_prior_df": d0,
        "cptac_n_significant": int(cptac["significant"].sum()),
        "cptac_n_tested": int(cptac["p_value"].notna().sum()),
        "sensitivity_all_peptides_n_significant": int(sensitivity["significant"].sum()),
        "pride_n_significant": int(pride["significant"].sum()),
        "pride_n_tested": int(pride["p_value"].notna().sum()),
        "pathway_concordance": conc,
        "pathway_n": n_terms,
        "pathway_n_significant_terms": n_sig_terms,
        "pathway_ci_low": ci_low,
        "pathway_ci_high": ci_high,
        "peptide_measurable_fraction": peptide_fraction,
        "peptide_n_measurable": peptide_ok_n,
        "peptide_n_discovery_significant": peptide_den,
        "r2_note": (
            "Fracción de proteínas treat-significativas en CPTAC cuyo UniProt tiene "
            f"Unique peptides ≥ {UNIQUE_PEPTIDES_MIN} en PXD034575. "
            "No identifica la misma secuencia en las dos plataformas."
            if peptide_fraction is not None
            else "Sin columna Unique peptides. R2 vuelve al proxy de log2FC finito."
        ),
        "r4_stroma": purity_stats or None,
        "r4_note": (
            "Fracción de proteínas certificadas que siguen pasando treat+BH al añadir el "
            "contraste estromal del par como covariable. No es el porcentaje de pureza del patólogo."
        ),
        "r5_empirical_p": empirical_p,
        "r5_note": (
            f"p empírico de {N_PERM} permutaciones de signo de los pares CPTAC "
            f"(semilla {SEED}). Observadas {observed_nsig} proteínas; mediana nula {null_median:.0f}."
        ),
        "gaps": gaps,
        "summary_metrics": summary.to_dict(orient="records"),
    }
    (out / "inference_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[certified] escrito en {out}")
