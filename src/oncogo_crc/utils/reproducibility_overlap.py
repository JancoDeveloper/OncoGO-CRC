"""Lógica compartida para cruce CPTAC ∩ PRIDE (puente UniProt + mapa gen↔UniProt)."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from oncogo_crc.utils.io import read_tsv


_UNIPROT_RE = re.compile(
    # 6 caracteres estándar (p. ej. P04637, O00560) + variante 10 caracteres (segundo grupo).
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$",
    re.IGNORECASE,
)


def uniprot_base(raw: str) -> str | None:
    if not raw or str(raw).lower() in ("nan", "none", ""):
        return None
    s = str(raw).strip().split(";")[0].strip()
    if not s:
        return None
    s = s.split("-")[0].strip().upper()
    if len(s) < 6:
        return None
    return s if _UNIPROT_RE.match(s) else None


def gene_to_uniprot_dict(mapping: pd.DataFrame) -> dict[str, str]:
    out: dict[str, str] = {}
    if mapping is None or mapping.empty:
        return out
    if "gene_symbol" not in mapping.columns or "uniprot_accession" not in mapping.columns:
        return out
    for _, r in mapping.iterrows():
        g = str(r.get("gene_symbol", "")).strip().upper()
        u = uniprot_base(str(r.get("uniprot_accession", "")))
        if g and u and g != "NAN":
            out[g] = u
    return out


def ambiguous_gene_symbols(mapping: pd.DataFrame) -> set[str]:
    """Genes (símbolo HGNC) que en el mapa se asocian a más de un UniProt canónico distinto."""
    amb: set[str] = set()
    if mapping is None or mapping.empty:
        return amb
    if "gene_symbol" not in mapping.columns or "uniprot_accession" not in mapping.columns:
        return amb
    tmp = mapping.copy()
    tmp["_g"] = tmp["gene_symbol"].astype(str).str.strip().str.upper()
    tmp["_u"] = tmp["uniprot_accession"].map(lambda x: uniprot_base(str(x)))
    tmp = tmp[tmp["_g"].notna() & (tmp["_g"] != "") & (tmp["_g"] != "NAN")]
    tmp = tmp[tmp["_u"].notna()]
    for g, grp in tmp.groupby("_g"):
        if len(set(grp["_u"].tolist())) > 1:
            amb.add(str(g))
    return amb


def uniprot_to_gene_dict(mapping: pd.DataFrame) -> dict[str, str]:
    """Primera aparición de cada accesión UniProt canónica → símbolo HGNC (para matrices PRIDE en UniProt)."""
    out: dict[str, str] = {}
    if mapping is None or mapping.empty:
        return out
    if "gene_symbol" not in mapping.columns or "uniprot_accession" not in mapping.columns:
        return out
    for _, r in mapping.iterrows():
        g = str(r.get("gene_symbol", "")).strip().upper()
        u = uniprot_base(str(r.get("uniprot_accession", "")))
        if g and u and g != "NAN" and u not in out:
            out[u] = g
    return out


def join_key_cptac(pid: str, gene_to_uni: dict[str, str]) -> str | None:
    s = str(pid).strip()
    if not s or s.upper() == "NAN":
        return None
    u = gene_to_uni.get(s.upper())
    if u:
        return u
    return uniprot_base(s)


def join_key_pride(pid: str) -> str | None:
    return uniprot_base(str(pid).strip())


def cptac_pride_overlap_tables(
    c: pd.DataFrame,
    p: pd.DataFrame,
    *,
    gene_to_uni: dict[str, str],
    alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Devuelve (merged_wide, summary_metrics)."""
    c = c.copy()
    p = p.copy()
    c["_join_uni"] = c["protein_id"].map(lambda x: join_key_cptac(str(x), gene_to_uni))
    p["_join_uni"] = p["protein_id"].map(lambda x: join_key_pride(str(x)))

    def _suffix(df: pd.DataFrame, tag: str) -> pd.DataFrame:
        skip = {"protein_id", "gene_symbol", "uniprot_accession", "_join_uni"}
        cols = [x for x in df.columns if x not in skip]
        base = df[["_join_uni", "protein_id"] + cols].copy()
        ren = {col: f"{col}_{tag}" for col in cols}
        ren["protein_id"] = f"protein_id_{tag}"
        return base.rename(columns=ren)

    # Nunca hacer merge por _join_uni = NaN: pandas empareja todos los NaN entre sí → producto cartesiano
    # (millones de filas falsas) cuando el mapa gen↔UniProt está vacío o casi no cruza CPTAC.
    c_suf = _suffix(c, "cptac")
    p_suf = _suffix(p, "pride")
    c_ok = c_suf[c_suf["_join_uni"].notna()].copy()
    p_ok = p_suf[p_suf["_join_uni"].notna()].copy()
    merged = c_ok.merge(p_ok, on="_join_uni", how="outer", indicator=True)
    merged = merged.copy()
    merged["_merge"] = merged["_merge"].astype(str)

    def _sig(merged_df: pd.DataFrame, prefix: str) -> pd.Series:
        col = f"padj_{prefix}"
        if col not in merged_df.columns:
            return pd.Series(False, index=merged_df.index)
        return merged_df[col].astype(float) < alpha

    in_both = merged["_merge"] == "both"
    sig_c = _sig(merged, "cptac")
    sig_p = _sig(merged, "pride")
    n_core = int((in_both & sig_c & sig_p).sum())
    raw_merged = c[["protein_id"]].merge(p[["protein_id"]], on="protein_id", how="inner")
    n_raw = int(len(raw_merged))
    n_c_key = int(c["_join_uni"].notna().sum())
    n_p_key = int(p["_join_uni"].notna().sum())
    n_both_key = int(in_both.sum())

    rows_sum: list[dict[str, float | int]] = [
        {"metric": "n_proteins_cptac", "value": int(c["protein_id"].nunique())},
        {"metric": "n_proteins_pride", "value": int(p["protein_id"].nunique())},
        {"metric": "n_intersection_raw_protein_id", "value": n_raw},
        {"metric": "n_cptac_with_join_uni", "value": n_c_key},
        {"metric": "n_pride_with_join_uni", "value": n_p_key},
        {"metric": "n_intersection_uniprot_bridge", "value": n_both_key},
        {"metric": f"n_both_significant_padj_lt_{alpha}", "value": n_core},
        {"metric": "n_cptac_sig_only", "value": int((in_both & sig_c & ~sig_p).sum())},
        {"metric": "n_pride_sig_only", "value": int((in_both & ~sig_c & sig_p).sum())},
    ]

    # Concordancia de dirección del efecto entre cohortes (misma clave UniProt).
    lc, lp = "log2fc_cptac", "log2fc_pride"
    if lc in merged.columns and lp in merged.columns:
        sub = merged.loc[in_both, [lc, lp]].dropna()
        sub = sub[np.isfinite(sub[lc].astype(float)) & np.isfinite(sub[lp].astype(float))]
        n_pair = int(len(sub))
        rows_sum.append({"metric": "n_pairs_log2fc_cptac_pride", "value": n_pair})
        if n_pair >= 5:
            pr, _ = stats.pearsonr(sub[lc].astype(float), sub[lp].astype(float))
            sr, _ = stats.spearmanr(sub[lc].astype(float), sub[lp].astype(float))
            rows_sum.append({"metric": "log2fc_pearson_r_cptac_pride", "value": float(pr)})
            rows_sum.append({"metric": "log2fc_spearman_r_cptac_pride", "value": float(sr)})
            same_sign = int((np.sign(sub[lc].astype(float)) == np.sign(sub[lp].astype(float))).sum())
            rows_sum.append({"metric": "log2fc_same_sign_fraction", "value": same_sign / n_pair})
        else:
            rows_sum.append({"metric": "log2fc_pearson_r_cptac_pride", "value": float("nan")})
            rows_sum.append({"metric": "log2fc_spearman_r_cptac_pride", "value": float("nan")})
            rows_sum.append({"metric": "log2fc_same_sign_fraction", "value": float("nan")})

    summary = pd.DataFrame(rows_sum)
    return merged, summary


def load_gene_uniprot_from_mapping(
    m: pd.DataFrame, *, exclude_ambiguous: bool = True
) -> dict[str, str]:
    """Construye gen→UniProt canónico desde un DataFrame de mapa (p. ej. maestro + suplementos)."""
    if m is None or m.empty:
        return {}
    full = gene_to_uniprot_dict(m)
    if not exclude_ambiguous:
        return full
    amb = ambiguous_gene_symbols(m)
    return {g: u for g, u in full.items() if g not in amb}


def load_gene_uniprot_from_settings_map(
    map_path: Path, *, exclude_ambiguous: bool = True
) -> dict[str, str]:
    if not map_path.exists():
        return {}
    return load_gene_uniprot_from_mapping(read_tsv(map_path), exclude_ambiguous=exclude_ambiguous)
