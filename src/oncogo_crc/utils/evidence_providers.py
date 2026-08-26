"""Proveedores pluggables de evidencia externa para la Capa C (priorización convergente).

Cada proveedor lee un **snapshot local** (TSV versionado en ``metadata/evidence/``) y devuelve una
tabla estandarizada indexada por ``gene_symbol`` (mayúsculas). Si el snapshot no existe, el proveedor
se marca ``available=False`` y **no penaliza** el score (diseño no-promisorio).

Contrato de cada snapshot (columnas mínimas; los nombres alternativos se detectan por heurística):

- ``depmap``      → ``gene_symbol``, ``chronos_mean`` (más negativo = más esencial); opc. ``essential_fraction``
- ``opentargets`` → ``gene_symbol``, ``tractability`` (0–1)
- ``dgidb``       → ``gene_symbol``, ``n_drugs`` (entero); opc. ``drugs``
- ``hpa``         → ``gene_symbol``, ``tumor_specificity`` (0–1 o score continuo)
- ``gtex``        → ``gene_symbol``, ``max_normal_tpm`` (menor = más seguro)
- ``cbioportal``  → ``gene_symbol``, ``alteration_frequency`` (0–1)

Cómo obtener cada snapshot: ver ``metadata/evidence/README.md``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class ProviderResult:
    """Resultado estandarizado de un proveedor de evidencia."""

    name: str
    available: bool
    df: pd.DataFrame  # index=gene_symbol (upper); columnas normalizadas del proveedor
    version: str = "unknown"
    path: str | None = None
    n_genes: int = 0
    columns: list[str] = field(default_factory=list)
    status: str = "ok"


def _pick_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def _read_version(path: Path) -> str:
    """Versión desde un sidecar ``<file>.version`` o mtime del snapshot."""
    sidecar = path.with_suffix(path.suffix + ".version")
    if sidecar.exists():
        return sidecar.read_text(encoding="utf-8").strip()
    try:
        import datetime as _dt

        return "snapshot_" + _dt.datetime.utcfromtimestamp(path.stat().st_mtime).strftime("%Y%m%d")
    except OSError:
        return "unknown"


def _gene_index(df: pd.DataFrame, gene_col: str) -> pd.DataFrame:
    df = df.copy()
    df["gene_symbol"] = df[gene_col].astype(str).str.strip().str.upper()
    df = df[df["gene_symbol"].notna() & (df["gene_symbol"] != "") & (df["gene_symbol"] != "NAN")]
    df = df.drop_duplicates("gene_symbol").set_index("gene_symbol")
    return df


# --- Especificación declarativa de proveedores -----------------------------------------------
# name: (filename, gene_col_candidates, {out_col: (value_col_candidates, dtype)})
_SPEC: dict[str, tuple[str, list[str], dict[str, tuple[list[str], str]]]] = {
    "depmap": (
        "depmap_crc_gene_effect.tsv",
        ["gene_symbol", "gene", "Gene"],
        {
            "depmap_chronos": (["chronos_mean", "gene_effect", "chronos", "effect"], "float"),
            "depmap_essential_fraction": (["essential_fraction", "frac_essential"], "float"),
        },
    ),
    "opentargets": (
        "opentargets_tractability.tsv",
        ["gene_symbol", "symbol", "target", "gene"],
        {
            "opentargets_tractability": (
                ["tractability", "overall_tractability", "tractability_score", "max_tractability"],
                "float",
            ),
            "opentargets_association": (["association_score", "overall_association"], "float"),
        },
    ),
    "dgidb": (
        "dgidb_gene_drugs.tsv",
        ["gene_symbol", "gene", "gene_name"],
        {
            "dgidb_n_drugs": (["n_drugs", "num_drugs", "interaction_count"], "float"),
        },
    ),
    "hpa": (
        "hpa_tumor_specificity.tsv",
        ["gene_symbol", "Gene", "gene"],
        {
            "hpa_tumor_specificity": (
                ["tumor_specificity", "specificity_score", "tumor_specificity_score"],
                "float",
            ),
        },
    ),
    "gtex": (
        "gtex_normal_max.tsv",
        ["gene_symbol", "Description", "gene"],
        {
            "gtex_max_normal_tpm": (["max_normal_tpm", "max_tpm", "max_tissue_tpm"], "float"),
        },
    ),
    "cbioportal": (
        "cbioportal_alteration_freq.tsv",
        ["gene_symbol", "Hugo_Symbol", "gene"],
        {
            "cbioportal_alteration_frequency": (
                ["alteration_frequency", "alt_freq", "frequency"],
                "float",
            ),
        },
    ),
}


def load_provider(name: str, evidence_dir: Path) -> ProviderResult:
    """Carga un proveedor por nombre desde ``evidence_dir``; tolera ausencia del snapshot."""
    if name not in _SPEC:
        return ProviderResult(name=name, available=False, df=pd.DataFrame(), status="unknown_provider")
    filename, gene_cands, out_map = _SPEC[name]
    path = evidence_dir / filename
    if not path.exists():
        return ProviderResult(name=name, available=False, df=pd.DataFrame(), status="snapshot_not_found")
    try:
        raw = pd.read_csv(path, sep="\t")
    except (OSError, pd.errors.ParserError, ValueError) as exc:  # pragma: no cover - defensivo
        return ProviderResult(name=name, available=False, df=pd.DataFrame(), status=f"read_error:{exc}")

    gene_col = _pick_col(raw, gene_cands)
    if gene_col is None:
        return ProviderResult(name=name, available=False, df=pd.DataFrame(), status="missing_gene_column")

    idx = _gene_index(raw, gene_col)
    out = pd.DataFrame(index=idx.index)
    used_cols: list[str] = []
    for out_col, (cands, dtype) in out_map.items():
        src = _pick_col(raw, cands)
        if src is None:
            continue
        series = idx[src]
        if dtype == "float":
            series = pd.to_numeric(series, errors="coerce")
        out[out_col] = series
        used_cols.append(out_col)

    if not used_cols:
        return ProviderResult(name=name, available=False, df=pd.DataFrame(), status="no_value_columns")

    return ProviderResult(
        name=name,
        available=True,
        df=out,
        version=_read_version(path),
        path=str(path),
        n_genes=int(len(out)),
        columns=used_cols,
        status="ok",
    )


def available_provider_names() -> list[str]:
    return list(_SPEC.keys())


def minmax_favorability(series: pd.Series, *, invert: bool = False) -> pd.Series:
    """Normaliza a [0,1] ignorando NaN. ``invert=True`` para señales donde *menor = mejor*.

    Devuelve NaN donde el valor original es NaN (evidencia ausente no cuenta en el promedio).
    """
    s = pd.to_numeric(series, errors="coerce").astype(float)
    valid = s.dropna()
    if valid.empty:
        return pd.Series(np.nan, index=series.index)
    lo, hi = float(valid.min()), float(valid.max())
    if hi <= lo:
        norm = pd.Series(np.where(s.notna(), 1.0, np.nan), index=series.index)
    else:
        norm = (s - lo) / (hi - lo)
    if invert:
        norm = 1.0 - norm
    return norm
