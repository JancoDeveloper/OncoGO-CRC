"""Acceso a datos para la app web (Nivel 2).

Lee artefactos reales del pipeline (``results/`` + ``data/processed/``) cuando existen y, si no,
genera un **conjunto de demostración sintético** para que la app sea explorable sin ejecutar el pipeline.
Todo lo sintético queda claramente etiquetado (``is_demo=True``).
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


def default_results_dir() -> Path:
    """Directorio de resultados por defecto (env > workspace conocido > repo)."""
    env = os.environ.get("ONCOGO_CRC_RESULTS")
    if env:
        return Path(env).expanduser()
    repo_root = Path(__file__).resolve().parents[2]
    # Si hay data_workspace en config, respétalo.
    cfg = repo_root / "config" / "config.yaml"
    if cfg.exists():
        try:
            import yaml

            paths = (yaml.safe_load(cfg.read_text()) or {}).get("paths", {})
            ws = paths.get("data_workspace")
            if ws and Path(ws).expanduser().exists():
                return Path(ws).expanduser() / "results"
        except Exception:  # noqa: BLE001 - fallback silencioso
            pass
    return repo_root / "results"


@dataclass
class Dataset:
    hypotheses: pd.DataFrame
    evidence_matrix: pd.DataFrame
    differentials: dict[str, pd.DataFrame]
    provenance: dict
    processed_matrix: pd.DataFrame | None
    is_demo: bool
    results_dir: Path
    notes: list[str]


# --------------------------------------------------------------------------------------
# Lectura de artefactos reales
# --------------------------------------------------------------------------------------
def _read_tsv(path: Path) -> pd.DataFrame | None:
    try:
        return pd.read_csv(path, sep="\t")
    except (OSError, pd.errors.ParserError, ValueError):
        return None


def _load_real(results_dir: Path) -> Dataset | None:
    hyp_path = results_dir / "prioritization" / "hypothesis_priority.tsv"
    hyp = _read_tsv(hyp_path) if hyp_path.exists() else None
    if hyp is None or hyp.empty:
        return None

    notes: list[str] = [f"Datos reales desde {results_dir}"]
    evm_path = results_dir / "prioritization" / "convergence_evidence_matrix.tsv"
    evm = _read_tsv(evm_path) if evm_path.exists() else pd.DataFrame()

    prov_path = results_dir / "prioritization" / "prioritization_provenance.json"
    prov: dict = {}
    if prov_path.exists():
        try:
            prov = json.loads(prov_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            prov = {}

    diffs: dict[str, pd.DataFrame] = {}
    diff_dir = results_dir / "differential"
    if diff_dir.is_dir():
        for p in sorted(diff_dir.glob("*_diff.tsv")):
            df = _read_tsv(p)
            if df is not None and "protein_id" in df.columns:
                diffs[p.name.replace("_diff.tsv", "")] = df

    # matriz procesada (para explorador de red por correlación real)
    processed = None
    primary = str(prov.get("primary_comparison", "cptac_tumor_vs_normal"))
    # deducir dataset primario del nombre; fallback cptac_colorectal
    proc_dir = results_dir.parent / "data" / "processed"
    for cand in ("cptac_colorectal", primary.split("_")[0]):
        pth = proc_dir / f"{cand}_zscore.tsv"
        if pth.exists():
            processed = _read_tsv(pth)
            break

    return Dataset(
        hypotheses=hyp,
        evidence_matrix=evm if evm is not None else pd.DataFrame(),
        differentials=diffs,
        provenance=prov,
        processed_matrix=processed,
        is_demo=False,
        results_dir=results_dir,
        notes=notes,
    )


# --------------------------------------------------------------------------------------
# Conjunto de demostración sintético
# --------------------------------------------------------------------------------------
def _demo_dataset(results_dir: Path, seed: int = 42) -> Dataset:
    rng = np.random.default_rng(seed)
    genes = [f"GENE{ i:03d}" for i in range(1, 121)]
    # nombres reconocibles para hacer la demo más ilustrativa
    known = ["EPCAM", "CEACAM5", "MKI67", "CDH1", "VIM", "COL1A1", "TIMP1", "MMP9", "FN1", "SPP1"]
    genes[: len(known)] = known
    n = len(genes)

    log2fc = rng.normal(0, 1.4, n)
    padj = np.clip(np.abs(rng.normal(0, 0.06, n)), 1e-8, 1)
    significant = padj < 0.05
    diff = pd.DataFrame(
        {
            "protein_id": genes,
            "gene_symbol": genes,
            "log2fc": log2fc,
            "p_value": padj * rng.uniform(0.3, 1.0, n),
            "padj": padj,
            "significant": significant,
            "comparison_id": "cptac_tumor_vs_normal",
        }
    )

    signals = [
        "diff_effect",
        "diff_significance",
        "reproducible",
        "network_hub",
        "depmap_essentiality",
        "opentargets_tractability",
        "dgidb_druggability",
        "hpa_tumor_specificity",
        "gtex_safety",
        "cbioportal_genomic",
    ]
    favor = pd.DataFrame({s: rng.uniform(0, 1, n) for s in signals}, index=genes)
    # introducir NaN (evidencia ausente) para realismo
    for s in signals[4:]:
        miss = rng.uniform(0, 1, n) < 0.3
        favor.loc[favor.index[miss], s] = np.nan

    weights = {
        "diff_effect": 1.0, "diff_significance": 1.0, "reproducible": 1.5, "network_hub": 0.5,
        "depmap_essentiality": 1.5, "opentargets_tractability": 1.5, "dgidb_druggability": 1.0,
        "hpa_tumor_specificity": 1.0, "gtex_safety": 1.0, "cbioportal_genomic": 0.75,
    }
    w = np.array([weights[s] for s in signals])
    F = favor[signals].to_numpy()
    mask = ~np.isnan(F)
    num = np.nansum(np.where(mask, F * w[None, :], 0.0), axis=1)
    den = np.where(mask, w[None, :], 0.0).sum(axis=1)
    score = np.divide(num, den, out=np.zeros_like(num), where=den > 0)

    hyp = pd.DataFrame(
        {
            "gene_symbol": genes,
            "log2fc": log2fc,
            "padj": padj,
            "significant": significant,
            "convergence_score": score,
            "evidence_count": mask.sum(axis=1),
            "external_evidence_count": (~np.isnan(favor[signals[4:]].to_numpy())).sum(axis=1),
        }
    )
    for s in signals[4:]:
        hyp[s] = favor[s].values
    hyp = hyp.sort_values("convergence_score", ascending=False).reset_index(drop=True)
    hyp.insert(0, "rank", np.arange(1, len(hyp) + 1))
    q = pd.qcut(hyp["convergence_score"].rank(method="first"), 4,
                labels=["tier_D", "tier_C", "tier_B", "tier_A"]).astype(str)
    hyp.insert(4, "evidence_tier", q.values)

    long_rows = []
    for gene in genes:
        for s in signals:
            v = favor.loc[gene, s]
            long_rows.append({
                "gene_symbol": gene, "signal": s,
                "favorability": (float(v) if pd.notna(v) else np.nan),
                "weight": weights[s], "available": bool(pd.notna(v)),
            })
    evm = pd.DataFrame(long_rows)

    # matriz procesada sintética (para explorador de red por correlación)
    n_samples = 40
    base = rng.normal(0, 1, (n, n_samples))
    # inducir módulos correlacionados
    for start in range(0, n, 15):
        factor = rng.normal(0, 1, n_samples)
        for g in range(start, min(start + 15, n)):
            base[g] += factor * rng.uniform(0.5, 1.2)
    proc = pd.DataFrame(base, columns=[f"S{i:02d}" for i in range(n_samples)])
    proc.insert(0, "protein_id", genes)

    prov = {
        "generated_at": "DEMO",
        "primary_comparison": "cptac_tumor_vs_normal",
        "candidate_source": "demo",
        "n_candidates": n,
        "weights": weights,
        "signals_used": signals,
        "providers": {s: {"available": True, "version": "DEMO"} for s in signals[4:]},
        "note": "CONJUNTO DE DEMOSTRACIÓN SINTÉTICO — no son resultados reales.",
    }

    return Dataset(
        hypotheses=hyp,
        evidence_matrix=evm,
        differentials={"cptac_tumor_vs_normal": diff},
        provenance=prov,
        processed_matrix=proc,
        is_demo=True,
        results_dir=results_dir,
        notes=["No se encontró hypothesis_priority.tsv: mostrando datos de DEMOSTRACIÓN sintéticos."],
    )


def load_dataset(results_dir: Path | None = None) -> Dataset:
    rd = Path(results_dir) if results_dir else default_results_dir()
    real = _load_real(rd)
    if real is not None:
        return real
    return _demo_dataset(rd)
