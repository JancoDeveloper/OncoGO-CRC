"""Capa C — priorización de candidatos por **evidencia convergente** (no-promisorio).

Toma la firma diferencial primaria (y opcionalmente la firma reproducible CPTAC∩PRIDE), la anota con
evidencia pública objetiva (esencialidad DepMap, tractabilidad Open Targets, druggability DGIdb,
especificidad tumoral HPA, seguridad GTEx, contexto genómico cBioPortal) y calcula un **score de
convergencia transparente** (media ponderada de señales normalizadas disponibles).

Salidas en ``results/prioritization/``:
- ``hypothesis_priority.tsv``            — tabla ancha priorizada (una fila por gen candidato).
- ``convergence_evidence_matrix.tsv``   — formato largo (gen × señal): valor crudo, favorabilidad, peso.
- ``prioritization_provenance.json``    — pesos, versiones/rutas de cada fuente, parámetros y timestamp.

Diseño no-promisorio: la evidencia **ausente no penaliza** (se excluye del promedio); se reporta
``evidence_count`` y la procedencia de cada fuente. El score prioriza **hipótesis**, no dianas validadas.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from oncogo_crc.certification.contract import drop_forbidden_evidence
from oncogo_crc.settings import Settings
from oncogo_crc.utils import run_manifest
from oncogo_crc.utils.evidence_providers import (
    available_provider_names,
    load_provider,
    minmax_favorability,
)
from oncogo_crc.utils.io import read_tsv, write_tsv

# Señales internas (siempre calculables desde artefactos del pipeline) y su favorabilidad.
# Las señales de proveedores externos se añaden dinámicamente según snapshots disponibles.
_DEFAULT_WEIGHTS: dict[str, float] = {
    "diff_effect": 1.0,
    "diff_significance": 1.0,
    "reproducible": 1.5,
    "network_hub": 0.5,
    "depmap_essentiality": 1.5,
    "opentargets_tractability": 1.5,
    "dgidb_druggability": 1.0,
    "hpa_tumor_specificity": 1.0,
    "gtex_safety": 1.0,
    "cbioportal_genomic": 0.75,
}


def _evidence_dir(settings: Settings, rel: str) -> Path:
    """Resuelve el directorio de snapshots en workspace y, si no existe, en el repo."""
    p_ws = settings.resolve(rel)
    if p_ws.exists():
        return p_ws
    p_repo = settings.root / rel
    return p_repo if p_repo.exists() else p_ws


def _load_reproducible_sets(overlap_path: Path, alpha: float) -> tuple[set[str], set[str]]:
    """Devuelve (genes_reproducibles, genes_reproducibles_mismo_signo) desde el overlap estricto."""
    if not overlap_path.exists():
        return set(), set()
    df = read_tsv(overlap_path)
    gcol = "protein_id_cptac" if "protein_id_cptac" in df.columns else None
    if gcol is None or "padj_cptac" not in df.columns or "padj_pride" not in df.columns:
        return set(), set()
    df = df.copy()
    df["_g"] = df[gcol].astype(str).str.strip().str.upper()
    both_sig = (pd.to_numeric(df["padj_cptac"], errors="coerce") < alpha) & (
        pd.to_numeric(df["padj_pride"], errors="coerce") < alpha
    )
    repro = set(df.loc[both_sig, "_g"].tolist())
    same_sign: set[str] = set()
    if "log2fc_cptac" in df.columns and "log2fc_pride" in df.columns:
        sign_ok = np.sign(pd.to_numeric(df["log2fc_cptac"], errors="coerce")) == np.sign(
            pd.to_numeric(df["log2fc_pride"], errors="coerce")
        )
        same_sign = set(df.loc[both_sig & sign_ok, "_g"].tolist())
    return repro, same_sign


def _load_hub_degree(hub_path: Path) -> dict[str, float]:
    if not hub_path.exists():
        return {}
    df = read_tsv(hub_path)
    if "protein_id" not in df.columns or "degree" not in df.columns:
        return {}
    df = df.copy()
    df["_g"] = df["protein_id"].astype(str).str.strip().str.upper()
    return dict(zip(df["_g"], pd.to_numeric(df["degree"], errors="coerce")))


def _tier_from_score(scores: pd.Series) -> pd.Series:
    """Cuartiles de score → tier interpretativo A/B/C/D (A = mayor convergencia)."""
    valid = scores.dropna()
    if valid.empty:
        return pd.Series(["tier_D"] * len(scores), index=scores.index)
    try:
        q = pd.qcut(scores.rank(method="first"), 4, labels=["tier_D", "tier_C", "tier_B", "tier_A"])
        return q.astype(str)
    except ValueError:
        return pd.Series(["tier_C"] * len(scores), index=scores.index)


def run(settings: Settings) -> None:
    cfg = settings.config.get("prioritization", {}) or {}
    if not bool(cfg.get("enabled", True)):
        print("[prioritization] deshabilitado (prioritization.enabled=false).")
        run_manifest.record_stage_meta(settings, "prioritization", executed=False, reason="disabled")
        return

    res = settings.resolve(settings.config["paths"]["results_dir"])
    out_dir = res / "prioritization"
    out_dir.mkdir(parents=True, exist_ok=True)

    comp = str(cfg.get("primary_comparison", "cptac_tumor_vs_normal"))
    alpha = float(cfg.get("alpha", settings.config.get("differential", {}).get("alpha", 0.05)))
    diff_path = res / "differential" / f"{comp}_diff.tsv"
    if not diff_path.exists():
        print(f"[prioritization] omitido (falta {diff_path.name}); ejecuta antes `differential`.")
        run_manifest.record_stage_meta(
            settings, "prioritization", executed=False, reason="missing_differential"
        )
        # Limpia salidas viejas para no publicar hipótesis obsoletas.
        for name in ("hypothesis_priority.tsv", "convergence_evidence_matrix.tsv"):
            (out_dir / name).unlink(missing_ok=True)
        return

    diff = read_tsv(diff_path)
    if "protein_id" not in diff.columns:
        print("[prioritization] diff sin columna protein_id; omitiendo.")
        run_manifest.record_stage_meta(settings, "prioritization", executed=False, reason="bad_diff")
        return

    cand = diff.copy()
    cand["gene_symbol"] = cand["protein_id"].astype(str).str.strip().str.upper()
    cand = cand[cand["gene_symbol"].notna() & (cand["gene_symbol"] != "") & (cand["gene_symbol"] != "NAN")]
    cand = cand.drop_duplicates("gene_symbol")

    # --- Selección de candidatos ---
    source = str(cfg.get("candidate_source", "significant"))
    max_candidates = int(cfg.get("max_candidates", 500) or 0)
    if source == "significant" and "significant" in cand.columns and bool(cand["significant"].any()):
        cand = cand[cand["significant"].astype(bool)]
    elif source == "all":
        pass
    else:
        # Fallback: sin significativas, prioriza por |log2fc|.
        if "log2fc" in cand.columns:
            cand = cand.reindex(cand["log2fc"].abs().sort_values(ascending=False).index)
    if max_candidates > 0 and len(cand) > max_candidates:
        if "log2fc" in cand.columns:
            cand = cand.reindex(cand["log2fc"].abs().sort_values(ascending=False).index).head(max_candidates)
        else:
            cand = cand.head(max_candidates)
    cand = cand.reset_index(drop=True)
    n_candidates = int(len(cand))
    if n_candidates == 0:
        print("[prioritization] 0 candidatos tras el filtro; omitiendo.")
        run_manifest.record_stage_meta(settings, "prioritization", executed=False, reason="no_candidates")
        return

    # --- Anotaciones internas ---
    repro_set, same_sign_set = _load_reproducible_sets(
        res / "reproducibility" / "cptac_pride_strict_overlap.tsv", alpha
    )
    have_repro = bool(repro_set) or (res / "reproducibility" / "cptac_pride_strict_overlap.tsv").exists()
    if have_repro:
        cand["reproducible"] = cand["gene_symbol"].isin(repro_set)
        cand["reproducible_same_sign"] = cand["gene_symbol"].isin(same_sign_set)

    primary_ds = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")
    hub_degree = _load_hub_degree(res / "networks" / f"{primary_ds}_hub_proteins.tsv")
    if hub_degree:
        cand["network_degree"] = cand["gene_symbol"].map(hub_degree)

    # --- Proveedores externos (snapshots locales pluggables) ---
    evidence_dir = _evidence_dir(settings, str(cfg.get("evidence_dir", "metadata/evidence")))
    providers_meta: dict[str, Any] = {}
    for name in available_provider_names():
        pr = load_provider(name, evidence_dir)
        providers_meta[name] = {
            "available": pr.available,
            "status": pr.status,
            "version": pr.version,
            "path": pr.path,
            "n_genes": pr.n_genes,
            "columns": pr.columns,
        }
        if pr.available:
            cand = cand.merge(pr.df, left_on="gene_symbol", right_index=True, how="left")

    dropped_forbidden = drop_forbidden_evidence(cand)

    # --- Cálculo de favorabilidades (0–1); NaN = evidencia ausente (no cuenta) ---
    favor = pd.DataFrame(index=cand.index)
    raw_values: dict[str, pd.Series] = {}

    if "log2fc" in cand.columns:
        raw_values["diff_effect"] = cand["log2fc"].abs()
        favor["diff_effect"] = minmax_favorability(cand["log2fc"].abs())
    if "padj" in cand.columns:
        neglog = -np.log10(pd.to_numeric(cand["padj"], errors="coerce").clip(lower=1e-300))
        raw_values["diff_significance"] = neglog
        favor["diff_significance"] = minmax_favorability(neglog)
    if "reproducible" in cand.columns:
        val = cand["reproducible"].astype(float) + 0.5 * cand.get(
            "reproducible_same_sign", pd.Series(False, index=cand.index)
        ).astype(float)
        raw_values["reproducible"] = val
        favor["reproducible"] = (val / val.max()) if val.max() > 0 else val
    if "network_degree" in cand.columns:
        raw_values["network_hub"] = cand["network_degree"]
        favor["network_hub"] = minmax_favorability(cand["network_degree"])
    if "depmap_chronos" in cand.columns:
        raw_values["depmap_essentiality"] = cand["depmap_chronos"]
        # Chronos más negativo = más esencial = más favorable → invertir.
        favor["depmap_essentiality"] = minmax_favorability(cand["depmap_chronos"], invert=True)
    if "opentargets_tractability" in cand.columns:
        v = pd.to_numeric(cand["opentargets_tractability"], errors="coerce").clip(0, 1)
        raw_values["opentargets_tractability"] = v
        favor["opentargets_tractability"] = v
    if "dgidb_n_drugs" in cand.columns:
        v = np.log1p(pd.to_numeric(cand["dgidb_n_drugs"], errors="coerce"))
        raw_values["dgidb_druggability"] = cand["dgidb_n_drugs"]
        favor["dgidb_druggability"] = minmax_favorability(v)
    if "hpa_tumor_specificity" in cand.columns:
        raw_values["hpa_tumor_specificity"] = cand["hpa_tumor_specificity"]
        favor["hpa_tumor_specificity"] = minmax_favorability(cand["hpa_tumor_specificity"])
    if "gtex_max_normal_tpm" in cand.columns:
        raw_values["gtex_safety"] = cand["gtex_max_normal_tpm"]
        favor["gtex_safety"] = minmax_favorability(cand["gtex_max_normal_tpm"], invert=True)
    if "cbioportal_alteration_frequency" in cand.columns:
        raw_values["cbioportal_genomic"] = cand["cbioportal_alteration_frequency"]
        favor["cbioportal_genomic"] = minmax_favorability(cand["cbioportal_alteration_frequency"])

    # --- Score de convergencia = media ponderada de señales disponibles ---
    weights = {**_DEFAULT_WEIGHTS, **(cfg.get("weights") or {})}
    signals = [s for s in favor.columns]
    w = np.array([float(weights.get(s, 0.0)) for s in signals], dtype=float)
    F = favor[signals].to_numpy(dtype=float)  # (n_cand, n_signals), con NaN
    mask = ~np.isnan(F)
    W = np.where(mask, w[None, :], 0.0)
    num = np.nansum(np.where(mask, F * w[None, :], 0.0), axis=1)
    den = W.sum(axis=1)
    score = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
    cand["convergence_score"] = score
    cand["evidence_count"] = mask.sum(axis=1).astype(int)
    n_external = [s for s in signals if s not in ("diff_effect", "diff_significance")]
    cand["external_evidence_count"] = (
        (~np.isnan(favor[n_external].to_numpy(dtype=float))).sum(axis=1).astype(int)
        if n_external
        else 0
    )
    cand["evidence_tier"] = _tier_from_score(cand["convergence_score"]).values
    cand = cand.sort_values("convergence_score", ascending=False, na_position="last").reset_index(drop=True)
    cand.insert(0, "rank", np.arange(1, len(cand) + 1))

    # --- Guardar tabla ancha (hipótesis priorizadas) ---
    front = ["rank", "gene_symbol", "convergence_score", "evidence_tier", "evidence_count", "external_evidence_count"]
    ordered = front + [c for c in cand.columns if c not in front]
    hyp = cand[ordered]
    write_tsv(hyp, out_dir / "hypothesis_priority.tsv")

    # --- Guardar matriz de evidencia larga (transparencia por señal) ---
    long_rows: list[dict[str, Any]] = []
    for i, gene in enumerate(cand["gene_symbol"].tolist()):
        for s in signals:
            long_rows.append(
                {
                    "gene_symbol": gene,
                    "signal": s,
                    "raw_value": (float(raw_values[s].iloc[i]) if s in raw_values and pd.notna(raw_values[s].iloc[i]) else np.nan),
                    "favorability": (float(favor[s].iloc[i]) if pd.notna(favor[s].iloc[i]) else np.nan),
                    "weight": float(weights.get(s, 0.0)),
                    "available": bool(pd.notna(favor[s].iloc[i])),
                }
            )
    write_tsv(pd.DataFrame(long_rows), out_dir / "convergence_evidence_matrix.tsv")

    # --- Procedencia (versiones de fuentes + parámetros) ---
    provenance = {
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "primary_comparison": comp,
        "candidate_source": source,
        "n_candidates": n_candidates,
        "alpha": alpha,
        "weights": weights,
        "signals_used": signals,
        "reproducibility_available": have_repro,
        "network_hub_available": bool(hub_degree),
        "evidence_dir": str(evidence_dir),
        "providers": providers_meta,
        "forbidden_evidence_dropped": dropped_forbidden,
        "note": (
            "Score de convergencia = media ponderada de señales disponibles. "
            "Evidencia ausente no penaliza. El association score de Open Targets está prohibido. "
            "Hipótesis exploratorias, no dianas validadas."
        ),
    }
    (out_dir / "prioritization_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")

    n_ext_available = sum(1 for m in providers_meta.values() if m["available"])
    run_manifest.record_stage_meta(
        settings,
        "prioritization",
        executed=True,
        primary_comparison=comp,
        n_candidates=n_candidates,
        n_signals=len(signals),
        n_external_providers_available=n_ext_available,
        outputs=["hypothesis_priority.tsv", "convergence_evidence_matrix.tsv", "prioritization_provenance.json"],
    )
    print(
        f"[prioritization] {n_candidates} candidatos priorizados con {len(signals)} señales "
        f"({n_ext_available}/{len(providers_meta)} proveedores externos disponibles) "
        f"-> {settings.relative_to_project(out_dir / 'hypothesis_priority.tsv')}"
    )
