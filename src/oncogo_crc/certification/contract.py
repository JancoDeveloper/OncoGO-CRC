"""Contrato estadístico y certificado R calculable con el cruce CPTAC ∩ PRIDE.

Los cinco subíndices del diseño 3.0 no se colapsan en un número. Con los
artefactos actuales solo R1 (dirección) y R2 (detectabilidad entre las dos
cohortes ya cruzadas) tienen datos. R3–R5 quedan explícitamente no calculados.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import binomtest

# El association score de Open Targets ya mezcla literatura, expresión
# diferencial y genética. Usarlo junto al diferencial propio cuenta dos veces
# la misma evidencia. Tractabilidad sí está permitida: es otro canal.
FORBIDDEN_EVIDENCE_COLUMNS = frozenset(
    {
        "opentargets_association",
        "association_score",
        "overall_association",
        "overallassociationscore",
    }
)

# Umbrales fijados en el diseño antes del análisis. No se ajustan para
# «mejorar» un resultado.
H1_PROTEIN_CONCORDANCE_MAX = 0.70
H1_PATHWAY_CONCORDANCE_MIN = 0.80
H4_MEASURABLE_FRACTION_REDUNDANT = 0.90
TREAT_LFC = 0.585  # log2(1.5); el diferencial actual aún no aplica treat()


@dataclass
class Subindex:
    code: str
    name: str
    value: float | None
    status: str  # calculado | no_calculable
    definition: str
    note: str = ""


@dataclass
class CertificationReport:
    n_proteins_both: int = 0
    n_discovery_significant: int = 0
    n_sign_evaluable: int = 0
    n_same_sign: int = 0
    sign_concordance: float | None = None
    sign_ci_low: float | None = None
    sign_ci_high: float | None = None
    sign_p_vs_half: float | None = None
    measurable_fraction: float | None = None
    measurable_fraction_log2fc: float | None = None
    h1_protein_clause: str = "no_evaluable"
    h1_pathway_clause: str = "no_evaluable"
    h1_status: str = "no_evaluable"
    pathway_concordance: float | None = None
    pathway_n: int = 0
    pathway_ci_low: float | None = None
    pathway_ci_high: float | None = None
    inference_method: str = ""
    h4_proxy: str = "no_evaluable"
    subindices: list[Subindex] = field(default_factory=list)
    gaps: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        return payload


def _finite_pair(a: pd.Series, b: pd.Series) -> pd.DataFrame:
    frame = pd.DataFrame({"a": pd.to_numeric(a, errors="coerce"), "b": pd.to_numeric(b, errors="coerce")})
    ok = np.isfinite(frame["a"]) & np.isfinite(frame["b"])
    # Un log2FC de exactamente 0 no tiene signo interpretable.
    ok &= (frame["a"] != 0) & (frame["b"] != 0)
    return frame.loc[ok]


def sign_concordance_discovery(
    overlap: pd.DataFrame,
    *,
    alpha: float = 0.05,
    discovery_padj: str = "padj_cptac",
    discovery_lfc: str = "log2fc_cptac",
    replica_lfc: str = "log2fc_pride",
) -> dict[str, Any]:
    """Concordancia de signo en proteínas significativas en la cohorte de descubrimiento.

    La cohorte de descubrimiento es CPTAC (PDC000116). La réplica es PRIDE.
    El intervalo es el exacto de Clopper–Pearson. El test compara contra p = 0,5
    (azar de signo).
    """
    empty = {
        "n_discovery_significant": 0,
        "n_sign_evaluable": 0,
        "n_same_sign": 0,
        "concordance": None,
        "ci_low": None,
        "ci_high": None,
        "p_vs_half": None,
    }
    needed = {discovery_padj, discovery_lfc, replica_lfc}
    if overlap is None or overlap.empty or not needed.issubset(overlap.columns):
        return empty

    disc = overlap.copy()
    if "_merge" in disc.columns:
        disc = disc.loc[disc["_merge"].astype(str) == "both"]
    padj = pd.to_numeric(disc[discovery_padj], errors="coerce")
    sig = disc.loc[padj < alpha]
    pairs = _finite_pair(sig[discovery_lfc], sig[replica_lfc])
    n_sig = int((padj < alpha).sum())
    n = int(len(pairs))
    if n == 0:
        empty["n_discovery_significant"] = n_sig
        return empty
    same = int((np.sign(pairs["a"]) == np.sign(pairs["b"])).sum())
    test = binomtest(same, n, p=0.5, alternative="two-sided")
    ci = test.proportion_ci(confidence_level=0.95, method="exact")
    return {
        "n_discovery_significant": n_sig,
        "n_sign_evaluable": n,
        "n_same_sign": same,
        "concordance": same / n,
        "ci_low": float(ci.low),
        "ci_high": float(ci.high),
        "p_vs_half": float(test.pvalue),
    }


def measurable_fraction_two_cohorts(
    overlap: pd.DataFrame,
    *,
    alpha: float = 0.05,
    discovery_padj: str = "padj_cptac",
    replica_lfc: str = "log2fc_pride",
) -> dict[str, Any]:
    """Proxy de medibilidad con el cruce outer de las dos cohortes.

    El denominador son todas las proteínas con q < alpha en descubrimiento,
    incluidas las que el cruce marca como ``left_only``. Cuentan como medibles
    solo las que además tienen log2FC finito en la réplica. Restringir el
    denominador a ``_merge == both`` daría fracción 1 por construcción.
    No sustituye el criterio de péptido proteotípico en al menos dos plataformas.
    """
    out = {"n_discovery_significant": 0, "n_measurable": 0, "fraction": None}
    if overlap is None or overlap.empty or discovery_padj not in overlap.columns:
        return out
    disc = overlap
    padj = pd.to_numeric(disc[discovery_padj], errors="coerce")
    sig = disc.loc[padj < alpha]
    n_sig = int(len(sig))
    out["n_discovery_significant"] = n_sig
    if n_sig == 0 or replica_lfc not in sig.columns:
        return out
    lfc = pd.to_numeric(sig[replica_lfc], errors="coerce")
    n_ok = int(np.isfinite(lfc).sum())
    out["n_measurable"] = n_ok
    out["fraction"] = n_ok / n_sig
    return out


def h1_protein_clause(concordance: float | None, *, threshold: float = H1_PROTEIN_CONCORDANCE_MAX) -> str:
    """H1 predice concordancia de proteína por debajo del umbral.

    ``sostiene_clausula_proteina`` no cierra H1: falta la concordancia de vía.
    """
    if concordance is None:
        return "no_evaluable"
    if concordance < threshold:
        return "sostiene_clausula_proteina"
    return "rechaza_clausula_proteina"


def h1_pathway_clause(concordance: float | None, *, threshold: float = H1_PATHWAY_CONCORDANCE_MIN) -> str:
    """H1 predice concordancia de vía por encima del umbral. El umbral no se incluye."""
    if concordance is None:
        return "no_evaluable"
    if concordance > threshold:
        return "sostiene_clausula_via"
    return "rechaza_clausula_via"


def h1_status(protein_clause: str, pathway_clause: str) -> str:
    """H1 es la conjunción. Se rechaza si falla cualquiera de las dos cláusulas."""
    if protein_clause == "rechaza_clausula_proteina" or pathway_clause == "rechaza_clausula_via":
        return "rechaza_h1"
    if protein_clause == "sostiene_clausula_proteina" and pathway_clause == "sostiene_clausula_via":
        return "sostiene_h1"
    if protein_clause == "sostiene_clausula_proteina":
        return "parcial_solo_proteina"
    return "no_evaluable"


def h4_proxy_verdict(
    fraction: float | None, *, redundant_at: float = H4_MEASURABLE_FRACTION_REDUNDANT
) -> str:
    """H4 se falsifica si al menos el 90 % de los candidatos naive ya es medible."""
    if fraction is None:
        return "no_evaluable"
    if fraction >= redundant_at:
        return "filtro_redundante_en_el_proxy"
    return "filtro_no_redundante_en_el_proxy"


def default_gaps() -> list[str]:
    return [
        "El diferencial sigue siendo un t de Welch o un t pareado sobre la matriz imputada "
        "(log2norm). El diseño 3.0 exige limma/treat (lfc = 0,585) dentro de cohorte y no "
        "imputar ausencias para inferencia. La matriz sin imputar se guarda ahora como "
        "*_log2_observed.tsv; el diferencial todavía no la lee.",
        "No hay meta-análisis de efectos aleatorios (REML). Fusionar CPTAC y PRIDE en un "
        "solo modelo, o corregir lote con ComBat para inferir, sigue prohibido.",
        "La pureza tumoral no entra como covariable. Sin ella, parte de la señal puede ser estroma.",
        "R3 (homogeneidad I²), R4 (robustez a pureza) y R5 (permutaciones de extremo a extremo) "
        "no se calculan: faltan el meta-análisis, el ajuste por pureza y las permutaciones.",
        "La concordancia de vía (umbral 0,80 de H1) no se certifica en esta etapa.",
        "La ventana de medibilidad de péptido proteotípico en al menos dos plataformas no está "
        "implementada. El proxy de dos cohortes no la reemplaza.",
        "El diseño factorial parcial (TMT, LFQ fresco, LFQ FFPE, iTRAQ) no está congelado en "
        "config/datasets.yaml. El cruce actual es CPTAC frente a un PRIDE agregado.",
        "H0 nivel A (checksum entre PXD019504 y PXD031705) está documentado como hallazgo "
        "exploratorio. Hay que repetirlo con TLS verificado antes de escribirlo como resultado "
        "de tesis. El script scripts/pride_independence_audit.py hace esa repetición.",
    ]


def build_certification_report(
    overlap: pd.DataFrame | None,
    *,
    alpha: float = 0.05,
    gaps: list[str] | None = None,
    pathway_concordance: float | None = None,
    pathway_n: int = 0,
    pathway_ci_low: float | None = None,
    pathway_ci_high: float | None = None,
    inference_method: str = "",
    measurable_override: float | None = None,
    r2_note: str | None = None,
) -> CertificationReport:
    sign = sign_concordance_discovery(overlap if overlap is not None else pd.DataFrame(), alpha=alpha)
    meas = measurable_fraction_two_cohorts(overlap if overlap is not None else pd.DataFrame(), alpha=alpha)
    conc = sign["concordance"]
    frac_log2 = meas["fraction"]
    frac = frac_log2 if measurable_override is None else measurable_override
    n_both = 0
    if overlap is not None and not overlap.empty and "_merge" in overlap.columns:
        n_both = int((overlap["_merge"].astype(str) == "both").sum())
    elif overlap is not None:
        n_both = int(len(overlap))

    r1 = Subindex(
        code="R1_direccion",
        name="Concordancia de signo",
        value=None if conc is None else round(float(conc), 4),
        status="calculado" if conc is not None else "no_calculable",
        definition=(
            "Fracción de proteínas con q < alpha en CPTAC cuyo log2FC tiene el mismo signo en PRIDE. "
            "Intervalo de Clopper–Pearson."
        ),
        note=(
            "Cláusula de proteína de H1. La de vía se informa aparte y no se pliega en este número."
        ),
    )
    r2_value = None if frac is None else round(float(frac), 4)
    r2 = Subindex(
        code="R2_detectabilidad",
        name="Detectabilidad entre las dos cohortes cruzadas",
        value=r2_value,
        status="calculado" if frac is not None else "no_calculable",
        definition=(
            "Fracción de proteínas significativas en CPTAC que también se miden en la réplica."
        ),
        note=r2_note
        or "Proxy por log2FC finito. No es el criterio de péptido proteotípico en ≥2 plataformas.",
    )
    pending = [
        Subindex(
            "R3_homogeneidad",
            "1 − I² del meta-análisis",
            None,
            "no_calculable",
            "Requiere meta-análisis de efectos aleatorios por proteína.",
        ),
        Subindex(
            "R4_pureza",
            "Fracción del efecto que sobrevive al ajuste por pureza",
            None,
            "no_calculable",
            "La pureza tumoral aún no es covariable del diferencial.",
        ),
        Subindex(
            "R5_nulo",
            "Resistencia a permutaciones de etiqueta",
            None,
            "no_calculable",
            "Requieren ≥100 permutaciones de extremo a extremo, todavía no ejecutadas.",
        ),
    ]
    return CertificationReport(
        n_proteins_both=n_both,
        n_discovery_significant=int(sign["n_discovery_significant"]),
        n_sign_evaluable=int(sign["n_sign_evaluable"]),
        n_same_sign=int(sign["n_same_sign"]),
        sign_concordance=None if conc is None else float(conc),
        sign_ci_low=sign["ci_low"],
        sign_ci_high=sign["ci_high"],
        sign_p_vs_half=sign["p_vs_half"],
        measurable_fraction=None if frac is None else float(frac),
        measurable_fraction_log2fc=None if frac_log2 is None else float(frac_log2),
        h1_protein_clause=h1_protein_clause(None if conc is None else float(conc)),
        h1_pathway_clause=h1_pathway_clause(pathway_concordance),
        h1_status=h1_status(
            h1_protein_clause(None if conc is None else float(conc)),
            h1_pathway_clause(pathway_concordance),
        ),
        pathway_concordance=pathway_concordance,
        pathway_n=int(pathway_n),
        pathway_ci_low=pathway_ci_low,
        pathway_ci_high=pathway_ci_high,
        inference_method=inference_method,
        h4_proxy=h4_proxy_verdict(None if frac is None else float(frac)),
        subindices=[r1, r2, *pending],
        gaps=default_gaps() if gaps is None else gaps,
    )


def drop_forbidden_evidence(frame: pd.DataFrame) -> list[str]:
    """Quita columnas de association score si un snapshot las arrastró al candidato."""
    dropped = [c for c in frame.columns if c.lower() in FORBIDDEN_EVIDENCE_COLUMNS or c in FORBIDDEN_EVIDENCE_COLUMNS]
    if dropped:
        frame.drop(columns=dropped, inplace=True)
    return dropped
