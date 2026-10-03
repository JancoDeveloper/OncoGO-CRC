"""Regresión de la corrida certificada v1 (1 de octubre de 2026).

Dos niveles.

El primero lee los artefactos congelados y comprueba que los números del
certificado siguen siendo los mismos. Protege contra sobrescrituras: si alguien
corre una etapa v2 apuntando a ``results/certified``, este test falla.

El segundo vuelve a ejecutar el contraste en un workspace temporal, con las
entradas enlazadas en solo lectura, y comprueba que salen los mismos números.
Eso es lo que demuestra que el resultado es determinista. Es lento, así que
corre solo con ``ONCOGO_REGRESION_V1=1``.

Nunca se escribe en ``results/certified`` de v1.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

WORKSPACE = Path("/Volumes/JANCODATA/OncoGO-CRC")

# Números de la corrida certificada v1. No se tocan.
V1 = {
    "cptac_n_significant": 688,
    "pride_n_significant": 140,
    "n_both_significant_padj_lt_0.05": 28,
    "n_same_sign": 104,
    "n_sign_evaluable": 119,
    "peptide_n_measurable": 226,
    "n_certified_surviving": 307,
}

requires_v1 = pytest.mark.skipif(
    not (WORKSPACE / "results/certified/inference_audit.json").is_file(),
    reason="workspace v1 no montado",
)


def _audit(base: Path) -> dict:
    return json.loads((base / "results/certified/inference_audit.json").read_text(encoding="utf-8"))


def _cert(base: Path) -> dict:
    return json.loads(
        (base / "results/certification/certificado_resumen.json").read_text(encoding="utf-8")
    )


def _summary(audit: dict) -> dict[str, float]:
    return {row["metric"]: row["value"] for row in audit.get("summary_metrics", [])}


def _observed(base: Path) -> dict[str, int]:
    audit = _audit(base)
    cert = _cert(base)
    summary = _summary(audit)
    return {
        "cptac_n_significant": int(audit["cptac_n_significant"]),
        "pride_n_significant": int(audit["pride_n_significant"]),
        "n_both_significant_padj_lt_0.05": int(summary["n_both_significant_padj_lt_0.05"]),
        "n_same_sign": int(cert["n_same_sign"]),
        "n_sign_evaluable": int(cert["n_sign_evaluable"]),
        "peptide_n_measurable": int(audit["peptide_n_measurable"]),
        "n_certified_surviving": int(audit["r4_stroma"]["n_certified_surviving"]),
    }


@requires_v1
def test_v1_congelada_mantiene_sus_numeros():
    assert _observed(WORKSPACE) == V1


@requires_v1
def test_v1_fracciones_son_numerador_sobre_denominador():
    audit = _audit(WORKSPACE)
    cert = _cert(WORKSPACE)
    assert cert["sign_concordance"] == pytest.approx(
        cert["n_same_sign"] / cert["n_sign_evaluable"], abs=1e-9
    )
    assert audit["peptide_measurable_fraction"] == pytest.approx(
        audit["peptide_n_measurable"] / audit["peptide_n_discovery_significant"], abs=1e-9
    )
    stroma = audit["r4_stroma"]
    assert stroma["fraction_surviving"] == pytest.approx(
        stroma["n_certified_surviving"] / stroma["n_certified"], abs=1e-9
    )


@requires_v1
def test_manifiesto_v1_existe_y_coincide():
    import subprocess
    import sys

    repo = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [sys.executable, str(repo / "scripts/v2/manifest_hashes.py"), "verify"],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


@requires_v1
@pytest.mark.skipif(
    os.environ.get("ONCOGO_REGRESION_V1") != "1",
    reason="re-ejecución completa: activar con ONCOGO_REGRESION_V1=1",
)
def test_reejecutar_v1_reproduce_los_numeros(tmp_path: Path):
    """Vuelve a correr el contraste en un workspace temporal y compara."""
    from oncogo_crc.settings import load_settings
    from oncogo_crc.stages import certification, certified_inference

    # Entradas enlazadas: se leen, no se copian ni se editan.
    for sub in ("data", "metadata"):
        (tmp_path / sub).symlink_to(WORKSPACE / sub)
    (tmp_path / "results").mkdir()

    settings = load_settings(data_workspace=tmp_path)
    assert settings.data_workspace == tmp_path.resolve()
    certified_inference.run(settings)
    certification.run(settings)

    assert _observed(tmp_path) == V1
