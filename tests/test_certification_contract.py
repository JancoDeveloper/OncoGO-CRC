"""Certificado R, canales prohibidos y colisión de ficheros."""

from __future__ import annotations

import pandas as pd

from oncogo_crc.certification.contract import (
    FORBIDDEN_EVIDENCE_COLUMNS,
    build_certification_report,
    drop_forbidden_evidence,
    h1_protein_clause,
    h4_proxy_verdict,
)
from oncogo_crc.certification.independence import FileRecord, detect_file_collisions
from oncogo_crc.utils.evidence_providers import load_provider


def _overlap() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "_merge": ["both", "both", "both", "both"],
            "padj_cptac": [0.01, 0.01, 0.01, 0.20],
            "log2fc_cptac": [1.2, -0.8, 0.5, 2.0],
            "log2fc_pride": [0.4, 0.3, -0.2, -1.0],
        }
    )


def test_sign_concordance_uses_only_discovery_hits():
    report = build_certification_report(_overlap(), alpha=0.05)
    # Tres significativas; dos mismo signo (1.2/0.4 y la tercera -0.2 no; -0.8/0.3 no).
    # 1.2 y 0.4 mismo; -0.8 y 0.3 distinto; 0.5 y -0.2 distinto. Concordancia 1/3.
    assert report.n_sign_evaluable == 3
    assert report.n_same_sign == 1
    assert report.sign_concordance == 1 / 3
    assert report.h1_protein_clause == "sostiene_clausula_proteina"
    assert report.sign_ci_low is not None and report.sign_ci_high is not None
    codes = {s.code: s.status for s in report.subindices}
    assert codes["R1_direccion"] == "calculado"
    assert codes["R3_homogeneidad"] == "no_calculable"
    assert codes["R5_nulo"] == "no_calculable"


def test_h1_rejects_protein_clause_at_threshold():
    assert h1_protein_clause(0.70) == "rechaza_clausula_proteina"
    assert h1_protein_clause(0.69) == "sostiene_clausula_proteina"
    assert h4_proxy_verdict(0.90) == "filtro_redundante_en_el_proxy"
    assert h4_proxy_verdict(0.40) == "filtro_no_redundante_en_el_proxy"


def test_association_score_is_forbidden_and_not_loaded(tmp_path):
    snap = tmp_path / "opentargets_tractability.tsv"
    snap.write_text(
        "gene_symbol\ttractability\tassociation_score\nTP53\t0.8\t0.99\n",
        encoding="utf-8",
    )
    loaded = load_provider("opentargets", tmp_path)
    assert loaded.available
    assert "opentargets_tractability" in loaded.columns
    assert "opentargets_association" not in loaded.columns
    frame = loaded.df.reset_index()
    frame["association_score"] = 0.99
    dropped = drop_forbidden_evidence(frame)
    assert "association_score" in dropped
    assert "association_score" not in frame.columns
    assert "opentargets_association" in FORBIDDEN_EVIDENCE_COLUMNS


def test_checksum_collision_is_literal_redeposit():
    shared = "abc123"
    records = [
        FileRecord("PXD019504", "1417AT_S1_12R1_CRC6T_i39_MichaelRoehrl.raw", 100, shared),
        FileRecord("PXD031705", "1417AT_S1_12R1_CRC6T_i39_MichaelRoehrl.raw", 100, shared),
        FileRecord("PXD034575", "CRC_6N.raw", 50, "other"),
        FileRecord("PXD019504", "proteinGroups.txt", 10, "gg"),
    ]
    hits = detect_file_collisions(records)
    checksums = [h for h in hits if h.kind == "checksum"]
    assert len(checksums) == 1
    assert checksums[0].severity == "prueba"
    assert {checksums[0].accession_a, checksums[0].accession_b} == {"PXD019504", "PXD031705"}
    assert all("proteinGroups" not in h.file_a for h in hits)
