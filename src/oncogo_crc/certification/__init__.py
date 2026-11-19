"""Contrato estadístico de la corrida certificada v1.

La detección de colisiones de ficheros de versiones anteriores no forma parte
de la corrida v3.7 y no se importa aquí.
"""

from oncogo_crc.certification.contract import (
    FORBIDDEN_EVIDENCE_COLUMNS,
    CertificationReport,
    build_certification_report,
)

__all__ = [
    "FORBIDDEN_EVIDENCE_COLUMNS",
    "CertificationReport",
    "build_certification_report",
]
