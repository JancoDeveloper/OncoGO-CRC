"""Capa de certificación de reproducibilidad (OncoGO-CRC 3.0).

No sustituye el pipeline de descubrimiento. Añade el contrato estadístico,
el índice R calculable con los artefactos que ya existen, y el detector de
colisión de ficheros (nivel A de la auditoría de independencia).
"""

from oncogo_crc.certification.contract import (
    FORBIDDEN_EVIDENCE_COLUMNS,
    CertificationReport,
    build_certification_report,
)
from oncogo_crc.certification.independence import Collision, FileRecord, detect_file_collisions

__all__ = [
    "FORBIDDEN_EVIDENCE_COLUMNS",
    "CertificationReport",
    "Collision",
    "FileRecord",
    "build_certification_report",
    "detect_file_collisions",
]
