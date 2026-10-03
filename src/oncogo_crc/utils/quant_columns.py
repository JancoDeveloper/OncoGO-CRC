"""Selección de columnas de cuantificación sin mezclar mediciones distintas.

La matriz PDC de colon trae, para el mismo espécimen, ``Log Ratio`` (todos los
péptidos) y ``Unshared Log Ratio`` (péptidos no compartidos con otra proteína).
Tratarlas como dos muestras duplica el n y mezcla dos cuantificaciones.
"""

from __future__ import annotations

SUMMARY_PROTEIN_IDS = frozenset({"MEAN", "MEDIAN", "STDDEV", "STD"})


def is_summary_protein(protein_id: object) -> bool:
    return str(protein_id).strip().upper() in SUMMARY_PROTEIN_IDS


def quant_column_kind(column: object) -> str | None:
    name = str(column)
    if name.endswith(" Unshared Log Ratio"):
        return "unshared_log_ratio"
    if name.endswith(" Log Ratio"):
        return "log_ratio"
    return None


def columns_for_family(columns: list[str], family: str | None) -> list[str]:
    """Devuelve las columnas de una sola familia. Si ninguna coincide, la lista original."""
    if not family:
        return list(columns)
    picked = [c for c in columns if quant_column_kind(c) == family]
    return picked if picked else list(columns)
