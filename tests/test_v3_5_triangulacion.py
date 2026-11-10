"""La triangulación no deja celdas vacías: lo que falta dice no_medida."""

from __future__ import annotations

from pathlib import Path

from oncogo_crc.v3.prereg import PreregistroV3
from oncogo_crc.v3.triangulacion import NO_MEDIDA, tabla_vacia

RAIZ = Path(__file__).resolve().parents[1]


def test_catorce_genes_y_ninguna_celda_vacia():
    prereg = PreregistroV3(RAIZ)
    genes = list(prereg.get("enmienda_v3_4.triangulacion.genes"))
    columnas = list(prereg.get("enmienda_v3_4.triangulacion.columnas"))
    tabla = tabla_vacia(genes, columnas)
    assert len(tabla) == 14
    assert set(columnas).issubset(tabla.columns)
    assert (tabla[columnas] == NO_MEDIDA).all().all()
    assert tabla.isna().sum().sum() == 0
