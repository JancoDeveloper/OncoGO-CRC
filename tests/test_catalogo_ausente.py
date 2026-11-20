"""El catálogo de literatura no es una entrada de la corrida sellada."""

from pathlib import Path

import pytest

from oncogo_crc.tesis.cierre import CatalogoNoDisponible, cargar_catalogo


def test_sin_archivo_el_catalogo_explica_la_ausencia(tmp_path: Path):
    with pytest.raises(CatalogoNoDisponible):
        cargar_catalogo(tmp_path / "candidatas_literatura.tsv")


def test_la_ruta_por_defecto_no_esta_en_el_arbol():
    with pytest.raises(CatalogoNoDisponible):
        cargar_catalogo()
