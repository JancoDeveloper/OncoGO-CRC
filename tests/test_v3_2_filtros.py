"""Filtros de v3.2 que no necesitan el report real."""

import pandas as pd

from oncogo_crc.stages.v3_2_tejido import ids_descartados
from oncogo_crc.v3.suero_v32 import ids_fuera_de_principal, semilla_v32
from oncogo_crc.v3.prereg import PreregistroV3


def test_el_tejido_saca_con_rev_y_la_marca_de_maxquant():
    ids = pd.Series(["P12345", "CON__Q9C075", "REV__F8VYV2"])
    grupos = pd.DataFrame({
        "Protein IDs": ["P99999", "Q11111"],
        "Reverse": ["", "+"],
        "Potential contaminant": ["+", ""],
        "Only identified by site": ["", ""],
    })
    malos = ids_descartados(ids, grupos)
    assert "CON__Q9C075" in malos
    assert "REV__F8VYV2" in malos
    assert "P99999" in malos
    assert "Q11111" in malos
    assert "P12345" not in malos


def test_el_contraste_principal_saca_contaminante_piel_y_panel():
    traza = pd.DataFrame({
        "protein_group": ["A", "B", "C", "D"],
        "contaminante": [True, False, False, False],
        "contaminante_piel": [False, True, False, False],
        "panel": ["", "", "coagulacion", ""],
    })
    assert ids_fuera_de_principal(traza) == {"A", "B", "C"}


def test_la_semilla_de_v32_sale_de_la_enmienda(tmp_path=None):
    prereg = PreregistroV3(__import__("pathlib").Path(__file__).resolve().parents[1])
    assert semilla_v32(prereg, "k4_alfa_efectivo") == 20260914 + 710
    assert semilla_v32(prereg, "h5_submuestreo") == 20260914 + 730
