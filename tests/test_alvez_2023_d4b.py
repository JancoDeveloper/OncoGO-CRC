"""Tests de D4b: comprobación externa de las proteínas de plasma de v1.

Lo que estos tests protegen es una distinción que es fácil de perder: que una
proteína no esté en el panel de Olink no es que no replique. Si el código las
confundiera, el informe diría que dos resultados de v1 han caído cuando en
realidad nadie los ha mirado.

Las matrices son sintéticas y viven solo aquí.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

RAIZ = Path(__file__).resolve().parents[1]


def _modulo():
    """Carga el script de descarga como módulo. No está en el paquete."""
    ruta = RAIZ / "scripts" / "download" / "alvez_2023.py"
    spec = importlib.util.spec_from_file_location("alvez_2023", ruta)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["alvez_2023"] = mod
    spec.loader.exec_module(mod)
    return mod


A = _modulo()


def _por_proteina(filas):
    return pd.DataFrame(
        filas,
        columns=[
            "proteina", "olink_id", "diferencia_npx", "p_ajustado",
            "significativa", "direccion", "etiqueta_del_articulo",
        ],
    )


# -------------------------------------------- la distinción que importa
def test_una_proteina_fuera_del_panel_no_cuenta_como_no_replica():
    d = _por_proteina([["LBP", "OID1", 0.14, 0.08, False, "sube", "not significant"]])
    res = A._comprobar_plasma_v1(d)
    estados = {p["gen"]: p["estado"] for p in res["proteinas"]}
    assert estados["SAA1"] == "no_en_el_panel"
    assert estados["CHMP4A"] == "no_en_el_panel"
    assert estados["LBP"] == "no_replica"
    # Solo una es comprobable: el denominador no son las tres.
    assert res["n_comprobables"] == 1
    assert res["n_declaradas_en_v1"] == 3
    assert res["n_que_replican"] == 0


def test_la_lectura_dice_cuantas_no_se_midieron():
    d = _por_proteina([["LBP", "OID1", 0.14, 0.08, False, "sube", "not significant"]])
    res = A._comprobar_plasma_v1(d)
    assert "no se midieron" in res["lectura"]
    assert "no es lo mismo que no replicar" in res["lectura"]


def test_una_proteina_fuera_del_panel_no_trae_numeros_inventados():
    res = A._comprobar_plasma_v1(_por_proteina([]))
    for p in res["proteinas"]:
        assert p["diferencia_npx"] is None
        assert p["p_ajustado"] is None
        assert p["direccion"] is None


def test_replica_exige_significativa_y_en_la_misma_direccion():
    """Significativa pero bajando no es replicar: v1 las encontró subiendo."""
    sube = _por_proteina([["LBP", "O", 1.0, 1e-9, True, "sube", "significant up"]])
    baja = _por_proteina([["LBP", "O", -1.0, 1e-9, True, "baja", "significant down"]])
    assert A._comprobar_plasma_v1(sube)["n_que_replican"] == 1
    assert A._comprobar_plasma_v1(baja)["n_que_replican"] == 0
    lbp = next(
        p for p in A._comprobar_plasma_v1(baja)["proteinas"] if p["gen"] == "LBP"
    )
    assert lbp["estado"] == "no_replica"


def test_no_significativa_pero_subiendo_tampoco_replica():
    d = _por_proteina([["LBP", "O", 0.14, 0.08, False, "sube", "not significant"]])
    assert A._comprobar_plasma_v1(d)["n_que_replican"] == 0


# -------------------------------------------- recuentos
def test_el_recuento_separa_las_que_suben_de_las_que_bajan():
    d = _por_proteina(
        [
            ["A", "O1", 1.0, 1e-9, True, "sube", "significant up"],
            ["B", "O2", -1.0, 1e-9, True, "baja", "significant down"],
            ["C", "O3", -2.0, 1e-9, True, "baja", "significant down"],
            ["D", "O4", 0.1, 0.5, False, "sube", "not significant"],
        ]
    )
    res = A._fase_aguda(d)
    assert res["n_medidas"] == 4
    assert res["n_significativas"] == 3
    assert res["n_suben"] == 1
    assert res["n_bajan"] == 2


def test_no_se_clasifica_la_fase_aguda_a_mano_y_se_dice_por_que():
    res = A._fase_aguda(_por_proteina([]))
    assert "no se clasifica" in res["nota"].lower()
    assert "preregistradas" in res["nota"]


# -------------------------------------------- parseo del suplementario
def _escribir_sintetico(tmp_path: Path) -> Path:
    """Un suplementario sintético con la estructura real del artículo."""
    path = tmp_path / "SINTETICO_alvez_moesm11.xlsx"
    hoja2a = pd.DataFrame(
        {
            "Cancer": ["Colorectal", "Colorectal", "Colorectal", "Glioma"],
            "Assay": ["ceacam5", " LBP ", "ZZZ", "CEACAM5"],
            "OlinkID": ["OID1", "OID2", "OID3", "OID4"],
            "NPX_difference": [0.52, 0.14, -2.0, 9.9],
            "p.adjusted": [4.7e-5, 0.082, 1e-12, 1e-3],
            "Significance": [
                "significant up", "not significant", "significant down",
                "significant up",
            ],
        }
    )
    with pd.ExcelWriter(path) as w:
        hoja2a.to_excel(w, sheet_name="2a", index=False)
    return path


def test_el_parser_se_queda_solo_con_el_colorrectal(tmp_path):
    d = A._por_proteina(_escribir_sintetico(tmp_path))
    assert len(d) == 3
    assert "Glioma" not in str(d.values)


def test_el_parser_normaliza_el_nombre_del_ensayo(tmp_path):
    """Olink trae espacios y mayúsculas irregulares. Sin normalizar, LBP no cruza."""
    d = A._por_proteina(_escribir_sintetico(tmp_path))
    assert set(d["proteina"]) == {"CEACAM5", "LBP", "ZZZ"}


def test_la_direccion_sale_del_signo_y_no_de_la_etiqueta(tmp_path):
    d = A._por_proteina(_escribir_sintetico(tmp_path))
    zzz = d[d["proteina"] == "ZZZ"].iloc[0]
    assert zzz["direccion"] == "baja"
    assert zzz["diferencia_npx"] == pytest.approx(-2.0)
    # La etiqueta del artículo se conserva aparte, para poder auditar discrepancias.
    assert zzz["etiqueta_del_articulo"] == "significant down"


def test_el_parser_ordena_por_p_ajustado(tmp_path):
    d = A._por_proteina(_escribir_sintetico(tmp_path))
    assert list(d["p_ajustado"]) == sorted(d["p_ajustado"])


def test_la_significacion_sale_de_la_etiqueta_del_articulo(tmp_path):
    """No se reaplica un umbral: el artículo ya corrigió por multiplicidad."""
    d = A._por_proteina(_escribir_sintetico(tmp_path))
    assert list(d.set_index("proteina")["significativa"][["CEACAM5", "LBP", "ZZZ"]]) == [
        True,
        False,
        True,
    ]
