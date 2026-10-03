"""Tests de D6: cáncer frente a enfermedad benigna.

Lo que estos tests protegen es el extractor de PDF. Un extractor que desplaza una
columna produce una tabla con pinta impecable, y nadie lo nota leyéndola. Las
comprobaciones de coherencia son la única defensa automática, así que se verifica
que de verdad fallan cuando tienen que fallar.

Los textos son sintéticos y viven solo aquí.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

RAIZ = Path(__file__).resolve().parents[1]


def _modulo():
    ruta = RAIZ / "scripts" / "download" / "wannberg_2026.py"
    spec = importlib.util.spec_from_file_location("wannberg_2026", ruta)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["wannberg_2026"] = mod
    spec.loader.exec_module(mod)
    return mod


W = _modulo()

TEXTO = """Supplementary Table S5
Cohort Assay OR SE Statistic P value Adjusted P value CI low CI high
Discovery RRM2 2.39 0.12 7.23 4.74e-13 4.63e-10 1.91 3.07
Discovery S100A12 1.27 0.06 3.72 1.97e-04 5.13e-04 1.12 1.43
Replication RRM2 2.86 0.22 4.75 2.02e-06 5.11e-04 1.93 4.61
Replication S100A12 1.45 0.11 3.29 1.01e-03 4.05e-03 1.16 1.81
"""


# --------------------------------------------------- extracción
def test_se_extraen_las_filas_de_las_dos_cohortes():
    d = W.extraer_s5(TEXTO)
    assert len(d) == 4
    assert set(d["cohorte"]) == {"descubrimiento", "replica"}
    assert set(d["proteina"]) == {"RRM2", "S100A12"}


def test_los_numeros_se_leen_en_el_orden_de_las_columnas():
    d = W.extraer_s5(TEXTO)
    r = d[(d["proteina"] == "RRM2") & (d["cohorte"] == "descubrimiento")].iloc[0]
    assert r["odds_ratio"] == pytest.approx(2.39)
    assert r["se"] == pytest.approx(0.12)
    assert r["p"] == pytest.approx(4.74e-13)
    assert r["p_ajustado"] == pytest.approx(4.63e-10)
    assert r["ic95_bajo"] == pytest.approx(1.91)
    assert r["ic95_alto"] == pytest.approx(3.07)


def test_una_fila_partida_a_medias_no_entra():
    """Mejor extraer menos filas que colar una con el intervalo de la siguiente."""
    roto = "Discovery RRM2 2.39 0.12 7.23 4.74e-13\nDiscovery PARP1 3.06 0.15 7.26 3.81e-13 4.46e-10 2.3 4.21"
    d = W.extraer_s5(roto)
    assert list(d["proteina"]) == ["PARP1"]


def test_el_salto_de_linea_del_pdf_no_rompe_la_extraccion():
    """El extractor de texto parte donde quiere; el patrón se aplica al texto plano."""
    pegado = TEXTO.replace("\n", " \n ").replace("2.39 0.12", "2.39\n0.12")
    assert len(W.extraer_s5(pegado)) == 4


def test_cada_fila_dice_de_que_tabla_y_pagina_sale():
    d = W.extraer_s5(TEXTO, pagina=6)
    assert set(d["fuente_tabla"]) == {"Supplementary Table S5"}
    assert set(d["fuente_pagina"]) == {6}


def test_la_columna_verificado_por_la_deja_vacia_el_script():
    d = W.extraer_s5(TEXTO)
    assert (d["verificado_por"] == "").all()


def test_la_localizacion_se_anota_y_lo_no_anotado_se_dice():
    d = W.extraer_s5(TEXTO + "Discovery ZZZZ 2.0 0.1 5.0 1e-05 1e-04 1.5 2.5\n")
    loc = dict(zip(d["proteina"], d["localizacion_uniprot"], strict=True))
    assert loc["S100A12"] == "inflamatoria_secretada"
    assert loc["ZZZZ"] == "sin_anotar"


def test_un_pdf_con_otra_maquetacion_devuelve_la_tabla_vacia():
    """Si el patrón no encuentra nada, no se inventa nada: la tabla sale vacía."""
    d = W.extraer_s5("Supplementary Table S5\nno hay ninguna fila con este formato")
    assert d.empty
    assert list(d.columns)[:2] == ["cohorte", "proteina"]


# --------------------------------------------------- comprobaciones
def test_la_extraccion_coherente_se_declara_coherente():
    res = W.comprobar_la_extraccion(W.extraer_s5(TEXTO))
    assert res["coherente"] is True
    assert res["problemas"] == []
    assert res["n_proteinas"] == 2


def test_un_or_fuera_de_su_intervalo_se_detecta():
    """Es la firma típica de una columna desplazada."""
    d = W.extraer_s5(TEXTO)
    d.loc[0, "odds_ratio"] = 99.0
    res = W.comprobar_la_extraccion(d)
    assert res["coherente"] is False
    assert any("fuera de su intervalo" in p for p in res["problemas"])


def test_un_p_ajustado_menor_que_el_crudo_se_detecta():
    d = W.extraer_s5(TEXTO)
    d.loc[0, "p_ajustado"] = d.loc[0, "p"] / 10
    res = W.comprobar_la_extraccion(d)
    assert res["coherente"] is False
    assert any("p ajustado menor" in p for p in res["problemas"])


def test_cohortes_con_proteinas_distintas_se_detecta():
    d = W.extraer_s5(TEXTO)
    d.loc[d.index[-1], "proteina"] = "OTRA"
    res = W.comprobar_la_extraccion(d)
    assert res["coherente"] is False
    assert any("no traen las mismas proteínas" in p for p in res["problemas"])


def test_se_cuenta_cuantas_filas_esperan_verificacion_humana():
    res = W.comprobar_la_extraccion(W.extraer_s5(TEXTO))
    assert res["pendiente_de_verificacion_humana"] == 4
    assert "una persona" in res["regla"]


# --------------------------------------------------- perfil y comparación
def test_el_perfil_cuenta_las_inflamatorias_aparte():
    res = W.perfil_de_las_proteinas(W.extraer_s5(TEXTO))
    assert res["n_inflamatorias"] == 1
    assert res["inflamatorias"] == ["S100A12"]
    assert res["por_localizacion"]["intracelular"] == 1


def test_el_perfil_dice_que_las_proteinas_de_v1_no_estan():
    res = W.perfil_de_las_proteinas(W.extraer_s5(TEXTO))
    assert res["plasma_de_v1_presente"] == {
        "SAA1": False,
        "LBP": False,
        "CHMP4A": False,
    }


def test_sin_d4b_no_se_inventa_el_auc_contra_sanos(tmp_path):
    res = W.comparacion_del_control(tmp_path / "no_existe.json")
    assert res["contra_personas_sanas"]["estado"] == "no_disponible"
    assert "lectura" not in res
    assert "diferencia_de_auc" not in res


def test_con_d4b_la_comparacion_sale_de_las_dos_fuentes(tmp_path):
    import json

    path = tmp_path / "estado.json"
    path.write_text(
        json.dumps(
            {
                "auc_segun_cuantas_proteinas": {
                    "con_el_panel_completo": {"n_proteinas": 1463, "auc": 0.994}
                },
                "cohorte": {"que_son_los_controles": "sanos del programa Wellness"},
            }
        ),
        encoding="utf-8",
    )
    res = W.comparacion_del_control(path)
    assert res["contra_personas_sanas"]["auc"] == 0.994
    assert res["diferencia_de_auc"] == pytest.approx(0.194, abs=1e-6)
    assert "contra quién se compara" in res["lectura"]
    # Los dos paneles tienen que tener el mismo tamaño para que comparar valga.
    assert (
        res["contra_personas_sanas"]["n_proteinas_del_panel"]
        == res["contra_pacientes_con_sintomas"]["n_proteinas_del_panel"]
    )


def test_el_perfil_no_falla_con_la_tabla_vacia():
    res = W.perfil_de_las_proteinas(
        pd.DataFrame(
            columns=["cohorte", "proteina", "odds_ratio", "localizacion_uniprot"]
        )
    )
    assert res["n_proteinas"] == 0
    assert res["or_minimo_entre_las_inflamatorias"] is None
