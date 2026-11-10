"""El cierre v3.7 se niega a firmar lo que no puede leer."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.v3.cierre_v37 import (
    auditar_rai,
    bh,
    emitir,
    exigir_carpeta_inicial,
    filtrar_terminos,
    sellos_coinciden,
    shap_de_un_ajuste,
)


def test_carpeta_con_restos_no_empieza(tmp_path):
    (tmp_path / "preanalitica_intentos.tsv").write_text("\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="restos"):
        exigir_carpeta_inicial(tmp_path)


def test_la_carpeta_de_crudos_si_puede_empezar(tmp_path):
    (tmp_path / "crudos_integridad.tsv").write_text("crudo\n", encoding="utf-8")
    (tmp_path / "crudos_integridad.tsv.sello.json").write_text("{}\n", encoding="utf-8")
    exigir_carpeta_inicial(tmp_path)


def test_fuente_ausente_es_falla(tmp_path):
    certificado = emitir(
        [{"nombre": "suero", "fuente": "results/v3_7/no_esta.tsv", "ruta": tmp_path / "no_esta.tsv",
          "criterio": "el fichero existe", "valor": None, "estado": "pasa"}],
        prohibidas=[],
    )
    assert certificado["casillas"][0]["estado"] == "falla"
    assert certificado["estado"] == "no_publicable"


def test_lecturas_prohibidas_dejan_el_certificado_no_publicable(tmp_path):
    fuente = tmp_path / "celda.tsv"
    fuente.write_text("n\t1\n", encoding="utf-8")
    certificado = emitir(
        [{"nombre": "suero", "fuente": "results/v3_7/celda.tsv", "ruta": fuente,
          "criterio": "el fichero existe", "valor": 1, "estado": "pasa"}],
        prohibidas=["results/v3_6/k1.tsv"],
    )
    assert certificado["estado"] == "no_publicable"
    assert certificado["casillas"][0]["sha256_fuente"]


def test_sellos_que_no_coinciden_fallan(tmp_path):
    comun = {"commit": "a", "arbol_limpio": True, "sha256_preregistro": "b",
             "sha256_enmiendas": "c", "sha256_manifest_entradas": "d", "etapa": "v3_7"}
    (tmp_path / "uno.sello.json").write_text(json.dumps(comun), encoding="utf-8")
    otro = dict(comun, sha256_manifest_entradas=None)
    (tmp_path / "dos.sello.json").write_text(json.dumps(otro), encoding="utf-8")
    assert sellos_coinciden(tmp_path)["estado"] == "falla"


def test_bh_de_los_p_conocidos():
    q = bh([
        ("CN:log2_n_precursores", 0.0027),
        ("CN:piel", 0.0112),
        ("relleno", 0.5),
        ("relleno2", 0.5),
        ("relleno3", 0.5),
        ("relleno4", 0.5),
        ("relleno5", 0.5),
        ("relleno6", 0.5),
        ("relleno7", 0.5),
        ("relleno8", 0.5),
    ])
    assert q["CN:log2_n_precursores"] == pytest.approx(0.027)
    assert q["CN:piel"] == pytest.approx(0.056)


def test_baiap3_sale_en_el_filtro_global_y_ceacam5_no_esta():
    traza = pd.DataFrame([{
        "protein_group": "O94812;O94812-2",
        "genes": "BAIAP3",
        "global_pg_q": 0.0115942,
        "pg_q_min": 0.00328947,
        "n_precursores": 55,
        "fraccion_proteotipicos": 0.0,
        "humano": True,
    }])
    tabla = auditar_rai(traza, matriz=pd.DataFrame({"protein_id": []}), celda=pd.DataFrame({"protein_id": []}))
    por_gen = dict(zip(tabla["gen"], tabla["paso"], strict=True))
    assert por_gen["BAIAP3"] == "filtro_global_1pct"
    assert por_gen["CEACAM5"] == "no_esta_en_el_report"
    assert float(tabla.loc[tabla["gen"].eq("BAIAP3"), "global_pg_q"].iloc[0]) == pytest.approx(0.0115942)


def test_un_termino_con_menos_de_15_genes_no_entra():
    universo = {f"G{i}" for i in range(501)}
    terminos = {
        "GO:pequeno": {f"G{i}" for i in range(14)},
        "GO:entra": {f"G{i}" for i in range(15)},
        "GO:grande": {f"G{i}" for i in range(501)},
    }
    filtrados = filtrar_terminos(terminos, universo, minimo=15, maximo=500)
    assert set(filtrados) == {"GO:entra"}


def test_los_pasos_del_embudo_son_subconjuntos():
    from oncogo_crc.v3.cierre_v37 import embudo_es_encajado
    tabla = pd.DataFrame({"gen": ["A", "B", "C"], "ultimo_paso": [6, 2, 0]})
    assert embudo_es_encajado(tabla, n_pasos=7)


def test_shap_lineal_suma_la_diferencia_del_logit():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(40, 6))
    y = np.array([1] * 20 + [0] * 20)
    shap, diferencia = shap_de_un_ajuste(x, y, semilla=0)
    assert np.allclose(shap.sum(axis=1), diferencia, atol=1e-8)
    from oncogo_crc.v3.cierre_v37 import comprobar_con_shap
    assert comprobar_con_shap(x, y, semilla=0) < 1e-8


def test_gsea_guarda_semilla_y_permutaciones(monkeypatch):
    from oncogo_crc.v3 import cierre_v37

    visto = {}

    def falso_prerank(**kwargs):
        visto.update(kwargs)
        class Resultado:
            res2d = pd.DataFrame(columns=["Term", "NES", "FDR q-val", "NOM p-val"])
        return Resultado()

    monkeypatch.setattr(cierre_v37, "prerank_gsea", falso_prerank)
    cierre_v37.gsea_prerank(pd.Series({"A": 1.0, "B": -1.0}), {"GO:1": ["A", "B"]}, semilla=42, permutaciones=200, minimo=15, maximo=500)
    assert visto["seed"] == 42
    assert visto["permutation_num"] == 200
