"""Cierre de tesis: etiquetas clínicas, especificidad y límite de los péptidos."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from oncogo_crc.settings import load_settings
from oncogo_crc.stages.tesis_cierre import run
from oncogo_crc.tesis.cierre import (
    cargar_catalogo,
    clasificar_corrida,
    cruzar_catalogo,
    etiqueta_tanaka,
    etiquetas_son_clinicas,
    grupos_de_muestras,
    lectura_con_ruido,
    puntuacion_enriquecimiento,
    veredicto_especificidad,
)
from oncogo_crc.utils.dia_nn_plasma import build_plasma_metadata_for_pxd055218_runs


def test_pxd055218_conserva_el_grupo_de_enfermedad():
    meta = build_plasma_metadata_for_pxd055218_runs(
        ["CP_01", "CN_01", "Normal_01", "otra_corrida"]
    )
    assert set(meta["condition"]) == {"crc", "colonoscopy_negative", "healthy"}
    assert clasificar_corrida("ruta/CP_03.raw") == "crc"
    assert clasificar_corrida("CN_3") == "colonoscopy_negative"


def test_el_corte_alfabetico_no_es_un_diagnostico():
    meta = pd.DataFrame(
        {
            "sample_id": ["a", "b"],
            "condition": ["sin_etiqueta_clinica", "sin_etiqueta_clinica"],
            "batch": ["PXD036594_sin_etiqueta_clinica", "PXD036594_sin_etiqueta_clinica"],
        }
    )
    assert etiquetas_son_clinicas(meta) is False
    grupos = grupos_de_muestras(["a", "b"], meta)
    assert grupos == {"crc": [], "colonoscopy_negative": [], "healthy": []}


def test_la_regla_de_especificidad_separa_los_tres_casos():
    assert veredicto_especificidad(0.82, 0.61, 0.55) == "especifica"
    assert veredicto_especificidad(0.60, 0.40, 0.90) == "compartida_con_enfermedad"
    assert veredicto_especificidad(0.62, 0.51, 0.58) == "inespecifica"
    assert veredicto_especificidad(None, None, 0.9) == "no_evaluable"


def test_el_cruce_acepta_accesion_uniprot():
    catalogo = cargar_catalogo()
    contraste = pd.DataFrame(
        {
            "protein_id": ["O00429"],
            "log2fc": [1.1],
            "padj": [0.01],
            "significant": [True],
        }
    )
    cruce = cruzar_catalogo(
        catalogo,
        contraste,
        cohorte="suero",
        mapa_uniprot={"BAIAP3": "O00429"},
    )
    fila = cruce.loc[cruce["gene"] == "BAIAP3"].iloc[0]
    assert fila["detectada"]
    assert fila["veredicto"] == "replicada"
    assert fila["origen"] == "literatura"


def test_el_catalogo_no_presenta_literatura_como_hallazgo_propio():
    catalogo = cargar_catalogo()
    assert (catalogo["origen"] == "literatura").all()
    assert "BAIAP3" in set(catalogo["gene"])
    contraste = pd.DataFrame(
        {
            "protein_id": ["BAIAP3", "TOP2A"],
            "log2fc": [1.2, -0.4],
            "padj": [0.01, 0.2],
            "significant": [True, False],
        }
    )
    cruce = cruzar_catalogo(catalogo, contraste, cohorte="sintetica")
    bai = cruce.loc[cruce["gene"] == "BAIAP3"].iloc[0]
    assert bai["origen"] == "literatura"
    assert bai["veredicto"] == "replicada"
    assert cruce.loc[cruce["gene"] == "CA1", "veredicto"].iloc[0] == "no_evaluable"


def test_la_etapa_corre_tres_grupos_y_no_confirma_peptidos(tmp_path):
    settings = load_settings(data_workspace=tmp_path)
    matriz = tmp_path / "data" / "raw" / "plasma"
    matriz.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    columnas = [f"CP_{i:02d}" for i in range(8)] + [f"CN_{i:02d}" for i in range(8)] + [
        f"Normal_{i:02d}" for i in range(8)
    ]
    proteinas = ["BAIAP3", "TFF3", "TOP2A"] + [f"P{i:03d}" for i in range(12)]
    filas = []
    for gen in proteinas:
        base = rng.normal(8, 0.3, size=24)
        if gen == "BAIAP3":
            base[:8] += 4
        if gen == "TFF3":
            base[:16] += 4
        filas.append([gen, *base])
    pd.DataFrame(filas, columns=["protein_id", *columnas]).to_csv(
        matriz / "plasma_crc_matrix.tsv", sep="\t", index=False
    )
    run(settings)
    salida = tmp_path / "results" / "tesis_cierre"
    auditoria = json.loads((salida / "cierre_audit.json").read_text())
    assert auditoria["estado_suero"] == "tres_grupos"
    assert auditoria["contrastes"]["crc_vs_enfermedad"]["estado"] == "calculado"
    especificidad = pd.read_csv(salida / "especificidad_candidatas.tsv", sep="\t")
    assert especificidad.loc[especificidad["gene"] == "BAIAP3", "veredicto"].iloc[0] == "especifica"
    assert especificidad.loc[especificidad["gene"] == "TFF3", "veredicto"].iloc[0] == "compartida_con_enfermedad"
    peptidos = pd.read_csv(salida / "peptidos_neoantigenos.tsv", sep="\t")
    assert (peptidos["origen"] == "literatura").all()
    assert "no_evaluable_en_matriz_proteica" in set(peptidos["nivel_en_esta_tesis"])


def test_la_columna_p_de_tanaka_es_el_tumor_primario():
    assert etiqueta_tanaka("LFQ intensity CRC02P") == ("02", "tumor_primario")
    assert etiqueta_tanaka("LFQ intensity CRC02M") == ("02", "metastasis_hepatica")
    assert etiqueta_tanaka("LFQ intensity CRC02N") == ("02", "mucosa_normal")
    assert etiqueta_tanaka("PRIDE-T-01") is None


def test_la_bandera_de_ruido_no_relaja_el_auc():
    assert veredicto_especificidad(0.82, 0.61, 0.55) == "especifica"
    assert lectura_con_ruido("especifica", "fase_aguda") == "auc_pasa_pero_fase_aguda"
    assert lectura_con_ruido("especifica", "") == "especifica"
    assert lectura_con_ruido("inespecifica", "hemolisis") == "inespecifica"


def test_el_conjunto_en_la_cabeza_del_ranking_tiene_enriquecimiento_positivo():
    puntajes = np.linspace(6, -6, 12)
    mascara = np.array([True] * 10 + [False, False])
    assert puntuacion_enriquecimiento(puntajes, mascara) > 0

