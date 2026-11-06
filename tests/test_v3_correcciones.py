"""Correcciones de v3 sobre datos sintéticos. Nada de esto es una cohorte real."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.v3.auc import auc_firmado, delong_ic, nulo_del_maximo_plegado
from oncogo_crc.v3.concordancia import concordancia_de_lin, concordancia_de_signo
from oncogo_crc.v3.etiquetas import CONTRASTE_PRINCIPAL, contraste_pareado, filas_desde_columnas, matriz_log2_lfq
from oncogo_crc.v3.faltantes import fisher_presencia
from oncogo_crc.v3.potencia import clasificar, p_tost, umbral_bh
from oncogo_crc.stages.v3_correcciones import _matriz_lineal_para_fraccion
from oncogo_crc.v3.rutas import exigir_salida_v3, rechazar_si_es_pipeline_antiguo

PROHIBIDA = "polyp"


def test_el_contraste_principal_es_primario_frente_a_mucosa():
    assert CONTRASTE_PRINCIPAL == ("primary", "mucosa")


def test_el_codigo_v3_no_usa_la_etiqueta_historica_de_pxd034575():
    raiz = Path(__file__).resolve().parents[1]
    for carpeta in (raiz / "src/oncogo_crc/v3",):
        for path in carpeta.rglob("*.py"):
            assert PROHIBIDA not in path.read_text(encoding="utf-8").lower(), path
    etapa = raiz / "src/oncogo_crc/stages/v3_correcciones.py"
    assert PROHIBIDA not in etapa.read_text(encoding="utf-8").lower()


def test_el_sufijo_p_es_primario_y_no_otra_lesion():
    columnas = [
        "LFQ intensity CRC02P",
        "LFQ intensity CRC02N",
        "LFQ intensity CRC02M",
    ]
    meta = filas_desde_columnas(columnas)
    tejidos = dict(zip(meta["suffix"], meta["tissue"], strict=True))
    assert tejidos == {"P": "primary", "N": "mucosa", "M": "liver_metastasis"}
    assert PROHIBIDA not in set(meta["tissue"])


def test_sintetico_el_contraste_pareado_ve_el_primario_y_no_la_metastasis():
    columnas = []
    filas = {"protein_id": ["SINTETICO_PROT"]}
    for i in (1, 2, 3, 4):
        for sufijo, valor in (("P", 32.0), ("N", 2.0), ("M", 4.0)):
            nombre = f"LFQ intensity CRC0{i}{sufijo}"
            columnas.append(nombre)
            filas[nombre] = [valor]
    cruda = pd.DataFrame(filas)
    meta = filas_desde_columnas(columnas)
    log2mat = matriz_log2_lfq(cruda)
    principal, n_pares = contraste_pareado(
        log2mat, meta, grupo_a="primary", grupo_b="mucosa",
        nombre="primario_frente_a_mucosa", lfc=0.585, alpha=0.05, min_pares=3,
    )
    assert n_pares == 4
    assert float(principal["log2fc"].iloc[0]) == pytest.approx(np.log2(32) - np.log2(2))


def test_sintetico_tost_declara_equivalente_un_efecto_nulo_preciso():
    coef = np.array([0.0, 3.0])
    se = np.array([0.05, 0.2])
    df = np.array([40.0, 40.0])
    p = p_tost(coef, se, df, 0.585)
    clases = clasificar(np.array([False, True]), np.array([p[0] < 0.05, False]))
    assert clases[0] == "equivalente"
    assert clases[1] == "diferente"
    assert "equivalente" not in clases[1]


def test_sintetico_auc_con_signo_no_se_pliega():
    rng = np.random.default_rng(0)
    positivo = rng.normal(0, 1, size=30)
    negativo = rng.normal(1, 1, size=30)
    auc = auc_firmado(positivo, negativo)
    assert auc < 0.5
    assert auc != max(auc, 1 - auc) or auc == 0.5
    punto, bajo, alto = delong_ic(positivo, negativo)
    assert bajo <= punto <= alto


def test_sintetico_el_nulo_del_maximo_no_da_p_cero():
    rng = np.random.default_rng(1)
    matriz = rng.normal(size=(8, 20))
    etiquetas = np.array([1] * 10 + [0] * 10)
    nulo = nulo_del_maximo_plegado(matriz, etiquetas, n_perm=50, min_n=5, rng=rng)
    assert nulo["p"] >= 1 / 51
    assert nulo["n_permutaciones"] == 50


def test_sintetico_fisher_separa_presencia_de_abundancia():
    matriz = np.array([
        [1.0, 1.0, 1.0, np.nan, np.nan, np.nan],
        [1.0, np.nan, 1.0, 1.0, np.nan, 1.0],
    ])
    etiquetas = np.array(["crc", "crc", "crc", "sano", "sano", "sano"])
    tabla = fisher_presencia(matriz, etiquetas, grupo_a="crc", grupo_b="sano", alpha=0.05)
    assert tabla.loc[0, "presentes_crc"] == 3
    assert tabla.loc[0, "presentes_sano"] == 0
    assert tabla.loc[0, "p_fisher"] < tabla.loc[1, "p_fisher"]


def test_v3_no_lee_el_pipeline_antiguo():
    with pytest.raises(RuntimeError):
        rechazar_si_es_pipeline_antiguo(Path("/tmp/results/differential/tabla.tsv"))
    with pytest.raises(RuntimeError):
        exigir_salida_v3(Path("/tmp/results/go/tabla.tsv"))
    exigir_salida_v3(Path("/tmp/results/v3/tabla.tsv"))


def test_sintetico_lin_vale_uno_si_los_efectos_son_iguales():
    x = np.array([0.2, -0.4, 1.1, 0.7])
    assert concordancia_de_lin(x, x)["ccc"] == pytest.approx(1.0)
    fraccion = concordancia_de_signo(x, x, nombre="SINTETICO", universo=4, universo_nombre="sintetico")
    assert fraccion.numerador == 4
    assert fraccion.denominador == 4


def test_bh_sin_rechazos_cae_en_alpha_sobre_m():
    assert umbral_bh(np.array([0.4, 0.5, 0.6, 0.7]), 0.05) == pytest.approx(0.05 / 4)


def test_sintetico_la_fraccion_de_alta_abundancia_suma_intensidad_lineal():
    mat = pd.DataFrame({"SINTETICO_s": [10.0, 20.0]}, index=["peq", "alb"])
    lineal, escala = _matriz_lineal_para_fraccion(mat, "log2_de_intensidad_lineal")
    assert escala == "intensidad_lineal_recuperada_con_exp2"
    frac = float(lineal.loc["alb", "SINTETICO_s"] / lineal["SINTETICO_s"].sum())
    assert frac == pytest.approx(2**20 / (2**10 + 2**20))
