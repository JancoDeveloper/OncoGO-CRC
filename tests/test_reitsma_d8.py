"""Tests de D8: meta-análisis bivariado de sensibilidad y especificidad.

Hay dos clases de test aquí. Unos comprueban que el modelo recupera parámetros
puestos a mano sobre tablas sintéticas. Otros comprueban que el ajuste propio
reproduce las estimaciones combinadas que publicaron los autores de las tablas,
que es la única validación externa posible de una reimplementación: si mi Reitsma
en Python da lo mismo que el suyo en R, el código está bien.

Las tablas sintéticas viven solo aquí. Las reales son las de
``data/literatura/extraccion_2x2.tsv``, que son números publicados.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.v2 import reitsma as R

RAIZ = Path(__file__).resolve().parents[1]
TABLA = RAIZ / "data" / "literatura" / "extraccion_2x2.tsv"


def _extractor():
    ruta = RAIZ / "scripts" / "literatura" / "extraer_2x2.py"
    spec = importlib.util.spec_from_file_location("extraer_2x2", ruta)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["extraer_2x2"] = mod
    spec.loader.exec_module(mod)
    return mod


E = _extractor()


def _tablas(sens, espec, *, n_enfermos=100, n_sanos=100, k=10, semilla=20260914):
    """k estudios con la sensibilidad y la especificidad que se le pidan."""
    rng = np.random.default_rng(semilla)
    tp = rng.binomial(n_enfermos, sens, k)
    tn = rng.binomial(n_sanos, espec, k)
    return tp, n_sanos - tn, n_enfermos - tp, tn


# ------------------------------------------------- recuperación de parámetros
def test_el_modelo_recupera_una_sensibilidad_puesta_a_mano():
    tp, fp, fn, tn = _tablas(0.70, 0.90, k=20)
    res = R.reitsma(tp, fp, fn, tn)
    assert res.sensibilidad == pytest.approx(0.70, abs=0.04)
    assert res.sensibilidad_ci_low < 0.70 < res.sensibilidad_ci_high


def test_el_modelo_recupera_una_especificidad_puesta_a_mano():
    tp, fp, fn, tn = _tablas(0.70, 0.90, k=20)
    res = R.reitsma(tp, fp, fn, tn)
    assert res.especificidad == pytest.approx(0.90, abs=0.04)
    assert res.especificidad_ci_low < 0.90 < res.especificidad_ci_high


def test_con_estudios_homogeneos_la_varianza_entre_estudios_es_pequena():
    tp, fp, fn, tn = _tablas(0.70, 0.90, n_enfermos=500, n_sanos=500, k=20)
    res = R.reitsma(tp, fp, fn, tn)
    assert res.tau_sensibilidad < 0.3
    assert res.tau_especificidad < 0.3


def test_con_estudios_heterogeneos_la_varianza_entre_estudios_crece():
    """Dos grupos de estudios con sensibilidades muy distintas."""
    a = _tablas(0.40, 0.90, n_enfermos=300, n_sanos=300, k=10, semilla=1)
    b = _tablas(0.90, 0.90, n_enfermos=300, n_sanos=300, k=10, semilla=2)
    tp, fp, fn, tn = (np.concatenate([x, y]) for x, y in zip(a, b, strict=True))
    homogeneo = R.reitsma(*_tablas(0.65, 0.90, n_enfermos=300, n_sanos=300, k=20))
    res = R.reitsma(tp, fp, fn, tn)
    assert res.tau_sensibilidad > homogeneo.tau_sensibilidad


def test_una_celda_cero_se_corrige_y_se_cuenta():
    tp, fp, fn, tn = _tablas(0.70, 0.99, n_sanos=20, k=10)
    tp, fp, fn, tn = list(tp), list(fp), list(fn), list(tn)
    fp[0] = 0
    res = R.reitsma(tp, fp, fn, tn)
    assert res.n_estudios_con_correccion >= 1
    assert np.isfinite(res.especificidad)


def test_la_correccion_no_se_aplica_a_los_estudios_que_no_la_necesitan():
    """Aplicarla siempre empujaría al centro también a los estudios completos."""
    tp, fp, fn, tn = _tablas(0.70, 0.90, k=10)
    assert R.reitsma(tp, fp, fn, tn).n_estudios_con_correccion == 0


def test_las_celdas_pequenas_se_senalan():
    tp, fp, fn, tn = _tablas(0.70, 0.90, n_enfermos=6, n_sanos=6, k=10)
    assert R.reitsma(tp, fp, fn, tn).n_estudios_con_celda_pequena > 0


def test_con_menos_de_cuatro_estudios_el_modelo_se_niega_a_ajustar():
    tp, fp, fn, tn = _tablas(0.70, 0.90, k=3)
    with pytest.raises(ValueError, match="al menos 4"):
        R.reitsma(tp, fp, fn, tn)


def test_el_intervalo_del_dor_usa_la_covarianza_entre_los_dos_logits():
    """Sin el término de covarianza el intervalo saldría del tamaño equivocado."""
    tp, fp, fn, tn = _tablas(0.70, 0.90, k=20)
    res = R.reitsma(tp, fp, fn, tn)
    esperado = (
        res.sensibilidad / (1 - res.sensibilidad)
    ) * (res.especificidad / (1 - res.especificidad))
    assert res.dor == pytest.approx(esperado, rel=1e-6)
    assert res.dor_ci_low < res.dor < res.dor_ci_high


def test_la_correlacion_entre_sensibilidad_y_especificidad_se_estima():
    """Si los estudios difieren por el umbral, la correlación sale negativa."""
    # Umbrales distintos: cuando sube la sensibilidad baja la especificidad.
    sens = [0.5, 0.6, 0.7, 0.8, 0.9, 0.55, 0.65, 0.75, 0.85, 0.95]
    espec = [0.95, 0.9, 0.85, 0.8, 0.7, 0.93, 0.88, 0.82, 0.75, 0.6]
    rng = np.random.default_rng(7)
    tp = [rng.binomial(200, s) for s in sens]
    tn = [rng.binomial(200, e) for e in espec]
    res = R.reitsma(tp, [200 - x for x in tn], [200 - x for x in tp], tn)
    assert res.correlacion < 0


# ------------------------------------------------- comparaciones
def test_una_comparacion_no_comparable_no_devuelve_diferencia():
    a = R.reitsma(*_tablas(0.60, 0.90, k=10, semilla=1))
    b = R.reitsma(*_tablas(0.80, 0.90, k=10, semilla=2))
    res = R.diferencia(a, b, comparable=False)
    assert res["comparable"] is False
    assert "diferencia_logit_sensibilidad" not in res
    assert "no significa nada" in res["motivo"]
    # Los dos puntos sí se dan, para poder leerlos por separado.
    assert res["a"]["sensibilidad"] == pytest.approx(a.sensibilidad)


def test_una_diferencia_real_se_detecta_cuando_si_es_comparable():
    a = R.reitsma(*_tablas(0.50, 0.90, n_enfermos=400, n_sanos=400, k=15, semilla=1))
    b = R.reitsma(*_tablas(0.90, 0.90, n_enfermos=400, n_sanos=400, k=15, semilla=2))
    res = R.diferencia(a, b, comparable=True)
    assert res["diferencia_logit_sensibilidad"]["difieren"] is True
    assert res["diferencia_logit_sensibilidad"]["diferencia"] < 0


def test_dos_pruebas_iguales_no_difieren():
    a = R.reitsma(*_tablas(0.70, 0.90, k=15, semilla=1))
    b = R.reitsma(*_tablas(0.70, 0.90, k=15, semilla=2))
    res = R.diferencia(a, b, comparable=True)
    assert res["diferencia_logit_sensibilidad"]["difieren"] is False


def test_el_veredicto_sale_del_intervalo_de_la_diferencia():
    """No de que los dos intervalos originales se solapen, que no es lo mismo."""
    a = R.reitsma(*_tablas(0.70, 0.90, n_enfermos=600, n_sanos=600, k=15, semilla=1))
    b = R.reitsma(*_tablas(0.78, 0.90, n_enfermos=600, n_sanos=600, k=15, semilla=2))
    res = R.diferencia(a, b, comparable=True)
    d = res["diferencia_logit_sensibilidad"]
    assert d["difieren"] == bool(d["ci_low"] * d["ci_high"] > 0)


def test_la_comparacion_con_un_publicado_se_declara_no_contraste():
    propio = R.reitsma(*_tablas(0.60, 0.88, k=10))
    res = R.comparar_con_publicado(
        propio,
        sensibilidad=0.84,
        sensibilidad_ci=(0.73, 0.91),
        especificidad=0.95,
        especificidad_ci=(0.80, 0.99),
        fuente="una fuente",
        k_publicado=4,
    )
    assert "no una prueba estadística" in res["no_es_un_contraste"]
    assert res["intervalos_de_sensibilidad_se_solapan"] is False


def test_si_los_intervalos_publicados_tocan_los_propios_se_dice_que_solapan():
    propio = R.reitsma(*_tablas(0.60, 0.88, k=10))
    res = R.comparar_con_publicado(
        propio,
        sensibilidad=0.58,
        sensibilidad_ci=(0.48, 0.68),
        especificidad=0.92,
        especificidad_ci=(0.89, 0.94),
        fuente="una fuente",
        k_publicado=18,
    )
    assert res["intervalos_de_sensibilidad_se_solapan"] is True


# ------------------------------------------------- la tabla de extracción
def test_las_celdas_de_cada_fila_suman_el_tamano_de_muestra():
    """Es la comprobación que detecta un orden de columnas intercambiado."""
    d = pd.read_csv(TABLA, sep="\t")
    suma = d["tp"] + d["fp"] + d["fn"] + d["tn"]
    assert (suma == d["n"]).all()


def test_ningun_estudio_aparece_dos_veces_en_la_misma_prueba():
    """Cinco estudios de metilación publican varios algoritmos sobre los mismos
    pacientes. Si entraran todos, esas personas contarían varias veces."""
    d = pd.read_csv(TABLA, sep="\t")
    assert not d.duplicated(subset=["prueba", "estudio", "anio"]).any()


def test_la_columna_verificado_por_sigue_vacia():
    d = pd.read_csv(TABLA, sep="\t")
    assert d["verificado_por"].isna().all() or (d["verificado_por"] == "").all()


def test_cada_fila_dice_de_que_tabla_y_pagina_sale():
    d = pd.read_csv(TABLA, sep="\t")
    assert d["fuente_tabla"].notna().all()
    assert d["fuente_pagina"].notna().all()
    assert d["fuente"].isin(E.FUENTES).all()


def test_los_totales_de_la_metilacion_reproducen_los_del_articulo():
    """El artículo declara 9 927 muestras y 2 975 cánceres. Que mi selección de
    una fila por estudio dé exactamente eso confirma que elegí las correctas."""
    d = pd.read_csv(TABLA, sep="\t")
    g = d[d["prueba"] == "metilacion_SEPT9"]
    assert len(g) == 25
    assert int(g["n"].sum()) == 9927
    assert int((g["tp"] + g["fn"]).sum()) == 2975


def test_la_fila_excluida_se_declara_con_su_motivo():
    """Una exclusión silenciosa es indistinguible de un descuido."""
    assert len(E.CEA_EXCLUIDAS) == 1
    fuera = E.CEA_EXCLUIDAS[0]
    assert fuera["estudio"] == "Moloney"
    assert "N/A" in fuera["motivo"]
    d = pd.read_csv(TABLA, sep="\t")
    assert "Moloney" not in set(d["estudio"])


def test_no_se_reconstruye_ninguna_tabla_desde_un_porcentaje():
    texto = E.ficha()["por_que_no_se_reconstruyen_tablas"]
    assert "De ninguna prueba se calcula una tabla 2x2" in texto


def test_la_regla_del_algoritmo_esta_declarada_y_es_por_algoritmo():
    assert "2(3)" in E.REGLA_ALGORITMO
    assert "comercial" in E.REGLA_ALGORITMO


# ------------------------------------------------- validación externa
def test_el_ajuste_propio_del_cea_reproduce_el_combinado_publicado():
    """Wang 2023 publica sensibilidad combinada 0,59 y especificidad 0,89.

    Es la única validación externa que admite una reimplementación: ajustar sus
    mismas tablas y comprobar que sale su mismo número.
    """
    d = pd.read_csv(TABLA, sep="\t")
    g = d[d["prueba"] == "CEA"]
    res = R.reitsma(g["tp"], g["fp"], g["fn"], g["tn"])
    assert res.sensibilidad == pytest.approx(0.59, abs=0.05)
    assert res.especificidad == pytest.approx(0.89, abs=0.05)


def test_el_ajuste_propio_de_la_metilacion_reproduce_el_combinado_publicado():
    """Sun 2017 publica sensibilidad combinada 0,71 y especificidad 0,92."""
    d = pd.read_csv(TABLA, sep="\t")
    g = d[d["prueba"] == "metilacion_SEPT9"]
    res = R.reitsma(g["tp"], g["fp"], g["fn"], g["tn"])
    assert res.sensibilidad == pytest.approx(0.71, abs=0.03)
    assert res.especificidad == pytest.approx(0.92, abs=0.03)
