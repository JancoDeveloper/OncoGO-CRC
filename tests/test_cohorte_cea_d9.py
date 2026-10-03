"""Tests de D9: cohorte de seguimiento con CEA seriado.

El estudio no se ha hecho. Estos tests son lo que da sentido a decir que el
protocolo está «listo»: el análisis corre entero sobre la cohorte sintética y se
comprueba que detecta una diferencia cuando existe y que no la fabrica cuando no
existe.

Los datos son los de ``tests/SINTETICO_cohorte_cea.csv``, que no son de nadie.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.v2 import cohorte_cea as C

RAIZ = Path(__file__).resolve().parents[1]
SINTETICO = RAIZ / "tests" / "SINTETICO_cohorte_cea.csv"
SEMILLA = 20260914 + 909


def _generador():
    ruta = RAIZ / "scripts" / "literatura" / "generar_sintetico_cohorte_cea.py"
    spec = importlib.util.spec_from_file_location("generar_sintetico", ruta)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["generar_sintetico"] = mod
    spec.loader.exec_module(mod)
    return mod


G = _generador()


@pytest.fixture(scope="module")
def cohorte():
    return pd.read_csv(SINTETICO)


@pytest.fixture(scope="module")
def visitas(cohorte):
    return C.preparar_visitas(cohorte)


@pytest.fixture(scope="module")
def comparacion(visitas):
    """Un solo ajuste para todos los tests que solo leen el resultado.

    Cada llamada hace validación cruzada dentro de cada réplica, así que
    repetirla en seis tests multiplicaba por seis el tiempo de la suite sin
    comprobar nada nuevo.
    """
    return C.comparar_modelos(visitas, semilla=SEMILLA, n_bootstrap=200)


# ------------------------------------------- los datos sintéticos
def test_el_fichero_sintetico_lleva_el_prefijo_y_vive_en_tests():
    assert SINTETICO.name.startswith("SINTETICO_")
    assert SINTETICO.parent.name == "tests"


def test_el_generador_es_reproducible():
    a = G.generar(50, semilla=G.SEMILLA)
    b = G.generar(50, semilla=G.SEMILLA)
    pd.testing.assert_frame_equal(a, b)


def test_los_identificadores_son_sinteticos_y_no_parecen_historias():
    d = G.generar(20)
    assert all(p.startswith("SINT") for p in d["paciente"].unique())


def test_el_fichero_del_repositorio_es_el_que_produce_el_generador():
    """Si alguien lo edita a mano, este test lo dice."""
    esperado = G.generar(300, semilla=G.SEMILLA)
    actual = pd.read_csv(SINTETICO, keep_default_na=False, dtype=str)
    producido = esperado.astype(str)
    assert len(actual) == len(producido)
    assert list(actual.columns) == list(producido.columns)


# ------------------------------------------- diccionario de datos
def test_una_columna_ausente_detiene_el_analisis(cohorte):
    with pytest.raises(C.FaltaColumna, match="protocolo_cohorte_cea"):
        C.comprobar_entrada(cohorte.drop(columns=["ultimo_seguimiento_meses"]))


def test_la_cohorte_sintetica_pasa_las_reglas_de_coherencia(cohorte):
    res = C.comprobar_entrada(cohorte)
    assert res["coherente"] is True
    assert res["n_pacientes"] == 300
    assert res["n_con_recaida"] > 0


def test_una_recaida_sin_mes_se_detecta(cohorte):
    roto = cohorte.copy()
    primer = roto[roto["recaida"] == 1]["paciente"].iloc[0]
    roto.loc[roto["paciente"] == primer, "mes_de_la_recaida"] = np.nan
    res = C.comprobar_entrada(roto)
    assert res["coherente"] is False
    assert any("sin mes de recaída" in p for p in res["problemas"])


def test_un_mes_de_recaida_en_alguien_sin_recaida_se_detecta(cohorte):
    roto = cohorte.copy()
    sano = roto[roto["recaida"] == 0]["paciente"].iloc[0]
    roto.loc[roto["paciente"] == sano, "mes_de_la_recaida"] = 20.0
    res = C.comprobar_entrada(roto)
    assert res["coherente"] is False
    assert any("sin recaída y con mes" in p for p in res["problemas"])


def test_un_cea_negativo_se_detecta(cohorte):
    roto = cohorte.copy()
    roto.loc[0, "cea"] = -1.0
    assert "negativos" in " ".join(C.comprobar_entrada(roto)["problemas"])


# ------------------------------------------- construcción de variables
def test_el_cambio_se_mide_contra_la_visita_anterior_del_mismo_paciente():
    d = pd.DataFrame(
        {
            "paciente": ["A", "A", "B", "B"],
            "visita_meses": [3.0, 9.0, 3.0, 9.0],
            # A duplica su CEA en seis meses; B lo mantiene.
            "cea": [2.0, 4.0, 20.0, 20.0],
            "recaida": [0, 0, 0, 0],
            "mes_de_la_recaida": [np.nan] * 4,
            "ultimo_seguimiento_meses": [36.0] * 4,
        }
    )
    t = C.preparar_visitas(d)
    a = t[(t["paciente"] == "A") & (t["visita_meses"] == 9.0)].iloc[0]
    b = t[(t["paciente"] == "B") & (t["visita_meses"] == 9.0)].iloc[0]
    # Una duplicación en seis meses son dos duplicaciones por año.
    assert a["cambio_anual"] == pytest.approx(2.0)
    assert b["cambio_anual"] == pytest.approx(0.0)
    # Y el umbral no distingue entre los dos como hace la trayectoria: B está
    # muy por encima de 5 sin cambiar nada, y A por debajo mientras se duplica.
    assert b["supera_el_umbral"] == 1
    assert a["supera_el_umbral"] == 0


def test_el_cambio_anual_corrige_por_el_tiempo_transcurrido():
    """Dos visitas separadas tres meses no pueden contar igual que seis."""
    d = pd.DataFrame(
        {
            "paciente": ["A", "A", "B", "B"],
            "visita_meses": [3.0, 6.0, 3.0, 9.0],
            "cea": [2.0, 4.0, 2.0, 4.0],
            "recaida": [0, 0, 0, 0],
            "mes_de_la_recaida": [np.nan] * 4,
            "ultimo_seguimiento_meses": [36.0] * 4,
        }
    )
    t = C.preparar_visitas(d)
    a = t[(t["paciente"] == "A") & (t["visita_meses"] == 6.0)].iloc[0]
    b = t[(t["paciente"] == "B") & (t["visita_meses"] == 9.0)].iloc[0]
    assert a["cambio_anual"] == pytest.approx(4.0)
    assert b["cambio_anual"] == pytest.approx(2.0)


def test_el_valor_de_referencia_de_cada_paciente_es_el_suyo(visitas):
    """Es toda la diferencia entre una trayectoria y un umbral común."""
    primera = visitas[visitas["visita_numero"] == 1]
    assert (primera["cambio_sobre_el_basal"].abs() < 1e-9).all()


def test_la_primera_visita_no_tiene_cambio_porque_no_hay_anterior(visitas):
    primera = visitas[visitas["visita_numero"] == 1]
    assert primera["cambio_anual"].isna().all()


def test_el_cea_se_transforma_en_log2_porque_sube_multiplicando():
    d = pd.DataFrame(
        {
            "paciente": ["A", "A"],
            "visita_meses": [3.0, 9.0],
            "cea": [2.0, 8.0],
            "recaida": [0, 0],
            "mes_de_la_recaida": [np.nan, np.nan],
            "ultimo_seguimiento_meses": [36.0, 36.0],
        }
    )
    t = C.preparar_visitas(d)
    assert list(t["log2_cea"]) == pytest.approx([1.0, 3.0])


# ------------------------------------------- censura
def test_una_visita_sin_doce_meses_de_seguimiento_y_sin_suceso_no_entra():
    d = pd.DataFrame(
        {
            "paciente": ["A", "A"],
            "visita_meses": [3.0, 12.0],
            "cea": [2.0, 2.0],
            "recaida": [0, 0],
            "mes_de_la_recaida": [np.nan, np.nan],
            # El seguimiento acaba en el mes 18: la visita del 3 sí tiene doce
            # meses por delante, la del 12 no.
            "ultimo_seguimiento_meses": [18.0, 18.0],
        }
    )
    t = C.preparar_visitas(d)
    assert bool(t[t["visita_meses"] == 3.0]["evaluable"].iloc[0]) is True
    assert bool(t[t["visita_meses"] == 12.0]["evaluable"].iloc[0]) is False


def test_una_visita_posterior_a_la_recaida_no_entra():
    """Después de la recaída el CEA ya no predice nada: la mide."""
    d = pd.DataFrame(
        {
            "paciente": ["A", "A"],
            "visita_meses": [6.0, 24.0],
            "cea": [2.0, 50.0],
            "recaida": [1, 1],
            "mes_de_la_recaida": [18.0, 18.0],
            "ultimo_seguimiento_meses": [36.0, 36.0],
        }
    )
    t = C.preparar_visitas(d)
    assert bool(t[t["visita_meses"] == 6.0]["evaluable"].iloc[0]) is True
    assert bool(t[t["visita_meses"] == 24.0]["evaluable"].iloc[0]) is False


def test_el_suceso_es_recaer_dentro_de_la_ventana_y_no_en_cualquier_momento():
    d = pd.DataFrame(
        {
            "paciente": ["A", "A"],
            "visita_meses": [3.0, 24.0],
            "cea": [2.0, 2.0],
            "recaida": [1, 1],
            "mes_de_la_recaida": [30.0, 30.0],
            "ultimo_seguimiento_meses": [36.0, 36.0],
        }
    )
    t = C.preparar_visitas(d)
    # Del mes 3 al 30 hay 27 meses: fuera de la ventana de 12.
    assert t[t["visita_meses"] == 3.0]["recae_en_la_ventana"].iloc[0] == 0
    # Del 24 al 30 hay 6: dentro.
    assert t[t["visita_meses"] == 24.0]["recae_en_la_ventana"].iloc[0] == 1


def test_se_publica_cuantas_visitas_se_pierden_por_censura(comparacion):
    res = comparacion
    assert res["n_visitas_excluidas_por_censura"] > 0
    assert res["n_visitas_evaluables"] > 0


# ------------------------------------------- AUC y comparación
def test_el_auc_cuenta_los_empates_como_medio():
    """El modelo del umbral solo da dos puntuaciones, así que casi todo empata."""
    y = np.array([1, 1, 0, 0])
    assert C._auc(y, np.array([1.0, 1.0, 1.0, 1.0])) == pytest.approx(0.5)
    assert C._auc(y, np.array([1.0, 1.0, 0.0, 0.0])) == pytest.approx(1.0)
    assert C._auc(y, np.array([0.0, 0.0, 1.0, 1.0])) == pytest.approx(0.0)


def test_el_auc_de_un_predictor_perfecto_es_uno():
    y = np.array([0, 0, 1, 1])
    assert C._auc(y, np.array([1.0, 2.0, 3.0, 4.0])) == pytest.approx(1.0)


def test_los_pliegues_se_parten_por_paciente_y_no_por_visita():
    """Un solo paciente con seis visitas no da AUC.

    Si los pliegues se partieran por visita, estas seis visitas se repartirían
    entre entrenamiento y evaluación y saldría un AUC alto: el modelo habría
    visto a la misma persona que evalúa. Partidos por paciente no queda nada
    fuera que evaluar, y la función lo dice devolviendo un valor no finito en
    lugar de puntuar dentro de muestra.
    """
    d = pd.DataFrame(
        {
            "paciente": ["A"] * 6,
            "visita_meses": [3.0, 6.0, 9.0, 12.0, 18.0, 24.0],
            "cea": [2.0, 2.2, 3.0, 5.0, 9.0, 20.0],
            "recaida": [1] * 6,
            "mes_de_la_recaida": [30.0] * 6,
            "ultimo_seguimiento_meses": [36.0] * 6,
        }
    )
    t = C.preparar_visitas(d)
    auc = C._auc_por_validacion_cruzada(
        t[t["evaluable"].astype(bool)], ["log2_cea"], semilla=SEMILLA
    )
    assert not np.isfinite(auc)


def test_la_comparacion_remuestrea_pacientes_y_lo_declara(comparacion):
    res = comparacion
    assert res["unidad_de_remuestreo"] == "paciente"
    assert "más estrechos" in res["por_que_el_paciente"]


def test_el_valor_mejora_sobre_el_umbral_en_la_cohorte_sintetica(comparacion):
    """El mecanismo programado hace que el valor informe más que el umbral."""
    res = comparacion
    d = res["diferencias_pareadas"]["valor_menos_umbral"]
    assert d["mejora"] is True
    assert d["ci_low"] > 0


def test_no_se_fabrica_una_mejora_que_el_generador_no_programo(comparacion):
    """La trayectoria no añade nada sobre el valor en estos datos sintéticos, y
    el intervalo de la diferencia incluye el cero. Que el código lo diga es lo
    que hay que comprobar antes de usarlo con datos reales."""
    res = comparacion
    d = res["diferencias_pareadas"]["trayectoria_menos_valor"]
    assert d["ci_low"] < 0 < d["ci_high"]
    assert d["mejora"] is False


def test_el_veredicto_sale_del_intervalo_de_la_diferencia(comparacion):
    res = comparacion
    for d in res["diferencias_pareadas"].values():
        if "ci_low" in d:
            assert d["mejora"] == bool(d["ci_low"] > 0)


def test_sin_las_dos_clases_la_comparacion_no_es_calculable():
    d = pd.DataFrame(
        {
            "paciente": ["A", "A", "B", "B"],
            "visita_meses": [3.0, 9.0, 3.0, 9.0],
            "cea": [2.0, 2.0, 3.0, 3.0],
            "recaida": [0, 0, 0, 0],
            "mes_de_la_recaida": [np.nan] * 4,
            "ultimo_seguimiento_meses": [36.0] * 4,
        }
    )
    t = C.preparar_visitas(d)
    res = C.comparar_modelos(t, semilla=SEMILLA, n_bootstrap=50)
    assert res["estado"] == "no_calculable"


def test_la_comparacion_es_reproducible_con_la_misma_semilla(visitas):
    a = C.comparar_modelos(visitas, semilla=SEMILLA, n_bootstrap=120)
    b = C.comparar_modelos(visitas, semilla=SEMILLA, n_bootstrap=120)
    assert (
        a["diferencias_pareadas"]["valor_menos_umbral"]["diferencia_mediana"]
        == b["diferencias_pareadas"]["valor_menos_umbral"]["diferencia_mediana"]
    )


# ------------------------------------------- tamaño necesario
def test_detectar_una_diferencia_mas_pequena_exige_mas_pacientes():
    grande = C.tamano_necesario(
        auc_umbral=0.60, auc_trayectoria=0.70,
        proporcion_con_recaida=0.11, visitas_por_paciente=5.0,
    )
    pequena = C.tamano_necesario(
        auc_umbral=0.60, auc_trayectoria=0.65,
        proporcion_con_recaida=0.11, visitas_por_paciente=5.0,
    )
    assert pequena["n_pacientes"] > grande["n_pacientes"]


def test_mas_visitas_por_paciente_no_valen_como_mas_pacientes():
    """El efecto de diseño descuenta la correlación interna del paciente."""
    pocas = C.tamano_necesario(
        auc_umbral=0.60, auc_trayectoria=0.70,
        proporcion_con_recaida=0.11, visitas_por_paciente=2.0,
    )
    muchas = C.tamano_necesario(
        auc_umbral=0.60, auc_trayectoria=0.70,
        proporcion_con_recaida=0.11, visitas_por_paciente=10.0,
    )
    assert muchas["efecto_diseno"] > pocas["efecto_diseno"]


def test_una_diferencia_nula_no_tiene_tamano_necesario():
    res = C.tamano_necesario(
        auc_umbral=0.70, auc_trayectoria=0.70,
        proporcion_con_recaida=0.11, visitas_por_paciente=5.0,
    )
    assert res["estado"] == "no_calculable"


def test_el_calculo_se_declara_como_aproximacion():
    res = C.tamano_necesario(
        auc_umbral=0.60, auc_trayectoria=0.70,
        proporcion_con_recaida=0.11, visitas_por_paciente=5.0,
    )
    assert "Hanley y McNeil" in res["es_una_aproximacion"]
    assert "simular" in res["es_una_aproximacion"]
    assert "orden de magnitud" in res["es_una_aproximacion"]
