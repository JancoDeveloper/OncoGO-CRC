"""Tests de D3: variación dentro de un mismo tumor.

Lo que se comprueba aquí es que las dos trampas del diseño no se cuelen: que el
número de observaciones esté igualado al comparar desviaciones, y que un
porcentaje de concordancia no se publique sin la referencia del azar. Con cuatro
pares, el azar ya pone de acuerdo a una proteína de cada ocho.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.v2 import heterogeneidad as H

SEMILLA = 20260914 + 303
REGIONES = [f"T4_{i}" for i in range(1, 9)]
OTROS = ["T1", "T2", "T3"]


def _matriz(sd_dentro: float, sd_entre: float, n_genes: int = 300, semilla: int = 1):
    """Matriz con la variación interna y la externa fijadas a mano."""
    rng = np.random.default_rng(semilla)
    base = rng.normal(0.0, sd_entre, size=(n_genes, 1 + len(OTROS)))
    datos = {}
    for i, c in enumerate(OTROS):
        datos[c] = base[:, i]
    for c in REGIONES:
        datos[c] = base[:, len(OTROS)] + rng.normal(0.0, sd_dentro, size=n_genes)
    return pd.DataFrame(datos, index=[f"G{i}" for i in range(n_genes)])


# ------------------------------------------- el diseño
def test_una_columna_del_diseno_que_falta_detiene_el_analisis():
    d = _matriz(0.5, 0.5).drop(columns=["T4_8"])
    with pytest.raises(H.FaltaColumna, match="T4_8"):
        H.variacion_dentro_y_entre(
            d,
            regiones_del_mismo_tumor=REGIONES,
            tumores_de_otros_pacientes=OTROS,
            semilla=SEMILLA,
        )


def test_sin_regiones_de_sobra_no_se_puede_igualar_el_numero_de_observaciones():
    """Con cuatro regiones no quedan cuatro distintas que elegir al azar."""
    cuatro = REGIONES[:4]
    d = _matriz(0.5, 0.5)[cuatro + OTROS]
    out = H.variacion_dentro_y_entre(
        d,
        regiones_del_mismo_tumor=cuatro,
        tumores_de_otros_pacientes=OTROS,
        semilla=SEMILLA,
    )
    assert out["estado"] == "no_calculable"
    assert "igualar el número de observaciones" in out["motivo"]


# ------------------------------------------- la comparación igualada
def test_el_numero_de_observaciones_queda_igualado():
    """Es la primera trampa: ocho contra cuatro no es una comparación justa."""
    d = _matriz(0.5, 0.5)
    out = H.variacion_dentro_y_entre(
        d,
        regiones_del_mismo_tumor=REGIONES,
        tumores_de_otros_pacientes=OTROS,
        semilla=SEMILLA,
    )
    assert out["n_pacientes_comparados"] == H.N_IGUALADO == 4


def test_cuando_la_variacion_interna_es_mayor_el_cociente_lo_dice():
    d = _matriz(sd_dentro=1.5, sd_entre=0.4, semilla=2)
    out = H.variacion_dentro_y_entre(
        d,
        regiones_del_mismo_tumor=REGIONES,
        tumores_de_otros_pacientes=OTROS,
        semilla=SEMILLA,
    )
    assert out["cociente_mediano"] > 1.5
    assert out["cociente_ci_low"] > 1.0
    assert out["el_intervalo_excluye_el_uno"] is True
    assert "no son del mismo tamaño" in out["veredicto"]


def test_cuando_la_variacion_interna_es_menor_el_cociente_tambien_lo_dice():
    d = _matriz(sd_dentro=0.2, sd_entre=1.5, semilla=3)
    out = H.variacion_dentro_y_entre(
        d,
        regiones_del_mismo_tumor=REGIONES,
        tumores_de_otros_pacientes=OTROS,
        semilla=SEMILLA,
    )
    assert out["cociente_mediano"] < 0.7
    assert out["cociente_ci_high"] < 1.0


def test_un_cociente_cuyo_intervalo_incluye_el_uno_no_se_declara_resultado():
    """Es la lectura honesta cuando hay un tumor y cuatro pacientes.

    Aquí cada proteína lleva su propia desviación, como en los datos de verdad.
    Eso ensancha el intervalo del cociente lo suficiente para que incluya el uno,
    que es exactamente la situación en la que no se puede declarar un resultado.
    """
    rng = np.random.default_rng(31)
    n_genes = 400
    # Cada proteína tiene su escala, y la interna y la externa son iguales.
    escala = rng.lognormal(mean=-0.5, sigma=1.2, size=n_genes)
    centro = rng.normal(0.0, 1.0, size=n_genes)
    datos = {}
    for c in OTROS:
        datos[c] = centro + rng.normal(0.0, 1.0, size=n_genes) * escala
    base4 = centro + rng.normal(0.0, 1.0, size=n_genes) * escala
    for c in REGIONES:
        datos[c] = base4 + rng.normal(0.0, 1.0, size=n_genes) * escala
    d = pd.DataFrame(datos, index=[f"G{i}" for i in range(n_genes)])
    out = H.variacion_dentro_y_entre(
        d,
        regiones_del_mismo_tumor=REGIONES,
        tumores_de_otros_pacientes=OTROS,
        semilla=SEMILLA,
    )
    assert out["cociente_ci_low"] < 1.0 < out["cociente_ci_high"]
    assert out["el_intervalo_excluye_el_uno"] is False
    assert "una señal, no un resultado establecido" in out["veredicto"]


def test_el_veredicto_concuerda_siempre_con_el_indicador_del_intervalo():
    """La frase y el indicador no pueden decir cosas distintas."""
    for sd_d in (0.2, 0.6, 1.5):
        out = H.variacion_dentro_y_entre(
            _matriz(sd_dentro=sd_d, sd_entre=0.6, semilla=9),
            regiones_del_mismo_tumor=REGIONES,
            tumores_de_otros_pacientes=OTROS,
            semilla=SEMILLA,
        )
        excluye = out["el_intervalo_excluye_el_uno"]
        assert excluye == (
            out["cociente_ci_low"] > 1.0 or out["cociente_ci_high"] < 1.0
        )
        if excluye:
            assert "no son del mismo tamaño" in out["veredicto"]
        else:
            assert "una señal, no un resultado establecido" in out["veredicto"]


def test_el_cociente_se_calcula_emparejado_dentro_de_cada_remuestreo():
    d = _matriz(0.8, 0.5)
    out = H.variacion_dentro_y_entre(
        d,
        regiones_del_mismo_tumor=REGIONES,
        tumores_de_otros_pacientes=OTROS,
        semilla=SEMILLA,
    )
    assert "emparejado" in out["el_cociente_es_pareado"]
    assert "no es el solapamiento" in out["el_cociente_es_pareado"]


def test_la_comparacion_es_reproducible_con_la_misma_semilla():
    d = _matriz(0.8, 0.5)
    kw = dict(
        regiones_del_mismo_tumor=REGIONES, tumores_de_otros_pacientes=OTROS
    )
    a = H.variacion_dentro_y_entre(d, semilla=SEMILLA, **kw)
    b = H.variacion_dentro_y_entre(d, semilla=SEMILLA, **kw)
    assert a["cociente_mediano"] == b["cociente_mediano"]


def test_una_semilla_distinta_no_cambia_la_conclusion():
    d = _matriz(sd_dentro=1.5, sd_entre=0.4, semilla=5)
    kw = dict(
        regiones_del_mismo_tumor=REGIONES, tumores_de_otros_pacientes=OTROS
    )
    a = H.variacion_dentro_y_entre(d, semilla=SEMILLA, **kw)
    b = H.variacion_dentro_y_entre(d, semilla=SEMILLA + 1, **kw)
    assert a["el_intervalo_excluye_el_uno"] == b["el_intervalo_excluye_el_uno"]


def test_se_declara_que_es_un_solo_tumor_de_un_solo_paciente():
    d = _matriz(0.8, 0.5)
    out = H.variacion_dentro_y_entre(
        d,
        regiones_del_mismo_tumor=REGIONES,
        tumores_de_otros_pacientes=OTROS,
        semilla=SEMILLA,
    )
    assert "no cuántos tumores son así" in out["lo_que_no_dice"]


# ------------------------------------------- por proteína
def test_una_proteina_medida_en_muy_pocas_columnas_no_recibe_desviacion():
    d = _matriz(0.5, 0.5, n_genes=50)
    d.loc["G0", REGIONES[2:]] = np.nan
    sd = H._sd(d, REGIONES)
    assert np.isnan(sd["G0"])
    assert np.isfinite(sd["G1"])


def test_las_certificadas_se_describen_aparte_solo_si_hay_suficientes():
    d = _matriz(0.8, 0.5, n_genes=200)
    pocas = H.proteinas_mas_heterogeneas(
        d,
        regiones_del_mismo_tumor=REGIONES,
        tumores_de_otros_pacientes=OTROS,
        certificadas={"G0", "G1"},
    )
    assert pocas["certificadas"]["estado"] == "no_calculable"
    muchas = H.proteinas_mas_heterogeneas(
        d,
        regiones_del_mismo_tumor=REGIONES,
        tumores_de_otros_pacientes=OTROS,
        certificadas={f"G{i}" for i in range(60)},
    )
    assert muchas["certificadas"]["n_presentes"] == 60
    assert "cociente_mediano" in muchas["certificadas"]


def test_describir_dos_medianas_no_es_contrastarlas():
    d = _matriz(0.8, 0.5, n_genes=200)
    out = H.proteinas_mas_heterogeneas(
        d,
        regiones_del_mismo_tumor=REGIONES,
        tumores_de_otros_pacientes=OTROS,
        certificadas={f"G{i}" for i in range(60)},
    )
    assert "no una prueba" in out["certificadas"]["no_es_un_contraste"]


def test_la_tabla_por_proteina_sale_ordenada_y_marca_las_certificadas():
    d = _matriz(0.8, 0.5, n_genes=120)
    out = H.proteinas_mas_heterogeneas(
        d,
        regiones_del_mismo_tumor=REGIONES,
        tumores_de_otros_pacientes=OTROS,
        certificadas={f"G{i}" for i in range(40)},
    )
    t = out["tabla"]
    assert list(t["cociente"]) == sorted(t["cociente"], reverse=True)
    assert int(t["certificada"].sum()) == 40


def test_la_tabla_descriptiva_declara_que_no_iguala_las_observaciones():
    d = _matriz(0.8, 0.5, n_genes=120)
    out = H.proteinas_mas_heterogeneas(
        d,
        regiones_del_mismo_tumor=REGIONES,
        tumores_de_otros_pacientes=OTROS,
        certificadas=set(),
    )
    assert "no está igualado" in out["nota_de_la_region_fija"]


# ------------------------------------------- el contraste pareado
def test_con_cuatro_pares_no_se_calcula_ningun_p():
    d = _matriz(0.5, 0.5, n_genes=100)
    d["NM1"] = d["T1"] - 2.0
    d["NM2"] = d["T2"] - 2.0
    out = H.contraste_pareado(
        d, pares=[("T1", "NM1"), ("T2", "NM2")], nombre="prueba"
    )
    assert "p" not in out
    assert "no se calcula ningún p" in out["sin_inferencia"]


def test_la_concordancia_se_publica_junto_a_lo_que_daria_el_azar():
    """Es la segunda trampa: un 45 % no se puede leer sin saber que el azar da 12 %."""
    d = _matriz(0.5, 0.5, n_genes=200)
    for i, t in enumerate(OTROS + ["T4_1"], start=1):
        d[f"NM{i}"] = d[t] - 3.0
    pares = [(t, f"NM{i}") for i, t in enumerate(OTROS + ["T4_1"], start=1)]
    out = H.contraste_pareado(d, pares=pares, nombre="tumor_vs_mucosa")
    assert out["fraccion_esperada_por_azar"] == pytest.approx(0.125)
    assert out["fraccion_con_todos_los_pares_en_la_misma_direccion"] == 1.0
    assert out["veces_por_encima_del_azar"] == pytest.approx(8.0)


def test_sin_efecto_la_concordancia_se_queda_en_lo_que_da_el_azar():
    rng = np.random.default_rng(21)
    n = 4000
    d = pd.DataFrame(
        {c: rng.normal(size=n) for c in OTROS + ["T4_1", "NM1", "NM2", "NM3", "NM4"]},
        index=[f"G{i}" for i in range(n)],
    )
    pares = [(t, f"NM{i}") for i, t in enumerate(OTROS + ["T4_1"], start=1)]
    out = H.contraste_pareado(d, pares=pares, nombre="nulo")
    assert out["fraccion_con_todos_los_pares_en_la_misma_direccion"] == pytest.approx(
        0.125, abs=0.03
    )
    assert out["veces_por_encima_del_azar"] == pytest.approx(1.0, abs=0.25)


def test_el_contraste_declara_que_no_entra_en_el_meta_analisis():
    d = _matriz(0.5, 0.5, n_genes=100)
    d["NM1"] = d["T1"] - 2.0
    out = H.contraste_pareado(d, pares=[("T1", "NM1")], nombre="prueba")
    assert "10 pacientes" in out["no_entra_en_el_meta_analisis"]


def test_un_par_con_una_columna_inexistente_se_detecta():
    d = _matriz(0.5, 0.5, n_genes=50)
    with pytest.raises(H.FaltaColumna, match="NO_EXISTE"):
        H.contraste_pareado(d, pares=[("T1", "NO_EXISTE")], nombre="prueba")
