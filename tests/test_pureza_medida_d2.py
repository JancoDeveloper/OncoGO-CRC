"""Tests de D2: el índice estromal de R4 contra una pureza medida.

Lo que hay que comprobar aquí no es que las correlaciones salgan: es que no se
pueda confundir «no hay relación» con «la medida no da para saberlo», y que los
dos trucos que esta fuente tiene preparados —el código de paciente escrito de
otra forma y las alícuotas repetidas— no pasen desapercibidos.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.v2 import pureza_medida as P


def _matriz(genes: list[str], alicuotas: list[str], semilla: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(semilla)
    return pd.DataFrame(
        rng.normal(25.0, 2.0, size=(len(genes), len(alicuotas))),
        index=genes,
        columns=alicuotas,
    )


def _alicuotas(cols: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "alicuota": cols,
            "paciente": ["TCGA-" + "-".join(c.split("-")[:2]) for c in cols],
        }
    )


# ------------------------------------------- el código de paciente
def test_una_alicuota_sin_su_paciente_detiene_el_cruce():
    """Sin este guarda el cruce con la clínica saldría incompleto en silencio."""
    cols = ["A6-3807-01A-22", "A6-3808-01A-22"]
    m = _matriz(["COL1A1", "FN1", "LUM"], cols)
    al = _alicuotas(cols).iloc[:1]
    with pytest.raises(P.FaltaColumna, match="no están en la tabla de"):
        P.por_paciente(m, al)


def test_el_prefijo_tcga_se_reconstruye_al_promediar():
    cols = ["A6-3807-01A-22", "A6-3808-01A-22"]
    m = _matriz(["COL1A1", "FN1", "LUM"], cols)
    out = P.por_paciente(m, _alicuotas(cols))
    assert list(out.columns) == ["TCGA-A6-3807", "TCGA-A6-3808"]


def test_las_dos_alicuotas_de_un_paciente_se_promedian_y_no_se_cuentan_dos_veces():
    """Dejarlas como dos columnas contaría a esa persona dos veces."""
    cols = ["AA-A00A-01A-22", "AA-A00A-01A-41", "A6-3807-01A-22"]
    m = pd.DataFrame(
        [[10.0, 20.0, 5.0]], index=["COL1A1"], columns=cols
    )
    out = P.por_paciente(m, _alicuotas(cols))
    assert out.shape[1] == 2
    assert out.loc["COL1A1", "TCGA-AA-A00A"] == pytest.approx(15.0)


# ------------------------------------------- la correlación
def test_la_correlacion_no_se_calcula_con_pocos_pacientes():
    x = pd.Series(np.arange(10.0))
    y = pd.Series(np.arange(10.0))
    r = P._correlacion(x, y)
    assert r["estado"] == "no_calculable"
    assert str(P.MIN_PARA_CORRELAR) in r["motivo"]


def test_una_relacion_perfecta_da_rho_uno():
    x = pd.Series(np.arange(40.0))
    r = P._correlacion(x, x * 3.0 + 1.0)
    assert r["rho"] == pytest.approx(1.0)
    assert r["varianza_compartida"] == pytest.approx(1.0)


def test_la_varianza_compartida_es_el_cuadrado_de_rho():
    rng = np.random.default_rng(7)
    x = pd.Series(rng.normal(size=60))
    y = pd.Series(x * 0.5 + rng.normal(size=60))
    r = P._correlacion(x, y)
    assert r["varianza_compartida"] == pytest.approx(r["rho"] ** 2)


def test_se_publican_las_dos_correlaciones_para_que_se_vea_la_discrepancia():
    """Si Spearman y Pearson discrepan es que hay extremos mandando."""
    x = pd.Series(list(np.arange(39.0)) + [1000.0])
    y = pd.Series(list(np.arange(39.0)) + [-1000.0])
    r = P._correlacion(x, y)
    assert r["rho"] > 0.8
    assert r["r_pearson"] < 0
    assert "aproximación" in r["intervalo_es_aproximado"]


def test_los_ausentes_no_entran_en_el_recuento():
    x = pd.Series([np.nan] * 10 + list(np.arange(30.0)))
    y = pd.Series(list(np.arange(40.0)))
    assert P._correlacion(x, y)["n"] == 30


# ------------------------------------------- la firma de GO
def test_una_tabla_de_go_sin_las_columnas_necesarias_detiene_el_analisis():
    with pytest.raises(P.FaltaColumna, match="go_id"):
        P.genes_del_termino(pd.DataFrame({"otra": [1]}), P.ECM_GO_TERM)


def test_la_firma_se_lee_en_mayusculas_y_solo_del_termino_pedido():
    go = pd.DataFrame(
        {
            "go_id": ["GO:0031012", "GO:0031012", "GO:0002376"],
            "gene_symbol": ["col1a1", " FN1 ", "PTPRC"],
        }
    )
    assert P.genes_del_termino(go, "GO:0031012") == {"COL1A1", "FN1"}


def test_el_termino_que_se_valida_es_el_que_usa_r4():
    """Validar otra firma no validaría R4."""
    assert P.ECM_GO_TERM == "GO:0031012"


# ------------------------------------------- la validación del índice
def test_sin_pureza_medida_la_etapa_no_tiene_nada_que_validar():
    cols = ["A6-3807-01A-22"] + [f"A6-38{i:02d}-01A-22" for i in range(10, 40)]
    m = _matriz(["COL1A1", "FN1", "LUM", "DCN"], cols)
    al = _alicuotas(cols)
    clin = pd.DataFrame({"otra": [1] * len(cols)}, index=al["paciente"])
    with pytest.raises(P.FaltaColumna, match="pureza_absolute"):
        P.validar_el_indice(m, al, clin, {"ecm": {"COL1A1", "FN1", "LUM"}})


def test_una_firma_corta_no_produce_un_indice():
    cols = [f"A6-38{i:02d}-01A-22" for i in range(10, 50)]
    m = _matriz(["COL1A1", "FN1", "LUM", "DCN"], cols)
    al = _alicuotas(cols)
    clin = pd.DataFrame(
        {"pureza_absolute": np.linspace(0.3, 0.9, len(cols))}, index=al["paciente"]
    )
    out = P.validar_el_indice(m, al, clin, {"ecm": {"COL1A1", "SOLO_UNO"}})
    assert out["ecm"]["estado"] == "no_calculable"
    assert str(P.MIN_PROTEINAS_DE_FIRMA) in out["ecm"]["motivo"]


def test_el_signo_esperado_se_declara_y_se_comprueba():
    """Un signo positivo no es un detalle: es que el índice no mide estroma."""
    cols = [f"A6-38{i:02d}-01A-22" for i in range(10, 60)]
    n = len(cols)
    pureza = np.linspace(0.3, 0.9, n)
    al = _alicuotas(cols)
    clin = pd.DataFrame({"pureza_absolute": pureza}, index=al["paciente"])
    # Un índice que sube con la pureza: el signo sale del revés a propósito.
    m = pd.DataFrame(
        np.vstack([pureza * 10.0 + i for i in range(4)]),
        index=["COL1A1", "FN1", "LUM", "DCN"],
        columns=cols,
    )
    out = P.validar_el_indice(m, al, clin, {"ecm": {"COL1A1", "FN1", "LUM", "DCN"}})
    assert out["ecm"]["signo_esperado"] == "negativo"
    assert out["ecm"]["signo_es_el_esperado"] is False


def test_un_indice_que_recupera_la_pureza_sale_con_rho_alto_y_signo_correcto():
    cols = [f"A6-38{i:02d}-01A-22" for i in range(10, 60)]
    n = len(cols)
    pureza = np.linspace(0.3, 0.9, n)
    al = _alicuotas(cols)
    clin = pd.DataFrame({"pureza_absolute": pureza}, index=al["paciente"])
    m = pd.DataFrame(
        np.vstack([-pureza * 10.0 + i for i in range(4)]),
        index=["COL1A1", "FN1", "LUM", "DCN"],
        columns=cols,
    )
    out = P.validar_el_indice(m, al, clin, {"ecm": {"COL1A1", "FN1", "LUM", "DCN"}})
    assert out["ecm"]["rho"] < -0.9
    assert out["ecm"]["signo_es_el_esperado"] is True
    assert out["ecm"]["n_proteinas_de_la_firma"] == 4


# ------------------------------------------- la reproducibilidad
def test_sin_replicas_se_dice_que_no_se_puede_interpretar_nada():
    cols = [f"A6-38{i:02d}-01A-22" for i in range(10, 40)]
    m = _matriz(["COL1A1", "FN1"], cols)
    out = P.reproducibilidad(m, _alicuotas(cols))
    assert out["estado"] == "no_calculable"
    assert "se puede interpretar" in out["consecuencia"]


def test_dos_medidas_identicas_del_mismo_tumor_dan_reproducibilidad_perfecta():
    cols = ["AA-A00A-01A-22", "AA-A00A-01A-41"]
    rng = np.random.default_rng(3)
    v = rng.normal(25.0, 2.0, size=40)
    m = pd.DataFrame({cols[0]: v, cols[1]: v}, index=[f"G{i}" for i in range(40)])
    out = P.reproducibilidad(m, _alicuotas(cols))
    assert out["r_mediana"] == pytest.approx(1.0)
    assert out["n_pares_de_replicas"] == 1


def test_una_proteina_no_medida_no_recibe_un_numero():
    cols = ["AA-A00A-01A-22", "AA-A00A-01A-41"]
    m = _matriz(["COL1A1"], cols)
    out = P.reproducibilidad_de_un_gen(m, _alicuotas(cols), "CEACAM5")
    assert out["estado"] == "no_medida"
    assert "rho" not in out


def test_el_ruido_frente_a_la_senal_distingue_una_medida_buena_de_una_mala():
    """Es lo que decide si un cero significa «no hay relación» o «no se sabe»."""
    cols = ["AA-A00A-01A-22", "AA-A00A-01A-41"] + [
        f"A6-38{i:02d}-01A-22" for i in range(10, 40)
    ]
    al = _alicuotas(cols)
    rng = np.random.default_rng(11)
    # Entre pacientes varía mucho; entre réplicas, casi nada: medida buena.
    buena = np.concatenate([[20.0, 20.1], rng.normal(25.0, 5.0, size=30)])
    # Entre réplicas varía tanto como entre pacientes: medida inservible.
    mala = np.concatenate([[20.0, 28.0], rng.normal(25.0, 4.0, size=30)])
    m = pd.DataFrame([buena, mala], index=["BUENA", "MALA"], columns=cols)
    rb = P.reproducibilidad_de_un_gen(m, al, "BUENA")
    rm = P.reproducibilidad_de_un_gen(m, al, "MALA")
    assert rb["ruido_frente_a_senal"] < 0.1
    assert rm["ruido_frente_a_senal"] > 1.0
    assert "no se puede saber" in rm["como_se_lee"]


def test_la_diferencia_se_traduce_a_veces_porque_el_log2_no_se_lee_solo():
    cols = ["AA-A00A-01A-22", "AA-A00A-01A-41"]
    m = pd.DataFrame([[20.0, 22.0]], index=["G"], columns=cols)
    out = P.reproducibilidad_de_un_gen(m, _alicuotas(cols), "G")
    assert out["diferencia_mediana_log2"] == pytest.approx(2.0)
    assert out["veces_de_diferencia_mediana"] == pytest.approx(4.0)


# ------------------------------------------- el CEA
def test_un_cea_de_cero_no_se_transforma():
    cols = [f"A6-38{i:02d}-01A-22" for i in range(10, 40)]
    al = _alicuotas(cols)
    rng = np.random.default_rng(5)
    cea = rng.uniform(1.0, 30.0, size=len(cols))
    cea[0] = 0.0
    clin = pd.DataFrame(
        {"cea_pretratamiento_ng_ml": cea, "pureza_absolute": 0.6}, index=al["paciente"]
    )
    m = _matriz(["CEACAM5", "CEACAM6", "CEACAM1"], cols)
    out = P.cea_de_proteina(m, al, clin)
    assert out["ceros_de_cea_no_transformados"] == 1
    assert out["n_con_cea_positivo"] == len(cols) - 1


def test_una_proteina_ausente_de_la_matriz_no_recibe_correlacion():
    cols = [f"A6-38{i:02d}-01A-22" for i in range(10, 40)]
    al = _alicuotas(cols)
    clin = pd.DataFrame(
        {"cea_pretratamiento_ng_ml": np.linspace(1.0, 20.0, len(cols))},
        index=al["paciente"],
    )
    m = _matriz(["CEACAM5"], cols)
    out = P.cea_de_proteina(m, al, clin)
    assert out["proteina"]["CEACAM6"]["estado"] == "no_medida"
    assert "no medida no es lo mismo" in out["proteina"]["CEACAM6"]["nota"]
    assert "rho" not in out["proteina"]["CEACAM6"]


def test_la_proteina_y_el_transcrito_se_comparan_sobre_los_mismos_pacientes():
    """Comparados sobre subconjuntos distintos no dirían cuál informa más."""
    cols = [f"A6-38{i:02d}-01A-22" for i in range(10, 60)]
    al = _alicuotas(cols)
    n = len(cols)
    rng = np.random.default_rng(13)
    cea = rng.uniform(1.0, 40.0, size=n)
    rna = rng.uniform(10.0, 2000.0, size=n)
    # El transcrito falta en diez pacientes: sin el cruce, los dos coeficientes
    # se calcularían sobre grupos distintos.
    rna[:10] = np.nan
    clin = pd.DataFrame(
        {"cea_pretratamiento_ng_ml": cea, "rnaseq_CEACAM5": rna}, index=al["paciente"]
    )
    m = _matriz(["CEACAM5"], cols)
    out = P.cea_de_proteina(m, al, clin)
    c = out["proteina_frente_a_transcrito"]
    assert c["n_con_las_tres_medidas"] == n - 10
    assert c["proteina_vs_cea"]["n"] == c["transcrito_vs_cea"]["n"]
    assert c["la_comparacion_es_sobre_los_mismos_pacientes"] is True


def test_comparar_dos_correlaciones_no_es_contrastarlas():
    cols = [f"A6-38{i:02d}-01A-22" for i in range(10, 60)]
    al = _alicuotas(cols)
    rng = np.random.default_rng(17)
    clin = pd.DataFrame(
        {
            "cea_pretratamiento_ng_ml": rng.uniform(1.0, 40.0, size=len(cols)),
            "rnaseq_CEACAM5": rng.uniform(10.0, 2000.0, size=len(cols)),
        },
        index=al["paciente"],
    )
    out = P.cea_de_proteina(_matriz(["CEACAM5"], cols), al, clin)
    assert "no_es_un_contraste" in out["proteina_frente_a_transcrito"]
    assert "potencia" in out["proteina_frente_a_transcrito"]["no_es_un_contraste"]


def test_sin_la_columna_del_cea_el_analisis_se_detiene():
    cols = [f"A6-38{i:02d}-01A-22" for i in range(10, 40)]
    al = _alicuotas(cols)
    clin = pd.DataFrame({"otra": [1] * len(cols)}, index=al["paciente"])
    with pytest.raises(P.FaltaColumna, match="cea_pretratamiento_ng_ml"):
        P.cea_de_proteina(_matriz(["CEACAM5"], cols), al, clin)
