"""Pruebas de D5: deconvolución con una referencia de composición externa.

Lo que hay que comprobar aquí no es que los mínimos cuadrados no negativos
funcionen —eso lo garantiza scipy— sino que el código se niegue a dar por bueno
un resultado que no lo es. La primera versión del diagnóstico miraba solo el
máximo por tipo celular y declaraba publicables unas proporciones con un residuo
de 0,9 y linfocitos T ausentes en 34 de 44 muestras. Varios de estos tests
existen para que eso no vuelva a pasar.

Los datos de este fichero son sintéticos y se quedan aquí: ninguna cifra suya
entra en los resultados.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.v2 import deconvolucion as D

SEMILLA = 20260914 + 415


@pytest.fixture
def firma_limpia() -> pd.DataFrame:
    """Tres tipos celulares con marcadores que no comparten.

    Es el caso favorable: cada tipo tiene sus genes y nadie más los expresa. Si
    el método no recupera las proporciones aquí, no las recupera en ningún sitio.
    """
    genes = [f"G{i:04d}" for i in range(300)]
    a = np.zeros((300, 3))
    for j in range(3):
        a[j * 100 : (j + 1) * 100, j] = 100.0
    return pd.DataFrame(a, index=genes, columns=["epitelio", "estroma", "inmune"])


@pytest.fixture
def firma_colineal() -> pd.DataFrame:
    """Dos tipos con perfiles casi idénticos, más uno distinto.

    Es el caso de los fibroblastos y el músculo liso: comparten las proteínas de
    matriz y la deconvolución no puede separarlos.
    """
    rng = np.random.default_rng(SEMILLA)
    genes = [f"G{i:04d}" for i in range(300)]
    base = rng.lognormal(2.0, 1.0, size=300)
    estroma = base * 10.0
    return pd.DataFrame(
        {
            "epitelio": np.r_[rng.lognormal(4.0, 0.5, 150), np.full(150, 0.1)],
            "fibroblastos": estroma,
            "musculo_liso": estroma * rng.lognormal(0.0, 0.05, size=300),
        },
        index=genes,
    )


def _mezclar(firma: pd.DataFrame, proporciones: dict[str, float], n: int, semilla: int):
    """Construye muestras mezclando los perfiles de la referencia."""
    rng = np.random.default_rng(semilla)
    a = firma.to_numpy(dtype=float)
    a = a / a.sum(axis=0, keepdims=True)
    x = np.array([proporciones[c] for c in firma.columns], dtype=float)
    x = x / x.sum()
    cols = {}
    for i in range(n):
        b = a @ x
        cols[f"M{i:02d}"] = b * rng.lognormal(0.0, 0.1, size=b.shape)
    return pd.DataFrame(cols, index=firma.index)


# --- la regla de marcadores -------------------------------------------------


def test_un_gen_expresado_igual_en_todos_los_tipos_no_entra():
    """No distingue nada, así que solo añadiría ruido al sistema."""
    from scripts.download.hpa_single_cell import seleccionar_marcadores

    firma = pd.DataFrame(
        {
            "a": [100.0] * 200 + [1.0] * 200,
            "b": [100.0] * 200 + [500.0] * 200,
            "c": [100.0] * 200 + [1.0] * 200,
        },
        index=[f"G{i}" for i in range(400)],
    )
    marcadores, info = seleccionar_marcadores(firma)
    assert "G0" not in marcadores.index
    assert "G200" in marcadores.index
    assert info["n_marcadores"] == 200


def test_un_gen_apenas_expresado_no_entra_por_tener_un_cociente_alto():
    """0,5 frente a 0,01 es un cociente de cincuenta y no significa nada.

    Sin el mínimo de expresión absoluta, el ruido de los genes casi no
    expresados entraría en la firma con la especificidad más alta de todas.
    """
    from scripts.download.hpa_single_cell import seleccionar_marcadores

    firma = pd.DataFrame(
        {
            "a": [0.5] * 200 + [1.0] * 200,
            "b": [0.01] * 200 + [500.0] * 200,
            "c": [0.01] * 200 + [1.0] * 200,
        },
        index=[f"G{i}" for i in range(400)],
    )
    marcadores, _ = seleccionar_marcadores(firma)
    assert "G0" not in marcadores.index
    assert "G200" in marcadores.index


def test_una_firma_sin_marcadores_suficientes_se_rechaza():
    from scripts.download.hpa_single_cell import seleccionar_marcadores

    firma = pd.DataFrame(
        {"a": [100.0, 1.0], "b": [1.0, 100.0]}, index=["G0", "G1"]
    )
    with pytest.raises(RuntimeError, match="menos de cien"):
        seleccionar_marcadores(firma)


# --- la mezcla --------------------------------------------------------------


def test_recupera_proporciones_conocidas_cuando_los_tipos_son_distinguibles(
    firma_limpia,
):
    res = D.recuperacion_en_mezclas_conocidas(
        firma_limpia,
        np.array([[0.6, 0.3, 0.1], [0.2, 0.2, 0.6]]),
        semilla=SEMILLA,
    )
    assert res["error_total_maximo"] < 0.01
    assert res["datos_sinteticos"] is True


def test_no_recupera_proporciones_cuando_dos_tipos_son_casi_iguales(firma_colineal):
    """El error se concentra en el reparto entre los dos tipos confundidos.

    Esto es lo que ocurre con fibroblastos y músculo liso en los datos reales, y
    es la razón de que las proporciones tipo a tipo no se publiquen.
    """
    res = D.recuperacion_en_mezclas_conocidas(
        firma_colineal,
        np.array([[0.5, 0.25, 0.25]]),
        semilla=SEMILLA,
        ruido_log=0.1,
    )
    assert res["error_total_mediano"] > 0.05


def test_un_tipo_con_perfil_de_mayor_magnitud_no_gana_por_la_escala(firma_limpia):
    """Sin normalizar las columnas, el tipo con genes más abundantes se
    llevaría la proporción por un motivo que no es biológico."""
    firma = firma_limpia.copy()
    firma["estroma"] = firma["estroma"] * 1000.0
    mezcla = _mezclar(
        firma_limpia, {"epitelio": 0.5, "estroma": 0.3, "inmune": 0.2}, 25, SEMILLA
    )
    fr = D.fracciones(mezcla, list(mezcla.columns), firma)
    assert abs(float(fr.tabla["estroma"].median()) - 0.3) < 0.05


def test_sin_marcadores_medidos_no_se_devuelve_ninguna_proporcion(firma_limpia):
    mezcla = _mezclar(
        firma_limpia, {"epitelio": 0.5, "estroma": 0.3, "inmune": 0.2}, 25, SEMILLA
    )
    mezcla.index = [f"OTRO{i}" for i in range(len(mezcla))]
    with pytest.raises(ValueError, match="marcadores"):
        D.fracciones(mezcla, list(mezcla.columns), firma_limpia)


def test_una_columna_que_no_existe_se_dice_en_vez_de_saltarse(firma_limpia):
    mezcla = _mezclar(firma_limpia, {"epitelio": 1.0, "estroma": 1.0, "inmune": 1.0}, 25, SEMILLA)
    with pytest.raises(D.FaltaColumna):
        D.fracciones(mezcla, [*mezcla.columns, "MUESTRA_QUE_NO_ESTA"], firma_limpia)


# --- el diagnóstico ---------------------------------------------------------


def test_un_residuo_alto_impide_publicar_las_proporciones(firma_limpia):
    """Si la mezcla no reproduce el perfil, las proporciones no describen nada.

    La muestra se construye con perfiles que no están en la referencia, así que
    ningún reparto puede explicarla.
    """
    rng = np.random.default_rng(SEMILLA)
    ajeno = pd.DataFrame(
        rng.lognormal(0.0, 2.0, size=(300, 25)),
        index=firma_limpia.index,
        columns=[f"M{i:02d}" for i in range(25)],
    )
    fr = D.fracciones(ajeno, list(ajeno.columns), firma_limpia)
    diag = D.diagnosticar(fr, firma_limpia)
    assert diag["residuo_relativo_mediano"] > D.RESIDUO_ALTO
    assert diag["la_mezcla_no_explica_el_perfil"] is True
    assert diag["proporciones_por_tipo_publicables"] is False
    assert "residuo" in diag["veredicto"].lower()


def test_un_tipo_ausente_en_casi_todas_las_muestras_se_detecta_aunque_el_maximo_no_sea_cero(
    firma_limpia,
):
    """Es el fallo que tenía la primera versión del diagnóstico.

    Mirar solo el máximo por tipo no detecta nada: basta una muestra con señal
    para que el máximo deje de ser cero, y el tipo seguía declarándose estimado
    cuando en realidad estaba descartado en el resto.
    """
    mezcla = _mezclar(
        firma_limpia, {"epitelio": 0.6, "estroma": 0.4, "inmune": 0.0}, 20, SEMILLA
    )
    inmune_solo = _mezclar(
        firma_limpia, {"epitelio": 0.0, "estroma": 0.0, "inmune": 1.0}, 1, SEMILLA
    )
    mezcla["M99"] = inmune_solo.iloc[:, 0]
    fr = D.fracciones(mezcla, list(mezcla.columns), firma_limpia)
    assert float(fr.tabla["inmune"].max()) > 0.9  # el máximo no delata nada
    diag = D.diagnosticar(fr, firma_limpia)
    assert "inmune" in diag["tipos_ausentes_en_demasiadas_muestras"]
    assert diag["proporciones_por_tipo_publicables"] is False


def test_dos_tipos_con_perfiles_casi_iguales_se_declaran_no_separables(
    firma_colineal,
):
    mezcla = _mezclar(
        firma_colineal,
        {"epitelio": 0.5, "fibroblastos": 0.25, "musculo_liso": 0.25},
        25,
        SEMILLA,
    )
    fr = D.fracciones(mezcla, list(mezcla.columns), firma_colineal)
    diag = D.diagnosticar(fr, firma_colineal)
    pares = {
        tuple(sorted((p["tipo_a"], p["tipo_b"]))) for p in diag["pares_no_separables"]
    }
    assert ("fibroblastos", "musculo_liso") in pares
    assert diag["proporciones_por_tipo_publicables"] is False


def test_las_proporciones_se_pueden_publicar_cuando_las_tres_comprobaciones_pasan(
    firma_limpia,
):
    """La dirección contraria. Un diagnóstico que siempre dijera «no publicable»
    no estaría comprobando nada."""
    mezcla = _mezclar(
        firma_limpia, {"epitelio": 0.5, "estroma": 0.3, "inmune": 0.2}, 25, SEMILLA
    )
    diag = D.diagnosticar(
        D.fracciones(mezcla, list(mezcla.columns), firma_limpia), firma_limpia
    )
    assert diag["la_mezcla_no_explica_el_perfil"] is False
    assert diag["pares_no_separables"] == []
    assert diag["tipos_ausentes_en_demasiadas_muestras"] == {}
    assert diag["proporciones_por_tipo_publicables"] is True


@pytest.mark.parametrize(
    "proporciones",
    [
        {"epitelio": 0.5, "estroma": 0.3, "inmune": 0.2},
        {"epitelio": 0.9, "estroma": 0.1, "inmune": 0.0},
        {"epitelio": 0.34, "estroma": 0.33, "inmune": 0.33},
    ],
)
def test_el_texto_del_veredicto_siempre_coincide_con_la_bandera(
    firma_limpia, proporciones
):
    """Que el texto y la bandera puedan discrepar es el modo de fallo peligroso:
    el informe diría una cosa y la tabla otra."""
    mezcla = _mezclar(firma_limpia, proporciones, 25, SEMILLA)
    diag = D.diagnosticar(
        D.fracciones(mezcla, list(mezcla.columns), firma_limpia), firma_limpia
    )
    dice_que_no = "**no** se publican" in diag["veredicto"]
    assert dice_que_no is (not diag["proporciones_por_tipo_publicables"])


def test_los_umbrales_del_diagnostico_se_declaran_como_enmienda(firma_limpia):
    mezcla = _mezclar(firma_limpia, {"epitelio": 0.5, "estroma": 0.3, "inmune": 0.2}, 25, SEMILLA)
    diag = D.diagnosticar(
        D.fracciones(mezcla, list(mezcla.columns), firma_limpia), firma_limpia
    )
    assert "E9" in diag["los_umbrales_no_estaban_preregistrados"]


# --- el eje epitelial -------------------------------------------------------


def test_el_agregado_epitelial_no_cambia_si_se_intercambian_dos_tipos_estromales(
    firma_colineal,
):
    """Es la razón de ser del módulo.

    Que la deconvolución no sepa si una proteína de matriz viene de un
    fibroblasto o de una célula de músculo liso hace que las dos proporciones
    sean inestables. La suma epitelial no se mueve, porque ninguno de los dos es
    epitelio.
    """
    mezcla = _mezclar(
        firma_colineal,
        {"epitelio": 0.5, "fibroblastos": 0.3, "musculo_liso": 0.2},
        25,
        SEMILLA,
    )
    fr = D.fracciones(mezcla, list(mezcla.columns), firma_colineal)
    intercambiada = firma_colineal.rename(
        columns={"fibroblastos": "musculo_liso", "musculo_liso": "fibroblastos"}
    )[firma_colineal.columns]
    fr2 = D.fracciones(mezcla, list(mezcla.columns), intercambiada)
    epi1 = D.eje_epitelial(fr, ["epitelio"])
    epi2 = D.eje_epitelial(fr2, ["epitelio"])
    assert float((epi1 - epi2).abs().max()) < 1e-6


def test_sin_ningun_tipo_epitelial_no_se_devuelve_un_eje_vacio(firma_limpia):
    mezcla = _mezclar(firma_limpia, {"epitelio": 0.5, "estroma": 0.3, "inmune": 0.2}, 25, SEMILLA)
    fr = D.fracciones(mezcla, list(mezcla.columns), firma_limpia)
    with pytest.raises(D.FaltaColumna):
        D.eje_epitelial(fr, ["un_tipo_que_no_esta"])


def test_excluir_un_gen_de_la_firma_cambia_el_eje_solo_si_ese_gen_estaba_dentro(
    firma_limpia,
):
    """La exclusión es la misma cautela de C7, y tiene que costar lo que cuesta.

    Para una proteína que no está en la firma, recalcular daría el mismo número,
    y por eso la etapa se ahorra el cálculo. Este test comprueba que ese atajo
    es correcto y no una suposición.
    """
    mezcla = _mezclar(firma_limpia, {"epitelio": 0.5, "estroma": 0.3, "inmune": 0.2}, 25, SEMILLA)
    cols = list(mezcla.columns)
    base = D.eje_epitelial(D.fracciones(mezcla, cols, firma_limpia), ["epitelio"])
    dentro = D.eje_epitelial_sin_gen(mezcla, cols, firma_limpia, "G0000", ["epitelio"])
    fuera = D.eje_epitelial_sin_gen(
        mezcla, cols, firma_limpia, "UN_GEN_QUE_NO_ESTA", ["epitelio"]
    )
    assert float((base - fuera).abs().max()) == 0.0
    assert float((base - dentro).abs().max()) > 0.0


def test_el_eje_se_usa_como_covariable_si_ordena_los_pares_y_coincide_con_el_indice():
    epi = pd.Series({f"T{i}": 0.2 for i in range(10)} | {f"N{i}": 0.5 for i in range(10)})
    pares = [(f"T{i}", f"N{i}") for i in range(10)]
    v = D.veredicto_del_eje_epitelial(
        epi, pares, {"coincide_con_lo_esperado": True}, {"residuo_relativo_mediano": 0.9}
    )
    assert v["pares_con_menos_epitelio_en_el_tumor"] == 10
    assert v["usable_como_covariable"] is True
    assert "orden" in v["no_es_una_proporcion_de_celulas"]


def test_el_eje_no_se_usa_si_no_ordena_los_pares_en_el_sentido_esperado():
    """La dirección contraria, que es la que importa: sin este test el filtro
    podría estar aprobando siempre."""
    epi = pd.Series(
        {f"T{i}": (0.6 if i < 6 else 0.2) for i in range(10)}
        | {f"N{i}": 0.4 for i in range(10)}
    )
    pares = [(f"T{i}", f"N{i}") for i in range(10)]
    v = D.veredicto_del_eje_epitelial(
        epi, pares, {"coincide_con_lo_esperado": True}, {"residuo_relativo_mediano": 0.9}
    )
    assert v["usable_como_covariable"] is False
    assert "No se usa" in v["por_que_se_puede_usar_igual"]


def test_el_eje_no_se_usa_si_no_coincide_con_el_indice_estromal():
    epi = pd.Series({f"T{i}": 0.2 for i in range(10)} | {f"N{i}": 0.5 for i in range(10)})
    pares = [(f"T{i}", f"N{i}") for i in range(10)]
    v = D.veredicto_del_eje_epitelial(
        epi, pares, {"coincide_con_lo_esperado": False}, {"residuo_relativo_mediano": 0.9}
    )
    assert v["usable_como_covariable"] is False


def test_el_nivel_absoluto_del_eje_nunca_se_presenta_como_un_porcentaje_de_celulas():
    epi = pd.Series({f"T{i}": 0.2 for i in range(10)} | {f"N{i}": 0.5 for i in range(10)})
    pares = [(f"T{i}", f"N{i}") for i in range(10)]
    v = D.veredicto_del_eje_epitelial(
        epi, pares, {"coincide_con_lo_esperado": True}, {"residuo_relativo_mediano": 0.93}
    )
    assert "tanto por ciento" in v["no_es_una_proporcion_de_celulas"]
    assert "0.93" in v["no_es_una_proporcion_de_celulas"]


# --- la comparación con el índice del propio proteoma -----------------------


def test_con_pocas_muestras_no_se_correlaciona_nada():
    a = pd.Series(np.linspace(0.1, 0.5, 10), index=[f"M{i}" for i in range(10)])
    b = pd.Series(np.linspace(0.5, 0.1, 10), index=[f"M{i}" for i in range(10)])
    res = D.comparar_con_indice(a, b)
    assert res["estado"] == "no_calculable"
    assert res["n"] == 10


def test_la_coincidencia_con_el_indice_no_se_presenta_como_una_validacion():
    """Las dos medidas salen del mismo proteoma. Que coincidan no demuestra que
    ninguna sea correcta."""
    rng = np.random.default_rng(SEMILLA)
    x = rng.normal(size=40)
    a = pd.Series(x, index=[f"M{i}" for i in range(40)])
    b = pd.Series(-x + rng.normal(0, 0.2, 40), index=[f"M{i}" for i in range(40)])
    res = D.comparar_con_indice(a, b)
    assert res["coincide_con_lo_esperado"] is True
    assert "no demuestra" in res["no_es_una_validacion"]
    assert "ABSOLUTE" in res["no_es_una_validacion"]


# --- por qué no en una matriz TMT -------------------------------------------


def test_un_log_cociente_no_admite_deconvolucion_y_se_comprueba(firma_limpia):
    """Con valores negativos no existe ninguna proporción que reproduzca el
    perfil, porque A·x con pesos no negativos solo da valores no negativos."""
    rng = np.random.default_rng(SEMILLA)
    log_ratio = pd.DataFrame(
        rng.normal(0.0, 1.0, size=(300, 20)),
        index=firma_limpia.index,
        columns=[f"C{i:02d} Unshared Log Ratio" for i in range(20)],
    )
    res = D.por_que_no_en_tmt(log_ratio, list(log_ratio.columns), firma_limpia)
    assert res["n_valores_negativos"] > 0
    assert 0.4 < res["fraccion_negativa"] < 0.6
    assert "negativo" in res["primer_motivo"]
    assert "logaritmo no es lineal" in res["segundo_motivo"]
    assert "canal de referencia" in res["tercer_motivo"]


def test_la_alternativa_para_cptac_se_dice_y_es_el_lenguaje_del_preregistro(
    firma_limpia,
):
    rng = np.random.default_rng(SEMILLA)
    log_ratio = pd.DataFrame(
        rng.normal(0.0, 1.0, size=(300, 20)),
        index=firma_limpia.index,
        columns=[f"C{i:02d}" for i in range(20)],
    )
    res = D.por_que_no_en_tmt(log_ratio, list(log_ratio.columns), firma_limpia)
    assert "cota del efecto de composición" in res["que_se_hace_en_su_lugar"]
    assert "no una pureza medida" in res["que_se_hace_en_su_lugar"]


def test_sin_columnas_de_la_matriz_tmt_se_dice_en_vez_de_devolver_un_numero(
    firma_limpia,
):
    log_ratio = pd.DataFrame(
        np.zeros((300, 2)), index=firma_limpia.index, columns=["A", "B"]
    )
    with pytest.raises(D.FaltaColumna):
        D.por_que_no_en_tmt(log_ratio, ["C", "D"], firma_limpia)


# --- los patrones de columna de cada cohorte --------------------------------


def test_el_tumor_va_primero_en_cada_par_y_es_M_en_una_cohorte_y_T_en_la_otra(
    tmp_path,
):
    """Usar el patrón de una cohorte en la otra deja los pares vacíos o, peor,
    los invierte y el contraste sale con el signo del revés.

    En PXD034575 el tumor es M (maligno), el normal N y el pólipo P. En
    PXD019504 el tumor es T. La P de pólipo no forma par en ninguna de las dos:
    el contraste de pólipos es otro análisis.
    """
    from oncogo_crc.stages.v2_deconvolucion import _pxd019504, _pxd034575

    replica = tmp_path / "pride_crc_matrix.tsv"
    pd.DataFrame(
        {
            "protein_id": ["P04637"],
            "LFQ intensity CRC02M": [10.0],
            "LFQ intensity CRC02N": [20.0],
            "LFQ intensity CRC02P": [30.0],
            "LFQ intensity CRC03T": [40.0],
        }
    ).to_csv(replica, sep="\t", index=False)
    _mat, pares = _pxd034575(replica, {"P04637": "TP53"})
    assert pares == [("LFQ intensity CRC02M", "LFQ intensity CRC02N")]

    independiente = tmp_path / "proteinGroups.txt"
    pd.DataFrame(
        {
            "Gene names": ["TP53"],
            "LFQ intensity CRC10T": [10.0],
            "LFQ intensity CRC10N": [20.0],
            "LFQ intensity CRC11M": [30.0],
        }
    ).to_csv(independiente, sep="\t", index=False)
    _mat2, pares2 = _pxd019504(independiente)
    assert pares2 == [("LFQ intensity CRC10T", "LFQ intensity CRC10N")]


def test_una_accesion_sin_simbolo_de_gen_no_se_inventa(tmp_path):
    """La réplica publica accesiones de UniProt y la firma está indexada por
    gen. Las accesiones que el mapa no conoce se pierden, y perderlas es lo
    correcto: asignarlas a un gen cualquiera metería la proteína equivocada en
    la firma.
    """
    from oncogo_crc.stages.v2_deconvolucion import _pxd034575

    p = tmp_path / "pride_crc_matrix.tsv"
    pd.DataFrame(
        {
            "protein_id": ["P04637", "UNA_ACCESION_DESCONOCIDA"],
            "LFQ intensity CRC02M": [10.0, 99.0],
            "LFQ intensity CRC02N": [20.0, 99.0],
        }
    ).to_csv(p, sep="\t", index=False)
    mat, _ = _pxd034575(p, {"P04637": "TP53"})
    assert list(mat.index) == ["TP53"]


def test_sin_simbolo_de_gen_la_cohorte_independiente_no_se_usa(tmp_path):
    """Una firma de ARN está indexada por gen. Cruzarla con accesiones de
    UniProt devuelve cero coincidencias y parece que no hay marcadores medidos,
    cuando lo que no hay es la misma forma de nombrar las proteínas."""
    from oncogo_crc.stages.v2_deconvolucion import _pxd019504

    p = tmp_path / "proteinGroups.txt"
    pd.DataFrame({"Protein IDs": ["P1"], "LFQ intensity CRC01T": [1.0]}).to_csv(
        p, sep="\t", index=False
    )
    with pytest.raises(KeyError, match="Gene names"):
        _pxd019504(p)
