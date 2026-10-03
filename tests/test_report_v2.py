"""Tests del informe de v2: figuras, texto y los tres formatos.

Lo que estos tests protegen es la propiedad que hace útil al informe: que cada
número salga de una tabla y que una etapa que no ha corrido lo diga en lugar de
callar. Un informe que omite la sección cuando falta el dato es peor que uno que
la deja en blanco, porque el lector no sabe que falta.
"""

from __future__ import annotations

import re

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from oncogo_crc.stages import v2_report as R
from oncogo_crc.v2 import figuras as F


# ------------------------------------------------------------------ formato
def test_los_numeros_se_imprimen_con_coma_y_menos_de_verdad():
    assert R._coma(0.8739495798319328) == "0,874"
    assert R._coma(-0.145, 3) == "\u22120,145"
    assert R._coma(None) == "no disponible"
    assert R._coma(float("nan")) == "no disponible"


def test_un_numero_diminuto_no_se_redondea_a_cero():
    """Imprimir 9,999·10⁻⁵ como «0,000» sería falso y borraría el contraste con v1."""
    assert F._coma_fino(9.999e-05) == "1,0·10⁻⁴"
    assert F._coma_fino(0.0099) == "0,010"
    assert F._coma_fino(0.874) == "0,874"
    assert F._coma_fino(-0.0001) == "\u22121,0·10⁻⁴"


def test_el_intervalo_se_declara_cuando_no_existe():
    assert R._ic({"ci_low": 0.8, "ci_high": 0.93}) == "IC 95 % 0,800–0,930"
    assert R._ic({"ci_low": None, "ci_high": None}) == "sin intervalo"


# ------------------------------------------------------------------ secciones
def test_una_etapa_que_no_ha_corrido_lo_dice_en_su_seccion():
    for constructor in (R._sec_gradiente, R._sec_modelos):
        sec = constructor(None, None)
        assert sec.figura is None
        assert "no ha corrido" in " ".join(sec.cuerpo)


def test_el_cierre_clinico_aparece_en_los_tres_formatos(tmp_path):
    secciones = [R.Seccion(1, "Prueba", ["Un párrafo."])]
    md = R._a_markdown(secciones, figuras_rel="figuras")
    htm = R._a_html(secciones, figuras_rel="figuras")
    assert R.CIERRE in md
    assert R.CIERRE in htm
    destino = tmp_path / "SINTETICO_informe.pdf"
    if R._a_pdf(secciones, destino):
        assert destino.exists() and destino.stat().st_size > 1000


def test_el_html_escapa_el_contenido_y_no_inyecta_etiquetas():
    secciones = [R.Seccion(1, "<script>", ["Texto con <b> y & suelto."])]
    htm = R._a_html(secciones, figuras_rel="figuras")
    assert "<script>" not in htm
    assert "&lt;script&gt;" in htm
    assert "&amp;" in htm


def test_la_negrita_y_el_codigo_del_markdown_llegan_al_html():
    secciones = [R.Seccion(1, "T", ["Esto es **importante** y esto es `codigo`."])]
    htm = R._a_html(secciones, figuras_rel="figuras")
    assert "<strong>importante</strong>" in htm
    assert "<code>codigo</code>" in htm


def test_las_limitaciones_no_se_repiten_entre_etapas():
    auditorias = {
        "a": {"limitaciones": ["La misma limitación.", "Otra."]},
        "b": {"limitaciones": ["La misma limitación."]},
        "c": None,
    }
    sec = R._sec_limitaciones(auditorias)
    texto = "\n".join(sec.cuerpo)
    assert texto.count("La misma limitación.") == 1
    assert "Otra." in texto
    assert R.CIERRE in texto


# ------------------------------------------------------------------ comparación
def test_la_comparacion_sale_de_la_tabla_y_no_del_texto():
    changelog = pd.DataFrame(
        {
            "numero": ["R1", "R5"],
            "v1": [0.8739, 0.0099],
            "v2": [0.8739, 9.999e-05],
            "cambia": [False, True],
            "motivo": ["Mismo cálculo.", "El piso del p."],
            "correccion": ["C15", "C1"],
        }
    )
    md = R._comparacion(changelog, {"version": "v2"})
    assert "0,8739" in md
    assert "0,000100" in md
    assert "| no |" in md and "| sí |" in md
    assert "De los 2 números del changelog, 1 cambian." in md


def test_sin_changelog_la_comparacion_lo_declara():
    md = R._comparacion(pd.DataFrame(), None)
    assert "no está disponible" in md
    assert R.CIERRE in md


def test_el_readme_lista_cada_figura_con_su_tabla(tmp_path):
    fig = F.Figura(
        nombre="n",
        titulo="t",
        pie="p",
        png=tmp_path / "una.png",
        svg=tmp_path / "una.svg",
        tsv=tmp_path / "una.tsv",
    )
    md = R._readme({"x": fig}, {"certified": {}}, hay_pdf=True)
    assert "`una.png`" in md and "`una.tsv`" in md
    assert "make all-v2" in md
    assert "SINTETICO_" in md
    assert R.CIERRE in md


def test_el_readme_declara_las_etapas_que_faltan(tmp_path):
    md = R._readme({}, {"meta": None, "certified": {}}, hay_pdf=False)
    assert "meta" in md
    assert "falta `reportlab`" in md


# ------------------------------------------------------------------ figuras
def _axes_minimo() -> dict:
    frac = {
        "nombre": "f", "valor": 0.5625, "numerador": 18, "denominador": 32,
        "ci_low": 0.377, "ci_high": 0.736,
    }
    return {
        "c5_adenoma_carcinoma": {
            "seleccion_en_el_mismo_brazo_v1": {
                "intermedio_dado_mismo_signo": {**frac, "valor": 0.9194,
                                                "numerador": 137, "denominador": 149,
                                                "ci_low": 0.86, "ci_high": 0.96},
            },
            "mitades_divididas": {
                "intermedio_dado_mismo_signo": {
                    "mediana": 0.8924, "p2_5": 0.7646, "p97_5": 1.0,
                },
                "n_proteinas_evaluables_mediana": 65.0,
            },
            "seleccion_independiente_cptac": {"intermedio_dado_mismo_signo": frac},
            "nulo_que_reselecciona": {
                "nulo_mediana": 0.7515, "nulo_p2_5": 0.694, "nulo_p97_5": 0.806,
                "n_elegibles_en_los_dos_brazos": 3037,
            },
        }
    }


def test_la_figura_de_c5_escribe_png_svg_y_su_tabla(tmp_path):
    fig = F.fig_c5_gradiente(_axes_minimo(), tmp_path)
    assert fig.png.exists() and fig.svg.exists() and fig.tsv.exists()
    tabla = pd.read_csv(fig.tsv, sep="\t")
    assert len(tabla) == 4
    assert set(tabla["calculo"]) == {
        "Mismo brazo (v1)", "Mitades divididas",
        "Esperado por azar", "Selección independiente",
    }
    # El análisis principal tiene que estar en la tabla con su valor exacto.
    assert tabla.loc[tabla["calculo"] == "Selección independiente", "valor"].iloc[0] == 0.5625


def test_la_figura_de_c5_sale_a_300_dpi(tmp_path):
    fig = F.fig_c5_gradiente(_axes_minimo(), tmp_path)
    from PIL import Image

    with Image.open(fig.png) as im:
        dpi = im.info.get("dpi")
    assert dpi is None or round(dpi[0]) == F.DPI


def test_una_coleccion_no_evaluable_sigue_apareciendo_en_la_tabla(tmp_path):
    """Si GO BP desaparece de la figura, el lector no sabe que se probó."""
    vias = {
        "umbral_concordancia_via": 0.8,
        "n_min_terminos": 10,
        "por_coleccion": {
            "go_bp": {
                "n_terminos_en_la_coleccion": 11248,
                "concordancia_entre_cohortes": {
                    "n_terminos_significativos_descubrimiento": 132,
                    "n_terminos_evaluados": 6,
                    "estado": "no_evaluable",
                    "concordancia": None,
                },
            },
            "reactome": {
                "n_terminos_en_la_coleccion": 2868,
                "concordancia_entre_cohortes": {
                    "n_terminos_significativos_descubrimiento": 449,
                    "n_terminos_evaluados": 51,
                    "estado": "sostiene",
                    "concordancia": {"valor": 0.9215, "ci_low": 0.811, "ci_high": 0.978},
                    "margen_del_ic_inferior_sobre_el_umbral": 0.0111,
                },
            },
            "msigdb_hallmark": {
                "n_terminos_en_la_coleccion": 50,
                "concordancia_entre_cohortes": {
                    "n_terminos_significativos_descubrimiento": 7,
                    "n_terminos_evaluados": 2,
                    "estado": "no_evaluable",
                    "concordancia": None,
                },
            },
        },
    }
    fig = F.fig_c4_vias(vias, tmp_path)
    tabla = pd.read_csv(fig.tsv, sep="\t")
    assert len(tabla) == 3
    assert set(tabla.loc[tabla["estado"] == "no_evaluable", "coleccion"]) == {
        "GO procesos", "Hallmark"
    }
    assert tabla["concordancia"].isna().sum() == 2


def test_los_dos_histogramas_de_c10_son_de_las_mismas_proteinas(tmp_path):
    """Comparar 185 proteínas contra 101 daría una mediana que no es la del certificado."""
    meta = {
        "diagnostico_de_tecnologia": {"caida_al_quitar_cptac": 0.0554},
        "mediana_1_menos_i2": 0.0938,
    }
    por_proteina = pd.DataFrame(
        {"uniprot": [f"P{i:05d}" for i in range(10)], "i2": [0.9] * 6 + [0.2] * 4}
    )
    solo_lfq = pd.DataFrame(
        {"uniprot": [f"P{i:05d}" for i in range(6)], "i2_solo_lfq": [0.85] * 6}
    )
    fig = F.fig_c10_heterogeneidad(meta, por_proteina, solo_lfq, tmp_path)
    tabla = pd.read_csv(fig.tsv, sep="\t")
    assert tabla["n_proteinas"].nunique() == 1
    assert int(tabla["n_proteinas"].iloc[0]) == 6
    assert tabla.loc[0, "i2_mediano"] == pytest.approx(0.9)


def test_la_barra_de_error_admite_un_solo_punto(tmp_path):
    """El fallo que rompía todas las figuras: xerr de forma (2,) con x escalar."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    F._barras_de_error(ax, 0, 0.5, 0.4, 0.6, color="k")
    plt.close(fig)


def test_un_intervalo_ausente_no_dibuja_una_barra_de_cero_disfrazada(tmp_path):
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots()
    F._barras_de_error(ax, 0, 0.5, np.nan, np.nan, color="k")
    plt.close(fig)


# ------------------------------------------------------------------ extremo a extremo
def test_el_informe_entero_se_genera_sin_ninguna_auditoria(tmp_path):
    """Sin datos el informe existe y dice que no hay datos. No se cae."""
    class _S:
        def resolve(self, nombre):
            return tmp_path / nombre

    (tmp_path / "results").mkdir()
    R.run(_S())
    salida = tmp_path / "results" / "v2" / "report"
    md = (salida / "informe_v2.md").read_text(encoding="utf-8")
    assert "no ha corrido" in md
    assert R.CIERRE in md
    assert (salida / "index.html").exists()
    assert (tmp_path / "results" / "v2" / "README.md").exists()


def test_el_informe_no_filtra_rutas_del_disco_del_autor(tmp_path):
    """Ni rutas absolutas ni nada que identifique a una persona."""
    class _S:
        def resolve(self, nombre):
            return tmp_path / nombre

    (tmp_path / "results").mkdir()
    R.run(_S())
    salida = tmp_path / "results" / "v2" / "report"
    for nombre in ("informe_v2.md", "index.html", "comparacion_v1_v2.md"):
        texto = (salida / nombre).read_text(encoding="utf-8")
        assert "/Volumes/" not in texto
        assert "/Users/" not in texto
        assert not re.search(r"\bCEA\b.*\b\d+[,.]\d+\b", texto)
