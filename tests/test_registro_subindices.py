"""Tests de C2: el registro de subíndices es la única fuente de verdad.

En v1, `stages/certification.py` actualizaba el valor, el estado, el nombre y la
nota de un subíndice, pero nunca la definición. Resultado: R4 figuraba como
`calculado` y su definición decía «La pureza tumoral aún no es covariable del
diferencial», y R5 decía «todavía no ejecutadas» con un valor ya calculado.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from oncogo_crc.v2.prereg import FaltaEnPreregistro, cargar

ROOT = Path(__file__).resolve().parents[1]
CODIGOS = {
    "R1_direccion",
    "R2_detectabilidad",
    "R3_homogeneidad",
    "R4_composicion",
    "R5_nulo",
}


@pytest.fixture(scope="module")
def prereg():
    return cargar(ROOT)


def test_el_registro_cubre_los_cinco_subindices(prereg):
    assert {s["codigo"] for s in prereg.subindices()} == CODIGOS


def test_ninguna_definicion_dice_que_no_se_calculo(prereg):
    """La regla que v1 rompía."""
    prohibidas = prereg.get("certificado.frases_prohibidas_si_calculado")
    assert prohibidas, "el preregistro debe declarar las frases prohibidas"
    for sub in prereg.subindices():
        texto = f"{sub['definicion']} {sub.get('nota_obligatoria', '')}".lower()
        for frase in prohibidas:
            assert frase.lower() not in texto, (
                f"{sub['codigo']}: la definición contiene {frase!r}, "
                "que es lo que hacía que R4 dijera 'calculado' y 'aún no es covariable' a la vez"
            )


def test_cada_subindice_nombra_su_numerador_y_su_denominador(prereg):
    for sub in prereg.subindices():
        assert "numerador" in sub, f"{sub['codigo']} no declara numerador"
        assert sub.get("denominador"), f"{sub['codigo']} no declara denominador"
        assert sub.get("unidad_denominador"), f"{sub['codigo']} no declara unidad"


def test_el_denominador_no_se_mide_en_filas(prereg):
    """702 filas del cruce vienen de 687 proteínas. La unidad tiene que decirlo."""
    permitidas = {"proteinas_unicas", "permutaciones", "terminos", "pacientes", "pares"}
    for sub in prereg.subindices():
        assert sub["unidad_denominador"] in permitidas
        assert "fila" not in sub["unidad_denominador"]


def test_cada_subindice_lleva_nota_obligatoria_y_metodo(prereg):
    for sub in prereg.subindices():
        assert sub.get("nota_obligatoria", "").strip(), f"{sub['codigo']} sin nota"
        assert sub.get("metodo", "").strip(), f"{sub['codigo']} sin método"


def test_r4_dice_que_es_una_cota_y_no_una_pureza(prereg):
    r4 = prereg.subindice("R4_composicion")
    texto = f"{r4['definicion']} {r4['nota_obligatoria']}".lower()
    assert "cota" in texto
    assert "no una pureza medida" in texto or "no es una pureza" in texto


def test_r5_declara_el_piso_del_p(prereg):
    r5 = prereg.subindice("R5_nulo")
    assert "1 / (b + 1)".replace(" ", "") in r5["nota_obligatoria"].lower().replace(" ", "")


def test_r1_obliga_a_declarar_la_cobertura(prereg):
    r1 = prereg.subindice("R1_direccion")
    assert "cobertura" in r1["nota_obligatoria"].lower()
    assert prereg.get("c15_denominador_r1.declarar_cobertura") is True
    assert prereg.get("c15_denominador_r1.cobertura_minima_para_generalizar") == 0.50


# --------------------------------------------------- preregistro en general
def test_un_umbral_que_no_esta_en_el_preregistro_falla(prereg):
    with pytest.raises(FaltaEnPreregistro):
        prereg.get("c4_vias.umbral_que_nadie_declaro")


def test_las_semillas_se_derivan_de_la_base_declarada(prereg):
    base = prereg.get("semillas.base")
    assert base == 20260914
    # La permutación de CPTAC conserva la semilla de v1: desplazamiento cero.
    assert prereg.semilla("r5_permutacion_cptac") == base
    assert prereg.semilla("c5_mitades_divididas") == base + 205
    with pytest.raises(FaltaEnPreregistro):
        prereg.semilla("analisis_inventado")


def test_c1_sube_las_permutaciones_y_enumera_la_replica(prereg):
    assert prereg.get("c1_permutaciones.cptac_n_permutaciones") == 10000
    assert prereg.get("c1_permutaciones.pride_exacto") is True
    assert prereg.get("c1_permutaciones.pride_n_inversiones") == 2**16


def test_c4_exige_diez_terminos_para_evaluar_la_clausula_de_via(prereg):
    """v1 declaró 'sostiene' con 4 términos e IC de 0,40 a 1,0."""
    assert prereg.get("c4_vias.n_min_terminos") == 10
    assert prereg.get("c4_vias.estado_si_insuficiente") == "no_evaluable"
    assert prereg.get("hipotesis.h1.clausula_via.requiere_n_min_terminos") == 10


def test_h2_queda_no_evaluada_sin_nulo_de_grados(prereg):
    assert prereg.get("hipotesis.h2.estado") == "no_evaluada"


def test_el_preregistro_prohibe_datos_sinteticos_en_resultados(prereg):
    assert prereg.get("integridad.prohibido_sintetico_en_resultados") is True
    assert prereg.get("integridad.sintetico_permitido_solo_en") == "tests/"
    assert prereg.get("datos_personales.permitidos") is False
