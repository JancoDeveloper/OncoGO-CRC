"""Tests del verificador de números (C14).

Un verificador que nunca falla no verifica nada. Estos tests comprueban las dos
mitades: que acepta un número que sale de una tabla y que rechaza uno escrito a
mano. Los ficheros son sintéticos y viven solo en tmp_path.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

RUTA = Path(__file__).resolve().parents[1] / "scripts" / "v2" / "verificar_numeros.py"
_spec = importlib.util.spec_from_file_location("verificar_numeros", RUTA)
vn = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vn)


@pytest.fixture
def raiz(tmp_path):
    base = tmp_path / "results" / "v2" / "certification"
    base.mkdir(parents=True)
    (base / "certificado_v2.json").write_text(
        json.dumps(
            {
                "subindices": {
                    "R1_direccion": {"valor": 0.8739495798319328, "n": 238},
                    "R3_homogeneidad": {"valor": 0.09378802770504324, "n": 101},
                },
                "nota": "181 proteínas significativas en la tercera cohorte.",
            }
        ),
        encoding="utf-8",
    )
    (base / "indice_R_v2.tsv").write_text(
        "subindice\tvalor\tn\nR1\t0.8739495798319328\t238\n", encoding="utf-8"
    )
    return tmp_path / "results" / "v2"


def _fuentes(raiz):
    valores, usadas = vn.recoger_fuentes(raiz, [])
    assert usadas
    return valores


def test_acepta_un_numero_que_sale_de_una_tabla(raiz, tmp_path):
    informe = tmp_path / "SINTETICO_informe.md"
    informe.write_text(
        "R1 vale 0,874 sobre 238 proteínas. R3 vale 0,094 sobre 101.\n", encoding="utf-8"
    )
    assert vn.verificar(informe, _fuentes(raiz)) == []


def test_rechaza_un_numero_escrito_a_mano(raiz, tmp_path):
    informe = tmp_path / "SINTETICO_informe.md"
    informe.write_text("R1 vale 0,812 sobre 238 proteínas.\n", encoding="utf-8")
    fallos = vn.verificar(informe, _fuentes(raiz))
    assert len(fallos) == 1
    assert fallos[0]["numero"] == "0,812"


def test_el_redondeo_fija_la_tolerancia(raiz, tmp_path):
    """0,87 vale porque 0,8739 redondea ahí. 0,8740 con cuatro decimales no."""
    bien = tmp_path / "SINTETICO_bien.md"
    bien.write_text("Sale 0,87.\n", encoding="utf-8")
    assert vn.verificar(bien, _fuentes(raiz)) == []

    mal = tmp_path / "SINTETICO_mal.md"
    mal.write_text("Sale 0,8651.\n", encoding="utf-8")
    assert len(vn.verificar(mal, _fuentes(raiz))) == 1


def test_los_identificadores_no_se_verifican(raiz, tmp_path):
    informe = tmp_path / "SINTETICO_informe.md"
    informe.write_text(
        "PXD019504 y PDC000116 aportan datos. El término GO:0006915 sube.\n"
        "La corrección C11 y la enmienda E7 lo explican. H1 queda rechazada.\n",
        encoding="utf-8",
    )
    assert vn.verificar(informe, _fuentes(raiz)) == []


def test_un_numero_dentro_de_una_cadena_del_json_cuenta_como_fuente(raiz, tmp_path):
    """Los motivos y las notas los escribe el código desde los datos."""
    informe = tmp_path / "SINTETICO_informe.md"
    informe.write_text("La tercera cohorte aporta 181 proteínas.\n", encoding="utf-8")
    assert vn.verificar(informe, _fuentes(raiz)) == []


def test_las_tablas_por_proteina_no_son_fuente(tmp_path):
    """Con cientos de miles de valores cualquier número coincidiría por azar."""
    base = tmp_path / "results" / "v2" / "meta"
    base.mkdir(parents=True)
    (base / "meta_por_proteina.tsv").write_text(
        "uniprot\tvalor\n" + "".join(f"P{i:05d}\t0.{i:04d}\n" for i in range(2000)),
        encoding="utf-8",
    )
    valores, usadas = vn.recoger_fuentes(tmp_path / "results" / "v2", [])
    assert valores == set()
    assert usadas == []


def test_un_numero_dentro_de_un_campo_de_texto_del_tsv_cuenta_como_fuente(tmp_path):
    """Los motivos del changelog los escribe el código desde los datos."""
    base = tmp_path / "results" / "v2" / "certification"
    base.mkdir(parents=True)
    (base / "changelog_numeros.tsv").write_text(
        "numero\tv1\tv2\tmotivo\n"
        "R1\t0.5\t0.6\tv1 tenía tres denominadores: 688, 702 y 244.\n",
        encoding="utf-8",
    )
    valores, _ = vn.recoger_fuentes(tmp_path / "results" / "v2", [])
    assert {688.0, 702.0, 244.0} <= valores

    informe = tmp_path / "SINTETICO_informe.md"
    informe.write_text("Los denominadores eran 688, 702 y 244.\n", encoding="utf-8")
    assert vn.verificar(informe, valores) == []


def test_el_numero_de_filas_de_una_tabla_de_resumen_es_fuente(tmp_path):
    base = tmp_path / "results" / "v2" / "certification"
    base.mkdir(parents=True)
    filas = "".join(f"R{i}\t0.{i}\n" for i in range(1, 12))
    (base / "indice_R_v2.tsv").write_text("codigo\tvalor\n" + filas, encoding="utf-8")
    valores, _ = vn.recoger_fuentes(tmp_path / "results" / "v2", [])

    informe = tmp_path / "SINTETICO_informe.md"
    informe.write_text("La tabla tiene 11 filas.\n", encoding="utf-8")
    assert vn.verificar(informe, valores) == []

    mal = tmp_path / "SINTETICO_mal.md"
    mal.write_text("La tabla tiene 37 filas.\n", encoding="utf-8")
    assert len(vn.verificar(mal, valores)) == 1


def test_el_orden_de_los_pasos_de_una_tabla_no_se_verifica(raiz, tmp_path):
    informe = tmp_path / "SINTETICO_informe.md"
    informe.write_text(
        "| Orden | Comando |\n| --- | --- |\n"
        "| 1 | make uno |\n| 13 | make trece |\n",
        encoding="utf-8",
    )
    assert vn.verificar(informe, _fuentes(raiz)) == []
