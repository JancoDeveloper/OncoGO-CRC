"""El sello exige un árbol limpio y la enmienda va antes de results/v3_2."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from oncogo_crc.utils.sello import ArbolSucio, exigir_arbol_limpio, hacer_sello

RAIZ = Path(__file__).resolve().parents[1]


def test_un_arbol_sucio_no_firma():
    with pytest.raises(ArbolSucio):
        exigir_arbol_limpio({"arbol_limpio": False})


def test_el_sello_trae_las_claves():
    sello = hacer_sello(RAIZ, etapa="test")
    for clave in (
        "commit",
        "arbol_limpio",
        "fecha_utc",
        "sha256_preregistro",
        "sha256_enmiendas",
        "python",
        "paquetes",
        "semilla",
        "etapa",
    ):
        assert clave in sello
    assert sello["semilla"] == 20260914
    assert sello["etapa"] == "test"
    assert sello["sha256_preregistro"]


def test_la_enmienda_es_anterior_a_cualquier_resultado_v3_2():
    enmienda = subprocess.check_output(
        ["git", "log", "-1", "--format=%H", "--grep=^prereg: enmienda v3.2"],
        cwd=RAIZ,
        text=True,
    ).strip()
    assert enmienda
    posteriores = subprocess.check_output(
        ["git", "log", "--reverse", "--format=%H", "--", "results/v3_2"],
        cwd=RAIZ,
        text=True,
    ).split()
    if not posteriores:
        return
    comprobacion = subprocess.run(
        ["git", "merge-base", "--is-ancestor", enmienda, posteriores[0]],
        cwd=RAIZ,
    )
    assert comprobacion.returncode == 0
