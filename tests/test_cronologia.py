"""La enmienda queda escrita antes de los resultados que interpreta."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PREFIJOS = ("prereg", "fix", "feat", "docs", "test", "tesis", "sitio", "chore", "ci")


def _git(*args: str) -> str:
    proceso = subprocess.run(
        ["git", *args],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    if proceso.returncode != 0:
        pytest.skip(proceso.stderr.strip() or "este entorno no tiene la historia de git")
    return proceso.stdout


def _primero_que_anade(ruta: str) -> str:
    salida = _git("log", "--diff-filter=A", "--reverse", "--format=%H", "--", ruta)
    commits = [linea for linea in salida.splitlines() if linea]
    if not commits:
        return ""
    return commits[0]


def _introduce_clave(clave: str) -> str:
    salida = _git(
        "log", "-S", clave, "--reverse", "--format=%H", "--", "config/preregistro_v3.yaml"
    )
    commits = [linea for linea in salida.splitlines() if linea]
    return commits[0] if commits else ""


def _es_anterior(anterior: str, posterior: str) -> bool:
    if not anterior or not posterior:
        return False
    proceso = subprocess.run(
        ["git", "merge-base", "--is-ancestor", anterior, posterior],
        cwd=REPO,
        capture_output=True,
        check=False,
    )
    return proceso.returncode == 0 and anterior != posterior


def test_enmienda_v3_2_es_anterior_a_sus_resultados():
    enmienda = _introduce_clave("enmienda_v3_2:")
    resultados = _primero_que_anade("results/v3_2")
    assert enmienda, "no aparece enmienda_v3_2 en el preregistro"
    assert resultados, "no hay ficheros bajo results/v3_2"
    assert _es_anterior(enmienda, resultados)


def test_enmienda_v3_3_es_anterior_a_sus_resultados():
    enmienda = _introduce_clave("enmienda_v3_3:")
    resultados = _primero_que_anade("results/v3_3")
    if not enmienda and not resultados:
        return
    if resultados and not enmienda:
        raise AssertionError("hay results/v3_3 sin la enmienda v3.3")
    if enmienda and resultados:
        assert _es_anterior(enmienda, resultados)


def test_enmienda_v3_4_es_anterior_a_sus_resultados():
    enmienda = _introduce_clave("enmienda_v3_4:")
    resultados = _primero_que_anade("results/v3_4")
    if not enmienda and not resultados:
        return
    if resultados and not enmienda:
        raise AssertionError("hay results/v3_4 sin la enmienda v3.4")
    if enmienda and resultados:
        assert _es_anterior(enmienda, resultados)


def test_enmienda_v3_6_es_anterior_a_sus_resultados():
    enmienda = _introduce_clave("enmienda_v3_6:")
    resultados = _primero_que_anade("results/v3_6")
    if not enmienda and not resultados:
        return
    if resultados and not enmienda:
        raise AssertionError("hay results/v3_6 sin la enmienda v3.6")
    if enmienda and resultados:
        assert _es_anterior(enmienda, resultados)


def test_enmienda_v3_5_es_anterior_a_sus_resultados():
    enmienda = _introduce_clave("enmienda_v3_5:")
    resultados = _primero_que_anade("results/v3_5")
    if not enmienda and not resultados:
        return
    if resultados and not enmienda:
        raise AssertionError("hay results/v3_5 sin la enmienda v3.5")
    if enmienda and resultados:
        assert _es_anterior(enmienda, resultados)


def test_v3_3_no_cambia_despues_de_la_etiqueta():
    proceso = subprocess.run(
        ["git", "rev-parse", "--verify", "v3.3-final"],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    if proceso.returncode != 0:
        pytest.skip("la etiqueta v3.3-final no está en este clon")
    posteriores = _git("rev-list", "v3.3-final..HEAD", "--", "results/v3_3")
    assert posteriores.strip() == ""


def test_v3_2_no_cambia_despues_de_la_etiqueta():
    etiqueta = _git("rev-parse", "v3.2-final")
    assert etiqueta.strip()
    posteriores = _git("rev-list", "v3.2-final..HEAD", "--", "results/v3_2")
    assert posteriores.strip() == ""


def test_los_commits_sin_prefijo_siguen_anotados():
    """50caef8 y c2be8e4 no llevan prefijo. 35fad60 sí lleva feat:; el encargo lo contó sin él."""
    texto = (REPO / "docs/cronologia_v3.md").read_text(encoding="utf-8")
    for corto in ("50caef8", "c2be8e4"):
        assert corto in texto
        linea = next(fila for fila in texto.splitlines() if corto in fila)
        assert "sin prefijo" in linea
    linea_resultados = next(fila for fila in texto.splitlines() if "35fad60" in fila)
    assert "feat" in linea_resultados
    assert any(prefijo + ":" in linea_resultados or "feat" in linea_resultados for prefijo in PREFIJOS)
