"""La web v3.7 cita results/v3_7 y no sale a la red."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import pandas as pd
import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]
SITIO = REPO / "sitio"
AVISO = "Ningún resultado de esta tesis sirve para decidir sobre un paciente."


def paginas_v37() -> list[Path]:
    encontradas = list(SITIO.glob("*.html"))
    indice = SITIO / "archivo" / "index.html"
    if indice.is_file():
        encontradas.append(indice)
    return encontradas


def cargar_json(relativo: str):
    return json.loads((REPO / relativo).read_text(encoding="utf-8"))


def caminar(obj, ruta: str):
    actual = obj
    for parte in ruta.split("."):
        if parte == "":
            continue
        if isinstance(actual, list):
            actual = actual[int(parte)]
        else:
            actual = actual[parte]
    return actual


def valor_de(fuente: str):
    if "#" not in fuente:
        raise AssertionError(fuente)
    ruta, campo = fuente.split("#", 1)
    if ruta == "docs/cambio_esm2_b_d.md":
        for linea in (REPO / ruta).read_text(encoding="utf-8").splitlines():
            if linea.startswith("| B |"):
                return int([parte.strip() for parte in linea.strip("|").split("|")][-1].split()[0])
        raise AssertionError(ruta)
    if ruta.endswith(".yaml"):
        return caminar(yaml.safe_load((REPO / ruta).read_text(encoding="utf-8")), campo)
    if ruta.endswith("manifest_entradas.sha256"):
        for linea in (REPO / ruta).read_text(encoding="utf-8").splitlines():
            huella, nombre = linea.split(maxsplit=1)
            if nombre == campo:
                return huella
        raise AssertionError(campo)
    if ruta.endswith(".json"):
        return caminar(cargar_json(ruta), campo)
    if ruta.endswith(".tsv"):
        tabla = pd.read_csv(REPO / ruta, sep="\t")
        if campo.endswith("#cuenta"):
            selector = campo[: -len("#cuenta")]
            clave, valor = selector.split("=", 1)
            return int((tabla[clave].astype(str) == valor).sum())
        selector, columna = campo.rsplit(".", 1)
        mascara = pd.Series(True, index=tabla.index)
        for parte in selector.split("|"):
            clave, valor = parte.split("=", 1)
            if clave == "fila":
                return tabla.iloc[int(valor)][columna]
            mascara &= tabla[clave].astype(str) == valor
        hallados = tabla.loc[mascara, columna]
        if hallados.empty:
            raise AssertionError(fuente)
        return hallados.iloc[0]
    raise AssertionError(ruta)


def coincide(esperado, bruto: str) -> bool:
    if isinstance(esperado, bool):
        return bruto == str(esperado).lower()
    if _numero(esperado) and _numero(bruto):
        return math.isclose(float(esperado), float(bruto), rel_tol=0, abs_tol=1e-8)
    return str(esperado) == bruto


def _numero(texto) -> bool:
    try:
        float(texto)
    except (TypeError, ValueError):
        return False
    return True


@pytest.mark.parametrize("pagina", paginas_v37(), ids=lambda p: p.name)
def test_cada_data_fuente_coincide(pagina: Path):
    texto = pagina.read_text(encoding="utf-8")
    pares = re.findall(r'data-value="([^"]*)"[^>]*data-fuente="([^"]*)"', texto)
    pares += re.findall(r'data-value-y="([^"]*)"[^>]*data-fuente-y="([^"]*)"', texto)
    vistos = []
    for primero, segundo in pares:
        if primero.startswith("results/") or primero.startswith("config/") or primero.startswith("docs/"):
            fuente, bruto = primero, segundo
        else:
            bruto, fuente = primero, segundo
        vistos.append((bruto, fuente))
    pares_data = re.findall(r'<data value="([^"]*)" data-fuente="([^"]*)">', texto)
    for bruto, fuente in pares_data + vistos:
        assert coincide(valor_de(fuente), bruto), (pagina.name, fuente, bruto)


@pytest.mark.parametrize("pagina", paginas_v37(), ids=lambda p: p.name)
def test_franja_aviso_y_cn(pagina: Path):
    texto = pagina.read_text(encoding="utf-8")
    cert = cargar_json("results/v3_7/certificado.json")
    assert texto.count(AVISO) == 2
    assert cert["estado"] in texto
    assert "v3.7-final" in texto
    visible = re.sub(r"<[^>]+>", " ", texto)
    assert not re.search(r"\bCN\b.{0,40}\bsanos\b|\bsanos\b.{0,40}\bCN\b", visible, re.I)
    assert "results/v3_6" not in texto
    assert "<script src" not in texto
    assert "fonts.googleapis" not in texto
    assert "cdn." not in texto.lower()
    assert pagina.stat().st_size < 400_000


def test_svg_tabla_y_movimiento():
    css = (SITIO / "css" / "v37.css").read_text(encoding="utf-8")
    assert "prefers-reduced-motion" in css
    assert "1.2s" in css
    for pagina in paginas_v37():
        texto = pagina.read_text(encoding="utf-8")
        script = re.search(r"<script>(.*?)</script>", texto, re.S)
        assert script and len(script.group(1).encode()) < 20_000
        for svg in re.findall(r"<svg\b.*?</svg>", texto, re.S):
            assert 'role="img"' in svg
            assert "<title>" in svg and "<desc>" in svg
        if "<svg" in texto:
            assert "<table" in texto


def test_tarjetas_glosario_y_enlaces():
    corr = pd.read_csv(REPO / "config" / "correspondencia_propuesta.tsv", sep="\t")
    prometido = (SITIO / "lo-prometido.html").read_text(encoding="utf-8")
    assert prometido.count("tarjeta-prometida") == len(corr)
    glosario = (SITIO / "glosario.html").read_text(encoding="utf-8")
    for i in range(1, 21):
        assert f'id="g{i}"' in glosario
    assert len(re.findall(r"<dt ", glosario)) >= 40
    externos = set()
    for pagina in paginas_v37():
        for href in re.findall(r'href="([^"]+)"', pagina.read_text(encoding="utf-8")):
            if href.startswith("http"):
                externos.add(href.split("#")[0])
            elif href.startswith("/") and not href.startswith("//"):
                limpio = href.split("#")[0].split("?")[0]
                if limpio.startswith("/figuras/") or limpio.startswith("/css/") or limpio.startswith("/archivo/v3_6"):
                    assert (SITIO / limpio.lstrip("/")).exists() or (SITIO / (limpio.lstrip("/") + "/index.html")).exists()
                    continue
                candidato = SITIO / (limpio.lstrip("/") + ".html")
                indice = SITIO / limpio.lstrip("/") / "index.html"
                archivo = SITIO / limpio.lstrip("/")
                assert candidato.is_file() or indice.is_file() or archivo.is_file(), href
    tabla = pd.read_csv(REPO / "config" / "enlaces_externos.tsv", sep="\t")
    comprobados = pd.read_csv(REPO / "docs" / "web_enlaces.tsv", sep="\t")
    for url in externos:
        fila = tabla.loc[tabla["url"].map(lambda base: url == base or url.startswith(str(base).rstrip("/") + "/"))]
        assert not fila.empty, url
        assert (fila["estado_http"].astype(str) == "200").all()
        otra = comprobados.loc[comprobados["url"].map(lambda base: url == base or url.startswith(str(base).rstrip("/") + "/"))]
        assert not otra.empty and (otra["estado_http"].astype(str) == "200").all()
