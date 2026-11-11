"""La web documental cita results/v3_7 y no sale a la red."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
SITIO = REPO / "sitio"
AVISO = "Ningún resultado de esta tesis sirve para decidir sobre un paciente."


def paginas():
    encontradas = []
    for ruta in SITIO.rglob("*.html"):
        if "archivo" in ruta.parts:
            continue
        encontradas.append(ruta)
    return encontradas


def contraste(fondo: str, tinta: str) -> float:
    def canal(valor: str) -> float:
        numero = int(valor, 16) / 255
        return numero / 12.92 if numero <= 0.04045 else ((numero + 0.055) / 1.055) ** 2.4

    def luz(color: str) -> float:
        rojo, verde, azul = color[1:3], color[3:5], color[5:7]
        return 0.2126 * canal(rojo) + 0.7152 * canal(verde) + 0.0722 * canal(azul)

    clara, oscura = sorted((luz(fondo), luz(tinta)))
    return (oscura + 0.05) / (clara + 0.05)


def test_el_contraste_de_los_tokens_pasa_aa():
    css = (REPO / "sitio_fuente" / "estilos" / "tokens.css").read_text(encoding="utf-8")
    bloques = re.findall(r"(?::root|\[[^\]]+\])\s*\{([^}]+)\}", css)
    assert len(bloques) >= 2
    for bloque in bloques:
        colores = dict(re.findall(r"--([\w-]+):\s*(#[0-9a-fA-F]{6})", bloque))
        pares = [
            ("papel", "tinta"),
            ("papel", "tinta-suave"),
            ("marca", "sobre-marca"),
            ("aviso-fondo", "aviso-tinta"),
        ]
        for fondo, tinta in pares:
            assert contraste(colores[fondo], colores[tinta]) >= 4.5, (fondo, tinta)


def test_cifras_aviso_y_lenguaje():
    for pagina in paginas():
        texto = pagina.read_text(encoding="utf-8")
        assert texto.count(AVISO) == 2, pagina
        visible = re.sub(r"<[^>]+>", " ", texto)
        assert not re.search(r"\bCN\b.{0,40}\bsanos\b|\bsanos\b.{0,40}\bCN\b", visible, re.I)
        if pagina.name != "index.html" or pagina.parent.name == "bibliografia":
            pass
        sin_titulos = re.sub(r'<span class="cita-ventana">.*?</span>', " ", texto, flags=re.S)
        sin_titulos = re.sub(r'<p class="titulo-citado">.*?</p>', " ", sin_titulos, flags=re.S)
        sin_titulos = sin_titulos.replace(
            "Candidata (paso 5): hipótesis para cohortes futuras, no biomarcador", " "
        )
        if "bibliografia" not in pagina.parts:
            assert "biomarcador" not in re.sub(r"<[^>]+>", " ", sin_titulos).casefold()
            assert "Biomarker" not in sin_titulos
    bib = (SITIO / "bibliografia" / "index.html").read_text(encoding="utf-8")
    assert "titulo-citado" in bib


def test_las_citas_estan_verificadas():
    enlaces = (REPO / "docs" / "web_enlaces.tsv").read_text(encoding="utf-8")
    bib = (REPO / "tesis-ucaldas" / "latex" / "references.bib").read_text(encoding="utf-8")
    for pagina in paginas():
        for clave in re.findall(r'href="/bibliografia#([^"]+)"', pagina.read_text(encoding="utf-8")):
            assert "@" in bib and clave in bib
            bloque = bib.split("@" ) 
            entrada = next(parte for parte in bib.split("\n@") if parte.startswith(f"article{{{clave},") or parte.startswith(f"misc{{{clave},") or f"{{{clave}," in parte[:80])
            doi = re.search(r"doi\s*=\s*\{([^}]+)\}", entrada).group(1)
            linea = next(fila for fila in enlaces.splitlines() if doi in fila)
            assert linea.split("\t")[-1] == "200"


def test_el_css_lleva_su_hash_y_no_hay_terceros():
    hojas = list((SITIO / "css").glob("v4.*.css"))
    assert len(hojas) == 1
    digest = hashlib.sha256(hojas[0].read_bytes()).hexdigest()[:10]
    assert hojas[0].name == f"v4.{digest}.css"
    for pagina in paginas():
        texto = pagina.read_text(encoding="utf-8")
        assert f"/css/{hojas[0].name}" in texto
        assert "<script src" not in texto
        assert "fonts.googleapis" not in texto
        assert "cdn." not in texto.lower()


def test_cada_nodo_de_diagrama_tiene_fichero_o_cita():
    bib = (REPO / "tesis-ucaldas" / "latex" / "references.bib").read_text(encoding="utf-8")
    manifiesto = (REPO / "results" / "v3_7" / "manifest_entradas.sha256").read_text(encoding="utf-8")
    for ruta in (REPO / "sitio_fuente" / "diagramas").glob("D*.yaml"):
        datos = yaml.safe_load(ruta.read_text(encoding="utf-8"))
        piezas = list(datos["nodos"]) + list(datos.get("aristas", []))
        for pieza in piezas:
            fuente = pieza["fuente"]
            assert (
                f"{{{fuente}," in bib
                or (REPO / fuente).is_file()
                or fuente in manifiesto
            ), (ruta.name, fuente)
        assert (SITIO / "figuras" / f"{datos['id']}.svg").is_file()
        assert (SITIO / "figuras" / f"{datos['id']}.sello.json").is_file()


def test_los_graficos_de_resultados_salen_de_la_corrida():
    suero = (SITIO / "resultados" / "suero" / "index.html").read_text(encoding="utf-8")
    assert 'id="r01"' in suero and 'id="f02"' in suero
    bosque = (SITIO / "figuras" / "r06.tsv").read_text(encoding="utf-8").splitlines()
    assert len(bosque) == 21
    svg = (SITIO / "figuras" / "f02.svg").read_text(encoding="utf-8")
    assert "currentColor" in svg and "#0b4f55" not in svg
    sello = json.loads((SITIO / "figuras" / "f02.sello.json").read_text(encoding="utf-8"))
    assert sello["svg_sha256"] == hashlib.sha256(svg.encode()).hexdigest()


def test_las_284_fichas_son_la_celda():
    import pandas as pd
    celda = pd.read_csv(REPO / "results" / "v3_7" / "suero_k4_crc_cn.tsv", sep="\t")
    otras = [
        set(pd.read_csv(REPO / "results" / "v3_7" / nombre, sep="\t")["protein_id"].astype(str))
        for nombre in ("suero_k4_crc_sanos.tsv", "suero_k4_cn_sanos.tsv")
    ]
    assert set(celda["protein_id"].astype(str)) == otras[0] == otras[1]
    paginas_proteina = list((SITIO / "proteinas").glob("*/index.html"))
    assert len(paginas_proteina) == len(celda)
    meta = json.loads((REPO / "sitio_fuente" / "fuentes" / "uniprot_284.meta.json").read_text(encoding="utf-8"))
    tsv = (REPO / "sitio_fuente" / "fuentes" / "uniprot_284.tsv").read_bytes()
    assert meta["sha256"] == hashlib.sha256(tsv).hexdigest()
    assert meta["faltan"] == []
    assert meta["x_uniprot_release"]
    etiqueta = "Candidata (paso 5): hipótesis para cohortes futuras, no biomarcador"
    con_etiqueta = 0
    for pagina in paginas_proteina:
        if pagina.parent.name == "proteinas":
            continue
        texto = pagina.read_text(encoding="utf-8")
        if etiqueta in texto:
            con_etiqueta += 1
        if pagina.parent.name in {"P05109", "P06702"}:
            assert "calprotectina" in texto
    assert con_etiqueta == 21


def test_d1_tiene_tabla_y_sello():
    pagina = (SITIO / "index.html").read_text(encoding="utf-8")
    assert 'role="img"' in pagina and "<title>" in pagina and "<desc>" in pagina
    assert (SITIO / "figuras" / "D1.sello.json").is_file()
    assert (SITIO / "figuras" / "D1.tsv").is_file()


def test_los_terminos_usados_tienen_fuente():
    glosario = yaml.safe_load((REPO / "sitio_fuente" / "glosario.yaml").read_text(encoding="utf-8"))
    nombres = {item["termino"].casefold() for item in glosario}
    for item in glosario:
        for campo in ("definicion", "en_palabras_simples", "rol_en_bioinformatica", "como_se_uso_aqui", "fuentes"):
            assert item.get(campo)
        assert item["fuentes"]
    for ruta in (REPO / "sitio_fuente" / "contenido").rglob("*.md"):
        for nombre in re.findall(r"\{\{termino:([^}]+)\}\}", ruta.read_text(encoding="utf-8")):
            assert nombre.casefold() in nombres
