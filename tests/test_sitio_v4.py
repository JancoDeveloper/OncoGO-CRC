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


def texto_bib() -> str:
    partes = [(REPO / "tesis-ucaldas" / "latex" / "references.bib").read_text(encoding="utf-8")]
    extra = REPO / "sitio_fuente" / "referencias_sitio.bib"
    if extra.is_file():
        partes.append(extra.read_text(encoding="utf-8"))
    return "\n".join(partes)


def test_las_citas_estan_verificadas():
    enlaces = (REPO / "docs" / "web_enlaces.tsv").read_text(encoding="utf-8")
    bib = texto_bib()
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
    paginas_proteina = [
        ruta for ruta in (SITIO / "proteinas").glob("*/index.html")
        if ruta.parent.name not in {"comparar", "candidatas"}
    ]
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
    indice = (SITIO / "proteinas" / "index.html").read_text(encoding="utf-8")
    grupos = [href for href in re.findall(r'href="(/proteinas/[^"]+)"', indice) if ";" in href or "%3B" in href]
    assert len(grupos) == 65
    for href in grupos:
        assert "%3B" not in href
        assert (SITIO / "proteinas" / href.removeprefix("/proteinas/") / "index.html").is_file()


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


def test_la_anotacion_coincide_con_el_manifiesto():
    manifiesto = json.loads((REPO / "sitio_fuente" / "fuentes" / "anotacion_manifest.json").read_text(encoding="utf-8"))
    for entrada in manifiesto["entradas"]:
        if entrada.get("estado") != "ok" or not entrada.get("fichero"):
            continue
        ruta = REPO / "sitio_fuente" / "fuentes" / entrada["fichero"]
        assert ruta.is_file(), entrada["id"]
        assert hashlib.sha256(ruta.read_bytes()).hexdigest() == entrada["sha256"]
    fichas = json.loads((REPO / "sitio_fuente" / "fuentes" / "anotacion" / "fichas.json").read_text(encoding="utf-8"))
    for registro in fichas["por_acceso"].values():
        assert registro.get("estado") in {"ok", "no_disponible"}
        if registro.get("estado") != "ok":
            assert registro.get("motivo")
            continue
        assert len(registro.get("secuencia") or "") == registro.get("longitud")
        for termino in registro.get("go_directo") or []:
            assert not str(termino.get("evidencia", "")).startswith("IEA")


def test_el_jaccard_de_diez_parejas():
    fichas = json.loads((REPO / "sitio_fuente" / "fuentes" / "anotacion" / "fichas.json").read_text(encoding="utf-8"))
    muestra = json.loads((SITIO / "fuentes" / "jaccard_muestra.json").read_text(encoding="utf-8"))
    assert len(muestra) == 10

    def conjunto(acceso):
        return {item["id"] for item in fichas["por_acceso"][acceso].get("go_directo") or []}

    for par in muestra:
        a, b = conjunto(par["a"]), conjunto(par["b"])
        union = a | b
        esperado = 0 if not union else len(a & b) / len(union)
        assert abs(esperado - par["directo"]) < 1e-12


def test_el_alineamiento_js_coincide_con_biopython():
    import random
    import subprocess
    import tempfile
    from Bio.Align import PairwiseAligner, substitution_matrices

    fichas = json.loads((REPO / "sitio_fuente" / "fuentes" / "anotacion" / "fichas.json").read_text(encoding="utf-8"))
    alfabeto = set("ARNDCQEGHILKMFPSTWYVBZX*")
    secuencias = [
        registro["secuencia"][:80]
        for registro in fichas["por_acceso"].values()
        if registro.get("estado") == "ok"
        and len(registro.get("secuencia") or "") >= 80
        and set(registro["secuencia"][:80]) <= alfabeto
    ]
    rng = random.Random(20260914)
    pares = [rng.sample(secuencias, 2) for _ in range(30)]
    script = (
        "const fs=require('fs');"
        "eval(fs.readFileSync(process.argv[1],'utf8'));"
        "const pares=JSON.parse(fs.readFileSync(process.argv[2],'utf8'));"
        "process.stdout.write(JSON.stringify(pares.map(p => Alineamiento.alinear(p[0], p[1]))));"
    )
    ruta = REPO / "sitio_fuente" / "js" / "alineamiento.js"
    with tempfile.TemporaryDirectory() as carpeta:
        entrada = Path(carpeta) / "pares.json"
        entrada.write_text(json.dumps(pares), encoding="utf-8")
        proc = subprocess.run(
            ["node", "-e", script, str(ruta), str(entrada)],
            text=True, capture_output=True, check=True,
        )
    resultados = json.loads(proc.stdout)
    matriz = substitution_matrices.load("BLOSUM62")
    aligner = PairwiseAligner(mode="global", substitution_matrix=matriz, open_gap_score=-10, extend_gap_score=-0.5)
    for (a, b), resultado in zip(pares, resultados):
        aln = next(iter(aligner.align(a, b)))
        cuentas = aln.counts()
        assert abs(resultado["puntuacion"] - aligner.score(a, b)) < 1e-6
        assert resultado["identidades"] == cuentas.identities
        assert resultado["huecos"] == cuentas.gaps


def test_la_atribucion_esta_y_el_visor_no_se_carga_solo():
    for pagina in paginas():
        texto = pagina.read_text(encoding="utf-8")
        assert "<script src" not in texto
    p02749 = next(ruta for ruta in (SITIO / "proteinas").glob("P02749/index.html"))
    texto = p02749.read_text(encoding="utf-8")
    assert "CC BY 4.0" in texto and "CC0" in texto
    candidatas = (SITIO / "proteinas" / "candidatas" / "index.html").read_text(encoding="utf-8")
    assert "CC BY 4.0" in candidatas
    assert candidatas.count("Candidata (paso 5): hipótesis para cohortes futuras, no biomarcador") >= 1


def test_las_rutas_viejas_apuntan_a_paginas_que_existen():
    hosting = json.loads((REPO / "firebase.json").read_text(encoding="utf-8"))["hosting"]
    obligatorias = {
        "/lo-prometido": "/investigacion/historia",
        "/historial": "/investigacion/historia",
        "/tejido": "/resultados/tejido",
        "/suero": "/resultados/suero",
        "/aprendizaje-automatico": "/resultados/aprendizaje",
        "/preanalitica": "/resultados/preanalitica",
        "/hipotesis": "/resultados/hipotesis",
        "/priorizacion": "/resultados/priorizacion",
        "/desviaciones": "/entregables",
        "/reproducir": "/investigacion/reproducibilidad",
        "/como-se-hizo": "/investigacion/reproducibilidad",
        "/aprende": "/datos/bases-de-datos",
        "/preguntas": "/ayuda",
    }
    por_origen = {item["source"]: item for item in hosting["redirects"]}
    for origen, destino in obligatorias.items():
        regla = por_origen[origen]
        assert regla["destination"] == destino
        assert regla["type"] == 301
        carpeta = SITIO / destino.strip("/")
        assert (carpeta / "index.html").is_file() or carpeta.with_suffix(".html").is_file()
    for pagina in paginas():
        assert 'href="/lo-prometido"' not in pagina.read_text(encoding="utf-8")

