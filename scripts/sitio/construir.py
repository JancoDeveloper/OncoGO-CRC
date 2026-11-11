#!/usr/bin/env python3
"""Construye sitio/ desde sitio_fuente/ y results/v3_7/. No modifica la corrida."""

from __future__ import annotations

import hashlib
import html
import json
import re
import shutil
from pathlib import Path

import yaml
from jinja2 import Environment, FileSystemLoader, select_autoescape

from oncogo_crc.formato import formatear

REPO = Path(__file__).resolve().parents[2]
V37 = REPO / "results" / "v3_7"
FUENTE = REPO / "sitio_fuente"
SITIO = REPO / "sitio"
BIB = REPO / "tesis-ucaldas" / "latex" / "references.bib"
AVISO = "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
SHA = "595118ae09ea6006c744c86b28cb0a10d42e0442"
GITHUB = f"https://github.com/JancoDeveloper/OncoGO-CRC/blob/{SHA}"


def e(texto) -> str:
    return html.escape(str(texto), quote=True)


def caminar(obj, ruta: str):
    actual = obj
    for parte in ruta.split("."):
        if isinstance(actual, list):
            actual = actual[int(parte)]
        else:
            actual = actual[parte]
    return actual


def valor_de(fuente: str):
    if not fuente.startswith(("results/v3_7/", "config/")):
        raise SystemExit(f"fuente no permitida: {fuente}")
    ruta, campo = fuente.split("#", 1)
    path = REPO / ruta
    if not path.is_file():
        raise SystemExit(f"no está el fichero {ruta}")
    if ruta.endswith(".json"):
        return caminar(json.loads(path.read_text(encoding="utf-8")), campo)
    if ruta.endswith(".yaml"):
        return caminar(yaml.safe_load(path.read_text(encoding="utf-8")), campo)
    raise SystemExit(f"atajo sin lector: {ruta}")


def dato_html(fuente: str, contador: list[int]) -> str:
    valor = valor_de(fuente)
    if isinstance(valor, float):
        bruto = format(valor, ".16g")
        visible = formatear(valor, "auc" if abs(valor) <= 1 else "decimal")
    elif isinstance(valor, bool):
        bruto = str(valor).lower()
        visible = bruto
    else:
        bruto = str(valor)
        visible = formatear(valor, "entero") if isinstance(valor, int) else bruto
    contador[0] += 1
    ident = f"ev{contador[0]}"
    return (
        f'<data value="{e(bruto)}" data-fuente="{e(fuente)}">{e(visible)}</data>'
        f'<details class="evidencia"><summary>fuente</summary>'
        f'<p id="{ident}">{e(fuente)}</p></details>'
    )


def leer_bib() -> dict[str, dict]:
    texto = BIB.read_text(encoding="utf-8")
    entradas = {}
    for bloque in re.split(r"\n(?=@)", texto):
        cabeza = re.match(r"@\w+\{([^,]+),", bloque.strip())
        if not cabeza:
            continue
        campos = {}
        for campo, valor in re.findall(r"(?m)^\s*(\w+)\s*=\s*\{(.*)\},?\s*$", bloque):
            campos[campo] = valor.replace("{", "").replace("}", "")
        entradas[cabeza.group(1)] = campos
    return entradas


def apellido(autor: str) -> str:
    primero = autor.split(" and ")[0]
    if "," in primero:
        return primero.split(",")[0].strip()
    return primero.split()[-1]


def cargar_glosario() -> list[dict]:
    datos = yaml.safe_load((FUENTE / "glosario.yaml").read_text(encoding="utf-8"))
    return datos or []


def termino_html(nombre: str, glosario: list[dict]) -> str:
    for item in glosario:
        if item["termino"].casefold() == nombre.casefold() or item.get("sigla", "").casefold() == nombre.casefold():
            ancla = re.sub(r"[^a-z0-9]+", "-", item["termino"].casefold()).strip("-")
            return (
                f'<a class="termino" href="/glosario#{e(ancla)}">{e(item["termino"])}'
                f'<span class="definicion-corta">{e(item["en_palabras_simples"])}</span></a>'
            )
    raise SystemExit(f"término sin ficha: {nombre}")


def cita_html(clave: str, bib: dict[str, dict]) -> str:
    if clave not in bib:
        raise SystemExit(f"cita sin entrada: {clave}")
    campos = bib[clave]
    rotulo = f"{apellido(campos.get('author', clave))} {campos.get('year', '')}".strip()
    titulo = campos.get("title", "")
    doi = campos.get("doi", "")
    return (
        f'<span class="cita-wrap"><a class="cita" href="/bibliografia#{e(clave)}">{e(rotulo)}</a>'
        f'<span class="cita-ventana">{e(titulo)} {e(doi)}</span></span>'
    )


def ancla_de(texto: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", texto.casefold()).strip("-") or "seccion"


def md_inline(texto: str, bib, glosario, contador) -> str:
    patron = re.compile(
        r"\{\{dato:([^}]+)\}\}|\{\{termino:([^}]+)\}\}|\{\{diagrama:([^}]+)\}\}|\[@([^\]]+)\]"
    )
    piezas = []
    pos = 0
    for hallado in patron.finditer(texto):
        piezas.append(e(texto[pos:hallado.start()]))
        if hallado.group(1):
            piezas.append(dato_html(hallado.group(1), contador))
        elif hallado.group(2):
            piezas.append(termino_html(hallado.group(2), glosario))
        elif hallado.group(3):
            piezas.append(diagrama_html(hallado.group(3), bib))
        else:
            piezas.append(cita_html(hallado.group(4), bib))
        pos = hallado.end()
    piezas.append(e(texto[pos:]))
    return "".join(piezas)


def md_a_html(texto: str, bib, glosario, contador) -> tuple[str, list[dict]]:
    indice = []
    bloques = re.split(r"\n\s*\n", texto.strip())
    piezas = []
    for bloque in bloques:
        lineas = bloque.splitlines()
        if bloque.strip().startswith("{{diagrama:"):
            piezas.append(md_inline(bloque.strip(), bib, glosario, contador))
        elif bloque.startswith("# "):
            piezas.append(f"<h1>{md_inline(bloque[2:].strip(), bib, glosario, contador)}</h1>")
        elif bloque.startswith("## "):
            titulo = bloque[3:].strip()
            ident = ancla_de(re.sub(r"\{\{[^}]+\}\}|\[@[^\]]+\]", "", titulo))
            indice.append({"id": ident, "texto": re.sub(r"\{\{[^}]+\}\}|\[@[^\]]+\]", "", titulo).strip()})
            piezas.append(f'<h2 id="{e(ident)}">{md_inline(titulo, bib, glosario, contador)}</h2>')
        elif bloque.startswith("### "):
            piezas.append(f"<h3>{md_inline(bloque[4:].strip(), bib, glosario, contador)}</h3>")
        elif all(re.match(r"(\d+\. |- )", linea) for linea in lineas):
            items = []
            for linea in lineas:
                limpio = re.sub(r"^(\d+\. |- )", "", linea)
                items.append(f"<li>{md_inline(limpio, bib, glosario, contador)}</li>")
            piezas.append("<ul>" + "".join(items) + "</ul>")
        else:
            piezas.append(f"<p>{md_inline(' '.join(lineas), bib, glosario, contador)}</p>")
    return "\n".join(piezas), indice


def diagrama_html(ident: str, bib) -> str:
    ruta = FUENTE / "diagramas" / f"{ident}.yaml"
    if not ruta.is_file():
        raise SystemExit(f"diagrama sin definición: {ident}")
    datos = yaml.safe_load(ruta.read_text(encoding="utf-8"))
    nodos = datos["nodos"]
    ancho = 160
    cuerpo = []
    filas = []
    for i, nodo in enumerate(nodos):
        x = 24 + i * (ancho + 28)
        cuerpo.append(
            f'<a href="{e(nodo["href"])}"><rect x="{x}" y="70" width="{ancho}" height="64" rx="8" '
            f'fill="none" stroke="currentColor" stroke-width="2"/></a>'
            f'<text x="{x + 16}" y="108" fill="currentColor" font-size="16">{e(nodo["texto"])}</text>'
        )
        if i:
            cuerpo.append(
                f'<line x1="{x - 28}" y1="102" x2="{x}" y2="102" stroke="currentColor" stroke-width="2"/>'
            )
        filas.append({"nodo": nodo["texto"], "fuente": nodo["fuente"]})
    for arista in datos.get("aristas", []):
        filas.append({"nodo": f'{arista["origen"]} → {arista["destino"]}', "fuente": arista["fuente"]})
    alto = 180
    largo = 24 + len(nodos) * (ancho + 28)
    svg = (
        f'<svg viewBox="0 0 {largo} {alto}" role="img" class="lienzo-svg">'
        f'<title>{e(datos["titulo"])}</title><desc>{e(datos["desc"])}</desc>'
        + "".join(cuerpo) + "</svg>"
    )
    carpeta = SITIO / "figuras"
    carpeta.mkdir(parents=True, exist_ok=True)
    (carpeta / f"{ident}.svg").write_text(svg, encoding="utf-8")
    sello = {
        "fuentes": sorted({fila["fuente"] for fila in filas}),
        "svg_sha256": hashlib.sha256(svg.encode()).hexdigest(),
        "commit_de_la_etiqueta": SHA,
    }
    (carpeta / f"{ident}.sello.json").write_text(json.dumps(sello, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    columnas = ["nodo", "fuente"]
    lineas = ["\t".join(columnas)] + ["\t".join(fila[col] for col in columnas) for fila in filas]
    (carpeta / f"{ident}.tsv").write_text("\n".join(lineas) + "\n", encoding="utf-8")
    tabla = "".join(
        f"<tr><td>{e(fila['nodo'])}</td><td>{e(fila['fuente'])}</td></tr>" for fila in filas
    )
    return (
        f'<figure class="grafico" id="{e(ident)}"><h3>{e(datos["pregunta"])}</h3>{svg}'
        f'<details><summary>Tabla de los datos del gráfico</summary>'
        f'<table><thead><tr><th>nodo</th><th>fuente</th></tr></thead><tbody>{tabla}</tbody></table></details>'
        f'<p><a href="/figuras/{e(ident)}.svg" download>Descargar SVG</a> '
        f'<a href="/figuras/{e(ident)}.tsv" download>Descargar TSV</a></p></figure>'
    )


def leer_paginas() -> list[dict]:
    paginas = []
    for ruta in sorted((FUENTE / "contenido").rglob("*.md")):
        texto = ruta.read_text(encoding="utf-8")
        if not texto.startswith("---"):
            raise SystemExit(f"sin cabecera: {ruta}")
        _, frente, cuerpo = texto.split("---", 2)
        meta = yaml.safe_load(frente)
        meta["cuerpo"] = cuerpo.strip()
        meta["origen"] = ruta
        paginas.append(meta)
    paginas.sort(key=lambda item: (item.get("orden", 0), item["ruta"]))
    return paginas


def menu_de(paginas: list[dict]) -> list[dict]:
    vistos = []
    for pagina in paginas:
        if pagina["seccion"] == "Inicio":
            continue
        if not any(item["nombre"] == pagina["seccion"] for item in vistos):
            vistos.append({"nombre": pagina["seccion"], "ruta": pagina["ruta"]})
    vistos.append({"nombre": "Glosario", "ruta": "/glosario"})
    vistos.append({"nombre": "Bibliografía", "ruta": "/bibliografia"})
    return vistos


def fuentes_de(cuerpo: str, bib) -> list[dict]:
    salida = []
    vistos = set()
    for fuente in re.findall(r'data-fuente="([^"]+)"', cuerpo):
        if fuente not in vistos:
            vistos.add(fuente)
            salida.append({"href": f"{GITHUB}/{fuente.split('#', 1)[0]}", "texto": fuente})
    for clave in re.findall(r'href="/bibliografia#([^"]+)"', cuerpo):
        if clave not in vistos and clave in bib:
            vistos.add(clave)
            doi = bib[clave].get("doi", "")
            salida.append({"href": f"https://doi.org/{doi}", "texto": f"{clave} · {doi}"})
    return salida


def escribir_estaticos() -> tuple[str, str]:
    css = (FUENTE / "estilos" / "tokens.css").read_text(encoding="utf-8")
    css += "\n" + (FUENTE / "estilos" / "componentes.css").read_text(encoding="utf-8")
    js = (FUENTE / "js" / "sitio.js").read_text(encoding="utf-8")
    css_nombre = "v4." + hashlib.sha256(css.encode()).hexdigest()[:10] + ".css"
    js_nombre = "v4." + hashlib.sha256(js.encode()).hexdigest()[:10] + ".js"
    (SITIO / "css").mkdir(parents=True, exist_ok=True)
    for viejo in (SITIO / "css").glob("v4.*.css"):
        if viejo.name != css_nombre:
            viejo.unlink()
    (SITIO / "css" / css_nombre).write_text(css, encoding="utf-8")
    (SITIO / "js").mkdir(parents=True, exist_ok=True)
    (SITIO / "js" / js_nombre).write_text(js, encoding="utf-8")
    destino_fuentes = SITIO / "fuentes"
    destino_fuentes.mkdir(parents=True, exist_ok=True)
    for archivo in (FUENTE / "estilos" / "fuentes").glob("*"):
        if archivo.suffix in {".woff2", ".txt"}:
            shutil.copy2(archivo, destino_fuentes / archivo.name)
    return css_nombre, js


def archivar_v37() -> None:
    destino = SITIO / "archivo" / "v3_7"
    if (destino / "index.html").exists():
        return
    destino.mkdir(parents=True, exist_ok=True)
    for pagina in SITIO.glob("*.html"):
        shutil.copy2(pagina, destino / pagina.name)
    indice = SITIO / "archivo" / "index.html"
    if indice.is_file():
        (destino / "archivo").mkdir(exist_ok=True)
        shutil.copy2(indice, destino / "archivo" / "index.html")
    for carpeta in ("css", "figuras"):
        origen = SITIO / carpeta
        if origen.is_dir() and not (destino / carpeta).exists():
            shutil.copytree(origen, destino / carpeta)
    for pagina in list(destino.glob("*.html")) + list((destino / "archivo").glob("*.html")):
        texto = pagina.read_text(encoding="utf-8")
        texto = texto.replace('href="/css/', 'href="/archivo/v3_7/css/')
        texto = texto.replace('href="/figuras/', 'href="/archivo/v3_7/figuras/')
        pagina.write_text(texto, encoding="utf-8")
    for pagina in SITIO.glob("*.html"):
        pagina.unlink()
    if indice.is_file():
        indice.unlink()
    for carpeta in ("css", "figuras"):
        origen = SITIO / carpeta
        if origen.is_dir():
            shutil.rmtree(origen)


def render(env, plantilla, **datos) -> str:
    return env.get_template("base.html.j2").render(**datos)


def main() -> None:
    archivar_v37()
    bib = leer_bib()
    glosario = cargar_glosario()
    css_nombre, js = escribir_estaticos()
    cert = json.loads((V37 / "certificado.json").read_text(encoding="utf-8"))
    contador = [0]
    estado = dato_html("results/v3_7/certificado.json#estado", contador)
    sello = dato_html("results/v3_7/certificado.json#sello.commit", contador)
    paginas = leer_paginas()
    menu = [{"nombre": "Inicio", "ruta": "/"}] + menu_de(paginas)
    env = Environment(
        loader=FileSystemLoader(FUENTE / "plantillas"),
        autoescape=select_autoescape(["html", "j2"]),
    )
    comunes = {
        "css": css_nombre,
        "js": js,
        "aviso": AVISO,
        "consulta": "",
        "franja_estado": estado,
        "franja_sello": sello,
        "menu": menu,
    }
    for pagina in paginas:
        cuerpo, indice = md_a_html(pagina["cuerpo"], bib, glosario, contador)
        if pagina["ruta"] == "/":
            tarjetas = []
            for item in menu:
                if item["ruta"] == "/":
                    continue
                tarjetas.append(
                    f'<a class="tarjeta" href="{e(item["ruta"])}"><h2>{e(item["nombre"])}</h2></a>'
                )
            cuerpo = cuerpo.replace("</h1>", "</h1>" + '<div class="tarjetas">' + "".join(tarjetas) + "</div>", 1)
        migas = "Inicio" if pagina["ruta"] == "/" else f'Inicio / {pagina["seccion"]} / {pagina["titulo"]}'
        html_pagina = render(
            env, "base.html.j2",
            titulo=pagina["titulo"],
            cuerpo=cuerpo,
            indice=indice,
            migas=migas,
            fuentes=fuentes_de(cuerpo, bib),
            **comunes,
        )
        destino = SITIO / "index.html" if pagina["ruta"] == "/" else SITIO / pagina["ruta"].strip("/") / "index.html"
        destino.parent.mkdir(parents=True, exist_ok=True)
        destino.write_text(html_pagina, encoding="utf-8")

    entradas = []
    for item in glosario:
        ident = ancla_de(item["termino"])
        fuentes = " ".join(cita_html(clave, bib) for clave in item["fuentes"])
        entradas.append(
            f'<dt id="{e(ident)}">{e(item["termino"])}</dt>'
            f'<dd><p>{e(item["definicion"])}</p>'
            f'<p>{e(item["en_palabras_simples"])}</p>'
            f'<p>{e(item["rol_en_bioinformatica"])}</p>'
            f'<p><a href="{e(item["como_se_uso_aqui"])}">Dónde se usa</a></p>'
            f'<p>{fuentes}</p></dd>'
        )
    glosario_html = render(
        env, "base.html.j2",
        titulo="Glosario",
        cuerpo="<h1>Glosario</h1><dl>" + "".join(entradas) + "</dl>",
        indice=[],
        migas="Inicio / Glosario",
        fuentes=[{"href": "/bibliografia", "texto": "Bibliografía"}],
        **comunes,
    )
    (SITIO / "glosario").mkdir(parents=True, exist_ok=True)
    (SITIO / "glosario" / "index.html").write_text(glosario_html, encoding="utf-8")

    citadas = set()
    for pagina in SITIO.rglob("*.html"):
        if "archivo" in pagina.parts:
            continue
        citadas.update(re.findall(r'href="/bibliografia#([^"]+)"', pagina.read_text(encoding="utf-8")))
    bloques = []
    for clave in sorted(citadas):
        campos = bib[clave]
        bloques.append(
            f'<article id="{e(clave)}" class="entrada"><h2>{e(apellido(campos.get("author", clave)))} {e(campos.get("year", ""))}</h2>'
            f'<p class="titulo-citado">{e(campos.get("title", ""))}</p>'
            f'<p><a href="https://doi.org/{e(campos.get("doi", ""))}">https://doi.org/{e(campos.get("doi", ""))}</a></p></article>'
        )
    bibliografia = render(
        env, "base.html.j2",
        titulo="Bibliografía",
        cuerpo="<h1>Bibliografía</h1>" + "".join(bloques),
        indice=[],
        migas="Inicio / Bibliografía",
        fuentes=[],
        **comunes,
    )
    (SITIO / "bibliografia").mkdir(parents=True, exist_ok=True)
    (SITIO / "bibliografia" / "index.html").write_text(bibliografia, encoding="utf-8")

    grupos = []
    for pagina in paginas:
        grupos.append(
            f'<li data-buscar="{e(pagina["titulo"])} {e(pagina["seccion"])}">'
            f'<a href="{e(pagina["ruta"])}">{e(pagina["titulo"])}</a></li>'
        )
    for item in glosario:
        grupos.append(
            f'<li data-buscar="{e(item["termino"])} {e(item.get("sigla", ""))}">'
            f'<a href="/glosario#{e(ancla_de(item["termino"]))}">{e(item["termino"])}</a></li>'
        )
    buscar = render(
        env, "base.html.j2",
        titulo="Buscar",
        cuerpo="<h1>Buscar</h1><p>Sin JavaScript se ve el índice completo.</p><ul>" + "".join(grupos) + "</ul>",
        indice=[],
        migas="Inicio / Buscar",
        fuentes=[],
        **comunes,
    )
    (SITIO / "buscar").mkdir(parents=True, exist_ok=True)
    (SITIO / "buscar" / "index.html").write_text(buscar, encoding="utf-8")
    faltante = render(
        env, "base.html.j2",
        titulo="No está",
        cuerpo='<h1>Esta página no está</h1><p><a href="/">Volver al inicio</a>.</p>',
        indice=[],
        migas="Inicio",
        fuentes=[],
        **comunes,
    )
    (SITIO / "404.html").write_text(faltante, encoding="utf-8")
    print("sitio documental", css_nombre, "paginas", len(paginas))


if __name__ == "__main__":
    main()
