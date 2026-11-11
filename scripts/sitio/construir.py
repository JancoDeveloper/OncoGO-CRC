#!/usr/bin/env python3
"""Construye sitio/ desde sitio_fuente/ y results/v3_7/. No modifica la corrida."""

from __future__ import annotations

import hashlib
import html
import importlib.util
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
FIGURAS: dict[str, str] = {}


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


def dato_html(fuente: str, contador: list[int], traza: bool = True) -> str:
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
    dato = f'<data value="{e(bruto)}" data-fuente="{e(fuente)}">{e(visible)}</data>'
    if not traza:
        return dato
    contador[0] += 1
    ident = f"ev{contador[0]}"
    return (
        dato
        + f'<details class="evidencia"><summary>Origen</summary>'
        + f'<p id="{ident}">{e(fuente)}</p></details>'
    )


def leer_bib() -> dict[str, dict]:
    textos = [BIB.read_text(encoding="utf-8")]
    extra = FUENTE / "referencias_sitio.bib"
    if extra.is_file():
        textos.append(extra.read_text(encoding="utf-8"))
    texto = "\n".join(textos)
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


def llamadas_html() -> str:
    return (
        '<p class="llamadas"><a class="boton" href="/investigacion/contexto">Explorar la investigación</a>'
        '<a class="boton boton-secundario" href="/resultados">Ver resultados</a></p>'
    )


CREDITO_SERVIER = (
    "Image provided by Servier Medical Art (https://smart.servier.com/), "
    "licensed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)."
)
NOTA_CONCEPTUAL = "Ilustración conceptual. No corresponde a una imagen experimental de esta investigación."


def _porcentaje(valor: int, total: int) -> str:
    return f"{100 * valor / total:.1f}".replace(".", ",") + " %"


def _lamina(archivo: str, respaldo: str, alt: str, ancho: int, alto: int, pie: str) -> str:
    return (
        '<figure class="lamina">'
        f'<img src="/ilustraciones/{archivo}" data-respaldo="/ilustraciones/{respaldo}" '
        f'alt="{e(alt)}" width="{ancho}" height="{alto}" loading="lazy" decoding="async">'
        f"<figcaption>{e(pie)} {e(CREDITO_SERVIER)} {e(NOTA_CONCEPTUAL)}</figcaption>"
        "</figure>"
    )


def celda_html(contador) -> str:
    return comparacion_html(contador)


def embudo_html(ident: str, contador) -> str:
    if ident != "priorizacion":
        raise SystemExit(f"embudo sin definición: {ident}")
    datos = json.loads((V37 / "priorizacion_embudo.json").read_text(encoding="utf-8"))
    pasos = datos["pasos"]
    mayor = max(int(paso["n"]) for paso in pasos) or 1
    filas = []
    for i, paso in enumerate(pasos):
        ancho = 100 * int(paso["n"]) / mayor
        fuente = f"results/v3_7/priorizacion_embudo.json#pasos.{i}.n"
        pieza = dato_html(fuente, contador)
        numero, separador, resto = pieza.partition("<details")
        filas.append(
            "<li><details>"
            "<summary>"
            f'<span class="embudo-nombre">{e(paso["nombre"])}</span>'
            f'<span class="embudo-riel"><span class="embudo-barra" style="--meta:{ancho / 100:.6f}"></span></span>'
            f"{numero}"
            "</summary>"
            f"{separador}{resto}"
            f"<p>Definición guardada en el JSON: {e(paso['nombre'])}. "
            f"El n es el de {e(fuente)}. No hay otra frase para este paso.</p>"
            "</details></li>"
        )
    return '<ol class="embudo">' + "".join(filas) + "</ol>"


def _hecho_estudio(contador) -> str:
    return (
        "<p>Lo que sigue es del estudio, no de este segmento. "
        "En tejido, el signo del tumor primario frente a la mucosa coincide en "
        f'{dato_html("results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.numerador", contador)}'
        " de "
        f'{dato_html("results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.denominador", contador)}'
        ". En suero, "
        f'{dato_html("results/v3_7/suero_embudo.json#celda_preregistrada.diferente", contador)}'
        " de "
        f'{dato_html("results/v3_7/suero_embudo.json#celda_preregistrada.n", contador)}'
        " separan el cáncer de la colonoscopia negativa. "
        "Eso no localiza un tumor en este segmento y no muestra una proteína saliendo del tumor.</p>"
    )


def anatomia_html(contador) -> str:
    nota = "Información no disponible en las fuentes del proyecto. Esta tesis no asigna muestras a un segmento del colon."
    regiones = [
        ("ciego", "Ciego"),
        ("ascendente", "Colon ascendente"),
        ("transverso", "Colon transverso"),
        ("descendente", "Colon descendente"),
        ("sigmoide", "Colon sigmoide"),
        ("recto", "Recto"),
    ]
    botones = []
    for ident, nombre in regiones:
        botones.append(
            "<li>"
            f'<button type="button" data-segmento aria-controls="seg-{ident}" aria-pressed="false">{e(nombre)}</button>'
            f'<div class="segmento" id="seg-{ident}"><p><strong>{e(nombre)}.</strong> {e(nota)}</p>'
            f"{_hecho_estudio(contador)}</div>"
            "</li>"
        )
    return (
        '<div class="colon">'
        + _lamina(
            "colon.png", "colon-respaldo.svg",
            "Dibujo del intestino grueso. No marca el lugar de las muestras.",
            661, 913,
            "Anatomía general del intestino grueso. No está a escala y no localiza las muestras.",
        )
        + "<ul class=\"segmentos\">" + "".join(botones) + "</ul></div>"
    )


def comparacion_html(contador) -> str:
    celda = json.loads((V37 / "suero_embudo.json").read_text(encoding="utf-8"))["celda_preregistrada"]
    total = int(celda["n"])
    clases = [
        ("equivalente", "Equivalentes", "TOST cabe en el margen de 1,5 veces."),
        ("indeterminada", "Indeterminadas", "Quedan fuera de las otras dos clases."),
        ("diferente", "Diferentes", "TREAT pasa el corte en el sentido observado."),
    ]
    if sum(int(celda[clave]) for clave, _, _ in clases) != total:
        raise SystemExit("la celda del suero no suma el total")
    base = "results/v3_7/suero_embudo.json#celda_preregistrada."
    tarjetas = []
    trozos = []
    leyenda = []
    etiqueta = []
    for clave, nombre, frase in clases:
        valor = int(celda[clave])
        porcentaje = _porcentaje(valor, total)
        tarjetas.append(
            f'<li class="clase clase-{clave}">'
            f'<p class="clase-nombre"><span class="muestra muestra-{clave}"></span> {e(nombre)}</p>'
            f'<p class="clase-valor">{dato_html(base + clave, contador, traza=False)}</p>'
            f'<p class="clase-pct">{e(porcentaje)}</p>'
            f"<p>{e(frase)}</p></li>"
        )
        if valor > 0:
            trozos.append(f'<span class="trozo trozo-{clave}" style="width:{100 * valor / total:.4f}%"></span>')
        leyenda.append(
            f'<li><span class="muestra muestra-{clave}"></span> {e(nombre)}: {e(porcentaje)}</li>'
        )
        etiqueta.append(f"{nombre} {valor}, {porcentaje}")
    return (
        '<figure class="comparacion">'
        "<figcaption>Distribución de los registros de la celda del suero según las categorías "
        "reportadas en la investigación. Las tres cuentas suman el total. "
        "No es un gráfico de una cohorte nueva.</figcaption>"
        f'<ul class="clases">{"".join(tarjetas)}</ul>'
        f'<p class="clase-total">Total de la celda: {dato_html(base + "n", contador, traza=False)} proteínas.</p>'
        f'<div class="apilada" role="img" aria-label="{e(". ".join(etiqueta))}">'
        + "".join(trozos) + "</div>"
        f'<ul class="leyenda">{"".join(leyenda)}</ul>'
        "<p>Diferentes no ocupa un tramo: el conteo es 0, y un ancho inventado lo haría parecer distinto de cero. "
        "El color nombra la clase. No dice que una clase sea mejor que otra.</p>"
        '<div class="lectura-clases">'
        "<h3>Equivalentes</h3>"
        "<p>Dos pruebas unilaterales que, juntas, declaran que el cambio cabe dentro de un margen. "
        "En esta celda, una proteína es equivalente si TOST cabe en el margen de 1,5 veces.</p>"
        "<h3>Indeterminadas</h3>"
        "<p>Quedan fuera de las otras dos clases. Un resultado negativo con poca potencia se queda "
        "en indeterminado, no en equivalente. No ver un cambio, con esta potencia, no las vuelve equivalentes.</p>"
        "<h3>Diferentes</h3>"
        "<p>TREAT exige que el cambio supere un margen, no solo que sea distinto de cero. "
        "Una proteína es diferente si TREAT pasa el corte en el sentido observado. "
        "El cero es el conteo de esa clase frente a CN. No dice que el tumor no cambie en el tejido.</p>"
        "</div>"
        '<details class="complemento"><summary>Definiciones, referencias y origen de los números</summary>'
        "<p>Los conteos salen de results/v3_7/suero_embudo.json, campo celda_preregistrada. "
        "El porcentaje de cada tarjeta es 100 × esa clase / el total, con un decimal.</p>"
        '<p>TOST está en el <a href="/glosario#tost">glosario</a>. '
        '<a href="/bibliografia#lakens2017">Lakens 2017</a>.</p>'
        '<p>TREAT está en el <a href="/glosario#treat">glosario</a>. '
        '<a href="/bibliografia#mccarthy2009">McCarthy 2009</a>.</p>'
        '<p>La potencia y la celda están en el <a href="/glosario#potencia">glosario</a> y en el de '
        '<a href="/glosario#celda">celda</a>. <a href="/bibliografia#rai2025">Rai 2025</a>.</p>'
        "</details></figure>"
    )


def matrices_html(contador) -> str:
    etiquetas = json.loads((V37 / "k1_etiquetas.json").read_text(encoding="utf-8"))
    contraste = etiquetas["contrastes"]["primario_frente_a_mucosa"]
    tejido_ancho = 100 * int(contraste["numerador"]) / int(contraste["denominador"])
    celda = json.loads((V37 / "suero_embudo.json").read_text(encoding="utf-8"))["celda_preregistrada"]
    suero_ancho = 100 * int(celda["diferente"]) / int(celda["n"])
    tejido = (
        "<h3>Tejido</h3>"
        "<p>Resultado calculado. Tumor y mucosa. El signo coincide en "
        f'{dato_html("results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.numerador", contador)}'
        " de "
        f'{dato_html("results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.denominador", contador)}'
        ". La barra es ese cociente.</p>"
        f'<span class="pista-riel"><span class="pista" style="--meta:{tejido_ancho / 100:.6f}"></span></span>'
    )
    suero = (
        "<h3>Suero</h3>"
        "<p>Resultado calculado. Otra matriz y otras personas. Frente a CN hay "
        f'{dato_html("results/v3_7/suero_embudo.json#celda_preregistrada.diferente", contador)}'
        " diferentes de "
        f'{dato_html("results/v3_7/suero_embudo.json#celda_preregistrada.n", contador)}'
        ". La barra es ese cociente.</p>"
        f'<span class="pista-riel"><span class="pista" style="--meta:{suero_ancho / 100:.6f}"></span></span>'
    )
    return (
        '<div class="matrices">'
        + _lamina(
            "tejido.png", "tejido-respaldo.svg",
            "Dibujo de epitelio cilíndrico con una célula caliciforme.",
            208, 258,
            "Epitelio en general. No es un corte del tumor ni de la mucosa de esta tesis.",
        )
        + _lamina(
            "suero.png", "suero-respaldo.svg",
            "Dibujo de un tubo de sangre.",
            42, 316,
            "Tubo de sangre. El suero medido es la fracción líquida; este dibujo no es la muestra.",
        )
        + f"<article>{tejido}</article><article>{suero}</article>"
        + "<p class=\"matrices-nota\">Las dos figuras no se unen. La tesis no fusiona las matrices "
        "y no muestra una proteína que salga del tumor hacia la sangre.</p>"
        "</div>"
    )


def explorar_html(paginas: list[dict]) -> str:
    queridas = [
        "/metodos/proteomica",
        "/resultados",
        "/proteinas",
        "/metodos/aprendizaje-automatico",
        "/conclusiones",
        "/bibliografia",
    ]
    por_ruta = {pagina["ruta"]: pagina for pagina in paginas}
    tarjetas = []
    for ruta in queridas:
        if ruta == "/bibliografia":
            titulo, resumen = "Bibliografía", "Referencias comprobadas contra su DOI."
        else:
            pagina = por_ruta[ruta]
            titulo, resumen = pagina["titulo"], pagina["resumen"]
        tarjetas.append(
            f'<a class="tarjeta" href="{e(ruta)}"><h2>{e(titulo)}</h2><p>{e(resumen)}</p></a>'
        )
    esquema = (
        '<figure class="esquema-proteina"><svg viewBox="0 0 240 110" role="img">'
        "<title>Esquema de una proteína</title>"
        "<desc>Cinta esquemática. No es la estructura de una proteína de esta tesis.</desc>"
        '<path d="M24 78 C48 18, 96 18, 118 58 S186 96, 214 36" fill="none" stroke="currentColor" '
        'stroke-width="8" stroke-linecap="round"/>'
        "</svg>"
        "<figcaption>Cinta esquemática. No es el modelo de una proteína de la celda. "
        "Cuando hay un modelo de AlphaFold, está en la ficha y se abre con el botón. "
        f"{e(NOTA_CONCEPTUAL)}</figcaption></figure>"
    )
    return esquema + '<div class="tarjetas">' + "".join(tarjetas) + "</div>"


def md_inline(texto: str, bib, glosario, contador) -> str:
    patron = re.compile(
        r"\{\{dato:([^}]+)\}\}|\{\{termino:([^}]+)\}\}|\{\{diagrama:([^}]+)\}\}|\{\{figura:([^}]+)\}\}|"
        r"\{\{embudo:([^}]+)\}\}|\{\{anatomia:([^}]+)\}\}|\{\{matrices\}\}|\{\{comparacion\}\}|"
        r"\{\{celda\}\}|\{\{llamadas\}\}|\{\{explorar\}\}|\[@([^\]]+)\]"
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
        elif hallado.group(4):
            if hallado.group(4) not in FIGURAS:
                raise SystemExit(f"figura sin datos: {hallado.group(4)}")
            piezas.append(FIGURAS[hallado.group(4)])
        elif hallado.group(5):
            piezas.append(embudo_html(hallado.group(5), contador))
        elif hallado.group(6):
            if hallado.group(6) != "colon":
                raise SystemExit(f"anatomía sin definición: {hallado.group(6)}")
            piezas.append(anatomia_html(contador))
        elif hallado.group(0) == "{{matrices}}":
            piezas.append(matrices_html(contador))
        elif hallado.group(0) == "{{comparacion}}":
            piezas.append(comparacion_html(contador))
        elif hallado.group(0) == "{{celda}}":
            piezas.append(celda_html(contador))
        elif hallado.group(0) == "{{llamadas}}":
            piezas.append(llamadas_html())
        elif hallado.group(0) == "{{explorar}}":
            piezas.append('<div class="tarjetas" id="explorar"></div>')
        else:
            piezas.append(cita_html(hallado.group(7), bib))
        pos = hallado.end()
    piezas.append(e(texto[pos:]))
    return "".join(piezas)


def md_a_html(texto: str, bib, glosario, contador) -> tuple[str, list[dict]]:
    indice = []
    bloques = re.split(r"\n\s*\n", texto.strip())
    piezas = []
    for bloque in bloques:
        lineas = bloque.splitlines()
        if bloque.strip().startswith(("{{diagrama:", "{{figura:", "{{embudo:", "{{anatomia:", "{{matrices}}", "{{comparacion}}", "{{celda}}", "{{llamadas}}", "{{explorar}}")):
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
    lugar = {}
    for i, nodo in enumerate(nodos):
        x = 24 + i * (ancho + 28)
        lugar[nodo["id"]] = x
        cuerpo.append(
            f'<g data-id="{e(nodo["id"])}"><title>{e(nodo["texto"])}</title>'
            f'<rect x="{x}" y="70" width="{ancho}" height="64" rx="8" '
            f'fill="none" stroke="currentColor" stroke-width="2"/></g>'
            f'<text x="{x + 16}" y="108" fill="currentColor" font-size="16">{e(nodo["texto"])}</text>'
        )
        filas.append({"nodo": nodo["texto"], "fuente": nodo["fuente"], "href": nodo["href"]})
    aristas = datos.get("aristas") or []
    if not aristas:
        aristas = [
            {"origen": nodos[i]["id"], "destino": nodos[i + 1]["id"], "fuente": nodos[i + 1]["fuente"]}
            for i in range(len(nodos) - 1)
        ]
    indice = {nodo["id"]: i for i, nodo in enumerate(nodos)}
    for arista in aristas:
        x1 = lugar[arista["origen"]] + ancho
        x2 = lugar[arista["destino"]]
        y = 102 if abs(indice[arista["origen"]] - indice[arista["destino"]]) == 1 else 36
        cuerpo.append(
            f'<line data-origen="{e(arista["origen"])}" data-destino="{e(arista["destino"])}" '
            f'x1="{x1}" y1="{y}" x2="{x2}" y2="{y}" stroke="currentColor" stroke-width="2"/>'
        )
        filas.append({"nodo": f'{arista["origen"]} → {arista["destino"]}', "fuente": arista["fuente"], "href": ""})
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
        "<tr><td>"
        + (f'<a href="{e(fila["href"])}">{e(fila["nodo"])}</a>' if fila.get("href") else e(fila["nodo"]))
        + f"</td><td>{e(fila['fuente'])}</td></tr>"
        for fila in filas
    )
    return (
        f'<figure class="grafico" id="{e(ident)}"><h3>{e(datos["pregunta"])}</h3>{svg}'
        f'<details><summary>Tabla de los datos del gráfico</summary>'
        f'<table><thead><tr><th>nodo</th><th>origen</th></tr></thead><tbody>{tabla}</tbody></table></details>'
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


def menu_principal() -> list[dict]:
    return [
        {"nombre": "Inicio", "ruta": "/"},
        {"nombre": "Investigación", "ruta": "/investigacion/contexto"},
        {"nombre": "Metodología", "ruta": "/metodos/proteomica"},
        {"nombre": "Resultados", "ruta": "/resultados"},
        {"nombre": "Explorar proteínas", "ruta": "/proteinas"},
    ]


def menu_secundario(paginas: list[dict]) -> list[dict]:
    cubiertas = {"Inicio", "Investigación", "Métodos", "Resultados", "Proteínas"}
    vistos = []
    for pagina in paginas:
        if pagina["seccion"] in cubiertas:
            continue
        if not any(item["nombre"] == pagina["seccion"] for item in vistos):
            vistos.append({"nombre": pagina["seccion"], "ruta": pagina["ruta"]})
    for nombre, ruta in (("Glosario", "/glosario"), ("Bibliografía", "/bibliografia")):
        if not any(item["ruta"] == ruta for item in vistos):
            vistos.append({"nombre": nombre, "ruta": ruta})
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


def _publicar_hash(origen: Path, carpeta: Path, prefijo: str) -> str:
    datos = origen.read_bytes()
    nombre = prefijo + hashlib.sha256(datos).hexdigest()[:10] + origen.suffix
    carpeta.mkdir(parents=True, exist_ok=True)
    for viejo in carpeta.glob(prefijo + "*" + origen.suffix):
        if viejo.name != nombre:
            viejo.unlink()
    (carpeta / nombre).write_bytes(datos)
    return nombre


def escribir_estaticos() -> tuple[str, str, str, str]:
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
    for viejo in (SITIO / "js").glob("v4.*.js"):
        if viejo.name != js_nombre:
            viejo.unlink()
    (SITIO / "js" / js_nombre).write_text(js, encoding="utf-8")
    molstar_js = _publicar_hash(FUENTE / "js" / "molstar" / "molstar.js", SITIO / "js", "molstar.")
    molstar_css = _publicar_hash(FUENTE / "js" / "molstar" / "molstar.css", SITIO / "css", "molstar.")
    shutil.copy2(FUENTE / "js" / "molstar" / "LICENSE", SITIO / "js" / "molstar-LICENSE.txt")
    destino_fuentes = SITIO / "fuentes"
    destino_fuentes.mkdir(parents=True, exist_ok=True)
    for archivo in (FUENTE / "estilos" / "fuentes").glob("*"):
        if archivo.suffix in {".woff2", ".txt"}:
            shutil.copy2(archivo, destino_fuentes / archivo.name)
    origen_ilustraciones = FUENTE / "ilustraciones"
    if origen_ilustraciones.is_dir():
        destino_ilustraciones = SITIO / "ilustraciones"
        if destino_ilustraciones.exists():
            shutil.rmtree(destino_ilustraciones)
        shutil.copytree(origen_ilustraciones, destino_ilustraciones)
    return css_nombre, js, molstar_js, molstar_css


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
    datos.setdefault("hermanas", [])
    datos.setdefault("portada", False)
    return env.get_template("base.html.j2").render(**datos)


def main() -> None:
    archivar_v37()
    bib = leer_bib()
    glosario = cargar_glosario()
    css_nombre, js, molstar_js, molstar_css = escribir_estaticos()
    cert = json.loads((V37 / "certificado.json").read_text(encoding="utf-8"))
    contador = [0]
    estado = dato_html("results/v3_7/certificado.json#estado", contador)
    sello = dato_html("results/v3_7/certificado.json#sello.commit", contador)
    global FIGURAS
    spec = importlib.util.spec_from_file_location("graficos_sitio", Path(__file__).with_name("graficos.py"))
    graficos = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(graficos)
    FIGURAS = graficos.generar()
    paginas = leer_paginas()
    menu = menu_principal()
    secundario = menu_secundario(paginas)
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
        "secundario": secundario,
        "js_extra": "",
        "molstar_js": molstar_js,
        "molstar_css": molstar_css,
    }
    for pagina in paginas:
        cuerpo, indice = md_a_html(pagina["cuerpo"], bib, glosario, contador)
        if pagina["ruta"] == "/" and 'id="explorar"' in cuerpo:
            cuerpo = cuerpo.replace(
                '<div class="tarjetas" id="explorar"></div>',
                explorar_html(paginas),
                1,
            )
        migas = "Inicio" if pagina["ruta"] == "/" else f'Inicio / {pagina["seccion"]} / {pagina["titulo"]}'
        hermanas = [
            {"nombre": otra["titulo"], "ruta": otra["ruta"]}
            for otra in paginas
            if otra["seccion"] == pagina["seccion"] and otra["ruta"] != pagina["ruta"] and pagina["seccion"] != "Inicio"
        ]
        html_pagina = render(
            env, "base.html.j2",
            titulo=pagina["titulo"],
            cuerpo=cuerpo,
            indice=indice,
            migas=migas,
            fuentes=fuentes_de(cuerpo, bib),
            hermanas=hermanas,
            portada=pagina["ruta"] == "/",
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
        cuerpo="<h1>Glosario</h1><dl>" + "".join(entradas) + "</dl>" + diagrama_html("D11", bib),
        indice=[],
        migas="Inicio / Glosario",
        fuentes=[{"href": "/bibliografia", "texto": "Bibliografía"}],
        **comunes,
    )
    (SITIO / "glosario").mkdir(parents=True, exist_ok=True)
    (SITIO / "glosario" / "index.html").write_text(glosario_html, encoding="utf-8")

    spec_proteinas = importlib.util.spec_from_file_location(
        "proteinas_sitio", Path(__file__).with_name("proteinas.py")
    )
    proteinas = importlib.util.module_from_spec(spec_proteinas)
    spec_proteinas.loader.exec_module(proteinas)
    fichas = proteinas.escribir(render, env, comunes, lambda clave: cita_html(clave, bib))
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
    for ficha in fichas:
        grupos.append(
            f'<li data-buscar="{e(ficha["id"])} {e(ficha["gen"])}">'
            f'<a href="{e(ficha["href"])}">{e(ficha["id"])}</a></li>'
        )
    for titulo, ruta in (("Comparador", "/proteinas/comparar"), ("Candidatas", "/proteinas/candidatas")):
        grupos.append(f'<li data-buscar="{e(titulo)}"><a href="{ruta}">{e(titulo)}</a></li>')
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
