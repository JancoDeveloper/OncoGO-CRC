"""Secciones de anotación, comparador y panorama de las candidatas."""

from __future__ import annotations

import hashlib
import html
import json
import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from oncogo_crc.formato import formatear

REPO = Path(__file__).resolve().parents[2]
SITIO = REPO / "sitio"
FUENTES = REPO / "sitio_fuente" / "fuentes"
V37 = REPO / "results" / "v3_7"
UMBRAL = 0.15
EJEMPLO = ("P02749", "P05155")
ETIQUETA = "Candidata (paso 5): hipótesis para cohortes futuras, no biomarcador"


def e(texto) -> str:
    return html.escape(str(texto), quote=True)


def cargar() -> dict:
    ruta = FUENTES / "anotacion" / "fichas.json"
    if not ruta.is_file():
        raise SystemExit("falta sitio_fuente/fuentes/anotacion/fichas.json")
    return json.loads(ruta.read_text(encoding="utf-8"))


def dato(valor, fuente: str, tipo: str | None = None) -> str:
    if valor is None or (isinstance(valor, float) and math.isnan(valor)):
        return "—"
    if isinstance(valor, bool):
        bruto, visible = ("true" if valor else "false"), ("true" if valor else "false")
    elif isinstance(valor, float):
        bruto = format(valor, ".16g")
        visible = formatear(valor, tipo) if tipo else bruto
    elif isinstance(valor, int):
        bruto = str(valor)
        visible = formatear(valor, tipo or "entero")
    else:
        bruto, visible = str(valor), str(valor)
    return (
        f'<data value="{e(bruto)}" data-fuente="{e(fuente)}">{e(visible)}</data>'
        f'<details class="evidencia"><summary>fuente</summary><p>{e(fuente)}</p></details>'
    )


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    union = a | b
    return len(a & b) / len(union) if union else 0.0


def ids_de(registro: dict, campo: str, aspecto: str = "") -> set[str]:
    salida = set()
    for item in registro.get(campo) or []:
        if aspecto and item.get("aspecto") != aspecto:
            continue
        if item.get("id"):
            salida.add(item["id"])
    return salida


def base_de(token: str) -> str:
    parte, _, sufijo = token.rpartition("-")
    if sufijo.isdigit() and parte:
        return parte
    return token


def registro_de(documento: dict, token: str) -> dict:
    tabla = documento["por_acceso"]
    if token in tabla:
        return tabla[token]
    return tabla.get(base_de(token), {"estado": "no_disponible", "motivo": "sin ficha", "acceso": token})


def atribucion(documento: dict, cita, con_string: bool = False) -> str:
    version = ""
    for registro in documento["por_acceso"].values():
        if registro.get("version_uniprot"):
            version = registro["version_uniprot"]
            break
    modelos = documento.get("modelos") or {}
    texto = (
        f"<p>UniProt {e(version)}, CC BY 4.0. {cita('uniprot2025')} "
        f"AlphaFold DB, CC BY 4.0. Suma de BinaryCIF alojados: {e(modelos.get('mb', '—'))} MB. "
        f"{cita('jumper2021')} {cita('varadi2024')} "
        f"InterPro y Pfam, por las referencias cruzadas de UniProt, CC BY 4.0. {cita('blum2025')} "
        f"Gene Ontology, obo {e(documento.get('obo_version', ''))}, que no es el de la tesis. "
        f"{cita('goconsortium2023')} Reactome, CC BY 4.0. {cita('milacic2024')} "
        "PDB, CC0: solo la lista, sin ficheros. "
        f"{e(documento.get('nota_go', ''))}</p>"
    )
    if con_string:
        texto += (
            "<p>STRING, humano 9606, puntuación combinada de al menos 0,7. "
            "Son asociaciones conocidas, con minería de texto, no un resultado de esta tesis. "
            f"CC BY 4.0. {cita('szklarczyk2023')}</p>"
        )
    return texto


def secuencia_html(secuencia: str, rangos: list[dict]) -> str:
    marcas = [""] * (len(secuencia) + 1)
    for rango in rangos:
        for i in range(int(rango["inicio"]), int(rango["fin"]) + 1):
            if 1 <= i <= len(secuencia) and rango["tipo"] == "signal":
                marcas[i] = "signal"
            elif 1 <= i <= len(secuencia) and not marcas[i]:
                marcas[i] = "domain"
    piezas = ['<p class="secuencia">']
    for inicio in range(1, len(secuencia) + 1, 10):
        fin = min(inicio + 9, len(secuencia))
        piezas.append(f'<span class="num">{inicio:4d}</span> ')
        i = inicio
        while i <= fin:
            clase = marcas[i]
            j = i
            while j <= fin and marcas[j] == clase:
                j += 1
            trozo = e(secuencia[i - 1:j - 1])
            piezas.append(f'<span class="{clase}">{trozo}</span>' if clase else trozo)
            i = j
        piezas.append("\n")
    piezas.append("</p>")
    return "".join(piezas)


def secciones(ficha: dict, documento: dict, cita, comunes: dict) -> str:
    token = ficha["id"].split(";")[0]
    registro = registro_de(documento, token)
    if registro.get("estado") != "ok":
        return (
            '<h2 id="secuencia">Secuencia</h2>'
            f'<p>Anotación no disponible. {e(registro.get("motivo", ""))}</p>'
        )
    fuente = f"sitio_fuente/fuentes/anotacion/fichas.json#por_acceso.{token}"
    secuencia = registro.get("secuencia") or ""
    partes = [
        '<h2 id="secuencia">Secuencia</h2>',
        f"<p>Longitud {dato(registro.get('longitud'), fuente + '.longitud', 'entero')}. "
        f"Masa {dato(registro.get('masa'), fuente + '.masa', 'entero')}. "
        f'<a href="/fuentes/fasta/{e(token)}.fasta">FASTA</a>. {cita("uniprot2025")}</p>',
    ]
    if ";" in ficha["id"]:
        partes.append(f"<p>El grupo de la celda es {e(ficha['id'])}. La secuencia es la de {e(token)}.</p>")
    partes.append("<p>El fondo marca el péptido señal. El subrayado marca un dominio de la entrada.</p>")
    partes.append(secuencia_html(secuencia, registro.get("rangos") or []))
    alfa = registro.get("alfafold") or {}
    base = base_de(token)
    fuente_af = f"sitio_fuente/fuentes/anotacion/alfafold/{base}.json#0"
    partes.append('<h2 id="estructura">Estructura</h2>')
    partes.append(
        "<p>Es una predicción. Un pLDDT bajo, por debajo de 70, indica una zona poco fiable: "
        "no se lee como un dominio bien resuelto. "
        f"{cita('jumper2021')} {cita('varadi2024')} {cita('sehnal2021')}</p>"
    )
    if alfa.get("estado") != "ok" or alfa.get("plddt_medio") is None:
        partes.append(f"<p>AlphaFold no tiene modelo para {e(base)}. {e(alfa.get('motivo', ''))}</p>")
    else:
        if alfa.get("modelo_es_canonico"):
            partes.append("<p>El modelo es el de la secuencia canónica, no el de esta isoforma.</p>")
        partes.append(
            "<p>pLDDT medio "
            + dato(alfa.get("plddt_medio"), fuente_af + ".globalMetricValue", "auc")
            + ". Muy alto "
            + dato(alfa.get("muy_alto"), fuente_af + ".fractionPlddtVeryHigh", "auc")
            + ". Fiable "
            + dato(alfa.get("fiable"), fuente_af + ".fractionPlddtConfident", "auc")
            + ". Bajo "
            + dato(alfa.get("bajo"), fuente_af + ".fractionPlddtLow", "auc")
            + ". Muy bajo "
            + dato(alfa.get("muy_bajo"), fuente_af + ".fractionPlddtVeryLow", "auc")
            + ".</p>"
        )
        if alfa.get("modelo_local") and comunes.get("molstar_js"):
            ancla = "visor-" + re.sub(r"[^A-Za-z0-9_-]", "-", token)
            partes.append(
                f'<p><button type="button" data-ver-estructura="{ancla}" '
                f'data-modelos="/{e(alfa["modelo_local"])}" '
                f'data-script="/js/{e(comunes["molstar_js"])}" '
                f'data-estilo="/css/{e(comunes["molstar_css"])}">Ver estructura</button></p>'
                f'<div id="{ancla}"></div>'
                "<p>Sin JavaScript queda la cifra de pLDDT. El visor no forma parte de la página hasta ese botón.</p>"
            )
        else:
            partes.append("<p>El BinaryCIF de este acceso no está alojado.</p>")
    if registro.get("pdb"):
        filas = []
        for pdb in registro["pdb"]:
            ident = pdb["id"]
            filas.append(
                "<tr><td>"
                f'<a href="https://www.ebi.ac.uk/pdbe/entry/pdb/{e(ident)}">{e(ident)}</a> '
                f'<a href="https://www.rcsb.org/structure/{e(ident)}">RCSB</a>'
                f"</td><td>{e(pdb.get('metodo', ''))}</td><td>{e(pdb.get('resolucion', ''))}</td>"
                f"<td>{e(pdb.get('cadenas', ''))}</td></tr>"
            )
        partes.append(
            "<table><thead><tr><th>PDB</th><th>Método</th><th>Resolución</th><th>Cadenas</th></tr></thead>"
            f"<tbody>{''.join(filas)}</tbody></table>"
        )
    else:
        partes.append("<p>UniProt no lista estructuras PDB para este acceso.</p>")
    partes.append('<h2 id="go">Términos GO</h2>')
    partes.append("<p>La vista directa exige una evidencia distinta de IEA. La propagación añade ancestros is_a del obo descargado, que no es el de la tesis.</p>")
    for aspecto, titulo in (("BP", "Proceso"), ("MF", "Función"), ("CC", "Componente")):
        filas = []
        for termino in registro.get("go_directo") or []:
            if termino.get("aspecto") != aspecto:
                continue
            filas.append(
                f'<tr><td><a href="https://www.ebi.ac.uk/QuickGO/term/{e(termino["id"])}">{e(termino["id"])}</a></td>'
                f'<td>{e(termino.get("nombre", ""))}</td><td>{e(termino.get("evidencia", ""))}</td></tr>'
            )
        cuerpo = "".join(filas) or "<tr><td colspan=\"3\">Sin términos directos en este aspecto.</td></tr>"
        partes.append(
            f"<h3>{titulo}</h3><table><thead><tr><th>Término</th><th>Nombre</th><th>Evidencia</th></tr></thead>"
            f"<tbody>{cuerpo}</tbody></table>"
        )
    n_prop = len(registro.get("go_propagado") or [])
    partes.append(f"<p>Términos en la vista propagada, contando los directos: {e(n_prop)}.</p>")
    partes.append('<h2 id="dominios">Dominios y vías</h2>')
    partes.append(f"<p>{cita('blum2025')} {cita('milacic2024')}</p>")
    for campo, titulo in (("interpro", "InterPro"), ("pfam", "Pfam"), ("reactome", "Reactome")):
        filas = []
        for item in registro.get(campo) or []:
            filas.append(f"<tr><td>{e(item.get('id', ''))}</td><td>{e(item.get('nombre', ''))}</td></tr>")
        cuerpo = "".join(filas) or "<tr><td colspan=\"2\">Sin entradas.</td></tr>"
        partes.append(f"<h3>{titulo}</h3><table><tbody>{cuerpo}</tbody></table>")
    partes.append(atribucion(documento, cita))
    return "\n".join(partes)


def veinte(fichas: list[dict]) -> list[dict]:
    vistos = set()
    salida = []
    for ficha in fichas:
        if not ficha.get("candidata") or ficha["gen"] in vistos:
            continue
        vistos.add(ficha["gen"])
        salida.append(ficha)
    return salida


def alinear_pares(pares: list[tuple[str, str]]) -> list[dict]:
    with tempfile.TemporaryDirectory() as carpeta:
        entrada = Path(carpeta) / "pares.json"
        salida = Path(carpeta) / "salida.json"
        entrada.write_text(json.dumps([{"a": a, "b": b} for a, b in pares]), encoding="utf-8")
        script = (
            "const fs=require('fs');"
            "eval(fs.readFileSync(process.argv[1],'utf8'));"
            "const pares=JSON.parse(fs.readFileSync(process.argv[2],'utf8'));"
            "fs.writeFileSync(process.argv[3], JSON.stringify(pares.map(p => Alineamiento.alinear(p.a, p.b))));"
        )
        subprocess.run(
            ["node", "-e", script, str(REPO / "sitio_fuente" / "js" / "alineamiento.js"), str(entrada), str(salida)],
            check=True,
        )
        return json.loads(salida.read_text(encoding="utf-8"))


def figura(fid: str, titulo: str, desc: str, cuerpo: str, filas: list[dict], pregunta: str, viewbox: str) -> str:
    svg = (
        f'<svg viewBox="{viewbox}" role="img" class="lienzo-svg">'
        f"<title>{e(titulo)}</title><desc>{e(desc)}</desc>{cuerpo}</svg>"
    )
    destino = SITIO / "figuras"
    destino.mkdir(parents=True, exist_ok=True)
    (destino / f"{fid}.svg").write_text(svg, encoding="utf-8")
    columnas = list(filas[0]) if filas else ["dato"]
    lineas = ["\t".join(columnas)]
    for fila in filas:
        lineas.append("\t".join(str(fila[col]) for col in columnas))
    (destino / f"{fid}.tsv").write_text("\n".join(lineas) + "\n", encoding="utf-8")
    sello = {
        "fuentes": ["sitio_fuente/fuentes/anotacion/fichas.json"],
        "svg_sha256": hashlib.sha256(svg.encode()).hexdigest(),
        "generador": "scripts/sitio/vistas_proteina.py",
        "umbral_jaccard": UMBRAL,
        "agrupamiento": "enlace promedio, distancia 1 menos Jaccard directo",
    }
    (destino / f"{fid}.sello.json").write_text(json.dumps(sello, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    cabecera = "".join(f"<th>{e(col)}</th>" for col in columnas)
    cuerpo_tabla = []
    for fila in filas:
        cuerpo_tabla.append("<tr>" + "".join(f"<td>{e(fila[col])}</td>" for col in columnas) + "</tr>")
    return (
        f'<figure class="grafico" id="{e(fid)}"><h3>{e(pregunta)}</h3>{svg}'
        f"<table><thead><tr>{cabecera}</tr></thead><tbody>{''.join(cuerpo_tabla)}</tbody></table>"
        f'<p><a href="/figuras/{e(fid)}.svg" download>Descargar SVG</a> '
        f'<a href="/figuras/{e(fid)}.tsv" download>Descargar TSV</a></p></figure>'
    )


def guardar_fasta(documento: dict) -> None:
    carpeta = SITIO / "fuentes" / "fasta"
    carpeta.mkdir(parents=True, exist_ok=True)
    for token, registro in documento["por_acceso"].items():
        secuencia = registro.get("secuencia") or ""
        if not secuencia:
            continue
        (carpeta / f"{token}.fasta").write_text(f">{token}\n{secuencia}\n", encoding="utf-8")


def copiar_modelos() -> None:
    origen = FUENTES / "anotacion" / "modelos"
    destino = SITIO / "fuentes" / "modelos"
    if destino.exists():
        shutil.rmtree(destino)
    if origen.is_dir() and any(origen.glob("*.bcif")):
        shutil.copytree(origen, destino)


def catalogo(documento: dict, fichas: list[dict], comunes: dict) -> dict:
    por_id = {ficha["id"]: ficha for ficha in fichas}
    suero = pd.read_csv(V37 / "suero_k4_crc_cn.tsv", sep="\t").set_index("protein_id")
    prio = pd.read_csv(V37 / "priorizacion.tsv", sep="\t").set_index("gen")
    proteinas = {}
    for token, registro in documento["por_acceso"].items():
        if registro.get("estado") != "ok":
            continue
        ficha = next((item for item in fichas if item["id"].split(";")[0] == token), None)
        protein_id = ficha["id"] if ficha else token
        clase, log2 = "", None
        if protein_id in suero.index:
            clase = str(suero.loc[protein_id, "clase"])
            log2 = float(suero.loc[protein_id, "log2fc"])
        gen = ficha["gen"] if ficha else registro.get("gen_uniprot") or ""
        paso = ""
        if gen and gen in prio.index:
            paso = str(int(prio.loc[gen, "ultimo_paso"]))
        alfa = registro.get("alfafold") or {}
        proteinas[token] = {
            "id": token,
            "gen": gen,
            "nombre": registro.get("nombre") or "",
            "secuencia": registro.get("secuencia") or "",
            "longitud": registro.get("longitud"),
            "masa": registro.get("masa"),
            "revisada": bool(registro.get("revisada")),
            "go": registro.get("go_directo") or [],
            "go_propagado": registro.get("go_propagado") or [],
            "interpro": [f"{item.get('id', '')} {item.get('nombre', '')}".strip() for item in registro.get("interpro") or []],
            "pfam": [f"{item.get('id', '')} {item.get('nombre', '')}".strip() for item in registro.get("pfam") or []],
            "reactome": [f"{item.get('id', '')} {item.get('nombre', '')}".strip() for item in registro.get("reactome") or []],
            "plddt": alfa.get("plddt_medio"),
            "modelo": ("/" + alfa["modelo_local"]) if alfa.get("modelo_local") else "",
            "clase": clase,
            "log2fc_cn": log2,
            "paso": paso,
            "candidata": bool(ficha and ficha.get("candidata")),
        }
    return {
        "proteinas": proteinas,
        "molstar_js": "/js/" + comunes.get("molstar_js", ""),
        "molstar_css": "/css/" + comunes.get("molstar_css", ""),
        "umbral": UMBRAL,
    }


def html_par(documento: dict, a: str, b: str, cita, alineacion: dict | None) -> str:
    ra, rb = registro_de(documento, a), registro_de(documento, b)
    ja = ids_de(ra, "go_directo")
    jb = ids_de(rb, "go_directo")
    pa = ids_de(ra, "go_propagado")
    pb = ids_de(rb, "go_propagado")
    directo = jaccard(ja, jb)
    propagado = jaccard(pa, pb)
    fuente = "sitio/fuentes/calculos_comparador.json"
    partes = [
        f"<p>Jaccard directo {dato(directo, fuente + f'#{a}|{b}.directo', 'auc')}. "
        f"Jaccard propagado {dato(propagado, fuente + f'#{a}|{b}.propagado', 'auc')}. "
        f"{cita('goconsortium2023')}</p>"
    ]
    for aspecto, titulo in (("BP", "Proceso"), ("MF", "Función"), ("CC", "Componente")):
        sa, sb = ids_de(ra, "go_directo", aspecto), ids_de(rb, "go_directo", aspecto)
        partes.append(
            f"<h3>{titulo}</h3><p>Solo {e(a)}: {e(len(sa - sb))}. "
            f"Compartidos: {e(len(sa & sb))}. Solo {e(b)}: {e(len(sb - sa))}.</p>"
        )
    if alineacion:
        columnas = alineacion["columnas"] or 1
        identidad = alineacion["identidades"] / columnas
        similitud = alineacion["similares"] / columnas
        partes.append(
            f"<p>Puntuación {dato(alineacion['puntuacion'], fuente + f'#{a}|{b}.puntuacion', 'auc')}. "
            f"Identidad {dato(100 * identidad, fuente + f'#{a}|{b}.identidad', 'porcentaje')}. "
            f"Similitud {dato(100 * similitud, fuente + f'#{a}|{b}.similitud', 'porcentaje')}. "
            f"Hueco de apertura 10 y de extensión 0,5. Matriz BLOSUM62. "
            f"{cita('needleman1970')} {cita('henikoff1992')} {cita('cock2009')}</p>"
            f'<p class="secuencia">{e(alineacion["arriba"])}\n{e(alineacion["abajo"])}</p>'
        )
    return "\n".join(partes)


def cerrar(render, env, comunes: dict, cita, fichas: list[dict]) -> None:
    documento = cargar()
    guardar_fasta(documento)
    copiar_modelos()
    datos = catalogo(documento, fichas, comunes)
    (SITIO / "fuentes" / "comparador.json").write_text(
        json.dumps(datos, ensure_ascii=False), encoding="utf-8"
    )
    ra = registro_de(documento, EJEMPLO[0])
    rb = registro_de(documento, EJEMPLO[1])
    pares = [(ra.get("secuencia") or "A", rb.get("secuencia") or "A")]
    grupo = veinte(fichas)
    secuencias = []
    for ficha in grupo:
        registro = registro_de(documento, ficha["id"].split(";")[0])
        secuencias.append(registro.get("secuencia") or "")
    for i in range(len(secuencias)):
        for j in range(i + 1, len(secuencias)):
            pares.append((secuencias[i], secuencias[j]))
    alineaciones = alinear_pares(pares)
    ejemplo = alineaciones[0]
    indice = 1
    identidad = [[None] * len(grupo) for _ in grupo]
    for i in range(len(grupo)):
        identidad[i][i] = 1.0
        for j in range(i + 1, len(grupo)):
            aln = alineaciones[indice]
            indice += 1
            valor = (aln["identidades"] / aln["columnas"]) if aln["columnas"] else 0.0
            identidad[i][j] = identidad[j][i] = valor
    calculos = {
        f"{EJEMPLO[0]}|{EJEMPLO[1]}": {
            "directo": jaccard(ids_de(ra, "go_directo"), ids_de(rb, "go_directo")),
            "propagado": jaccard(ids_de(ra, "go_propagado"), ids_de(rb, "go_propagado")),
            "puntuacion": ejemplo["puntuacion"],
            "identidad": (100 * ejemplo["identidades"] / ejemplo["columnas"]) if ejemplo["columnas"] else 0,
            "similitud": (100 * ejemplo["similares"] / ejemplo["columnas"]) if ejemplo["columnas"] else 0,
        }
    }
    (SITIO / "fuentes" / "calculos_comparador.json").write_text(
        json.dumps(calculos, ensure_ascii=False), encoding="utf-8"
    )
    tokens = sorted(
        token for token, registro in documento["por_acceso"].items()
        if registro.get("estado") == "ok"
    )
    muestra = []
    for i in range(10):
        a, b = tokens[i], tokens[i + 1]
        muestra.append({
            "a": a,
            "b": b,
            "directo": jaccard(ids_de(registro_de(documento, a), "go_directo"), ids_de(registro_de(documento, b), "go_directo")),
        })
    (SITIO / "fuentes" / "jaccard_muestra.json").write_text(
        json.dumps(muestra, ensure_ascii=False), encoding="utf-8"
    )
    extra = dict(comunes)
    extra["js_extra"] = (
        (REPO / "sitio_fuente" / "js" / "alineamiento.js").read_text(encoding="utf-8")
        + "\n"
        + (REPO / "sitio_fuente" / "js" / "comparar.js").read_text(encoding="utf-8")
    )
    cuerpo = (
        "<h1>Comparador</h1>"
        "<p>Describe conocimiento previo: secuencia, estructura, GO, dominios y vías. "
        "No es un resultado de la tesis. Compartir términos, dominios o secuencia no significa "
        "que dos proteínas se comporten igual en la sangre. Las proteínas más estudiadas tienen más anotaciones.</p>"
        '<p><a href="/proteinas/candidatas">Las 20 candidatas</a>. '
        f'<a href="/proteinas/comparar?p={EJEMPLO[0]},{EJEMPLO[1]}">Ejemplo {EJEMPLO[0]} y {EJEMPLO[1]}</a>. '
        '<a href="/proteinas/comparar?p=P05109,P06702">Calprotectina (S100A8 y S100A9)</a>.</p>'
        '<form action="/proteinas/comparar" method="get"><label for="p">De 2 a 4 accesos, separados por comas</label> '
        f'<input id="p" name="p" value="{EJEMPLO[0]},{EJEMPLO[1]}"> <button type="submit">Comparar</button></form>'
        "<p>Sin JavaScript se muestra el par de ejemplo. Otro par se calcula en el navegador a partir de la anotación guardada.</p>"
        '<nav class="pestanias" aria-label="Vistas de la comparación">'
        '<a href="#vista-resumen">Resumen</a><a href="#vista-go">GO</a>'
        '<a href="#vista-secuencia">Secuencia</a><a href="#vista-dominios">Dominios</a>'
        '<a href="#vista-estructura">Estructura</a></nav>'
        f'<div id="comparacion" data-par="{EJEMPLO[0]},{EJEMPLO[1]}">'
        '<h2 id="vista-resumen">Resumen</h2>'
        f"<p>{e(EJEMPLO[0])} y {e(EJEMPLO[1])}, lado a lado, en la anotación guardada.</p>"
        '<h2 id="vista-go">GO</h2>'
        + html_par(documento, EJEMPLO[0], EJEMPLO[1], cita, None)
        + '<h2 id="vista-secuencia">Secuencia</h2>'
        + html_par(documento, EJEMPLO[0], EJEMPLO[1], cita, ejemplo)
        + '<h2 id="vista-dominios">Dominios</h2>'
        + dominios_lado(documento, EJEMPLO)
        + '<h2 id="vista-estructura">Estructura</h2>'
        + estructura_lado(documento, EJEMPLO, comunes)
        + "</div>"
        + atribucion(documento, cita, con_string=False)
    )
    pagina = render(
        env, "base.html.j2",
        titulo="Comparador",
        cuerpo=cuerpo,
        indice=[
            {"id": "vista-resumen", "texto": "Resumen"},
            {"id": "vista-go", "texto": "GO"},
            {"id": "vista-secuencia", "texto": "Secuencia"},
            {"id": "vista-dominios", "texto": "Dominios"},
            {"id": "vista-estructura", "texto": "Estructura"},
        ],
        migas="Inicio / Proteínas / Comparador",
        fuentes=[{"href": "/fuentes/comparador.json", "texto": "Anotación para el comparador"}],
        **extra,
    )
    destino = SITIO / "proteinas" / "comparar" / "index.html"
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_text(pagina, encoding="utf-8")
    cuerpo_c = panorama(documento, grupo, identidad, cita)
    pagina_c = render(
        env, "base.html.j2",
        titulo="Candidatas",
        cuerpo=cuerpo_c,
        indice=[
            {"id": "lectura", "texto": "Lectura"},
            {"id": "mapa", "texto": "Mapa"},
            {"id": "red", "texto": "Red"},
            {"id": "identidad", "texto": "Identidad"},
            {"id": "compartidos", "texto": "Términos"},
        ],
        migas="Inicio / Proteínas / Candidatas",
        fuentes=[{"href": "/fuentes/comparador.json", "texto": "Anotación"}],
        **comunes,
    )
    destino_c = SITIO / "proteinas" / "candidatas" / "index.html"
    destino_c.parent.mkdir(parents=True, exist_ok=True)
    destino_c.write_text(pagina_c, encoding="utf-8")
    anotar_omitido(documento)


def dominios_lado(documento: dict, accesos: tuple[str, str]) -> str:
    filas = []
    for acceso in accesos:
        registro = registro_de(documento, acceso)
        for campo in ("interpro", "pfam", "reactome"):
            for item in registro.get(campo) or []:
                filas.append(f"<tr><td>{e(acceso)}</td><td>{e(campo)}</td><td>{e(item.get('id', ''))}</td><td>{e(item.get('nombre', ''))}</td></tr>")
    cuerpo = "".join(filas) or "<tr><td colspan=\"4\">Sin entradas.</td></tr>"
    return f"<table><thead><tr><th>Acceso</th><th>Base</th><th>Id</th><th>Nombre</th></tr></thead><tbody>{cuerpo}</tbody></table>"


def estructura_lado(documento: dict, accesos: tuple[str, str], comunes: dict) -> str:
    modelos = []
    partes = ["<table><tbody>"]
    for acceso in accesos:
        alfa = (registro_de(documento, acceso).get("alfafold") or {})
        partes.append(f"<tr><th>{e(acceso)}</th><td>{e(alfa.get('plddt_medio', 'sin modelo'))}</td></tr>")
        if alfa.get("modelo_local"):
            modelos.append("/" + alfa["modelo_local"])
    partes.append("</tbody></table>")
    if len(modelos) >= 2 and comunes.get("molstar_js"):
        partes.append(
            '<p><button type="button" data-ver-estructura="visor-ejemplo" data-modelos="'
            + e(",".join(modelos[:2]))
            + f'" data-script="/js/{e(comunes["molstar_js"])}" data-estilo="/css/{e(comunes["molstar_css"])}">'
            "Ver estructuras</button></p><div id=\"visor-ejemplo\"></div>"
        )
    partes.append("<p>No se superponen estructuras: no hay un método de superposición citado.</p>")
    return "\n".join(partes)


def panorama(documento: dict, grupo: list[dict], identidad: list[list], cita) -> str:
    n = len(grupo)
    matriz = []
    for i, ficha in enumerate(grupo):
        fila = []
        for j, otra in enumerate(grupo):
            a = ids_de(registro_de(documento, ficha["id"].split(";")[0]), "go_directo")
            b = ids_de(registro_de(documento, otra["id"].split(";")[0]), "go_directo")
            fila.append(jaccard(a, b))
        matriz.append(fila)
    dist = [[0.0 if i == j else 1 - matriz[i][j] for j in range(n)] for i in range(n)]
    if n >= 2:
        orden = list(leaves_list(linkage(squareform(dist, checks=False), method="average")))
    else:
        orden = list(range(n))
    lado = 28
    ancho = 160 + n * lado
    alto = 40 + n * lado
    cuerpo = []
    filas = []
    for y, i in enumerate(orden):
        cuerpo.append(f'<text x="4" y="{48 + y * lado}" fill="currentColor" font-size="10">{e(grupo[i]["gen"])}</text>')
        for x, j in enumerate(orden):
            valor = matriz[i][j]
            cuerpo.append(
                f'<rect x="{150 + x * lado}" y="{32 + y * lado}" width="{lado - 2}" height="{lado - 2}" '
                f'fill="none" stroke="currentColor" stroke-opacity="{max(0.15, valor):.3f}"/>'
                f'<text x="{152 + x * lado}" y="{46 + y * lado}" fill="currentColor" font-size="7">'
                f'{valor:.2f}</text>'
            )
            filas.append({"gen_a": grupo[i]["gen"], "gen_b": grupo[j]["gen"], "jaccard": f"{valor:.6f}"})
    mapa = figura(
        "r07", "Jaccard GO de las candidatas",
        "Cada celda es el Jaccard de la vista directa, sin IEA. El orden es enlace promedio.",
        "".join(cuerpo), filas, "¿Qué términos directos comparten?",
        f"0 0 {ancho} {max(alto, 80)}",
    )
    radio = 180
    cx, cy = 280, 230
    red = []
    filas_red = []
    posiciones = []
    for k, i in enumerate(orden):
        angulo = 2 * math.pi * k / max(n, 1)
        posiciones.append((cx + radio * math.cos(angulo), cy + radio * math.sin(angulo)))
    for a in range(n):
        for b in range(a + 1, n):
            if matriz[a][b] < UMBRAL:
                continue
            ia, ib = orden.index(a), orden.index(b)
            x1, y1 = posiciones[ia]
            x2, y2 = posiciones[ib]
            red.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="currentColor"/>')
            filas_red.append({"gen_a": grupo[a]["gen"], "gen_b": grupo[b]["gen"], "jaccard": f"{matriz[a][b]:.6f}"})
    suero = pd.read_csv(V37 / "suero_k4_crc_cn.tsv", sep="\t").set_index("protein_id")
    for k, i in enumerate(orden):
        x, y = posiciones[k]
        clase = "indeterminada"
        if grupo[i]["id"] in suero.index:
            clase = str(suero.loc[grupo[i]["id"], "clase"])
        red.append(
            f'<circle class="nodo-{e(clase)}" cx="{x:.1f}" cy="{y:.1f}" r="14" fill="none" stroke-width="3"/>'
            f'<text x="{x + 16:.1f}" y="{y:.1f}" fill="currentColor" font-size="11">{e(grupo[i]["gen"])}</text>'
        )
    red_svg = figura(
        "r08", "Red de las candidatas",
        "Una arista existe si el Jaccard directo es al menos 0,15. El trazo distingue la clase en CP frente a CN.",
        "".join(red), filas_red or [{"gen_a": "", "gen_b": "", "jaccard": ""}],
        "¿Qué pares superan el umbral fijado?",
        "0 0 640 460",
    )
    cabecera = "".join(f"<th>{e(ficha['gen'])}</th>" for ficha in grupo)
    filas_id = []
    for i, ficha in enumerate(grupo):
        celdas = []
        for j in range(n):
            celdas.append(f"<td>{e(f'{100 * identidad[i][j]:.1f}')}</td>")
        filas_id.append(f"<tr><th>{e(ficha['gen'])}</th>{''.join(celdas)}</tr>")
    conteo: dict[str, dict] = {}
    for ficha in grupo:
        for termino in registro_de(documento, ficha["id"].split(";")[0]).get("go_directo") or []:
            actual = conteo.setdefault(termino["id"], {"nombre": termino.get("nombre", ""), "n": 0, "aspecto": termino.get("aspecto", "")})
            actual["n"] += 1
    comunes = sorted(conteo.items(), key=lambda item: (-item[1]["n"], item[0]))[:25]
    filas_go = "".join(
        f"<tr><td>{e(ident)}</td><td>{e(info['aspecto'])}</td><td>{e(info['nombre'])}</td><td>{e(info['n'])}</td></tr>"
        for ident, info in comunes
    )
    string = documento.get("string") or {}
    bloque_string = ""
    if string.get("estado") == "ok":
        crudo = json.loads((FUENTES / "anotacion" / "string.json").read_text(encoding="utf-8"))
        filas_s = []
        for arista in crudo:
            filas_s.append(
                "<tr>"
                f"<td>{e(arista.get('preferredName_A', ''))}</td>"
                f"<td>{e(arista.get('preferredName_B', ''))}</td>"
                f"<td>{dato(arista.get('score'), 'sitio_fuente/fuentes/anotacion/string.json', 'auc')}</td>"
                f"<td>{dato(arista.get('tscore'), 'sitio_fuente/fuentes/anotacion/string.json', 'auc')}</td>"
                f"<td>{dato(arista.get('escore'), 'sitio_fuente/fuentes/anotacion/string.json', 'auc')}</td>"
                f"<td>{dato(arista.get('dscore'), 'sitio_fuente/fuentes/anotacion/string.json', 'auc')}</td>"
                "</tr>"
            )
        bloque_string = (
            "<h2>STRING</h2><p>Asociaciones conocidas entre estas veinte, con minería de texto. "
            "No es un resultado de la tesis.</p>"
            "<table><thead><tr><th>A</th><th>B</th><th>Combinada</th><th>Texto</th><th>Experimento</th><th>Bases</th></tr></thead>"
            f"<tbody>{''.join(filas_s)}</tbody></table>"
        )
    return (
        "<h1>Candidatas</h1>"
        f'<p id="lectura">La priorización deja {e(n)} genes en el paso 5 y ninguno en el paso 6. '
        f"{e(ETIQUETA)} {cita('tanaka2024')}</p>"
        "<p>La tesis lee ese paso como proteínas abundantes del plasma, no como proteínas que salgan del tumor. "
        "En los términos directos, lo que más se repite es localización extracelular: "
        "región extracelular, exosoma y micropartícula sanguínea. "
        f"Eso no añade un hallazgo tumoral. {cita('goconsortium2023')}</p>"
        f'<h2 id="mapa">Mapa</h2>{mapa}'
        f'<h2 id="red">Red</h2><p>El umbral 0,15 quedó fijado en docs/sitio_v4_diseno.md antes de dibujar la red.</p>{red_svg}'
        f'<h2 id="identidad">Identidad</h2><p>Porcentaje de columnas idénticas sobre el alineamiento global, huecos incluidos.</p>'
        f"<table><thead><tr><th>Gen</th>{cabecera}</tr></thead><tbody>{''.join(filas_id)}</tbody></table>"
        f'<h2 id="compartidos">Términos</h2><table><thead><tr><th>GO</th><th>Aspecto</th><th>Nombre</th><th>n</th></tr></thead>'
        f"<tbody>{filas_go}</tbody></table>"
        f"{bloque_string}{atribucion(documento, cita, con_string=string.get('estado') == 'ok')}"
        f"<p>Cada gen del paso 5 entra una vez. En esta página: {e(', '.join(item['gen'] for item in grupo))}.</p>"
    )


def anotar_omitido(documento: dict) -> None:
    ruta = REPO / "docs" / "sitio_v4_omitido.md"
    actual = ruta.read_text(encoding="utf-8")
    marca = "\n## Anotación y modelos\n"
    if marca in actual:
        actual = actual.split(marca)[0].rstrip() + "\n"
    modelos = documento.get("modelos") or {}
    string = documento.get("string") or {}
    texto = (
        f"\nLos BinaryCIF de AlphaFold suman {modelos.get('mb', '—')} MB "
        f"({modelos.get('n', 0)} ficheros). El corte fijado era 80 MB. "
        f"{'Quedaron alojados.' if modelos.get('alojados') else 'No se alojaron: ' + str(modelos.get('motivo', ''))}\n\n"
        "El mapa de las 284 no se dibuja: un SVG con una celda por par supera 1,5 MB y no se lee. "
        "Queda el panorama de las 20.\n\n"
        "No se calculó la similitud semántica de Wang.\n\n"
        "No se superponen estructuras.\n\n"
    )
    if string.get("estado") != "ok":
        texto += f"STRING no se incorporó: {string.get('motivo', 'sin respuesta')}.\n"
    ruta.write_text(actual + marca + texto, encoding="utf-8")
