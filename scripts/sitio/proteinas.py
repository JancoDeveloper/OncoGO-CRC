"""Fichas de las 284 proteínas de la celda. Lee la descarga de UniProt ya guardada."""

from __future__ import annotations

import importlib.util
import json
import math
import numbers
from pathlib import Path
from urllib.parse import quote

import pandas as pd

from oncogo_crc.formato import formatear

REPO = Path(__file__).resolve().parents[2]
V37 = REPO / "results" / "v3_7"
FUENTE = REPO / "sitio_fuente" / "fuentes"
SITIO = REPO / "sitio"
ETIQUETA = "Candidata (paso 5): hipótesis para cohortes futuras, no biomarcador"
NOTA_CALPROTECTINA = (
    "S100A8 y S100A9 forman la calprotectina, de neutrófilos. "
    "También sube en la enfermedad inflamatoria intestinal. "
    "Como CN incluye esos pacientes, no se espera que esta pareja separe CP de CN."
)
CONTRASTES = (
    ("CP frente a CN", "suero_k4_crc_cn.tsv"),
    ("CP frente a HC", "suero_k4_crc_sanos.tsv"),
    ("CN frente a HC", "suero_k4_cn_sanos.tsv"),
)
COLUMNAS = (
    "log2fc", "ic95_bajo", "ic95_alto", "p_value", "padj", "p_tost",
    "clase", "efecto_minimo_detectable", "n_obs",
)


def e(texto) -> str:
    import html
    return html.escape(str(texto), quote=True)


def genes_de(valor) -> list[str]:
    if valor is None or (isinstance(valor, float) and math.isnan(valor)):
        return []
    return [gen for gen in str(valor).replace(",", ";").split(";") if gen and gen != "nan"]


def base_de(token: str) -> str:
    parte, _, sufijo = token.rpartition("-")
    if sufijo.isdigit() and parte:
        return parte
    return token


TIPOS = {
    "log2fc": "log2fc",
    "ic95_bajo": "log2fc",
    "ic95_alto": "log2fc",
    "tanaka_log2fc": "log2fc",
    "cptac_coad_log2fc": "log2fc",
    "efecto_minimo_detectable": "log2fc",
    "p_value": "p",
    "padj": "p",
    "p_tost": "p",
    "pg_q_min": "p",
    "global_pg_q": "p",
    "n_obs": "entero",
    "n_precursores": "entero",
    "ultimo_paso": "entero",
    "fraccion_proteotipicos": "auc",
    "media_abs_shap": "auc_completo",
}


def cifra(valor, fuente: str, columna: str = "") -> str:
    if valor is None or (isinstance(valor, float) and math.isnan(valor)) or pd.isna(valor):
        return "—"
    tipo = TIPOS.get(columna)
    if tipo == "entero" or (tipo is None and isinstance(valor, numbers.Integral) and not isinstance(valor, bool)):
        bruto = str(int(valor))
        visible = formatear(int(valor), "entero")
    elif tipo:
        try:
            valor = float(valor)
        except (TypeError, ValueError):
            bruto = str(valor)
            visible = bruto
        else:
            bruto = format(valor, ".16g")
            visible = formatear(valor, tipo)
    else:
        bruto = str(valor)
        visible = bruto
    return (
        f'<data value="{e(bruto)}" data-fuente="{e(fuente)}">{e(visible)}</data>'
        f'<details class="evidencia"><summary>fuente</summary><p>{e(fuente)}</p></details>'
    )


def escribir(render, env, comunes, cita) -> list[dict]:
    meta = json.loads((FUENTE / "uniprot_284.meta.json").read_text(encoding="utf-8"))
    uniprot = pd.read_csv(FUENTE / "uniprot_284.tsv", sep="\t").drop_duplicates("Entry")
    uniprot = uniprot.set_index("Entry")
    celda = pd.read_csv(V37 / "suero_k4_crc_cn.tsv", sep="\t")
    tablas = {
        nombre: pd.read_csv(V37 / nombre, sep="\t").set_index("protein_id")
        for _, nombre in CONTRASTES
    }
    traza = pd.read_csv(V37 / "suero_trazabilidad_grupos.tsv", sep="\t").set_index("protein_group")
    prio = pd.read_csv(V37 / "priorizacion.tsv", sep="\t").set_index("gen")
    tri = pd.read_csv(V37 / "triangulacion.tsv", sep="\t").set_index("gen")
    shap = pd.read_csv(V37 / "ml_explicacion.tsv", sep="\t")
    shap = shap[shap["modelo"] == "M2_proteinas_principal"].set_index("gen")
    paso5 = set(prio.index[prio["ultimo_paso"] == 5].astype(str))
    fichas = []
    for protein_id in celda["protein_id"].astype(str):
        fila_traza = traza.loc[protein_id] if protein_id in traza.index else None
        genes = genes_de(None if fila_traza is None else fila_traza["genes"])
        gen = genes[0] if len(genes) == 1 else ""
        representate = base_de(protein_id.split(";")[0])
        nombre = ""
        if representate in uniprot.index:
            bruto = str(uniprot.loc[representate].get("Protein names", ""))
            if bruto not in {"nan", "None", ""}:
                nombre = bruto
        fichas.append({
            "id": protein_id,
            "href": "/proteinas/" + quote(protein_id, safe="-;"),
            "genes": genes,
            "gen": gen,
            "nombre": nombre,
            "candidata": len(genes) == 1 and gen in paso5,
        })
    indice = [
        {"id": "resumen", "texto": "Resumen"},
        {"id": "suero", "texto": "Suero"},
        {"id": "calidad", "texto": "Calidad"},
        {"id": "tejido", "texto": "Tejido"},
        {"id": "priorizacion", "texto": "Priorización"},
        {"id": "shap", "texto": "SHAP"},
        {"id": "secuencia", "texto": "Secuencia"},
        {"id": "estructura", "texto": "Estructura"},
        {"id": "go", "texto": "GO"},
        {"id": "dominios", "texto": "Dominios"},
    ]
    vistas = _vistas()
    documento = vistas.cargar()
    for ficha in fichas:
        cuerpo = cuerpo_de(ficha, uniprot, tablas, traza, prio, tri, shap, cita, vistas, documento, comunes)
        html_pagina = render(
            env, "base.html.j2",
            titulo=ficha["id"],
            cuerpo=cuerpo,
            indice=indice,
            migas=f'Inicio / Proteínas / {ficha["id"]}',
            fuentes=fuentes_de(meta),
            **comunes,
        )
        destino = SITIO / "proteinas" / ficha["id"] / "index.html"
        destino.parent.mkdir(parents=True, exist_ok=True)
        destino.write_text(html_pagina, encoding="utf-8")
    lista = render(
        env, "base.html.j2",
        titulo="Proteínas",
        cuerpo=lista_html(fichas, meta, cita),
        indice=[],
        migas="Inicio / Proteínas",
        fuentes=fuentes_de(meta),
        **comunes,
    )
    destino = SITIO / "proteinas" / "index.html"
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_text(lista, encoding="utf-8")
    salida = SITIO / "fuentes"
    salida.mkdir(parents=True, exist_ok=True)
    for nombre in ("uniprot_284.tsv", "uniprot_284.meta.json"):
        (salida / nombre).write_bytes((FUENTE / nombre).read_bytes())
    vistas.cerrar(render, env, comunes, cita, fichas)
    return fichas


def _vistas():
    spec = importlib.util.spec_from_file_location(
        "vistas_proteina", Path(__file__).with_name("vistas_proteina.py")
    )
    modulo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(modulo)
    return modulo


def fuentes_de(meta: dict) -> list[dict]:
    return [
        {"href": "/fuentes/uniprot_284.tsv", "texto": f'UniProt {meta["x_uniprot_release"]}'},
        {"href": "https://github.com/JancoDeveloper/OncoGO-CRC/blob/595118ae09ea6006c744c86b28cb0a10d42e0442/results/v3_7/suero_k4_crc_cn.tsv", "texto": "suero_k4_crc_cn.tsv"},
    ]


def lista_html(fichas: list[dict], meta: dict, cita) -> str:
    candidatas = [ficha for ficha in fichas if ficha["candidata"]]
    repetidos = sorted({
        ficha["gen"] for ficha in candidatas
        if sum(otra["gen"] == ficha["gen"] for otra in candidatas) > 1
    })
    nota_repetidos = (
        " Genes del paso 5 en más de una ficha: " + ", ".join(repetidos) + "."
        if repetidos else ""
    )
    n_celda = json.loads((V37 / "suero_embudo.json").read_text(encoding="utf-8"))["celda_preregistrada"]["n"]
    filas = []
    for ficha in sorted(fichas, key=lambda item: item["id"]):
        marca = " candidata" if ficha["candidata"] else ""
        nombre = ficha.get("nombre") or ""
        filas.append(
            f'<li data-buscar="{e(ficha["id"])} {e(ficha["gen"])} {e(nombre)}{marca}">'
            f'<a href="{e(ficha["href"])}">{e(ficha["id"])}</a>'
            f' <span>{e(ficha["gen"])}</span>'
            + (f' <span class="nombre-proteina">{e(nombre)}</span>' if nombre else "")
            + "</li>"
        )
    return (
        "<h1>Proteínas</h1>"
        f'<p>La lista es la celda, {cifra(n_celda, "results/v3_7/suero_embudo.json#celda_preregistrada.n")} proteínas, '
        f'la misma en los tres contrastes del suero. '
        f'UniProt {e(meta["x_uniprot_release"])}, descargado el {e(meta["fecha_utc"])}. '
        f'Ninguna consulta volvió inactiva. {cita("rai2025")}</p>'
        "<p>El nombre sale de UniProt cuando la descarga lo trajo. "
        "Una ficha lleva la etiqueta de candidata cuando su único gen tiene ultimo_paso 5 en priorizacion.tsv. "
        "Esa etiqueta es una hipótesis para cohortes futuras."
        f"{e(nota_repetidos)} {cita('tanaka2024')}</p>"
        '<p><a href="/proteinas/comparar">Comparador</a> · <a href="/proteinas/candidatas">Candidatas</a></p>'
        '<p class="filtro-proteinas"><label for="filtro-proteinas">Buscar en esta lista por identificador, gen o nombre</label>'
        '<input id="filtro-proteinas" type="search" autocomplete="off">'
        '<label><input id="solo-paso-5" type="checkbox"> Solo el paso 5, hipótesis para cohortes futuras</label></p>'
        '<ul class="explorador">' + "".join(filas) + "</ul>"
    )


def cuerpo_de(ficha, uniprot, tablas, traza, prio, tri, shap, cita, vistas, documento, comunes) -> str:
    protein_id = ficha["id"]
    representate = base_de(protein_id.split(";")[0])
    entrada = uniprot.loc[representate] if representate in uniprot.index else None
    partes = [f"<h1>{e(protein_id)}</h1>"]
    if ficha["candidata"]:
        partes.append(f'<p class="etiqueta-candidata">{e(ETIQUETA)}</p>')
    if ficha["gen"] in {"S100A8", "S100A9"}:
        partes.append(f"<p>{e(NOTA_CALPROTECTINA)} {cita('rai2025')}</p>")
    partes.append(resumen(ficha, protein_id, representate, entrada, cita))
    partes.append(suero(protein_id, tablas, cita))
    partes.append(calidad(protein_id, traza))
    partes.append(tejido(ficha, tri, cita))
    partes.append(priorizacion(ficha, prio))
    partes.append(explicacion(ficha, shap))
    partes.append(vistas.secciones(ficha, documento, cita, comunes))
    return "\n".join(partes)


def resumen(ficha, protein_id, representate, entrada, cita) -> str:
    if entrada is None:
        nombre = "UniProt no devolvió esta entrada en la descarga guardada."
        funcion = ""
        revisada = "—"
        gen_up = ""
    else:
        nombre = str(entrada.get("Protein names", ""))
        funcion = str(entrada.get("Function [CC]", ""))
        if funcion == "nan":
            funcion = ""
        revisada = str(entrada.get("Reviewed", ""))
        gen_up = str(entrada.get("Gene Names (primary)", ""))
        if gen_up == "nan":
            gen_up = ""
    grupo = ""
    if ";" in protein_id:
        grupo = f"<p>El grupo de la celda es {e(protein_id)}. La ficha de UniProt es {e(representate)}.</p>"
    funcion_html = f"<p>{e(funcion)}</p>" if funcion else "<p>El campo Function vino vacío en la descarga.</p>"
    return (
        '<h2 id="resumen">Resumen</h2>'
        f"<p>{e(nombre)} Gen en la corrida: {e(', '.join(ficha['genes']) or 'sin gen')}. "
        f"Gen en UniProt: {e(gen_up or 'sin gen')}. Revisada: {e(revisada)}. {cita('rai2025')}</p>"
        + grupo
        + "<h3>Función</h3>"
        + "<p>Texto de UniProt, no un resultado de la corrida. Si más abajo aparece la etiqueta de candidata, esa etiqueta es una hipótesis.</p>"
        + funcion_html
        + f'<p>Texto del campo Function, UniProt, acceso {e(representate)}. {cita("rai2025")}</p>'
    )


def suero(protein_id, tablas, cita) -> str:
    bloques = [
        '<h2 id="suero">Suero</h2>',
        f"<p>Resultado calculado de la corrida. Tres contrastes, sin mezclar las matrices. {cita('rai2025')}</p>",
    ]
    for titulo, nombre in CONTRASTES:
        tabla = tablas[nombre]
        if protein_id not in tabla.index:
            bloques.append(f"<h3>{e(titulo)}</h3><p>No está en {e(nombre)}.</p>")
            continue
        fila = tabla.loc[protein_id]
        celdas = []
        for columna in COLUMNAS:
            fuente = f"results/v3_7/{nombre}#protein_id={protein_id}.{columna}"
            celdas.append(f"<tr><th>{e(columna)}</th><td>{cifra(fila[columna], fuente, columna)}</td></tr>")
        bloques.append(f"<h3>{e(titulo)}</h3><table><tbody>{''.join(celdas)}</tbody></table>")
    return "\n".join(bloques)


def calidad(protein_id, traza) -> str:
    if protein_id not in traza.index:
        return '<h2 id="calidad">Calidad</h2><p>No está en suero_trazabilidad_grupos.tsv.</p>'
    fila = traza.loc[protein_id]
    columnas = (
        "pg_q_min", "global_pg_q", "n_precursores", "fraccion_proteotipicos",
        "contaminante", "contaminante_piel", "queratina_epitelial",
    )
    celdas = []
    for columna in columnas:
        fuente = f"results/v3_7/suero_trazabilidad_grupos.tsv#protein_group={protein_id}.{columna}"
        celdas.append(f"<tr><th>{e(columna)}</th><td>{cifra(fila[columna], fuente, columna)}</td></tr>")
    return '<h2 id="calidad">Calidad</h2><table><tbody>' + "".join(celdas) + "</tbody></table>"


def tejido(ficha, tri, cita) -> str:
    if len(ficha["genes"]) != 1:
        return '<h2 id="tejido">Tejido</h2><p>Unión ambigua: el grupo tiene más de un gen, y no se elige una fila.</p>'
    gen = ficha["gen"]
    if gen not in tri.index:
        return f'<h2 id="tejido">Tejido</h2><p>{e(gen)} no está en triangulacion.tsv. {cita("vasaikar2019")}</p>'
    fila = tri.loc[gen]
    columnas = ("tanaka_log2fc", "cptac_coad_log2fc", "clase_h5")
    celdas = []
    for columna in columnas:
        fuente = f"results/v3_7/triangulacion.tsv#gen={gen}.{columna}"
        celdas.append(f"<tr><th>{e(columna)}</th><td>{cifra(fila[columna], fuente, columna)}</td></tr>")
    return (
        f'<h2 id="tejido">Tejido</h2><p>Resultado calculado. Una fila, el gen {e(gen)}. {cita("vasaikar2019")} {cita("tanaka2024")}</p>'
        f"<table><tbody>{''.join(celdas)}</tbody></table>"
    )


def priorizacion(ficha, prio) -> str:
    if len(ficha["genes"]) != 1:
        return '<h2 id="priorizacion">Priorización</h2><p>Unión ambigua: no se asigna un paso.</p>'
    gen = ficha["gen"]
    if gen not in prio.index:
        return f'<h2 id="priorizacion">Priorización</h2><p>{e(gen)} no está en priorizacion.tsv.</p>'
    fila = prio.loc[gen]
    columnas = ("ultimo_paso", "motivo_no_pasa_el_siguiente", "clase_h5")
    celdas = []
    for columna in columnas:
        fuente = f"results/v3_7/priorizacion.tsv#gen={gen}.{columna}"
        celdas.append(f"<tr><th>{e(columna)}</th><td>{cifra(fila[columna], fuente, columna)}</td></tr>")
    return '<h2 id="priorizacion">Priorización</h2><table><tbody>' + "".join(celdas) + "</tbody></table>"


def explicacion(ficha, shap) -> str:
    if len(ficha["genes"]) != 1:
        return '<h2 id="shap">SHAP</h2><p>Unión ambigua: no se copia una explicación.</p>'
    gen = ficha["gen"]
    if gen not in shap.index:
        return f'<h2 id="shap">SHAP</h2><p>{e(gen)} no está en la explicación de M2.</p>'
    valor = shap.loc[gen, "media_abs_shap"]
    fuente = f"results/v3_7/ml_explicacion.tsv#modelo=M2_proteinas_principal|variable={shap.loc[gen, 'variable']}.media_abs_shap"
    return (
        '<h2 id="shap">SHAP</h2><p>Media del valor absoluto en M2: '
        + cifra(valor, fuente, "media_abs_shap") + "</p>"
    )
