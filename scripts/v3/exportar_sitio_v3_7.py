#!/usr/bin/env python3
"""Sitio v3.7 a partir de results/v3_7/, config/ y docs/.

Reutiliza el aviso y el formateo de exportar_sitio_v3_6.py sin modificarlo.
No escribe en results/v3_7/.
"""

from __future__ import annotations

import hashlib
import html
import importlib.util
import json
import shutil
from pathlib import Path

import pandas as pd
import yaml

from oncogo_crc.formato import formatear

REPO = Path(__file__).resolve().parents[2]
V37 = REPO / "results" / "v3_7"
SITIO = REPO / "sitio"
SHA = "595118ae09ea6006c744c86b28cb0a10d42e0442"
GITHUB = f"https://github.com/JancoDeveloper/OncoGO-CRC/blob/{SHA}"

_spec = importlib.util.spec_from_file_location(
    "exportar_sitio_v3_6", REPO / "scripts" / "v3" / "exportar_sitio_v3_6.py"
)
_v36 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_v36)
AVISO = _v36.AVISO

PAGINAS_V36 = [
    "index.html", "datos.html", "reproducir.html", "historial.html", "aprende.html",
    "resultados.html", "404.html", "preguntas.html", "como-se-hizo.html", "glosario.html",
]


def e(texto) -> str:
    return html.escape(str(texto), quote=True)


def dato(valor, tipo: str, fuente: str) -> str:
    if isinstance(valor, float):
        bruto = format(valor, ".16g")
    else:
        bruto = str(valor)
    visible = valor if tipo == "texto" else formatear(valor, tipo)
    return f'<data value="{e(bruto)}" data-fuente="{e(fuente)}">{e(visible)}</data>'


def cargar(nombre: str):
    ruta = V37 / nombre
    if not ruta.is_file():
        return None
    return json.loads(ruta.read_text(encoding="utf-8"))


def tabla_tsv(nombre: str) -> pd.DataFrame | None:
    ruta = V37 / nombre
    if not ruta.is_file():
        return None
    return pd.read_csv(ruta, sep="\t")


def falta(nombre: str) -> str:
    return (
        f'<p class="ausente">Falta <span class="mono">{e(nombre)}</span> en results/v3_7. '
        "Esta sección no inventa el dato.</p>"
    )


def caminar(obj, ruta: str):
    actual = obj
    for parte in ruta.split("."):
        if isinstance(actual, list):
            actual = actual[int(parte)]
        else:
            actual = actual[parte]
    return actual


def archivar_v36() -> None:
    destino = SITIO / "archivo" / "v3_6"
    if (destino / "index.html").exists():
        return
    destino.mkdir(parents=True, exist_ok=True)
    for nombre in PAGINAS_V36:
        origen = SITIO / nombre
        if origen.is_file():
            shutil.move(origen, destino / nombre)
    for carpeta in ("css", "js", "figuras", "datos"):
        origen = SITIO / carpeta
        if origen.is_dir():
            shutil.move(origen, destino / carpeta)
    for pagina in destino.glob("*.html"):
        texto = pagina.read_text(encoding="utf-8")
        for viejo, nuevo in (
            ("/css/", "/archivo/v3_6/css/"),
            ("/js/", "/archivo/v3_6/js/"),
            ("/figuras/", "/archivo/v3_6/figuras/"),
            ("/datos/", "/archivo/v3_6/datos/"),
        ):
            texto = texto.replace(viejo, nuevo)
        pagina.write_text(texto, encoding="utf-8")


def sello_de(svg: str, fuentes: list[str]) -> dict:
    return {
        "svg_sha256": hashlib.sha256(svg.encode()).hexdigest(),
        "fuentes": fuentes,
        "commit_de_la_etiqueta": SHA,
        "generador": "scripts/v3/exportar_sitio_v3_7.py",
    }


def guardar_figura(fid: str, svg: str, filas: list[dict], fuentes: list[str]) -> None:
    carpeta = SITIO / "figuras"
    carpeta.mkdir(parents=True, exist_ok=True)
    (carpeta / f"{fid}.svg").write_text(svg, encoding="utf-8")
    (carpeta / f"{fid}.sello.json").write_text(
        json.dumps(sello_de(svg, fuentes), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if not filas:
        (carpeta / f"{fid}.tsv").write_text("vacio\n", encoding="utf-8")
        return
    columnas = list(filas[0])
    lineas = ["\t".join(columnas)]
    for fila in filas:
        lineas.append("\t".join(str(fila[col]) for col in columnas))
    (carpeta / f"{fid}.tsv").write_text("\n".join(lineas) + "\n", encoding="utf-8")


def ruta_en_el_commit(fuente: str) -> str:
    ruta = fuente.split("#", 1)[0]
    if ruta.startswith(("results/", "config/", "docs/")):
        return ruta
    return f"results/v3_7/{ruta}"


def figura(fid: str, titulo: str, desc: str, cuerpo: str, filas: list[dict], fuentes: list[str], pregunta: str, viewbox: str = "0 0 640 360") -> str:
    svg = (
        f'<svg viewBox="{viewbox}" role="img" class="lienzo-svg">'
        f"<title>{e(titulo)}</title><desc>{e(desc)}</desc>{cuerpo}</svg>"
    )
    guardar_figura(fid, svg, filas, fuentes)
    cabecera = "".join(f"<th>{e(col)}</th>" for col in (filas[0] if filas else {"dato": ""}))
    cuerpo_tabla = []
    for fila in filas:
        celdas = "".join(f"<td>{e(fila[col])}</td>" for col in fila)
        cuerpo_tabla.append(f"<tr>{celdas}</tr>")
    if not cuerpo_tabla:
        cuerpo_tabla.append("<tr><td>Sin filas en el fichero de origen.</td></tr>")
    enlaces = []
    for fuente in fuentes:
        ruta = ruta_en_el_commit(fuente)
        enlaces.append(f'<a href="{GITHUB}/{e(ruta)}">origen {e(ruta)}</a>')
    return f"""<figure class="grafico" id="{e(fid)}">
<h3>{e(pregunta)}</h3>
{svg}
<details><summary>Tabla de los datos del gráfico</summary>
<table class="ordenable"><thead><tr>{cabecera}</tr></thead><tbody>{''.join(cuerpo_tabla)}</tbody></table>
</details>
<p class="descargas"><a href="/figuras/{e(fid)}.svg" download>Descargar SVG</a>
<a href="/figuras/{e(fid)}.tsv" download>Descargar TSV</a>
{' '.join(enlaces)}</p>
</figure>"""


def escala(valores: list[float], inicio: float, fin: float):
    bajo, alto = min(valores), max(valores)
    if alto == bajo:
        alto = bajo + 1
    def punto(valor: float) -> float:
        return inicio + (fin - inicio) * (valor - bajo) / (alto - bajo)
    return punto, bajo, alto


def menu() -> str:
    enlaces = [
        ("/", "Inicio"), ("/lo-prometido", "Lo prometido"), ("/datos", "Datos"),
        ("/tejido", "Tejido"), ("/suero", "Suero"), ("/aprendizaje-automatico", "Aprendizaje"),
        ("/preanalitica", "Preanalítica"), ("/hipotesis", "Hipótesis"),
        ("/priorizacion", "Priorización"), ("/desviaciones", "Desviaciones"),
        ("/reproducir", "Reproducir"), ("/como-se-hizo", "Cómo se hizo"),
        ("/aprende", "Aprende"), ("/glosario", "Glosario"), ("/preguntas", "Preguntas"),
        ("/historial", "Historial"), ("/archivo/", "Archivo"),
    ]
    return "".join(f'<a href="{href}">{nombre}</a>' for href, nombre in enlaces)


def pagina(titulo: str, cuerpo: str, cert: dict | None, css: str) -> str:
    if cert is None:
        franja = '<p class="franja">Falta results/v3_7/certificado.json. No hay franja de certificado.</p>'
    else:
        franja = (
            '<p class="franja" role="status">Certificado '
            + dato(cert["estado"], "estado", "results/v3_7/certificado.json#estado")
            + ". Commit del sello "
            + dato(cert["sello"]["commit"], "texto", "results/v3_7/certificado.json#sello.commit")
            + ". Etiqueta v3.7-final.</p>"
        )
    return f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{e(titulo)} · OncoGO-CRC</title>
<link rel="stylesheet" href="/css/{e(css)}">
</head>
<body>
<a class="saltar" href="#contenido">Saltar al contenido</a>
<header class="barra">
<a href="/">OncoGO-CRC</a>
<nav>{menu()}</nav>
<button type="button" id="tema">Tema</button>
</header>
{franja}
<p class="aviso">{AVISO}</p>
<main id="contenido" class="pagina">
{cuerpo}
</main>
<p class="pie">{AVISO}</p>
<script>{JS}</script>
</body>
</html>
"""


JS = r"""
(function () {
  var reducido = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  var boton = document.getElementById("tema");
  if (boton) boton.addEventListener("click", function () {
    var actual = document.documentElement.getAttribute("data-tema");
    document.documentElement.setAttribute("data-tema", actual === "oscuro" ? "claro" : "oscuro");
  });
  document.querySelectorAll("table.ordenable").forEach(function (tabla) {
    var filas = Array.prototype.slice.call(tabla.tBodies[0].rows);
    Array.prototype.forEach.call(tabla.tHead.rows[0].cells, function (th, indice) {
      th.tabIndex = 0;
      function ordenar() {
        var abajo = th.getAttribute("aria-sort") !== "ascending";
        filas.sort(function (a, b) {
          var x = a.cells[indice].textContent.trim();
          var y = b.cells[indice].textContent.trim();
          var nx = parseFloat(x.replace(",", "."));
          var ny = parseFloat(y.replace(",", "."));
          var cmp = (isFinite(nx) && isFinite(ny)) ? nx - ny : x.localeCompare(y, "es");
          return abajo ? cmp : -cmp;
        });
        filas.forEach(function (fila) { tabla.tBodies[0].appendChild(fila); });
        th.setAttribute("aria-sort", abajo ? "ascending" : "descending");
      }
      th.addEventListener("click", ordenar);
      th.addEventListener("keydown", function (ev) { if (ev.key === "Enter") ordenar(); });
    });
  });
  document.querySelectorAll("[data-filtro]").forEach(function (botonFiltro) {
    botonFiltro.addEventListener("click", function () {
      var campo = botonFiltro.getAttribute("data-filtro");
      var valor = botonFiltro.getAttribute("data-valor");
      document.querySelectorAll(".celda-multi").forEach(function (celda) {
        var ocultar = valor !== "todas" && celda.getAttribute(campo) !== valor;
        if (ocultar) celda.setAttribute("hidden", "");
        else celda.removeAttribute("hidden");
      });
    });
  });
  if (!reducido) {
    document.querySelectorAll("data.contador").forEach(function (nodo) {
      var fin = parseFloat(nodo.getAttribute("value"));
      if (!isFinite(fin)) return;
      var t0 = performance.now();
      function paso(ahora) {
        var t = Math.min(1, (ahora - t0) / 1200);
        nodo.textContent = String(Math.round(fin * t));
        if (t < 1) requestAnimationFrame(paso);
        else nodo.textContent = nodo.getAttribute("data-visible");
      }
      requestAnimationFrame(paso);
    });
  }
})();
"""


CSS = """
:root { color-scheme: light; --papel:#f7f4ef; --tinta:#14202b; --linea:#6d6256; --marca:#0b4f55; --aviso:#6b2d12; --ficha:#fffdf8; }
:root[data-tema="oscuro"] { color-scheme: dark; --papel:#14202b; --tinta:#f4efe6; --linea:#c8bba8; --marca:#8fd0d4; --aviso:#f0b59a; --ficha:#1c2c38; }
@media (prefers-color-scheme: dark) {
  :root:not([data-tema="claro"]) { color-scheme: dark; --papel:#14202b; --tinta:#f4efe6; --linea:#c8bba8; --marca:#8fd0d4; --aviso:#f0b59a; --ficha:#1c2c38; }
}
* { box-sizing: border-box; }
body { margin:0; background:var(--papel); color:var(--tinta); font:1.05rem/1.55 Georgia, "Iowan Old Style", serif; }
a { color:var(--marca); }
a:focus, button:focus, [tabindex]:focus { outline:3px solid var(--marca); outline-offset:2px; }
.saltar { position:absolute; left:-999px; }
.saltar:focus { left:1rem; top:1rem; background:var(--ficha); padding:.4rem; z-index:2; }
.barra, .franja, .aviso, .pie { padding:.7rem 1.2rem; }
.barra { display:flex; flex-wrap:wrap; gap:.6rem 1rem; align-items:center; border-bottom:1px solid var(--linea); }
nav { display:flex; flex-wrap:wrap; gap:.4rem .8rem; }
.franja { background:#d8efe8; color:#14202b; font-family:ui-sans-serif, sans-serif; }
:root[data-tema="oscuro"] .franja { background:#1e4a45; color:#f4efe6; }
.aviso, .pie { background:#f3e1d4; color:#14202b; }
:root[data-tema="oscuro"] .aviso, :root[data-tema="oscuro"] .pie { background:#4a2c22; color:#f4efe6; }
.pagina { max-width:68rem; margin:0 auto; padding:1rem 1.2rem 3rem; }
.grafico, .tarjeta, .ficha { background:var(--ficha); border:1px solid var(--linea); padding:1rem; margin:1rem 0; }
.lienzo-svg { width:100%; height:auto; background:var(--papel); }
.rejilla { display:grid; grid-template-columns:repeat(auto-fit, minmax(16rem, 1fr)); gap:1rem; }
table { border-collapse:collapse; width:100%; font-family:ui-sans-serif, sans-serif; font-size:.92rem; }
th, td { border-bottom:1px solid var(--linea); text-align:left; padding:.35rem; vertical-align:top; }
th { cursor:pointer; }
.mono { font-family:ui-monospace, monospace; font-size:.9rem; overflow-wrap:anywhere; }
.ausente { border-left:4px solid var(--aviso); padding-left:.6rem; }
[data-tip]:focus::after, [data-tip]:hover::after { content:attr(data-tip); display:block; background:var(--tinta); color:var(--papel); padding:.3rem .5rem; }
svg [hidden] { display:none; }
.trazo { stroke-dasharray:900; stroke-dashoffset:900; animation:dibuja 1.2s ease forwards; }
@keyframes dibuja { to { stroke-dashoffset:0; } }
@media (prefers-reduced-motion: reduce) {
  .trazo { animation:none; stroke-dashoffset:0; }
}
@media print {
  .barra, button, .saltar { display:none; }
  a { text-decoration:underline; }
  body { background:#fff; color:#000; }
}
"""


def mapa() -> str:
    cuerpo = """
<rect x="40" y="40" width="150" height="64" fill="none" stroke="#0b4f55" stroke-width="2"/>
<text x="115" y="78" text-anchor="middle" fill="currentColor">Tejido</text>
<rect x="245" y="40" width="150" height="64" fill="none" stroke="#0b4f55" stroke-width="2"/>
<text x="320" y="78" text-anchor="middle" fill="currentColor">Suero</text>
<rect x="450" y="40" width="150" height="64" fill="none" stroke="#0b4f55" stroke-width="2"/>
<text x="525" y="78" text-anchor="middle" fill="currentColor">Hipótesis</text>
<rect x="245" y="220" width="150" height="64" fill="none" stroke="#6b2d12" stroke-width="2" stroke-dasharray="4 3"/>
<text x="320" y="258" text-anchor="middle" fill="currentColor">Priorización</text>
<path class="trazo" d="M190 72 H245 M395 72 H450 M320 104 V220" fill="none" stroke="#0b4f55" stroke-width="3"/>
"""
    return figura(
        "f01", "Mapa del estudio",
        "Cuatro bloques unidos por un trazo: tejido, suero, hipótesis y priorización. El trazo se dibuja una vez.",
        cuerpo, [{"bloque": "Tejido"}, {"bloque": "Suero"}, {"bloque": "Hipótesis"}, {"bloque": "Priorización"}],
        ["certificado.json"], "Cómo se recorre el estudio",
    )


def embudo(emb) -> str:
    c, cel = emb["controles"], emb["celda_preregistrada"]
    tramos = [
        ("Grupos", c["grupos"], "results/v3_7/suero_embudo.json#controles.grupos"),
        ("q global", c["grupos_global_pg_q_0_01"], "results/v3_7/suero_embudo.json#controles.grupos_global_pg_q_0_01"),
        ("Humanas", c["humanos_pg_y_global"], "results/v3_7/suero_embudo.json#controles.humanos_pg_y_global"),
        ("Celda", cel["n"], "results/v3_7/suero_embudo.json#celda_preregistrada.n"),
    ]
    alto = max(n for _, n, _ in tramos)
    partes = []
    filas = []
    for i, (nombre, n, fuente) in enumerate(tramos):
        ancho = 40 + 520 * n / alto
        x = 40 + i * 150
        partes.append(
            f'<rect x="{x}" y="{300 - 240 * n / alto:.1f}" width="36" height="{240 * n / alto:.1f}" '
            f'fill="none" stroke="#0b4f55" stroke-width="3" data-value="{n}" data-fuente="{fuente}" '
            f'tabindex="0" data-tip="{e(nombre)} {n}"><title>{e(nombre)}: {n}</title></rect>'
            f'<text x="{x}" y="330" fill="currentColor">{e(nombre)}</text>'
        )
        filas.append({"tramo": nombre, "n": n})
    return figura("f02", "Embudo del suero", "Cuatro alturas: grupos, filtro global, proteínas humanas y celda.",
                  "".join(partes), filas, ["suero_embudo.json"], "¿Cuántas proteínas llegan a la celda?")


def celda_svg(emb, crc_hc: pd.DataFrame, cn_hc: pd.DataFrame) -> str:
    cel = emb["celda_preregistrada"]
    partes_celda = [
        ("Diferente", cel["diferente"], "results/v3_7/suero_embudo.json#celda_preregistrada.diferente", "4 3"),
        ("Equivalente", cel["equivalente"], "results/v3_7/suero_embudo.json#celda_preregistrada.equivalente", "0"),
        ("Indeterminada", cel["indeterminada"], "results/v3_7/suero_embudo.json#celda_preregistrada.indeterminada", "1 2"),
    ]
    total = sum(n for _, n, _, _ in partes_celda) or 1
    x = 20
    cuerpo = []
    filas = []
    for nombre, n, fuente, raya in partes_celda:
        ancho = 400 * n / total
        cuerpo.append(
            f'<rect x="{x:.1f}" y="30" width="{ancho:.1f}" height="36" fill="none" stroke="#0b4f55" '
            f'stroke-width="3" stroke-dasharray="{raya}" data-value="{n}" data-fuente="{fuente}" tabindex="0" '
            f'data-tip="{e(nombre)} {n}"><title>{e(nombre)}: {n}</title></rect>'
        )
        x += ancho
        filas.append({"clase": nombre, "n": n, "contraste": "CP frente a CN"})
    n_cp = int((crc_hc["clase"] == "diferente").sum())
    n_cn = int((cn_hc["clase"] == "diferente").sum())
    cuerpo.append(
        f'<text x="20" y="120" fill="currentColor">CP frente a HC</text>'
        f'<rect x="180" y="96" width="{20 + n_cp * 30}" height="28" fill="none" stroke="#6b2d12" stroke-width="3" '
        f'data-value="{n_cp}" data-fuente="results/v3_7/suero_k4_crc_sanos.tsv#clase=diferente#cuenta" tabindex="0" '
        f'data-tip="CP frente a HC {n_cp}"/>'
        f'<text x="20" y="170" fill="currentColor">CN frente a HC</text>'
        f'<rect x="180" y="146" width="{20 + n_cn * 30}" height="28" fill="none" stroke="#0b4f55" stroke-width="3" '
        f'stroke-dasharray="4 3" data-value="{n_cn}" data-fuente="results/v3_7/suero_k4_cn_sanos.tsv#clase=diferente#cuenta" '
        f'tabindex="0" data-tip="CN frente a HC {n_cn}"/>'
    )
    filas.append({"clase": "diferente", "n": n_cp, "contraste": "CP frente a HC"})
    filas.append({"clase": "diferente", "n": n_cn, "contraste": "CN frente a HC"})
    return figura("f03", "Celda y contrastes con HC",
                  "La barra superior es CP frente a CN. Las dos de abajo son los contrastes con HC, con otro trazo.",
                  "".join(cuerpo), filas, ["suero_embudo.json", "suero_k4_crc_sanos.tsv", "suero_k4_cn_sanos.tsv"],
                  "¿Qué clases salen?")


def multiverso(resumen, tabla: pd.DataFrame) -> str:
    columnas = ["cuantificacion", "contaminantes", "paneles", "presencia"]
    cuerpo = []
    filas = []
    for i, fila in tabla.reset_index(drop=True).iterrows():
        col = int(i) % 24
        ren = int(i) // 24
        diferente = int(fila["diferente"])
        raya = "4 2" if diferente else "0"
        attrs = " ".join(f'data-{nombre}="{e(fila[nombre])}"' for nombre in columnas)
        fuente = f"results/v3_7/multiverso_suero.tsv#fila={i}.diferente"
        cuerpo.append(
            f'<rect class="celda-multi" x="{20 + col * 25}" y="{20 + ren * 22}" width="22" height="18" '
            f'fill="none" stroke="#0b4f55" stroke-width="2" stroke-dasharray="{raya}" {attrs} '
            f'data-value="{diferente}" data-fuente="{fuente}" tabindex="0" data-tip="diferente {diferente}"/>'
        )
        filas.append({**{k: fila[k] for k in columnas}, "diferente": diferente, "fila": int(i)})
    botones = ['<p><button type="button" data-filtro="data-contaminantes" data-valor="todas">Todas</button>']
    for valor in sorted(tabla["contaminantes"].astype(str).unique()):
        botones.append(
            f'<button type="button" data-filtro="data-contaminantes" data-valor="{e(valor)}">{e(valor)}</button>'
        )
    botones.append("</p>")
    nota = (
        "<p>Celdas de suero: "
        + dato(resumen["n_celdas_suero"], "entero", "results/v3_7/multiverso_resumen.json#n_celdas_suero")
        + ". Con alguna proteína diferente: "
        + dato(resumen["n_con_diferente"], "entero", "results/v3_7/multiverso_resumen.json#n_con_diferente")
        + ". El trazo discontinuo marca esas celdas. Sin JavaScript se ven todas.</p>"
    )
    return nota + "".join(botones) + figura(
        "f04", "Rejilla del multiverso",
        "288 celdas. El trazo discontinuo marca las que tienen alguna proteína diferente.",
        "".join(cuerpo), filas, ["multiverso_suero.tsv", "multiverso_resumen.json"],
        "¿En cuántas celdas aparece una proteína diferente?",
    )


def r1(k1, h7: pd.DataFrame, tri: pd.DataFrame) -> str:
    cruz = h7[["gen", "log2fc_cptac"]].merge(tri[["gen", "tanaka_log2fc"]], on="gen", how="inner")
    cruz["tanaka"] = pd.to_numeric(cruz["tanaka_log2fc"], errors="coerce")
    cruz["cptac"] = pd.to_numeric(cruz["log2fc_cptac"], errors="coerce")
    cruz = cruz.dropna(subset=["tanaka", "cptac"])
    sx, *_ = escala(cruz["cptac"].tolist(), 50, 600)
    sy, *_ = escala(cruz["tanaka"].tolist(), 300, 40)
    puntos = []
    filas = []
    for _, fila in cruz.iterrows():
        gen = str(fila["gen"])
        x, y = float(fila["cptac"]), float(fila["tanaka"])
        puntos.append(
            f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3" fill="none" stroke="#0b4f55" '
            f'data-value="{format(x, ".16g")}" data-fuente="results/v3_7/h7_tabla_proteinas.tsv#gen={e(gen)}.log2fc_cptac" '
            f'data-value-y="{format(y, ".16g")}" data-fuente-y="results/v3_7/triangulacion.tsv#gen={e(gen)}.tanaka_log2fc" '
            f'tabindex="0" data-tip="{e(gen)}"/>'
        )
        filas.append({"gen": gen, "log2fc_cptac": x, "tanaka_log2fc": y})
    primario = k1["contrastes"]["primario_frente_a_mucosa"]
    meta = k1["contrastes"]["metastasis_hepatica_frente_a_mucosa"]
    cob = k1["cobertura"]
    texto = (
        "<p>Concordancia de signo del primario frente a la mucosa: "
        + dato(primario["numerador"], "entero", "results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.numerador")
        + " de "
        + dato(primario["denominador"], "entero", "results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.denominador")
        + ". Cobertura de R1: "
        + dato(cob["numerador"], "entero", "results/v3_7/k1_etiquetas.json#cobertura.numerador")
        + " de "
        + dato(cob["denominador"], "entero", "results/v3_7/k1_etiquetas.json#cobertura.denominador")
        + ", frente al azar p "
        + dato(cob["p_vs_azar"], "p", "results/v3_7/k1_etiquetas.json#cobertura.p_vs_azar")
        + ". Metástasis hepática frente a mucosa: "
        + dato(meta["numerador"], "entero", "results/v3_7/k1_etiquetas.json#contrastes.metastasis_hepatica_frente_a_mucosa.numerador")
        + " de "
        + dato(meta["denominador"], "entero", "results/v3_7/k1_etiquetas.json#contrastes.metastasis_hepatica_frente_a_mucosa.denominador")
        + ". Los puntos son los genes de la tabla de H7 con log2FC numérico de Tanaka en la triangulación. "
        "Ese recuento no es el denominador del contraste.</p>"
    )
    return texto + figura(
        "f05", "log2FC de CPTAC y de Tanaka",
        "Cada punto es un gen de la tabla de H7 medido también en Tanaka. El eje horizontal es CPTAC y el vertical es Tanaka.",
        "".join(puntos), filas, ["h7_tabla_proteinas.tsv", "triangulacion.tsv", "k1_etiquetas.json"],
        "¿Van en el mismo sentido CPTAC y Tanaka?",
    )


def controles(tabla: pd.DataFrame) -> str:
    queridos = ["AQP8", "GUCA2A", "SLC26A3"]
    filas_src = tabla[tabla["control"].isin(queridos)]
    cuerpo = []
    filas = []
    for i, (_, fila) in enumerate(filas_src.iterrows()):
        fuente = f"results/v3_7/controles.tsv#control={fila['control']}.pasa"
        cuerpo.append(
            f'<rect x="40" y="{40 + i * 70}" width="560" height="48" fill="none" stroke="#6b2d12" '
            f'stroke-width="3" stroke-dasharray="5 3" data-value="{e(fila["pasa"])}" data-fuente="{fuente}" '
            f'tabindex="0" data-tip="{e(fila["control"])} {e(fila["pasa"])}"/>'
            f'<text x="56" y="{70 + i * 70}" fill="currentColor">{e(fila["control"])}: {e(fila["observado"])}</text>'
        )
        filas.append({"control": fila["control"], "pasa": fila["pasa"], "observado": fila["observado"]})
    return figura("f06", "Controles de tejido no medibles",
                  "AQP8, GUCA2A y SLC26A3 quedan marcados con trazo discontinuo porque el fichero no los da por medibles.",
                  "".join(cuerpo), filas, ["controles.tsv"], "¿Pasaron los controles de tejido?")


def go(go_json, gsea: pd.DataFrame) -> str:
    fdr = float(go_json["parametros"]["fdr"])
    cptac = gsea[gsea["cohorte"] == "cptac_tumor_frente_a_normal"]
    tanaka = gsea[gsea["cohorte"] == "tanaka_primario_frente_a_mucosa"].set_index("Term")
    bajos = cptac[pd.to_numeric(cptac["FDR q-val"], errors="coerce") <= fdr]
    puntos_datos = []
    for _, fila in bajos.iterrows():
        termino = str(fila["Term"])
        if termino not in tanaka.index:
            continue
        otra = tanaka.loc[termino]
        if isinstance(otra, pd.DataFrame):
            otra = otra.iloc[0]
        a, b = float(fila["NES"]), float(otra["NES"])
        if a == 0 or b == 0:
            continue
        puntos_datos.append((termino, a, b))
    if puntos_datos:
        sx, *_ = escala([p[1] for p in puntos_datos], 50, 600)
        sy, *_ = escala([p[2] for p in puntos_datos], 300, 40)
    cuerpo, filas = [], []
    for termino, a, b in puntos_datos:
        cuerpo.append(
            f'<circle cx="{sx(a):.1f}" cy="{sy(b):.1f}" r="5" fill="none" stroke="#0b4f55" stroke-width="2" '
            f'data-value="{format(a, ".16g")}" data-fuente="results/v3_7/gsea_tejido.tsv#Term={e(termino)}|cohorte=cptac_tumor_frente_a_normal.NES" '
            f'data-value-y="{format(b, ".16g")}" data-fuente-y="results/v3_7/gsea_tejido.tsv#Term={e(termino)}|cohorte=tanaka_primario_frente_a_mucosa.NES" '
            f'tabindex="0" data-tip="{e(termino)}"/>'
        )
        filas.append({"termino": termino, "nes_cptac": a, "nes_tanaka": b})
    conc = go_json["concordancia_nes"]
    texto = (
        "<p>Términos de CPTAC con FDR dentro del umbral que también tienen NES en Tanaka, y el mismo signo: "
        + dato(conc["numerador"], "entero", "results/v3_7/go_tejido.json#concordancia_nes.numerador")
        + " de "
        + dato(conc["denominador"], "entero", "results/v3_7/go_tejido.json#concordancia_nes.denominador")
        + ". GO no emite el veredicto de H1: "
        + dato(go_json["veredicto_h1"], "estado", "results/v3_7/go_tejido.json#veredicto_h1")
        + ".</p>"
    )
    coh = go_json["cohortes"]
    for nombre in ("cptac_tumor_frente_a_normal", "tanaka_primario_frente_a_mucosa"):
        bloque = coh[nombre]
        texto += "<p>" + e(nombre) + ": términos "
        texto += dato(bloque["n_terminos_probados"], "entero", f"results/v3_7/go_tejido.json#cohortes.{nombre}.n_terminos_probados")
        texto += ", ORA suben "
        texto += dato(bloque["n_ora_significativos_suben"], "entero", f"results/v3_7/go_tejido.json#cohortes.{nombre}.n_ora_significativos_suben")
        texto += ", ORA bajan "
        texto += dato(bloque["n_ora_significativos_bajan"], "entero", f"results/v3_7/go_tejido.json#cohortes.{nombre}.n_ora_significativos_bajan")
        texto += ", GSEA "
        texto += dato(bloque["n_gsea_fdr"], "entero", f"results/v3_7/go_tejido.json#cohortes.{nombre}.n_gsea_fdr")
        texto += ".</p>"
    return texto + figura("f07", "NES de CPTAC frente al de Tanaka",
                          "Cada punto es un término con FDR de CPTAC dentro del umbral y con NES en Tanaka.",
                          "".join(cuerpo), filas, ["gsea_tejido.tsv", "go_tejido.json"],
                          "¿El signo del NES coincide?")


def bosque(ml, k5) -> str:
    versiones = ml["versiones"]
    nombres = list(versiones)
    valores = [float(versiones[n]["auc_oof_media"]) for n in nombres]
    nulo = float(k5["nulos"]["sin_contaminantes_ni_paneles"]["percentil95_nulo"])
    sx, *_ = escala(valores + [0.5, nulo], 160, 600)
    cuerpo = [
        f'<line x1="{sx(nulo):.1f}" y1="20" x2="{sx(nulo):.1f}" y2="300" stroke="#6b2d12" stroke-dasharray="4 3" '
        f'data-value="{format(nulo, ".16g")}" data-fuente="results/v3_7/k5_auc.json#nulos.sin_contaminantes_ni_paneles.percentil95_nulo" '
        f'tabindex="0" data-tip="percentil 95 del nulo"/>'
    ]
    filas = []
    for i, nombre in enumerate(nombres):
        auc = float(versiones[nombre]["auc_oof_media"])
        y = 36 + i * 48
        fuente = f"results/v3_7/ml_suero_permutacion.json#versiones.{nombre}.auc_oof_media"
        cuerpo.append(
            f'<line x1="{sx(0.5):.1f}" y1="{y}" x2="{sx(auc):.1f}" y2="{y}" stroke="#0b4f55" stroke-width="3"/>'
            f'<circle cx="{sx(auc):.1f}" cy="{y}" r="6" fill="none" stroke="#0b4f55" stroke-width="2" '
            f'data-value="{format(auc, ".16g")}" data-fuente="{fuente}" tabindex="0" data-tip="{e(nombre)}"/>'
            f'<text x="8" y="{y + 4}" fill="currentColor">{e(nombre.split("_")[0])}</text>'
        )
        filas.append({"modelo": nombre, "auc_oof_media": auc})
    filas.append({"modelo": "K5_percentil95_nulo_sin_contaminantes_ni_paneles", "auc_oof_media": nulo})
    return figura("f08", "AUC de M1 a M6 y nulo de K5",
                  "Cada punto es el AUC fuera de pliegue. La línea discontinua es el percentil 95 del nulo de K5 en la matriz sin contaminantes ni paneles.",
                  "".join(cuerpo), filas, ["ml_suero_permutacion.json", "k5_auc.json"],
                  "¿Qué modelos separan, y contra qué nulo?")


def shap(tabla: pd.DataFrame) -> str:
    cuerpo, filas = [], []
    y = 16
    for modelo in ("M2_proteinas_principal", "M3_solo_qc"):
        bloque = tabla[tabla["modelo"] == modelo].sort_values("media_abs_shap", ascending=False)
        valores = [float(v) for v in bloque["media_abs_shap"]]
        if not valores:
            continue
        ancho = max(valores)
        for _, fila in bloque.iterrows():
            valor = float(fila["media_abs_shap"])
            gen = str(fila["gen"]) if pd.notna(fila["gen"]) else str(fila["variable"])
            fuente = f"results/v3_7/ml_explicacion.tsv#modelo={modelo}|variable={fila['variable']}.media_abs_shap"
            cuerpo.append(
                f'<rect x="150" y="{y}" width="{400 * valor / ancho:.1f}" height="12" fill="none" stroke="#0b4f55" '
                f'data-value="{format(valor, ".16g")}" data-fuente="{e(fuente)}" tabindex="0" data-tip="{e(gen)}"/>'
                f'<text x="4" y="{y + 11}" fill="currentColor" font-size="10">{e(modelo[:2])} {e(gen)}</text>'
            )
            filas.append({"modelo": modelo, "variable": fila["variable"], "gen": gen, "media_abs_shap": valor})
            y += 16
    return figura("f09", "SHAP de M2 y de M3",
                  "Barras de la media del valor absoluto de SHAP. M3 tiene las variables de calidad que trae la tabla, no veinte proteínas.",
                  "".join(cuerpo), filas, ["ml_explicacion.tsv"], "¿Qué variables empujan el modelo?",
                  viewbox="0 0 640 460")


def muestras(orden: pd.DataFrame, calidad: pd.DataFrame) -> str:
    orden = orden.copy()
    orden["id"] = orden["crudo"].str.replace(".raw", "", regex=False)
    calidad = calidad.set_index("sample_id")
    cuerpo, filas = [], []
    for i, fila in orden.sort_values("orden").reset_index(drop=True).iterrows():
        grupo = str(fila["grupo"])
        raya = {"CP": "0", "CN": "4 2", "HC": "1 2"}.get(grupo, "2 2")
        identificador = str(fila["id"])
        if identificador in calidad.index:
            n = int(calidad.loc[identificador, "n_precursores"])
            fuente = f"results/v3_7/suero_calidad_muestra.tsv#sample_id={identificador}.n_precursores"
            alto = 8 + (n / 6000) * 200
            valor = str(n)
        else:
            n = None
            fuente = f"results/v3_7/preanalitica_orden.tsv#crudo={fila['crudo']}.orden"
            alto = 8
            valor = str(int(fila["orden"]))
        x = 16 + int(i) * 10
        cuerpo.append(
            f'<rect x="{x}" y="{300 - alto:.1f}" width="8" height="{alto:.1f}" fill="none" stroke="#0b4f55" '
            f'stroke-dasharray="{raya}" data-value="{valor}" data-fuente="{e(fuente)}" tabindex="0" '
            f'data-tip="{e(identificador)} {e(grupo)}"/>'
        )
        filas.append({"crudo": fila["crudo"], "grupo": grupo, "orden": int(fila["orden"]), "n_precursores": n if n is not None else "sin_union"})
    nota = (
        "<p>El orden sale de preanalitica_orden.tsv. El número de precursores se une solo cuando el identificador coincide "
        "con suero_calidad_muestra.tsv. Las filas HC_ no se igualan a las filas Normal_.</p>"
    )
    return nota + figura("n10", "Sesenta muestras por orden de adquisición",
                         "Cada barra es una muestra. El trazo distingue el grupo. La altura es el número de precursores cuando el identificador coincide.",
                         "".join(cuerpo), filas, ["preanalitica_orden.tsv", "suero_calidad_muestra.tsv"],
                         "¿En qué orden se adquirieron las muestras?")


def tablero(k1, h5, h6, h7, h8, h9, h10, pre) -> str:
    h2 = caminar(pre, "hipotesis.h2.estado_heredado")
    filas_src = [
        ("H1", k1["h1_veredicto"], "estado", "results/v3_7/k1_etiquetas.json#h1_veredicto"),
        ("H2", h2, "texto", "config/preregistro_v3.yaml#hipotesis.h2.estado_heredado"),
        ("H3", "el preregistro v3 no tiene h3", "texto", ""),
        ("H4", caminar(pre, "hipotesis.h4.umbral"), "auc", "config/preregistro_v3.yaml#hipotesis.h4.umbral"),
        ("H5", h5["estado"], "estado", "results/v3_7/h5_especificidad_organo.json#estado"),
        ("H6", h6["estado"], "estado", "results/v3_7/h6_origen_celular.json#estado"),
        ("H7", h7["estado"], "estado", "results/v3_7/h7_visibilidad_en_sangre.json#estado"),
        ("H8", h8["estado"], "estado", "results/v3_7/h8_carga_tumoral.json#estado"),
        ("H9", h9["estado"], "estado", "results/v3_7/h9_wald.json#estado"),
        ("H10", h10["estado"], "estado", "results/v3_7/h10_especificidad.json#estado"),
    ]
    cuerpo, filas = [], []
    for i, (nombre, valor, tipo, fuente) in enumerate(filas_src):
        col, ren = i % 5, i // 5
        visible = valor if tipo == "texto" and not isinstance(valor, (int, float)) else (
            formatear(valor, tipo) if tipo != "texto" else str(valor)
        )
        bruto = format(valor, ".16g") if isinstance(valor, float) else str(valor)
        attrs = f' data-value="{e(bruto)}" data-fuente="{e(fuente)}"' if fuente else ""
        cuerpo.append(
            f'<rect x="{20 + col * 122}" y="{40 + ren * 120}" width="112" height="90" fill="none" stroke="#0b4f55" '
            f'stroke-width="2"{attrs} tabindex="0" data-tip="{e(nombre)} {e(visible)}"/>'
            f'<text x="{28 + col * 122}" y="{78 + ren * 120}" fill="currentColor">{e(nombre)}</text>'
            f'<text x="{28 + col * 122}" y="{100 + ren * 120}" fill="currentColor" font-size="11">{e(visible)[:18]}</text>'
        )
        filas.append({"hipotesis": nombre, "valor": bruto, "fuente": fuente})
    return figura("n11", "Tablero de H1 a H10",
                  "Cada celda copia el estado del fichero de esa hipótesis. H3 no está en el preregistro v3. H4 muestra el umbral, no un veredicto.",
                  "".join(cuerpo), filas,
                  ["k1_etiquetas.json", "h5_especificidad_organo.json", "h6_origen_celular.json",
                   "h7_visibilidad_en_sangre.json", "h8_carga_tumoral.json", "h9_wald.json", "h10_especificidad.json"],
                  "¿Cuál es el veredicto de cada hipótesis?")


def embudo_priori(pri) -> str:
    pasos = pri["pasos"]
    alto = max(p["n"] for p in pasos) or 1
    cuerpo, filas = [], []
    for i, paso in enumerate(pasos):
        n = int(paso["n"])
        fuente = f"results/v3_7/priorizacion_embudo.json#pasos.{i}.n"
        cuerpo.append(
            f'<rect x="{20 + i * 88}" y="{320 - 260 * n / alto:.1f}" width="36" height="{260 * n / alto:.1f}" '
            f'fill="none" stroke="#0b4f55" stroke-width="3" data-value="{n}" data-fuente="{fuente}" tabindex="0" '
            f'data-tip="paso {paso["paso"]} {n}"/>'
            f'<text x="{20 + i * 88}" y="340" fill="currentColor">{paso["paso"]}</text>'
        )
        filas.append({"paso": paso["paso"], "nombre": paso["nombre"], "n": n})
    return figura("n12", "Embudo de la priorización", "Siete pasos, en el orden del fichero.",
                  "".join(cuerpo), filas, ["priorizacion_embudo.json"], "¿Cuántas proteínas quedan en cada paso?")


def calor(pri_tsv: pd.DataFrame, tri: pd.DataFrame) -> str:
    genes = pri_tsv.loc[pri_tsv["ultimo_paso"] == 5, "gen"].astype(str).tolist()
    sub = tri[tri["gen"].isin(genes)].set_index("gen")
    columnas = ["tanaka_log2fc", "cptac_coad_log2fc"]
    cuerpo, filas = [], []
    numeros = []
    for gen in genes:
        if gen not in sub.index:
            continue
        for col in columnas:
            valor = pd.to_numeric(pd.Series([sub.loc[gen, col]]), errors="coerce").iloc[0]
            if pd.notna(valor):
                numeros.append(float(valor))
    limite = max(abs(n) for n in numeros) if numeros else 1
    for i, gen in enumerate(genes):
        if gen not in sub.index:
            continue
        for j, col in enumerate(columnas):
            bruto = sub.loc[gen, col]
            valor = pd.to_numeric(pd.Series([bruto]), errors="coerce").iloc[0]
            if pd.isna(valor):
                continue
            tono = int(255 * (1 - min(1, abs(float(valor)) / limite)))
            raya = "0" if float(valor) >= 0 else "3 2"
            fuente = f"results/v3_7/triangulacion.tsv#gen={gen}.{col}"
            cuerpo.append(
                f'<rect x="{140 + j * 80}" y="{16 + i * 16}" width="70" height="14" fill="none" stroke="rgb({tono},{tono},{tono})" '
                f'stroke-width="4" stroke-dasharray="{raya}" data-value="{format(float(valor), ".16g")}" data-fuente="{e(fuente)}" '
                f'tabindex="0" data-tip="{e(gen)} {col}"/>'
            )
            filas.append({"gen": gen, "columna": col, "valor": float(valor)})
        cuerpo.append(f'<text x="4" y="{28 + i * 16}" fill="currentColor" font-size="10">{e(gen)}</text>')
    return figura("n13", "Triangulación de las 20 del paso 5",
                  "Dos columnas numéricas, Tanaka y CPTAC. El trazo discontinuo es un log2FC negativo. El gris más oscuro es la magnitud mayor.",
                  "".join(cuerpo), filas, ["triangulacion.tsv", "priorizacion.tsv"],
                  "¿Qué signo tienen en tejido las 20 del paso 5?")


def certificado_svg(cert) -> str:
    cuerpo, filas = [], []
    for i, casilla in enumerate(cert["casillas"]):
        col, ren = i % 4, i // 4
        raya = "0" if casilla["estado"] == "pasa" else "4 2"
        fuente = f"results/v3_7/certificado.json#casillas.{i}.estado"
        cuerpo.append(
            f'<rect x="{16 + col * 156}" y="{16 + ren * 52}" width="148" height="44" fill="none" stroke="#0b4f55" '
            f'stroke-dasharray="{raya}" data-value="{e(casilla["estado"])}" data-fuente="{fuente}" tabindex="0" '
            f'data-tip="{e(casilla["nombre"])} {e(casilla["estado"])}"/>'
            f'<text x="{24 + col * 156}" y="{42 + ren * 52}" fill="currentColor" font-size="11">{e(casilla["nombre"])}</text>'
        )
        filas.append({"casilla": casilla["nombre"], "estado": casilla["estado"]})
    return figura("n14", "Casillas del certificado",
                  "Cada casilla copia su estado. El trazo discontinuo no es pasa.",
                  "".join(cuerpo), filas, ["certificado.json"], "¿Qué casillas pasan?")


def n_b_desde_doc() -> int:
    for linea in (REPO / "docs" / "cambio_esm2_b_d.md").read_text(encoding="utf-8").splitlines():
        if linea.startswith("| B |"):
            celda = [parte.strip() for parte in linea.strip("|").split("|")][-1]
            return int(celda.split()[0])
    raise RuntimeError("docs/cambio_esm2_b_d.md no tiene la fila B")


def enlaces_publicables() -> list[dict]:
    ruta = REPO / "docs" / "web_enlaces.tsv"
    if not ruta.is_file():
        return []
    lineas = ruta.read_text(encoding="utf-8").splitlines()
    cabecera = lineas[0].split("\t")
    salen = []
    for linea in lineas[1:]:
        fila = dict(zip(cabecera, linea.split("\t")))
        if fila.get("estado_http") == "200":
            salen.append(fila)
    return salen


def glosario() -> str:
    terminos = [
        ("g1", "AUC", "Área bajo la curva de una clasificación. Resume una separación, no decide un tratamiento."),
        ("g2", "FDR", "Tasa esperada de falsos hallazgos cuando se miran muchas proteínas a la vez."),
        ("g3", "BH", "Corrección de Benjamini-Hochberg. Controla el FDR dentro del grupo que se declara."),
        ("g4", "CN", "Pacientes sintomáticos con otra enfermedad intestinal. No son el grupo de comparación sano."),
        ("g5", "CP", "Grupo de cáncer colorrectal en el suero de esta tesis."),
        ("g6", "HC", "Grupo de personas sanas del suero. Se adquirió aparte de CP y de CN."),
        ("g7", "Celda", "La matriz principal del suero, después de los filtros escritos antes de mirar."),
        ("g8", "Diferente", "La proteína pasa el contraste con el efecto mínimo y el FDR del preregistro."),
        ("g9", "Equivalente", "La proteína queda dentro del margen de equivalencia del contraste."),
        ("g10", "Indeterminada", "Ni diferente ni equivalente. El estudio no la resuelve."),
        ("g11", "log2FC", "Logaritmo en base dos del cambio entre dos grupos."),
        ("g12", "Treat", "Umbral de magnitud. Un cambio pequeño puede ser significativo y no pasar treat."),
        ("g13", "GO", "Gene Ontology. Vocabulario de funciones, usado aquí sin propagar y sin IEA."),
        ("g14", "ORA", "Análisis de sobrerrepresentación. Pregunta si un término acumula genes de una lista."),
        ("g15", "GSEA", "Enriquecimiento sobre una lista ordenada, no sobre un corte."),
        ("g16", "NES", "Puntuación de enriquecimiento normalizada de GSEA."),
        ("g17", "SHAP", "Reparto aditivo del valor de un modelo ya ajustado. No es una proteína nueva."),
        ("g18", "DIA", "Adquisición independiente de datos. El espectrómetro no elige un precursor por su intensidad."),
        ("g19", "Precursor", "El ion del péptido antes de fragmentarlo. Su recuento habla de la profundidad de la corrida."),
        ("g20", "Sello", "Bloque con el commit, la fecha y los hash de las entradas de una salida."),
        ("g21", "Certificado", "Lista de casillas. Publicable significa que las casillas exigidas pasan o están declaradas."),
        ("g22", "Multiverso", "Las mismas proteínas leídas con otras reglas de cuantificación, contaminantes y presencia."),
        ("g23", "Preregistro", "Lo que se escribió antes de calcular: umbrales, semillas y qué no se iba a ejecutar."),
        ("g24", "Enmienda", "Un cambio posterior, con fecha y con el valor anterior. No reescribe el preregistro en silencio."),
        ("g25", "Semilla", "El entero que fija el azar de permutaciones y pliegues."),
        ("g26", "Determinismo", "Dos procesos, con el mismo código y las mismas entradas, dejan los mismos números."),
        ("g27", "Contaminante", "Proteína que no viene de la muestra. La biblioteca se declara y se puede dejar fuera."),
        ("g28", "Cohorte", "Un conjunto de muestras medido con un diseño propio. Aquí no se fusionan cohortes."),
        ("g29", "Tejido", "Las comparaciones de tumor, mucosa y metástasis. No es el suero."),
        ("g30", "Suero", "La matriz de PXD055218. CP, CN y HC no se interpretan como la misma pregunta."),
        ("g31", "Mucosa", "El tejido no tumoral pareado en el diseño de Tanaka."),
        ("g32", "Metástasis", "La muestra hepática del mismo diseño. No sustituye al tumor primario."),
        ("g33", "Colonoscopia", "La exploración que define el grupo CN cuando no muestra este cáncer."),
        ("g34", "Fase aguda", "Respuesta del organismo que sube proteínas abundantes. No identifica el tumor."),
        ("g35", "Imputación", "Rellenar un hueco. Esta corrida no imputa la matriz principal."),
        ("g36", "ESM-2", "Modelo de lenguaje de proteínas. En esta tesis el vector se lee de un fichero ya calculado."),
        ("g37", "TMbed", "Predicción de hélices transmembrana y péptido señal. También se lee de un fichero."),
        ("g38", "Ventana", "Un trozo consecutivo de una secuencia larga. El vector final promedia las ventanas."),
        ("g39", "Veredicto", "La palabra que la regla escrita asigna a una hipótesis. El secundario de H7 no entra en el suyo."),
        ("g40", "Desviación", "Una pieza anunciada que no se ejecutó, y que se declara en vez de presentarse como resultado."),
        ("g41", "Réplica", "La segunda cohorte con la que se contrasta un hallazgo de la primera."),
        ("g42", "Cobertura", "Qué fracción de las proteínas de la primera cohorte llega a medirse en la segunda."),
    ]
    items = []
    for ident, nombre, texto in terminos:
        items.append(f'<dt id="{ident}">{e(nombre)}</dt><dd>{e(texto)}</dd>')
    return "<dl>" + "".join(items) + "</dl>"


def preguntas() -> str:
    return """
<details><summary>¿Esto sirve para un paciente?</summary><p>No. El aviso de cada página es la respuesta.</p></details>
<details><summary>¿CN son personas sanas?</summary><p>No. CN son pacientes sintomáticos con otra enfermedad intestinal.</p></details>
<details><summary>¿Por qué el resultado principal es negativo?</summary><p>En la celda preregistrada, CP frente a CN no deja proteínas diferentes.</p></details>
<details><summary>¿Qué es el archivo v1?</summary><p>Páginas anteriores. No forman parte del certificado de esta tesis.</p></details>
"""


def cuerpo_index(emb) -> str:
    cel = emb["celda_preregistrada"]
    return (
        "<h1>OncoGO-CRC v3.7</h1>"
        "<p>Lo primero: en la celda, CP frente a CN deja "
        + dato(cel["diferente"], "entero", "results/v3_7/suero_embudo.json#celda_preregistrada.diferente")
        + " proteínas diferentes, "
        + dato(cel["equivalente"], "entero", "results/v3_7/suero_embudo.json#celda_preregistrada.equivalente")
        + " equivalentes y "
        + dato(cel["indeterminada"], "entero", "results/v3_7/suero_embudo.json#celda_preregistrada.indeterminada")
        + " indeterminadas, de "
        + dato(cel["n"], "entero", "results/v3_7/suero_embudo.json#celda_preregistrada.n")
        + ". CN son pacientes sintomáticos con otra enfermedad intestinal.</p>"
        "<h2>En un minuto</h2>"
        "<p>El tejido se replica en parte. El suero, frente a CN, no separa proteínas de la celda. "
        "Las hipótesis conservan el veredicto de su fichero. La priorización llega a veinte proteínas en el paso cinco y a ninguna en el seis.</p>"
        + mapa()
    )


def cuerpo_prometido(corr: pd.DataFrame, figuras_por_punto: dict[str, str]) -> str:
    tarjetas = []
    for _, fila in corr.iterrows():
        punto = str(fila["punto"])
        grafico = figuras_por_punto.get(punto, "")
        fichero = str(fila["fichero"]) if pd.notna(fila["fichero"]) else ""
        fuente = f"config/correspondencia_propuesta.tsv#punto={punto}.estado"
        tarjetas.append(
            f'<article class="tarjeta tarjeta-prometida"><h2>{e(fila["prometido"])}</h2>'
            f'<p>Estado {dato(fila["estado"], "texto", fuente)}.</p>'
            f'<p>{e(fila["explicacion"])}</p>'
            + (f'<p class="mono">{e(fichero)}</p>' if fichero else "<p>Esta fila no tiene fichero de resultado.</p>")
            + grafico + "</article>"
        )
    return "<h1>Lo prometido</h1>" + "".join(tarjetas)


def cuerpo_datos(manifest: str) -> str:
    filas = []
    for linea in manifest.splitlines():
        if not linea.strip():
            continue
        huella, ruta = linea.split(maxsplit=1)
        filas.append(
            "<tr><td class=\"mono\">"
            + dato(huella, "texto", f"results/v3_7/manifest_entradas.sha256#{ruta}")
            + f"</td><td class=\"mono\">{e(ruta)}</td></tr>"
        )
    return (
        "<h1>Datos</h1>"
        "<p>CP es el grupo de cáncer. CN son pacientes sintomáticos con otra enfermedad intestinal. "
        "HC son personas sanas, medidas en otras corridas del mismo instrumento.</p>"
        "<h2>Entradas del manifiesto</h2><table class=\"ordenable\"><thead><tr><th>sha256</th><th>ruta</th></tr></thead><tbody>"
        + "".join(filas) + "</tbody></table>"
        "<h2>Cohortes</h2><ul>"
        "<li>Tejido de descubrimiento: PDC000116, a través de results/certified/cptac_unshared_treat.tsv.</li>"
        "<li>Tejido de réplica: PXD034575.</li>"
        "<li>Suero: PXD055218.</li></ul>"
    )


def cuerpo_desviaciones(h7, pre) -> str:
    sec = h7["modelo_secundario"]
    n_b = n_b_desde_doc()
    ausencias = caminar(pre, "enmienda_v3_7_cierre.no_se_ejecutan_en_esta_tesis")
    items = "".join(f"<li><strong>{e(clave)}</strong>: {e(texto)}</li>" for clave, texto in ausencias.items())
    return (
        "<h1>Desviaciones</h1><p>Lo que no se ejecutó no se presenta como resultado.</p><ul>"
        + items + "</ul><h2>ESM-2</h2><p>En la corrida B el documento de cambio cuenta "
        + dato(n_b, "entero", "docs/cambio_esm2_b_d.md#tabla.B.proteinas")
        + " proteínas con vector, de "
        + dato(sec["n_con_secuencia"], "entero", "results/v3_7/h7_visibilidad_en_sangre.json#modelo_secundario.n_con_secuencia")
        + " con secuencia, porque el cálculo en MPS omitía errores. La corrida certificada lee el fichero de CPU: "
        + dato(sec["n_con_esm2"], "entero", "results/v3_7/h7_visibilidad_en_sangre.json#modelo_secundario.n_con_esm2")
        + " con vector y "
        + dato(sec["n_omitidas_esm2"], "entero", "results/v3_7/h7_visibilidad_en_sangre.json#modelo_secundario.n_omitidas_esm2")
        + " omitidas, todas sin secuencia. El AUC secundario es "
        + dato(sec["auc_secundario"], "auc_completo", "results/v3_7/h7_visibilidad_en_sangre.json#modelo_secundario.auc_secundario")
        + " y el incremento sobre el principal es "
        + dato(sec["incremento"], "auc_completo", "results/v3_7/h7_visibilidad_en_sangre.json#modelo_secundario.incremento")
        + ". Ese secundario no entra en el veredicto, que sigue en "
        + dato(h7["estado"], "estado", "results/v3_7/h7_visibilidad_en_sangre.json#estado")
        + ".</p>"
    )


def cuerpo_reproducir(cert) -> str:
    return (
        "<h1>Reproducir</h1>"
        "<pre class=\"mono\">make v3-7</pre>"
        "<p>Esa orden repite la etapa. No recalcula ESM-2: lo lee de external/modelos/h7_esm2_embeddings.tsv.</p>"
        "<h2>Determinismo</h2><ul>"
        "<li>A quedó archivada en el volumen. docs/determinismo_dos_procesos.md la compara con la B que entonces ocupaba results/v3_7, no con esta E.</li>"
        "<li>B y C, con el mismo sello, no coinciden en el secundario de H7. El detalle está en docs/determinismo_b_c.md.</li>"
        "<li>La primera pareja D y E se descartó: secuencia.json quedó en no_evaluable. Está en el volumen y en docs/determinismo_d_e.md.</li>"
        "<li>D y E, ya con el estado calculado, son idénticas en los números. docs/determinismo_d_e.md.</li>"
        "<li>La comprobación estricta, que no da por iguales dos sha256 distintos, no encuentra otra diferencia. docs/determinismo_d_e_estricto.md.</li>"
        "</ul>"
        + certificado_svg(cert)
        + "<p>Los tests de la corrida están en results/v3_7/tests.txt.</p>"
    )


def cuerpo_aprende(enlaces: list[dict]) -> str:
    if not enlaces:
        return "<h1>Aprende</h1>" + falta("docs/web_enlaces.tsv")
    items = []
    for fila in enlaces:
        items.append(
            f'<li><a href="{e(fila["url"])}">{e(fila["id"])}</a>. {e(fila.get("titulo", ""))}</li>'
        )
    return (
        "<h1>Aprende</h1>"
        "<p>Estas lecciones apuntan fuera. Solo se enlaza lo que respondió 200 en docs/web_enlaces.tsv.</p>"
        "<ul>" + "".join(items) + "</ul>"
        "<p>Una cohorte no es una persona. Un AUC no es una decisión. Un término GO no es un mecanismo demostrado.</p>"
    )


def mini(fid: str, punto: str, titulo: str, valor, tipo: str, fuente: str) -> str:
    bruto = format(valor, ".16g") if isinstance(valor, float) else str(valor)
    visible = valor if tipo == "texto" else formatear(valor, tipo)
    cuerpo = (
        f'<text x="40" y="180" fill="currentColor" font-size="28">{e(visible)}</text>'
        f'<rect x="40" y="80" width="200" height="40" fill="none" stroke="#0b4f55" stroke-width="3" '
        f'data-value="{e(bruto)}" data-fuente="{e(fuente)}" tabindex="0" data-tip="{e(punto)}"/>'
    )
    return figura(fid, titulo, titulo, cuerpo, [{"punto": punto, "valor": bruto}], [fuente], titulo)


def escribir_css() -> str:
    contenido = CSS.encode()
    nombre = "v37." + hashlib.sha256(contenido).hexdigest()[:10] + ".css"
    carpeta = SITIO / "css"
    carpeta.mkdir(parents=True, exist_ok=True)
    for viejo in carpeta.glob("v37*.css"):
        if viejo.name != nombre:
            viejo.unlink()
    (carpeta / nombre).write_bytes(contenido)
    return nombre


def main() -> None:
    archivar_v36()
    css = escribir_css()
    cert = cargar("certificado.json")
    emb = cargar("suero_embudo.json")
    if cert is None or emb is None:
        raise SystemExit("Faltan certificado.json o suero_embudo.json en results/v3_7")
    k1 = cargar("k1_etiquetas.json")
    go_json = cargar("go_tejido.json")
    ml = cargar("ml_suero_permutacion.json")
    k5 = cargar("k5_auc.json")
    h5, h6 = cargar("h5_especificidad_organo.json"), cargar("h6_origen_celular.json")
    h7, h8 = cargar("h7_visibilidad_en_sangre.json"), cargar("h8_carga_tumoral.json")
    h9, h10 = cargar("h9_wald.json"), cargar("h10_especificidad.json")
    pre_json = cargar("preanalitica.json")
    pri = cargar("priorizacion_embudo.json")
    multi = cargar("multiverso_resumen.json")
    corr = pd.read_csv(REPO / "config" / "correspondencia_propuesta.tsv", sep="\t")
    pre_yaml = yaml.safe_load((REPO / "config" / "preregistro_v3.yaml").read_text(encoding="utf-8"))
    crc_hc = tabla_tsv("suero_k4_crc_sanos.tsv")
    cn_hc = tabla_tsv("suero_k4_cn_sanos.tsv")
    crc_cn = tabla_tsv("suero_k4_crc_cn.tsv")
    multi_tsv = tabla_tsv("multiverso_suero.tsv")
    h7_tsv = tabla_tsv("h7_tabla_proteinas.tsv")
    tri = tabla_tsv("triangulacion.tsv")
    gsea = tabla_tsv("gsea_tejido.tsv")
    cont = tabla_tsv("controles.tsv")
    expl = tabla_tsv("ml_explicacion.tsv")
    orden = tabla_tsv("preanalitica_orden.tsv")
    calidad = tabla_tsv("suero_calidad_muestra.tsv")
    pri_tsv = tabla_tsv("priorizacion.tsv")
    aud = tabla_tsv("auditoria_rai2025.tsv")
    manifest = (V37 / "manifest_entradas.sha256").read_text(encoding="utf-8")

    f_emb = embudo(emb)
    f_cel = celda_svg(emb, crc_hc, cn_hc)
    f_mul = multiverso(multi, multi_tsv)
    f_r1 = r1(k1, h7_tsv, tri)
    f_con = controles(cont)
    f_go = go(go_json, gsea)
    f_auc = bosque(ml, k5)
    f_sha = shap(expl)
    f_mue = muestras(orden, calidad)
    f_tab = tablero(k1, h5, h6, h7, h8, h9, h10, pre_yaml)
    f_pri = embudo_priori(pri)
    f_cal = calor(pri_tsv, tri)
    sec = h7["modelo_secundario"]
    figuras_punto = {
        "integrar_proteomica": mini("p01", "celda", "Proteínas de la celda", emb["celda_preregistrada"]["n"], "entero", "results/v3_7/suero_embudo.json#celda_preregistrada.n"),
        "proteinas_diferenciales": mini("p02", "diferente", "Diferentes frente a CN", emb["celda_preregistrada"]["diferente"], "entero", "results/v3_7/suero_embudo.json#celda_preregistrada.diferente"),
        "reproducibilidad": mini("p03", "cobertura", "Cobertura de R1", k1["cobertura"]["numerador"], "entero", "results/v3_7/k1_etiquetas.json#cobertura.numerador"),
        "priorizacion": mini("p04", "paso5", "Paso 5", pri["pasos"][4]["n"], "entero", "results/v3_7/priorizacion_embudo.json#pasos.4.n"),
        "go_ora_gsea": mini("p05", "nes", "Concordancia de NES", go_json["concordancia_nes"]["numerador"], "entero", "results/v3_7/go_tejido.json#concordancia_nes.numerador"),
        "redes": mini("p06", "redes", "Redes", "desviacion_declarada", "texto", "config/correspondencia_propuesta.tsv#punto=redes.estado"),
        "ml_interpretable": mini("p07", "m3", "AUC de M3", ml["versiones"]["M3_solo_qc"]["auc_oof_media"], "auc", "results/v3_7/ml_suero_permutacion.json#versiones.M3_solo_qc.auc_oof_media"),
        "score": mini("p08", "paso6", "Paso 6", pri["pasos"][5]["n"], "entero", "results/v3_7/priorizacion_embudo.json#pasos.5.n"),
        "ampliar_cohortes": mini("p09", "cohortes", "Cohortes nuevas", "no_cumplido", "texto", "config/correspondencia_propuesta.tsv#punto=ampliar_cohortes.estado"),
        "robustecer_validacion": mini("p10", "certificado", "Certificado", cert["estado"], "estado", "results/v3_7/certificado.json#estado"),
    }
    aud_html = "<h2>Auditoría de BAIAP3 y CEACAM5</h2><table class=\"ordenable\"><thead><tr><th>gen</th><th>paso</th><th>motivo</th></tr></thead><tbody>"
    for _, fila in aud.iterrows():
        aud_html += f"<tr><td>{e(fila['gen'])}</td><td>{e(fila['paso'])}</td><td>{e(fila['motivo'])}</td></tr>"
    aud_html += "</tbody></table>"
    cp_genes = ", ".join(crc_hc.loc[crc_hc["clase"] == "diferente", "protein_id"].astype(str))
    cn_genes = ", ".join(cn_hc.loc[cn_hc["clase"] == "diferente", "protein_id"].astype(str))
    tejido = (
        "<h1>Tejido</h1>" + f_r1 + f_con + f_go
        + "<h2>H5</h2><p>" + dato(h5["estado"], "estado", "results/v3_7/h5_especificidad_organo.json#estado") + ".</p>"
        + "<h2>H6</h2><p>" + dato(h6["estado"], "estado", "results/v3_7/h6_origen_celular.json#estado")
        + ". " + e(h6.get("motivo_degenerada") or h6.get("nota") or "") + "</p>"
    )
    suero = (
        "<h1>Suero</h1><p>CN son pacientes sintomáticos con otra enfermedad intestinal. "
        "La celda preregistrada compara CP con CN.</p>"
        + f_emb + f_cel + f_mul + f_cal + aud_html
        + "<p>Identificadores con clase diferente, CP frente a HC: " + e(cp_genes)
        + ". CN frente a HC: " + e(cn_genes) + ".</p>"
        + "<p>Frente a CN, la clase diferente suma "
        + dato(int((crc_cn["clase"] == "diferente").sum()), "entero", "results/v3_7/suero_k4_crc_cn.tsv#clase=diferente#cuenta")
        + ".</p>"
    )
    aprendizaje = (
        "<h1>Aprendizaje automático</h1>"
        "<p>M3 usa solo la calidad y ya separa los grupos. El AUC fuera de pliegue es "
        + dato(ml["versiones"]["M3_solo_qc"]["auc_oof_media"], "auc", "results/v3_7/ml_suero_permutacion.json#versiones.M3_solo_qc.auc_oof_media")
        + ".</p>" + f_auc + f_sha
        + "<p>Comprobación lineal de SHAP en M2: "
        + dato(cargar("ml_explicacion.json")["comprobacion_shap"]["M2_proteinas_principal"], "auc_completo", "results/v3_7/ml_explicacion.json#comprobacion_shap.M2_proteinas_principal")
        + ".</p>"
    )
    preanalitica = (
        "<h1>Preanalítica</h1><p>Estado "
        + dato(pre_json["estado"], "estado", "results/v3_7/preanalitica.json#estado")
        + ". Crudos válidos "
        + dato(pre_json["n_crudos_validos"], "entero", "results/v3_7/preanalitica.json#n_crudos_validos")
        + ". Bloques perfectos "
        + dato(str(pre_json["bloques_perfectos"]).lower(), "texto", "results/v3_7/preanalitica.json#bloques_perfectos")
        + ".</p>" + f_mue
    )
    hipotesis = "<h1>Hipótesis</h1>" + f_tab + "<p>H7 secundario, que no entra en el veredicto: " + dato(sec["auc_secundario"], "auc_completo", "results/v3_7/h7_visibilidad_en_sangre.json#modelo_secundario.auc_secundario") + ".</p>"
    priorizacion = (
        "<h1>Priorización</h1>" + f_pri
        + "<p>Nadie pasa el paso 6: "
        + dato(str(pri["nadie_pasa_el_paso_6"]).lower(), "texto", "results/v3_7/priorizacion_embudo.json#nadie_pasa_el_paso_6")
        + ".</p><h2>Las del paso 5</h2><ul>"
        + "".join(f"<li>{e(gen)}</li>" for gen in pri_tsv.loc[pri_tsv["ultimo_paso"] == 5, "gen"].astype(str))
        + "</ul>"
    )
    enlaces = enlaces_publicables()
    paginas = {
        "index.html": ("Inicio", cuerpo_index(emb)),
        "lo-prometido.html": ("Lo prometido", cuerpo_prometido(corr, figuras_punto)),
        "datos.html": ("Datos", cuerpo_datos(manifest)),
        "tejido.html": ("Tejido", tejido),
        "suero.html": ("Suero", suero),
        "aprendizaje-automatico.html": ("Aprendizaje automático", aprendizaje),
        "preanalitica.html": ("Preanalítica", preanalitica),
        "hipotesis.html": ("Hipótesis", hipotesis),
        "priorizacion.html": ("Priorización", priorizacion),
        "desviaciones.html": ("Desviaciones", cuerpo_desviaciones(h7, pre_yaml)),
        "reproducir.html": ("Reproducir", cuerpo_reproducir(cert)),
        "como-se-hizo.html": ("Cómo se hizo", "<h1>Cómo se hizo</h1><p>La corrida lee el sello, el preregistro y las entradas declaradas. ESM-2 no se calcula dentro de la corrida. La semilla del sello es " + dato(cert["sello"]["semilla"], "entero", "results/v3_7/certificado.json#sello.semilla") + ".</p>"),
        "aprende.html": ("Aprende", cuerpo_aprende(enlaces)),
        "glosario.html": ("Glosario", "<h1>Glosario</h1>" + glosario()),
        "preguntas.html": ("Preguntas", "<h1>Preguntas</h1>" + preguntas()),
        "historial.html": ("Historial", "<h1>Historial</h1><ol><li>v1: las páginas de docs/web, ahora en el archivo.</li><li>v3.3 a v3.6: el sitio que está en /archivo/v3_6/.</li><li>v3.7: esta corrida, con certificado publicable y etiqueta v3.7-final.</li></ol>"),
        "404.html": ("No está", "<h1>Esta página no está</h1><p><a href=\"/\">Volver al inicio</a>.</p>"),
    }
    for nombre, (titulo, cuerpo) in paginas.items():
        (SITIO / nombre).write_text(pagina(titulo, cuerpo, cert, css), encoding="utf-8")
    archivo = pagina(
        "Archivo",
        "<h1>Archivo</h1><ul><li><a href=\"/archivo/v3_6/\">Sitio anterior, v3.6</a></li>"
        "<li><a href=\"/archivo/v3_3\">v3.3</a></li><li><a href=\"/archivo/v3_4\">v3.4</a></li></ul>",
        cert,
        css,
    )
    (SITIO / "archivo" / "index.html").write_text(archivo, encoding="utf-8")
    (SITIO / "archivo" / "v3_6" / "index.html")  # existe por el traslado
    print("sitio v3.7", len(paginas), "paginas", css)


if __name__ == "__main__":
    main()
