"""Sitio estático de la corrida v3.5. Cada cifra sale de un JSON o un TSV."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[2]
SITIO = RAIZ / "sitio"
AVISO = "Ningún resultado de esta tesis sirve para decidir sobre un paciente."


def texto(valor: object) -> str:
    if isinstance(valor, bool):
        return "sí" if valor else "no"
    if isinstance(valor, float):
        return format(valor, ".12g")
    return str(valor)


def leer(ruta: str, campo: str) -> object:
    archivo = RAIZ / ruta
    if archivo.suffix == ".json":
        dato = json.loads(archivo.read_text(encoding="utf-8"))
        for parte in campo.split("."):
            dato = dato[parte]
        return dato
    with archivo.open(encoding="utf-8", newline="") as fh:
        filas = list(csv.DictReader(fh, delimiter="\t"))
    clave, columna = campo.split(".", 1)
    fila = next(item for item in filas if item[next(iter(item))] == clave)
    return fila[columna]


def cifra(ruta: str, campo: str) -> str:
    valor = texto(leer(ruta, campo))
    return f'<data data-fuente="{ruta}#{campo}">{valor}</data>'


def pagina(titulo: str, cuerpo: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{titulo}</title>
<link rel="stylesheet" href="css/v35.css">
</head>
<body>
<a class="saltar" href="#contenido">Saltar al contenido</a>
<p class="aviso" role="note">{AVISO}</p>
<header>
<p class="marca">OncoGO-CRC · corrida v3.5</p>
<nav aria-label="Secciones">
<a href="index.html">Inicio</a>
<a href="metodos.html">Métodos</a>
<a href="resultados.html">Resultados</a>
<a href="datos.html">Datos</a>
<a href="glosario.html">Glosario</a>
<a href="preguntas.html">Preguntas</a>
<a href="reproducir.html">Reproducir</a>
<a href="historial.html">Historial</a>
</nav>
</header>
<main id="contenido">
{cuerpo}
</main>
<footer><p>{AVISO}</p></footer>
<script src="js/v35.js"></script>
</body>
</html>
"""


def tarjeta(titulo: str, simple: str, hecho: str, salio: str, significa: str, no: str, enlace: str) -> str:
    return f"""<article class="tarjeta">
<h3>{titulo}</h3>
<h4>En palabras simples</h4><p>{simple}</p>
<h4>Qué hicimos</h4><p>{hecho}</p>
<h4>Qué salió</h4><p>{salio}</p>
<h4>Qué significa</h4><p>{significa}</p>
<h4>Qué no significa</h4><p>{no}</p>
<p><a href="{enlace}">Ficha técnica</a></p>
</article>"""


def escribir() -> None:
    v34 = SITIO / "v3_4.html"
    if not v34.exists():
        shutil.copyfile(SITIO / "index.html", v34)
    (SITIO / "css").mkdir(exist_ok=True)
    (SITIO / "js").mkdir(exist_ok=True)
    (SITIO / "css" / "v35.css").write_text(CSS, encoding="utf-8")
    (SITIO / "js" / "v35.js").write_text(JS, encoding="utf-8")
    suero = (
        f"Ninguna proteína separa el cáncer de la colonoscopia negativa: "
        f"{cifra('results/v3_5/multiverso_v3_4.json', 'celda_preregistrada.diferente')} de "
        f"{cifra('results/v3_5/multiverso_v3_4.json', 'celda_preregistrada.n')}."
    )
    paginas = {
        "index.html": pagina("OncoGO-CRC", inicio(suero)),
        "metodos.html": pagina("Métodos", metodos()),
        "resultados.html": pagina("Resultados", resultados(suero)),
        "datos.html": pagina("Datos", datos()),
        "glosario.html": pagina("Glosario", glosario()),
        "preguntas.html": pagina("Preguntas", preguntas(suero)),
        "reproducir.html": pagina("Reproducir", reproducir()),
        "historial.html": pagina("Historial", historial()),
    }
    for nombre, html in paginas.items():
        (SITIO / nombre).write_text(html, encoding="utf-8")
    print("[v3.5] sitio escrito")


def inicio(suero: str) -> str:
    return f"""
<h1>¿Las proteínas que cambian en un tumor de colon dejan una huella que se pueda ver en la sangre?</h1>
<p class="respuesta">No en este suero. {suero} La separación del modelo se asocia con la calidad de la muestra.</p>
<p>El certificado de esta corrida está en {cifra('results/v3_5/certificado_v3_5.json', 'estado')}.</p>
<section class="rejilla">
<article><h2>Tejido</h2><p>El signo del cambio tumoral se repite: {cifra('results/v3_3/k1_etiquetas.json', 'contrastes.primario_frente_a_mucosa.numerador')} de {cifra('results/v3_3/k1_etiquetas.json', 'contrastes.primario_frente_a_mucosa.denominador')}.</p></article>
<article><h2>Sangre</h2><p>{suero} {cifra('results/v3_5/multiverso_v3_4.json', 'celda_preregistrada.equivalente')} quedan equivalentes y {cifra('results/v3_5/multiverso_v3_4.json', 'celda_preregistrada.indeterminada')} indeterminadas.</p></article>
<article><h2>Predicción</h2><p>Todas las proteínas dan {cifra('results/v3_5/ml_celda_preregistrada.json', 'auc')}. La calidad sola da {cifra('results/v3_3/ml_suero_modelos.tsv', 'M3_solo_qc.auc_oof_media')}. Sin la calidad queda {cifra('results/v3_3/ml_suero_modelos.tsv', 'M5_principal_residualizada_qc.auc_oof_media')}.</p></article>
</section>
<h2>Las hipótesis</h2>
{tarjeta("H5. ¿El colon cambia distinto?", "El cambio del colon se parece al de otros carcinomas.", "Se comparó el signo del cambio en un universo común de genes.", f"El estado es {cifra('results/v3_5/h5_especificidad_organo.json', 'estado')}.", "La proliferación es común a los carcinomas.", "No demuestra una proteína exclusiva del colon.", "resultados.html#h5")}
{tarjeta("H6. ¿Importa el tipo de célula?", "La tabla no puede responder.", "Se cruzó el origen celular con el signo del cambio.", f"El estado es {cifra('results/v3_5/h6_origen_celular.json', 'estado')}. La tabla degenerada es {cifra('results/v3_5/h6_origen_celular.json', 'tabla_degenerada')}.", "No hay variación en el desenlace.", "No es evidencia ni a favor ni en contra.", "resultados.html#h6")}
{tarjeta("H7. ¿El tejido anticipa la sangre?", "Con las variables fijadas de antemano, no.", "Se entrenó fuera de cada recurso con suficientes casos.", f"El estado es {cifra('results/v3_5/h7_visibilidad_en_sangre.json', 'estado')}.", "Los recursos no miden lo mismo.", "No es una prueba diagnóstica.", "resultados.html#h7")}
{tarjeta("H8. ¿El estadio IV lleva más proteína del tumor?", "Esta corrida no alcanza a medirlo.", "directLFQ agregó las corridas y Fisher comparó las tumorales.", f"El estado es {cifra('results/v3_5/h8_carga_tumoral.json', 'estado')}.", "Falta profundidad en un DIA sin depleción.", "No demuestra que esa proteína no exista.", "resultados.html#h8")}
{tarjeta("H9. ¿La genética explica el nivel en sangre?", "El ensayo exploratorio no lo muestra.", "Se usó un instrumento por proteína.", f"El estado es {cifra('results/v3_5/h9_wald.json', 'estado')}.", "Es una desviación de procedimiento declarada.", "No decide la tesis.", "resultados.html#h9")}
{tarjeta("H10. ¿Alguna proteína es propia del CRC?", "El atlas marca unas pocas.", "Se contaron las específicas del atlas.", f"El estado es {cifra('results/v3_5/h10_especificidad.json', 'estado')} y el conteo es {cifra('results/v3_5/h10_especificidad.json', 'n_especificas_crc')}.", "CEACAM5 ilustra el seguimiento clínico.", "No es un hallazgo nuevo ni un cribado.", "resultados.html#h10")}
<h2>Por qué la sangre engaña</h2>
<p>La muestra pasa por el tubo, la espera, la centrífuga, el congelador y el espectrómetro. En cada paso puede entrar piel, eritrocito, plaqueta o coagulación. Esas señales pueden separar grupos sin que haya una proteína del tumor.</p>
<h2>Qué queda abierto</h2>
<p>La preanalítica queda en {cifra('results/v3_5/preanalitica_v3_5.json', 'estado')}: hay fecha en {cifra('results/v3_5/preanalitica_v3_5.json', 'n_con_fecha')} de {cifra('results/v3_5/preanalitica_v3_5.json', 'n_crudos_validos')} crudos. No se inventan los tiempos que faltan.</p>
"""


def metodos() -> str:
    return """
<h1>Cómo se hizo</h1>
<ol>
<li>Se fijó la pregunta y las reglas antes de mirar el resultado nuevo.</li>
<li>En el tejido se comparó el signo del cambio entre dos cohortes públicas.</li>
<li>En el suero se clasificó cada proteína como diferente, equivalente o indeterminada.</li>
<li>Un modelo mezcló las proteínas. Otro usó solo la calidad de la muestra.</li>
<li>Cada hipótesis tuvo su regla escrita. Si faltó un dato, el estado es no evaluable.</li>
</ol>
<p>El cambio que la regla llama grande es el que ya estaba escrito. No se movió después de ver el resultado. directLFQ agrega iones a proteína. ThermoRawFileParser lee la fecha del crudo. Un recurso con pocos casos no entrena a los demás.</p>
"""


def resultados(suero: str) -> str:
    return f"""
<h1>Resultados</h1>
<p>Cada número nombra su fichero. Un estado vacío no se rellena.</p>
<section id="suero"><h2>Suero</h2><p>{suero}</p>
<p>Equivalentes: {cifra('results/v3_5/multiverso_v3_4.json', 'celda_preregistrada.equivalente')}. Indeterminadas: {cifra('results/v3_5/multiverso_v3_4.json', 'celda_preregistrada.indeterminada')}.</p>
<svg role="img" aria-labelledby="embudo-titulo" viewBox="0 0 420 140">
<title id="embudo-titulo">Embudo del suero, de grupos del informe a la celda principal</title>
<rect x="10" y="20" width="90" height="36"></rect>
<rect x="120" y="28" width="80" height="28"></rect>
<rect x="220" y="34" width="70" height="24"></rect>
<rect x="310" y="40" width="60" height="20"></rect>
</svg>
<p>Grupos del informe: {cifra('results/v3_3/suero_resumen.json', 'n_grupos_report')}. Con el filtro global: {cifra('results/v3_3/suero_controles.json', 'controles.grupos_global_pg_q_0_01')}. Humanos: {cifra('results/v3_3/suero_resumen.json', 'n_matriz_humana_fdr')}. Celda principal: {cifra('results/v3_3/suero_resumen.json', 'n_principal')}.</p>
</section>
<section id="modelo"><h2>Modelo</h2>
<p>La calidad sola separa casi tan bien como todas las proteínas. Sin la calidad, la separación se pierde. Con contaminantes, la mediana sube a {cifra('results/v3_5/multiverso_resumen_v3_5.json', 'auc_ml_sin_quitar_mediana')}.</p>
<p>La celda de esta corrida reproduce {cifra('results/v3_5/ml_celda_preregistrada.json', 'auc')}. La permutación de la versión sin calidad da {cifra('results/v3_3/ml_suero_modelos.tsv', 'M5_principal_residualizada_qc.p_permutacion')}.</p>
</section>
<section id="pre"><h2>Orden de adquisición</h2>
<p>Estado {cifra('results/v3_5/preanalitica_v3_5.json', 'estado')}. El área del orden, solo con los crudos legibles, es {cifra('results/v3_5/preanalitica_v3_5.json', 'calculado_con_legibles.auc.valor')}. Los grupos legibles se midieron por bloques. Eso no cambia el titular.</p>
</section>
<section id="h5"><h2>H5</h2><p>Estado {cifra('results/v3_5/h5_especificidad_organo.json', 'estado')}. El universo común tiene {cifra('results/v3_5/h5_especificidad_organo.json', 'universo_principal_n')} genes. No demuestra que el colon tenga una dirección propia.</p></section>
<section id="h6"><h2>H6</h2><p>Estado {cifra('results/v3_5/h6_origen_celular.json', 'estado')}. Tabla degenerada: {cifra('results/v3_5/h6_origen_celular.json', 'tabla_degenerada')}. No es evidencia ni a favor ni en contra.</p></section>
<section id="h7"><h2>H7</h2><p>Estado {cifra('results/v3_5/h7_visibilidad_en_sangre.json', 'estado')}. Álvez: {cifra('results/v3_5/h7_visibilidad_en_sangre.json', 'por_recurso.alvez_2023.auc')}. HPA: {cifra('results/v3_5/h7_visibilidad_en_sangre.json', 'por_recurso.hpa_disease_blood_atlas.auc')}. Papier queda {cifra('results/v3_5/h7_visibilidad_en_sangre.json', 'por_recurso.papier_2024.estado')}.</p></section>
<section id="h8"><h2>H8</h2><p>Estado {cifra('results/v3_5/h8_carga_tumoral.json', 'estado')}. Corridas unidas: {cifra('results/v3_5/h8_carga_tumoral.json', 'mapa.n_unidas')}. Estándares fuera: {cifra('results/v3_5/h8_carga_tumoral.json', 'n_estandares_excluidos')}. Tumorales epiteliales: {cifra('results/v3_5/h8_carga_tumoral.json', 'n_tumorales_epiteliales')}. Versión de directLFQ: {cifra('results/v3_5/h8_carga_tumoral.json', 'directlfq_version')}.</p></section>
<section id="h9"><h2>H9</h2><p>Estado {cifra('results/v3_5/h9_wald.json', 'estado')}. No decide la tesis.</p></section>
<section id="h10"><h2>H10 y triangulación</h2><p>Estado {cifra('results/v3_5/h10_especificidad.json', 'estado')}. Específicas: {cifra('results/v3_5/h10_especificidad.json', 'n_especificas_crc')}. Filas de la triangulación: {cifra('results/v3_5/triangulacion_v3_5.json', 'n_filas')}. Celdas vacías: {cifra('results/v3_5/triangulacion_v3_5.json', 'celdas_vacias')}. No es una validación.</p></section>
"""


def datos() -> str:
    filas = [
        ("PXD055218", "Suero de esta tesis. Cáncer, colonoscopia negativa y personas sanas. No se fusiona con las demás."),
        ("PXD034575", "Tejido de Tanaka. Primario, metástasis hepática y mucosa."),
        ("PDC000116", "Tejido de CPTAC. Sirve para repetir el signo del cambio."),
        ("PXD064655", "Suero por estadio. Entra solo en H8."),
        ("PXD033236", "Vesículas. Entra como variable secundaria."),
        ("HPA y Álvez", "Atlas de sangre. No se llaman personas sanas de esta cohorte."),
        ("Papier", "Riesgo prospectivo. No entrena a los otros recursos."),
        ("GCST90255676", "Genética de proteínas. H9 exploratorio."),
    ]
    cuerpo = "".join(f"<tr><td>{a}</td><td>{b}</td></tr>" for a, b in filas)
    return f"<h1>De dónde salen los datos</h1><table><caption>Cada cohorte se lee sola</caption><thead><tr><th>Accesión</th><th>Para qué se usa</th></tr></thead><tbody>{cuerpo}</tbody></table>"


def glosario() -> str:
    terminos = {
        "proteina": "Molécula que la célula fabrica a partir de un gen.",
        "proteoma": "El conjunto de proteínas que se midió en una muestra.",
        "espectrometria": "Técnica que pesa los fragmentos de las proteínas.",
        "dia": "Modo de medición que parte los fragmentos en ventanas.",
        "peptido": "Trozo de proteína que el aparato mide.",
        "precursor": "El péptido con su carga, antes de fragmentarlo.",
        "suero": "La parte líquida de la sangre después de coagular.",
        "plasma": "La parte líquida cuando se impide la coagulación.",
        "colonoscopia-negativa": "Exploración del colon sin hallar ese cáncer.",
        "cohorte": "Un grupo de personas medido con el mismo diseño.",
        "contraste": "La comparación escrita entre dos grupos.",
        "log2fc": "Cuántas veces cambia una proteína, en escala de doble.",
        "cambio-regla": "El cambio que la regla, escrita antes, llama grande.",
        "fdr": "La fracción de hallazgos falsos que se acepta.",
        "bh": "El método de Benjamini y Hochberg para esa fracción.",
        "treat": "Prueba que pide un cambio mínimo, no solo uno distinto de cero.",
        "equivalencia": "El cambio grande queda descartado.",
        "indeterminada": "Ni se demuestra el cambio ni se descarta.",
        "auc": "Cuánto separa un modelo a dos grupos. El azar queda a la mitad.",
        "validacion-cruzada": "Se entrena en una parte y se prueba en la otra.",
        "permutacion": "Se mezclan las etiquetas para ver el azar.",
        "preregistro": "Las reglas escritas antes de mirar el resultado.",
        "enmienda": "Un cambio de reglas fechado, sin borrar lo anterior.",
        "multiverso": "La misma pregunta con varias reglas razonables.",
        "contaminante": "Proteína que no viene de la muestra biológica.",
        "hemolisis": "Rotura de glóbulos rojos dentro del tubo.",
        "plaquetas": "Células de la coagulación que pueden contaminar el plasma.",
        "coagulacion": "El proceso de formar coágulo. Puede ser del tubo o del cáncer.",
        "fase-aguda": "Proteínas que suben con la inflamación del huésped.",
        "tnm": "Sistema que describe el tamaño y la extensión del tumor.",
        "metastasis-hepatica": "Tumor del colon que creció en el hígado.",
        "cea": "CEACAM5. Se usa en el seguimiento, no como cribado.",
        "pqtl": "Una variante genética asociada al nivel de una proteína.",
        "mendeliana": "Usar esa variante como si fuera un experimento natural.",
        "pepe": "Fases para evaluar una señal. Esta tesis es la primera.",
        "directlfq": "Programa que resume los iones en una proteína.",
        "thermoparse": "Programa que lee la fecha dentro del crudo.",
    }
    items = "".join(f'<dt id="{k}">{k.replace("-", " ")}</dt><dd>{v}</dd>' for k, v in terminos.items())
    return f"<h1>Glosario</h1><dl>{items}</dl>"


def preguntas(suero: str) -> str:
    return f"""
<h1>Preguntas frecuentes</h1>
<details open><summary>¿Esto sirve para detectar el cáncer?</summary><p>No. {AVISO} {suero}</p></details>
<details><summary>¿Por qué un resultado negativo vale?</summary><p>Porque descarta el cambio que la regla había fijado. No dice que la proteína no exista.</p></details>
<details><summary>¿Qué es una colonoscopia negativa?</summary><p>Es una persona con síntomas a la que no se le encontró este cáncer en la exploración.</p></details>
<details><summary>¿Por qué el grupo es pequeño?</summary><p>Porque muchas proteínas quedan indeterminadas. Falta potencia para un cambio moderado.</p></details>
<details><summary>¿Qué significa que la calidad explique la diferencia?</summary><p>La piel y la coagulación separan los grupos. Esa señal puede nacer en el tubo.</p></details>
"""


def reproducir() -> str:
    return f"""
<h1>Cómo repetir la corrida</h1>
<p>El comando es <code>make all-v3-5</code>. Escribe solo <code>results/v3_5/</code>. El certificado queda {cifra('results/v3_5/certificado_v3_5.json', 'estado')} y se puede publicar: {cifra('results/v3_5/certificado_v3_5.json', 'puede_publicar')}.</p>
<p>Los sellos están junto a cada fichero. El manifiesto es <code>results/v3_5_manifest.sha256</code>.</p>
"""


def historial() -> str:
    return """
<h1>Qué cambió entre versiones</h1>
<p>v3.2 contó 205 proteínas y una concordancia de 89/90. Esas cifras pertenecen a esa versión. v3.3 rehízo el tejido y el suero. v3.4 añadió vesículas y la secuencia. En esa página se decía que directLFQ no está instalado y se contaban 120 crudos. v3.5 corrige las dos frases: directLFQ se ejecutó y los crudos válidos no incluyen los ficheros de AppleDouble. H5 de v3.3 decía refuta. La corrida v3.5 dice sostiene. Las etiquetas del tejido no quedan confirmadas del todo: solo CA1 se comprueba.</p>
<p>La tabla de cambios está en <code>results/v3_5/changelog_numeros_v3_5.tsv</code>.</p>
"""


CSS = """
:root { color-scheme: light dark; --tinta:#1c1915; --fondo:#f7f3ea; --aviso:#6b1d1d; --linea:#c8bba6; }
@media (prefers-color-scheme: dark) {
  :root { --tinta:#f4efe4; --fondo:#1c1915; --aviso:#ffb4a8; --linea:#4a4036; }
}
* { box-sizing: border-box; }
body { margin:0; font:1.05rem/1.5 Georgia, "Iowan Old Style", serif; color:var(--tinta); background:var(--fondo); }
a { color:inherit; }
header, main, footer { padding:1rem 1.2rem; max-width:46rem; margin:auto; }
nav { display:flex; flex-wrap:wrap; gap:.7rem; }
.aviso, footer { background:#f3d6d0; color:#3d1010; }
@media (prefers-color-scheme: dark) { .aviso, footer { background:#3d1010; color:#ffd9d3; } }
.saltar { position:absolute; left:-999px; }
.saltar:focus { left:1rem; background:var(--fondo); }
h1, h2, h3 { line-height:1.2; }
.respuesta { font-size:1.25rem; }
.rejilla { display:grid; gap:1rem; }
@media (min-width:40rem) { .rejilla { grid-template-columns:repeat(3,1fr); } }
article, section { border-top:1px solid var(--linea); padding-top:.6rem; }
data { font-variant-numeric:tabular-nums; font-weight:700; }
svg rect { fill:currentColor; opacity:.75; }
table { border-collapse:collapse; width:100%; }
td, th { border-bottom:1px solid var(--linea); text-align:left; padding:.3rem; }
@media (prefers-reduced-motion: reduce) { * { scroll-behavior:auto; transition:none; } }
@media print { nav, .saltar, script { display:none; } }
"""

JS = """
document.querySelectorAll("details").forEach((nodo) => {
  nodo.addEventListener("toggle", () => nodo.setAttribute("data-abierto", nodo.open ? "si" : "no"));
});
"""


if __name__ == "__main__":
    escribir()
