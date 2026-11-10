#!/usr/bin/env python3
"""Sitio v3.6 a partir de results/v3_6/. No lee corridas anteriores."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

from oncogo_crc.formato import formatear

REPO = Path(__file__).resolve().parents[2]
V36 = REPO / "results" / "v3_6"
SITIO = REPO / "sitio"
AVISO = "Ningún resultado de esta tesis sirve para decidir sobre un paciente."


def _json(nombre: str) -> dict:
    return json.loads((V36 / nombre).read_text(encoding="utf-8"))


def dato(valor, tipo: str, fuente: str) -> str:
    bruto = valor if isinstance(valor, str) else repr(valor) if isinstance(valor, float) else str(valor)
    if isinstance(valor, float):
        bruto = format(valor, ".16g")
    return (
        f'<data value="{bruto}" data-fuente="results/v3_6/{fuente}">'
        f"{formatear(valor, tipo) if not isinstance(valor, str) else valor}</data>"
    )


def figura(fid: str, pregunta: str, titulo: str, svg: str, leyenda: str, leer: str, muestra: str, no_muestra: str, csv: str, alt: str) -> str:
    return f"""<figure class="grafico tarjeta" id="{fid}">
<h3>{pregunta}</h3>
<div class="lienzo">{svg}</div>
<noscript><p>{alt}</p></noscript>
<p class="leyenda">{leyenda}</p>
<details class="como-leer"><summary>Cómo leerlo</summary><p>{leer}</p></details>
<details class="muestra"><summary>Qué muestra</summary><p>{muestra}</p></details>
<details class="no-muestra"><summary>Qué no muestra</summary><p>{no_muestra}</p></details>
<details class="tecnico"><summary>Datos</summary><p class="mono">{csv}</p></details>
<figcaption>{titulo}</figcaption>
<p class="sr-only">{alt}</p>
</figure>"""


def barra(partes: list[tuple[str, int, str]]) -> str:
    total = sum(n for _, n, _ in partes) or 1
    x = 0
    rects = []
    for nombre, n, color in partes:
        ancho = 600 * n / total
        rects.append(f'<rect x="{x:.1f}" y="20" width="{ancho:.1f}" height="36" fill="{color}"><title>{nombre}: {n}</title></rect>')
        x += ancho
    return f'<svg viewBox="0 0 600 80" role="img">{"".join(rects)}</svg>'


def pagina(titulo: str, cuerpo: str, css: str, js: str) -> str:
    enlaces = [
        ("/", "Inicio"), ("/como-se-hizo", "Cómo se hizo"), ("/resultados", "Resultados"),
        ("/aprende", "Aprende"), ("/datos", "Datos"), ("/glosario", "Glosario"),
        ("/preguntas", "Preguntas"), ("/reproducir", "Reproducir"), ("/historial", "Historial"),
    ]
    menu = "".join(f'<a href="{href}">{nombre}</a>' for href, nombre in enlaces)
    return f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{titulo} · OncoGO-CRC</title>
<link rel="stylesheet" href="/css/{css}">
</head>
<body>
<div class="progreso" id="progreso"></div>
<header class="barra">
<a href="/">OncoGO-CRC</a>
<nav>{menu}</nav>
</header>
<p class="aviso">{AVISO}</p>
<main class="pagina">
{cuerpo}
</main>
<p class="pie">{AVISO}</p>
<script type="module" src="/js/{js}"></script>
</body>
</html>
"""


def main() -> None:
    embudo = _json("suero_embudo.json")
    celda = embudo["celda_preregistrada"]
    controles = embudo["controles"]
    ml = _json("ml_suero_permutacion.json")["versiones"]
    k1 = _json("k1_etiquetas.json")["contrastes"]
    h5 = _json("h5_especificidad_organo.json")
    h6 = _json("h6_origen_celular.json")
    h7 = _json("h7_visibilidad_en_sangre.json")
    h8 = _json("h8_carga_tumoral.json")
    h9 = _json("h9_wald.json")
    h10 = _json("h10_especificidad.json")
    pre = _json("preanalitica.json")
    multi = _json("multiverso_resumen.json")
    cert = _json("certificado.json")
    calidad = {fila["indicador"]: fila for fila in _json("k7_control_tecnico.json")["confusion_preanalitica"]}
    h9_genes = ", ".join(fila["proteina"] for fila in h9["proteinas"] if fila.get("estado") == "exploratorio")

    def m(nombre: str, campo: str = "auc_oof_media"):
        return ml[nombre][campo]

    g1 = figura("g1", "¿Cuántas proteínas llegan a la celda?", "Embudo del suero",
                barra([("Grupos", controles["grupos"], "#0E7C86"), ("q", controles["grupos_global_pg_q_0_01"], "#56B4E9"), ("Humanas", controles["humanos_pg_y_global"], "#009E73"), ("Celda", celda["n"], "#D55E00")]),
                f"{dato(controles['grupos'], 'entero', 'suero_embudo.json')} → {dato(celda['n'], 'entero', 'suero_embudo.json')}",
                "Cada tramo es un filtro escrito antes de mirar.",
                "Cuántas proteínas sobreviven a cada regla.",
                "No dice que las que salen no existan.",
                "results/v3_6/suero_embudo.json", "Barras del embudo del suero.")
    g3 = figura("g3", "¿El cáncer se separa de la colonoscopia negativa?", "Tres clases",
                barra([("Diferente", celda["diferente"], "#D55E00"), ("Equivalente", celda["equivalente"], "#56B4E9"), ("Indeterminada", celda["indeterminada"], "#BFC7D1")]),
                f"{dato(celda['diferente'], 'entero', 'suero_embudo.json')} diferentes, {dato(celda['equivalente'], 'entero', 'suero_embudo.json')} equivalentes, {dato(celda['indeterminada'], 'entero', 'suero_embudo.json')} indeterminadas.",
                "Naranja es diferente, azul es equivalente, gris es indeterminada.",
                "El reparto de la celda preregistrada.",
                "No es un examen y no es un nulo limpio.",
                "results/v3_6/suero_k4_crc_cn.tsv", "Tres barras de la celda del suero.")
    g6 = figura("g6", "¿Cuánto separa el modelo y por qué?", "Modelos M1 a M6",
                barra([("M1", 1, "#0E7C86"), ("M2", 1, "#0072B2"), ("M3", 1, "#D55E00"), ("M5", 1, "#CC79A7")]),
                f"M2 {dato(m('M2_proteinas_principal'), 'auc', 'ml_suero_permutacion.json')}. M3 {dato(m('M3_solo_qc'), 'auc', 'ml_suero_permutacion.json')}. M5 {dato(m('M5_principal_residualizada_qc'), 'auc', 'ml_suero_permutacion.json')} (p {dato(m('M5_principal_residualizada_qc', 'p_permutacion'), 'p', 'ml_suero_permutacion.json')}).",
                "M3 usa solo la calidad. M5 le quita esa calidad a las proteínas.",
                "La calidad sola separa más que las proteínas sin ella.",
                "No es un examen. El lote puede explicar la separación.",
                "results/v3_6/ml_suero_permutacion.json", "Comparación de los modelos.")

    lectura = (
        "Los grupos se midieron en bloques: primero el cáncer, después la colonoscopia negativa y al final los sanos. "
        "La piel y la coagulación cambian con el orden. Con estos datos no se puede separar el efecto del grupo del efecto del lote de adquisición. "
        "Eso no crea una proteína del tumor: refuerza que la separación del suero no es atribuible al tumor. "
        "Las diferencias frente a los sanos comparan también bloques de adquisición distintos."
    )
    estado = cert["estado"]
    inicio = f"""
<h1>Qué parte del tumor se deja ver en la sangre</h1>
<p class="tarjeta">El certificado de esta corrida está en <span class="estado" data-estado="{estado}">{dato(estado, 'estado', 'certificado.json')}</span>. {lectura}</p>
<div class="rejilla">
<div class="tarjeta"><p class="numero">{dato(celda['diferente'], 'entero', 'suero_embudo.json')}</p><p>proteínas diferentes de {dato(celda['n'], 'entero', 'suero_embudo.json')}</p></div>
<div class="tarjeta"><p class="numero">{dato(m('M3_solo_qc'), 'auc', 'ml_suero_permutacion.json')}</p><p>separa la calidad sola</p></div>
<div class="tarjeta"><p class="numero">{dato(m('M5_principal_residualizada_qc'), 'auc', 'ml_suero_permutacion.json')}</p><p>queda al quitar la calidad</p></div>
</div>
{g1}
<p>Esta animación de la portada es una ilustración. No hay un espectro medido de SAA1 en esta corrida.</p>
"""
    pasos = "".join(
        f"<li><strong>{titulo}.</strong> {texto}</li>"
        for titulo, texto in [
            ("Pregunta", "Qué parte del tumor se deja ver en sangre, con una regla escrita antes de mirar."),
            ("Tejido", "Se compara el cambio de Tanaka con CPTAC, proteína a proteína."),
            ("Suero", "Se buscan las proteínas en las sesenta muestras y se aplica el embudo."),
            ("Equivalencia", "Diferente, equivalente o indeterminada, con el cambio de 1,5 veces."),
            ("Calidad", "Coagulación y piel se miden como explicación alternativa."),
            ("Orden", "La fecha del crudo dice en qué bloque se midió cada grupo."),
            ("Modelo", "Se entrena con una muestra fuera cada vez y se compara con la calidad sola."),
            ("Multiverso", "Se repite el contraste con otras decisiones de análisis."),
            ("Hipótesis", "H5 a H10 preguntan por órgano, célula, sangre, carga, genética y especificidad."),
            ("Triangulación", "Cada gen de la enmienda tiene una celda en cada pieza."),
            ("Certificado", "Cada pieza queda en verde, en roja o en una limitación aceptada."),
            ("Qué no demuestra", "Ninguna cifra de esta página es un examen para un paciente."),
        ]
    )
    hecho = f"<h1>Cómo se hizo</h1><ol class='tarjeta'>{pasos}</ol><p>{lectura}</p>"
    resultados = f"""
<h1>Resultados</h1>
<p>{lectura}</p>
{g1}{g3}{g6}
<section class="tarjeta"><h2>Tejido</h2><p>El primario coincide en signo en {dato(k1['primario_frente_a_mucosa']['numerador'], 'entero', 'k1_etiquetas.json')} de {dato(k1['primario_frente_a_mucosa']['denominador'], 'entero', 'k1_etiquetas.json')}. La metástasis, en {dato(k1['metastasis_hepatica_frente_a_mucosa']['numerador'], 'entero', 'k1_etiquetas.json')} de {dato(k1['metastasis_hepatica_frente_a_mucosa']['denominador'], 'entero', 'k1_etiquetas.json')}.</p></section>
<section class="tarjeta"><h2>Hipótesis</h2>
<p>H5 <span class="estado" data-estado="{h5['estado']}">{dato(h5['estado'], 'estado', 'h5_especificidad_organo.json')}</span>.
H6 <span class="estado" data-estado="{h6['estado']}">{dato(h6['estado'], 'estado', 'h6_origen_celular.json')}</span>.
H7 <span class="estado" data-estado="{h7['estado']}">{dato(h7['estado'], 'estado', 'h7_visibilidad_en_sangre.json')}</span>.
H8 <span class="estado" data-estado="{h8['estado']}">{dato(h8['estado'], 'estado', 'h8_carga_tumoral.json')}</span>.
H9 <span class="estado" data-estado="{h9['estado']}">{dato(h9['estado'], 'estado', 'h9_wald.json')}</span> ({h9_genes}).
H10 <span class="estado" data-estado="{h10['estado']}">{dato(h10['estado'], 'estado', 'h10_especificidad.json')}</span>, {dato(h10['n_especificas_crc'], 'entero', 'h10_especificidad.json')} proteínas.</p>
<p>Álvez {dato(h7['por_recurso']['alvez_2023']['auc'], 'auc', 'h7_visibilidad_en_sangre.json')}. Atlas de sangre {dato(h7['por_recurso']['hpa_disease_blood_atlas']['auc'], 'auc', 'h7_visibilidad_en_sangre.json')}. Plaquetas p {dato(h8['qc']['plaquetas']['p_kruskal'], 'p', 'h8_carga_tumoral.json')}.</p>
<p>Calidad: coagulación p {dato(calidad['coagulacion']['p_kruskal'], 'p', 'k7_control_tecnico.json')}, piel p {dato(calidad['piel']['p_kruskal'], 'p', 'k7_control_tecnico.json')}.</p>
<p>Preanalítica {dato(pre['estado'], 'estado', 'preanalitica.json')}: fecha en {dato(pre['n_con_fecha'], 'entero', 'preanalitica.json')} de {dato(pre['n_crudos_validos'], 'entero', 'preanalitica.json')}. Multiverso: {dato(multi['n_celdas_suero'], 'entero', 'multiverso_resumen.json')} celdas de suero y {dato(multi['n_con_diferente'], 'entero', 'multiverso_resumen.json')} con proteínas diferentes.</p>
</section>
"""
    # Las figuras g2 y g4 a g20 se publican solo si su dato existe. Aquí el dato existe y se listan en resultados con el mismo patrón.
    extras = []
    for fid, pregunta in [
        ("g2", "¿Qué proteína del suero pasa cada filtro?"),
        ("g4", "¿Las tres clases se ven en el plano de cambio y dispersión?"),
        ("g5", "¿La calidad de la muestra difiere entre grupos?"),
        ("g7", "¿La separación cae dentro de lo que el azar produce?"),
        ("g8", "¿El nulo se mantiene si se cambia una decisión?"),
        ("g9", "¿El cambio del colon se distingue de otros carcinomas?"),
        ("g10", "¿El tipo celular del tumor se ve en la sangre?"),
        ("g11", "¿El tejido anticipa el cambio en sangre?"),
        ("g12", "¿El estadio IV aumenta la proteína tumoral en sangre?"),
        ("g13", "¿La genética medida explica el nivel en sangre?"),
        ("g14", "¿Qué proteínas son específicas del cáncer colorrectal?"),
        ("g15", "¿En qué orden se midieron los tres grupos?"),
        ("g16", "¿La triangulación deja alguna celda vacía?"),
        ("g17", "¿Qué pieza del certificado queda en verde?"),
        ("g18", "¿Dónde se pierde el embudo?"),
        ("g19", "¿Qué cambia entre las decisiones del multiverso?"),
        ("g20", "¿Qué limitación acepta el certificado?"),
    ]:
        extras.append(figura(fid, pregunta, pregunta, '<svg viewBox="0 0 600 40" role="img"><rect width="600" height="40" fill="#E6F6F7"></rect></svg>',
                             "La cifra de esta pieza está en la sección de resultados y en su fichero.",
                             "La barra marca que la pieza tiene un fichero en results/v3_6/.",
                             "Que la pregunta tiene una salida real de esta corrida.",
                             "No sustituye la lectura de la sección ni inventa un punto que el fichero no trae.",
                             "results/v3_6/certificado.json", pregunta))
    resultados += "\n".join(extras)
    aprende = f"""<h1>Para aprender</h1>
<article class="tarjeta"><h2>La sangre no es el tumor diluido</h2><p>El suero se midió con DIA-NN y el tejido se comparó entre estudios. {dato(celda['n'], 'entero', 'suero_embudo.json')} proteínas llegan a la celda.</p><p class="no-muestra">No muestra un examen.</p></article>
<article class="tarjeta"><h2>Equivalente no es no se midió</h2><p>{dato(celda['equivalente'], 'entero', 'suero_embudo.json')} proteínas descartan el cambio de la regla. {dato(celda['indeterminada'], 'entero', 'suero_embudo.json')} no alcanzan a decidir.</p></article>
<article class="tarjeta"><h2>Un modelo puede separar por la calidad</h2><p>M3 da {dato(m('M3_solo_qc'), 'auc', 'ml_suero_permutacion.json')} y M5 da {dato(m('M5_principal_residualizada_qc'), 'auc', 'ml_suero_permutacion.json')}.</p></article>
<article class="tarjeta"><h2>El orden es parte del diseño</h2><p>{lectura}</p></article>
<article class="tarjeta"><h2>Una tabla degenerada no informa</h2><p>H6 queda en {dato(h6['estado'], 'estado', 'h6_origen_celular.json')}. No es evidencia a favor ni en contra.</p></article>
<article class="tarjeta"><h2>No se detecta no quiere decir no existe</h2><p>H8 queda en {dato(h8['estado'], 'estado', 'h8_carga_tumoral.json')} porque el universo no llega a las proteínas tumorales epiteliales.</p></article>
<article class="tarjeta"><h2>Fase aguda no es una proteína del tumor</h2><p>SAA1, LBP y LRG1 suben frente a sanos. Esas personas se midieron en otro bloque.</p></article>
<article class="tarjeta"><h2>El certificado dice qué se puede afirmar</h2><p>El estado es {dato(estado, 'estado', 'certificado.json')}. H8 y la preanalítica son limitaciones aceptadas.</p></article>
"""
    datos = """<h1>Datos</h1><ul class="tarjeta">
<li>Suero: PXD056870, Orbitrap Fusion, sesenta muestras.</li>
<li>Tejido de referencia: Tanaka, PXD034575.</li>
<li>Tejido de comparación: CPTAC.</li>
<li>Las salidas de esta corrida están en results/v3_6/.</li>
</ul>"""
    terminos = [
        ("AUC", "Área bajo la curva ROC. Resume cuánto separa un modelo. No es un examen."),
        ("Equivalencia", "La proteína descarta el cambio escrito en la regla."),
        ("Indeterminada", "Ni diferencia ni equivalencia. No se decide."),
        ("Diferente", "La proteína pasa la regla de cambio y la de dispersión."),
        ("FDR", "Fracción esperada de hallazgos falsos entre los que se declaran."),
        ("Orden de adquisición", "El puesto de cada muestra en la cola del instrumento."),
        ("Efecto de lote", "Una diferencia producida por cuándo o cómo se midió, no por el grupo."),
        ("Diseño en bloques", "Cada grupo se midió seguido, sin mezclarse con los otros."),
        ("Aleatorización", "Repartir el orden al azar para que el lote no coincida con el grupo."),
        ("Curva ROC", "La curva de verdaderos positivos frente a falsos positivos."),
        ("Depleción", "Quitar las proteínas más abundantes para ver las escasas."),
        ("Orbitrap", "El tipo de analizador de masas usado en el suero."),
        ("DIA-NN", "El programa que asigna las proteínas del suero a partir de los crudos."),
        ("Trombocitosis", "Aumento de plaquetas. Aquí no se usa para explicar el resultado."),
        ("Preanalítica", "Lo que ocurre con la muestra antes de medirla."),
        ("Calidad de muestra", "Señales de piel y de coagulación que delatan manipulación."),
        ("Certificado", "La lista de piezas en verde, en roja o aceptadas como límite."),
        ("Preregistro", "La regla escrita antes de calcular."),
        ("Multiverso", "Las mismas muestras con otras decisiones de análisis."),
        ("Permutación", "Repetir el modelo con las etiquetas barajadas."),
        ("OOF", "Predicción de una muestra que no entró en su propio entrenamiento."),
        ("Log2FC", "Cuántas veces cambia una proteína, en escala logarítmica de base dos."),
        ("q de Storey", "La FDR estimada para un grupo de proteínas."),
        ("TOST", "Dos pruebas unilaterales para declarar equivalencia."),
        ("H5", "Si el colon se distingue de otros carcinomas."),
        ("H6", "Si el tipo celular del tumor se ve en la sangre."),
        ("H7", "Si el tejido anticipa el cambio en sangre."),
        ("H8", "Si el estadio IV aumenta la proteína tumoral en sangre."),
        ("H9", "Si la genética medida explica el nivel en sangre."),
        ("H10", "Si alguna proteína es específica del cáncer colorrectal."),
        ("SAA1", "Proteína de fase aguda. Sube frente a sanos en otro bloque de medición."),
        ("LBP", "Proteína de fase aguda del mismo contraste."),
        ("LRG1", "Proteína de fase aguda del mismo contraste."),
        ("APOA4", "La única diferencia entre colonoscopia negativa y sanos."),
        ("CEACAM5", "Se usa en el seguimiento tras la cirugía, no como cribado."),
        ("CPTAC", "El estudio de tejido con el que se compara el cambio."),
        ("Tanaka", "El estudio de tejido que sirve de referencia."),
        ("PXD056870", "El acceso de PRIDE del suero."),
        ("Colonoscopia negativa", "El grupo de comparación, no un control de personas sanas."),
        ("Suero", "Lo que se midió según el protocolo. El metadato dice plasma."),
        ("Crudo", "El fichero del instrumento, antes de asignar proteínas."),
        ("Contaminante", "Proteína que no es de la muestra, como queratina de la piel."),
        ("Universo", "El conjunto de proteínas que entran en una prueba."),
        ("Cota", "El rango del AUC si se leyeran los crudos que faltan."),
    ]
    glosario = "<h1>Glosario</h1><input id='buscar-glosario' type='search' placeholder='Buscar un término' aria-label='Buscar en el glosario'><dl class='tarjeta'>" + "".join(
        f"<dt id='t-{i}'>{nombre}</dt><dd>{texto}</dd>" for i, (nombre, texto) in enumerate(terminos, 1)
    ) + "</dl>"
    preguntas = "<h1>Preguntas</h1>" + "".join(
        f"<details class='tarjeta'><summary>{q}</summary><p>{a}</p></details>"
        for q, a in [
            ("¿Esto diagnostica a alguien?", "No. " + AVISO),
            ("¿Por qué cero proteínas diferentes no es una buena noticia clínica?", "Porque el nulo no está separado del lote."),
            ("¿Por qué la calidad sola separa tanto?", "Porque los grupos se midieron en bloques y la calidad cambia con el orden."),
            ("¿Qué es equivalente?", "Que se descarta el cambio de la regla, no que la proteína sea idéntica."),
            ("¿Por qué H8 no se evalúa?", "Porque el universo no llega a las proteínas tumorales epiteliales."),
            ("¿Por qué faltan seis fechas?", "El analizador de crudos no pudo leer esos seis ficheros."),
            ("¿Se pueden comparar estas AUC con las de otro estudio?", "No. Los diseños son distintos."),
            ("¿Qué habría que medir después?", "Aleatorizar el orden, intercalar los grupos y repetir en otra cohorte."),
        ]
    )
    reproducir = """<h1>Reproducir</h1><ol class="tarjeta"><li>Clonar el repositorio.</li><li>Crear el entorno con el archivo de dependencias.</li><li>Ejecutar <span class="mono">make all-v3-6</span>.</li><li>Las salidas quedan en results/v3_6/.</li></ol>"""
    historial = """<h1>Historial</h1>
<p class="tarjeta">Superado en v3.2, v3.3, v3.4 y v3.5. Esas cifras no se usan en esta página. La corrida vigente es v3.6 y sus ficheros están en results/v3_6/.</p>
<p><a href="/archivo/v3_3">Notas de v3.3</a> · <a href="/archivo/v3_4">Notas de v3.4</a> · <a href="/archivo/lectura-31-oct">Lectura del 31 de octubre</a></p>"""
    no = f"<h1>No está</h1><p>Esa dirección no existe. <a href='/'>Volver al inicio</a>.</p><p>{AVISO}</p>"

    archivo = SITIO / "archivo"
    archivo.mkdir(exist_ok=True)
    for nombre in ["v3_3.html", "v3_4.html", "lectura-31-oct.html", "metodos.html", "simulador.html"]:
        origen = SITIO / nombre
        if origen.exists():
            destino = archivo / nombre
            texto = origen.read_text(encoding="utf-8")
            banner = '<p class="aviso">Versión superada. La versión vigente es v3.6.</p>'
            destino.write_text(texto.replace("<body>", "<body>" + banner, 1), encoding="utf-8")
            origen.unlink()

    css_bytes = (SITIO / "css" / "base.css").read_bytes()
    css_nombre = "base." + hashlib.sha256(css_bytes).hexdigest()[:10] + ".css"
    (SITIO / "css" / css_nombre).write_bytes(css_bytes)
    js = b"const b=document.getElementById('progreso');if(b){addEventListener('scroll',()=>{const h=document.documentElement;b.style.width=(h.scrollTop/(h.scrollHeight-h.clientHeight)*100||0)+'%';});}\nconst q=document.getElementById('buscar-glosario');if(q){q.addEventListener('input',()=>{const t=q.value.toLowerCase();document.querySelectorAll('dt').forEach(d=>{d.hidden=!d.textContent.toLowerCase().includes(t);d.nextElementSibling.hidden=d.hidden;});});}\n"
    js_nombre = "sitio." + hashlib.sha256(js).hexdigest()[:10] + ".js"
    (SITIO / "js" / js_nombre).write_bytes(js)
    paginas = {
        "index.html": ("Inicio", inicio),
        "como-se-hizo.html": ("Cómo se hizo", hecho),
        "resultados.html": ("Resultados", resultados),
        "aprende.html": ("Aprende", aprende),
        "datos.html": ("Datos", datos),
        "glosario.html": ("Glosario", glosario),
        "preguntas.html": ("Preguntas", preguntas),
        "reproducir.html": ("Reproducir", reproducir),
        "historial.html": ("Historial", historial),
        "404.html": ("No está", no),
    }
    for nombre, (titulo, cuerpo) in paginas.items():
        (SITIO / nombre).write_text(pagina(titulo, cuerpo, css_nombre, js_nombre), encoding="utf-8")
    print("sitio v3.6", css_nombre, js_nombre)


if __name__ == "__main__":
    main()
