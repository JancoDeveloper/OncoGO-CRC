"""Escribe sitio/datos/v3_4.json y la página desde results/v3_4.

No toca el sitio si el certificado no dice publicable. Las cifras salen de
los JSON de esa corrida. El texto explica la pregunta, el orden y la lectura.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from oncogo_crc.stages.v3_4_report import puede_publicar

RAIZ = Path(__file__).resolve().parents[2]
SITIO = RAIZ / "sitio"
RESULTADOS = Path("/Volumes/JANCODATA/OncoGO-CRC/results/v3_4")


def _oracion(texto: object) -> str:
    limpio = str(texto or "").strip()
    if limpio and limpio[0].islower():
        limpio = limpio[0].upper() + limpio[1:]
    return _html(limpio)


def _html(texto: object) -> str:
    return str(texto).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def _finito(valor: object) -> float | None:
    try:
        numero = float(valor)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numero):
        return None
    return numero


def _num(valor: object, decimales: int = 3) -> str:
    numero = _finito(valor)
    if numero is None:
        return "—"
    return f"{numero:.{decimales}f}".replace(".", ",")


def _si(valor: object) -> str:
    if valor is True:
        return "sí"
    if valor is False:
        return "no"
    if valor is None:
        return "—"
    return _html(valor)


def _conteo(bloque: object) -> str:
    if not isinstance(bloque, dict) or not bloque:
        return "—"
    return ", ".join(f"{_entero(cantidad)} {str(nombre).replace('_', ' ')}" for nombre, cantidad in bloque.items())


NOMBRES = {
    "magnitud_en_tejido": "magnitud en el tejido",
    "abundancia_en_tejido": "abundancia en el tejido",
    "go_extracelular": "componente extracelular",
    "go_membrana": "membrana",
    "go_exosoma": "exosoma",
    "hpa_secreted_to_blood": "secretada a la sangre",
    "concentracion_plasmatica_hpa_log10": "concentración plasmática",
    "especificidad_tisular": "especificidad tisular",
    "intestino_entre_los_elevados": "intestino entre los tejidos elevados",
    "clase_h5_pan_carcinoma": "clase pan-carcinoma",
    "clase_h5_especifica_colon": "clase específica de colon",
    "clase_h5_intermedia": "clase intermedia",
    "h6_epitelial": "origen epitelial",
    "h6_no_epitelial": "origen no epitelial",
    "h6_mixta": "origen mixto",
    "presencia_en_vesiculas_pxd033236": "presencia en vesículas de tumor",
}


def _variables(nombres: object) -> str:
    if not isinstance(nombres, list) or not nombres:
        return "—"
    return ", ".join(NOMBRES.get(str(nombre), str(nombre).replace("_", " ")) for nombre in nombres)


def _entero(valor: object) -> str:
    numero = _finito(valor)
    if numero is None:
        return "—"
    return str(int(round(numero)))


def _limpio(valor: object) -> object:
    if isinstance(valor, float) and not math.isfinite(valor):
        return None
    if isinstance(valor, dict):
        return {clave: _limpio(item) for clave, item in valor.items() if clave != "sello"}
    if isinstance(valor, list):
        return [_limpio(item) for item in valor]
    return valor


def _leer(base: Path, nombre: str) -> dict | None:
    ruta = base / nombre
    if not ruta.is_file():
        return None
    datos = json.loads(ruta.read_text(encoding="utf-8"))
    return _limpio(datos) if isinstance(datos, dict) else None


def _figuras(sitio: Path) -> list[dict]:
    carpeta = sitio / "figuras"
    salida = []
    for sello in sorted(carpeta.glob("F*.sello.json")):
        datos = json.loads(sello.read_text(encoding="utf-8"))
        ident = str(datos.get("figura") or sello.name.split(".")[0])
        titulo = str(datos.get("titulo") or ident)
        if titulo.endswith(": False."):
            titulo = titulo[:-len(": False.")] + ": no."
        elif titulo.endswith(": True."):
            titulo = titulo[:-len(": True.")] + ": sí."
        png = carpeta / f"{ident}.png"
        if not png.is_file():
            continue
        salida.append({
            "id": ident,
            "titulo": titulo,
            "etapa": datos.get("etapa"),
            "fecha_utc": datos.get("fecha_utc"),
        })
    return salida


def datos_publicos(base: Path, sitio: Path = SITIO) -> dict:
    certificado = _leer(base, "certificado_v3_4.json")
    if certificado is None:
        raise FileNotFoundError(base / "certificado_v3_4.json")
    if not puede_publicar(certificado):
        raise RuntimeError(f"el certificado está en {certificado.get('estado')}; no se escribe el sitio")
    sello = certificado.get("sello") or {}
    if not sello:
        crudo = json.loads((base / "certificado_v3_4.json").read_text(encoding="utf-8"))
        sello = crudo.get("sello") or {}
    h7 = _leer(base, "h7_visibilidad_en_sangre.json") or {}
    secundario = h7.get("modelo_secundario") or {}
    return {
        "certificado": {clave: valor for clave, valor in certificado.items() if clave != "sello"},
        "sello": {
            "commit": sello.get("commit"),
            "fecha_utc": sello.get("fecha_utc"),
            "etapa": sello.get("etapa"),
        },
        "h5": _leer(base, "h5_especificidad_organo.json") or {},
        "controles": _leer(base, "controles_v3_4.json") or {},
        "h6": _leer(base, "h6_origen_celular.json") or {},
        "h7": h7,
        "h8": _leer(base, "h8_carga_tumoral.json") or {},
        "h9": _leer(base, "h9_wald.json") or {},
        "h10": _leer(base, "h10_especificidad.json") or {},
        "multiverso": _leer(base, "multiverso_v3_4.json") or {},
        "preanalitica": _leer(base, "preanalitica_v3_4.json") or {},
        "secuencia": {
            "estado": secundario.get("estado") or (_leer(base, "secuencia_v3_4.json") or {}).get("estado"),
            "n_con_secuencia": secundario.get("n_con_secuencia"),
            "n_sin_secuencia": secundario.get("n_sin_secuencia"),
            "auc_principal": secundario.get("auc_principal"),
            "auc_secundario": secundario.get("auc_secundario"),
            "incremento": secundario.get("incremento"),
            "nota": secundario.get("nota"),
            "recurso": secundario.get("recurso"),
        },
        "triangulacion": _leer(base, "triangulacion_v3_4.json") or {},
        "figuras": _figuras(sitio),
    }


def _estado(valor: object) -> str:
    texto = _html(valor or "sin estado")
    return f'<span class="estado e-{_html(valor)}">{texto}</span>'


def _auc(nombre: str, bloque: dict) -> str:
    auc = _finito(bloque.get("auc"))
    bajo = _finito(bloque.get("ic95_bajo"))
    alto = _finito(bloque.get("ic95_alto"))
    if auc is None or bajo is None or alto is None:
        return ""

    def x(valor: float) -> float:
        return max(0.0, min(100.0, (valor - 0.30) / 0.40 * 100.0))

    return f"""
    <figure class="auc">
      <figcaption>{_html(nombre)} · AUC {_num(auc)} · IC 95 % {_num(bajo)} a {_num(alto)} · {_entero(bloque.get("n_pos"))} con la etiqueta y {_entero(bloque.get("n_neg"))} sin ella</figcaption>
      <svg viewBox="0 0 100 16" role="img" aria-label="Intervalo del AUC de {_html(nombre)}">
        <line x1="{x(0.5)}" y1="0" x2="{x(0.5)}" y2="16" class="ref"/>
        <line x1="{x(bajo)}" y1="8" x2="{x(alto)}" y2="8" class="ic"/>
        <circle cx="{x(auc)}" cy="8" r="2.2" class="punto"/>
      </svg>
      <p class="eje"><span>0,30</span><span>0,50 azar</span><span>0,70</span></p>
    </figure>"""


def _recorrido(datos: dict) -> str:
    h5, h6, h7 = datos["h5"], datos["h6"], datos["h7"]
    h8, h9, h10 = datos["h8"], datos["h9"], datos["h10"]
    multi = datos["multiverso"]
    celda = multi.get("celda_preregistrada") or {}
    pasos = [
        ("Tejido", "H5", h5.get("estado"), "Se mira si el cambio tumoral del colon es compartido con otros carcinomas, antes de hablar de sangre."),
        ("Controles", "Controles", (datos["controles"].get("veredicto_v3_4") or {}).get("estado") if isinstance(datos["controles"].get("veredicto_v3_4"), dict) else datos["controles"].get("veredicto_v3_4"), "Se comprueba que el contraste de tejido vea los marcadores que debería ver y que la mucosa pase el control de calidad."),
        ("Célula", "H6", h6.get("estado"), "Se pregunta si la concordancia entre cohortes de tejido depende del tipo celular. El sentido no está predicho."),
        ("Sangre", "Multiverso", "calculado" if celda else None, "En plasma se contrasta cáncer frente a colonoscopia negativa, y se repite el análisis en la rejilla escrita."),
        ("Predicción", "H7", h7.get("estado"), "Con lo que se sabe del tejido, antes de mirar la sangre, se intenta predecir qué proteínas tumorales cambian en sangre."),
        ("Estadio", "H8", h8.get("estado"), "En otra cohorte de plasma se pregunta si el estadio IV concentra proteínas epiteliales del tumor."),
        ("Causalidad", "H9", h9.get("estado"), "Un instrumento genético por proteína. Esta pieza es exploratoria y no decide la tesis."),
        ("Atlas", "H10", h10.get("estado"), "Se describe cuántas proteínas del atlas de sangre son específicas de cáncer colorrectal."),
    ]
    tarjetas = []
    for i, (titulo, etiqueta, estado, texto) in enumerate(pasos, start=1):
        tarjetas.append(
            f'<a class="paso" href="#r-{_html(etiqueta.lower())}"><span class="n">{i}</span>'
            f"<strong>{_html(titulo)}</strong><small>{_html(etiqueta)}</small>{_estado(estado)}<p>{_html(texto)}</p></a>"
        )
    return "".join(tarjetas)


def _fuentes() -> str:
    filas = [
        ("Tanaka, tumor frente a mucosa", "El contraste de tejido del que salen las proteínas que cambian en el colon.", "H5, H6, H7 y H8"),
        ("CPTAC y LinkedOmics, otros carcinomas", "Sirve para ver si ese cambio es de colon o aparece en otros órganos.", "H5"),
        ("PXD055218, plasma DIA-NN", "Cáncer colorrectal, colonoscopia negativa y personas sanas. Es la sangre de la tesis.", "Multiverso, aprendizaje y H7"),
        ("Álvez 2023", "El artículo compara colorrectal con otros cánceres. No es una comparación con personas sanas.", "H7"),
        ("Atlas de sangre del HPA, 2025", "Colorrectal frente a Healthy. Es el segundo recurso de sangre de H7.", "H7 y H10"),
        ("Papier 2024", "Riesgo prospectivo. Entra solo si el fichero permite etiquetar todas las proteínas.", "H7"),
        ("PXD033236, vesículas", "Una proteína cuenta como presente en vesículas de tumor si pasa el FDR escrito en el propio archivo y tiene abundancia en una columna de tumor.", "H7"),
        ("UniProt, ESM-2 y TMbed", "Secuencia, hélices y péptido señal. Entran en un modelo secundario. El principal no las usa.", "H7"),
        ("PXD064655", "Plasma con estadio. La enmienda pide agregar los precursores con directLFQ.", "H8"),
        ("GWAS Catalog GCST90255676", "El instrumento cis de cada proteína. Un solo instrumento no decide la tesis.", "H9"),
        ("Paneles de Geyer", "Eritrocitos, plaquetas y coagulación. Miden calidad de la muestra, no el tumor.", "Multiverso"),
    ]
    cuerpo = "".join(
        f"<tr><td>{_html(a)}</td><td>{_html(b)}</td><td>{_html(c)}</td></tr>" for a, b, c in filas
    )
    return f"<table><thead><tr><th>Fuente</th><th>Para qué entra</th><th>Dónde se usa</th></tr></thead><tbody>{cuerpo}</tbody></table>"


def _h5(datos: dict) -> str:
    h5 = datos["h5"]
    sub = h5.get("submuestreo") or {}
    universo = h5.get("concordancias_universo_principal") or {}
    colon = _finito(h5.get("concordancia_colon_tanaka"))
    otro = _finito(h5.get("maxima_concordancia_otro_organo"))
    brecha = None if colon is None or otro is None else colon - otro
    excluidos = h5.get("carcinomas_excluidos") or {}
    if isinstance(excluidos, dict) and excluidos:
        texto_excluido = "; ".join(f"{nombre}: {motivo}" for nombre, motivo in excluidos.items())
    else:
        texto_excluido = "ninguno"
    filas = []
    for organo, bloque in universo.items():
        if not isinstance(bloque, dict):
            continue
        filas.append(
            f"<tr><td>{_html(organo)}</td><td>{_entero(bloque.get('numerador'))}/{_entero(bloque.get('denominador'))}</td><td>{_num(bloque.get('estimacion'))}</td></tr>"
        )
    return f"""
    <article id="r-h5">
      <header><h3>H5. ¿El cambio del colon es común a otros carcinomas?</h3>{_estado(h5.get("estado"))}</header>
      <p>La predicción escrita es que una parte importante del cambio de tumor frente a normal en colon también aparece en otros carcinomas. La refutación pide, a la vez, una fracción pan-carcinoma menor que 0,20 y una concordancia con el otro órgano más alto al menos 0,15 por debajo de la de colon.</p>
      <p>De {_entero(h5.get("n_clasificadas"))} proteínas clasificadas, {_entero((h5.get("conteo") or {}).get("pan_carcinoma"))} quedan pan-carcinoma. La fracción es {_num(h5.get("fraccion_pan_carcinoma"))}. En el universo común de {_entero(h5.get("universo_principal_n"))} genes, colon tiene concordancia {_num(colon)} y el otro órgano más alto {_num(otro)}. La diferencia es {_num(brecha)}. Hace falta una diferencia de al menos 0,15, junto con la fracción bajo 0,20, para refutar. El estado que escribió la corrida es {_html(h5.get("estado"))}. {_html(h5.get("sostiene_significa"))}</p>
      <p>El submuestreo usa {_entero(sub.get("n_repeticiones"))} repeticiones y {_entero(sub.get("n_sub"))} pares. La mediana de la fracción pan-carcinoma es {_num(sub.get("mediana_fraccion_pan"))}, con percentiles 2,5 y 97,5 en {_num(sub.get("p2_5"))} y {_num(sub.get("p97_5"))}. La fracción de repeticiones al mismo lado de 0,20 es {_num(sub.get("fraccion_mismo_lado_0_20"))}. Carcinomas excluidos: {_html(texto_excluido)}.</p>
      <table><thead><tr><th>Carcinoma</th><th>Mismo signo</th><th>Concordancia</th></tr></thead><tbody>{"".join(filas)}</tbody></table>
    </article>"""


def _controles(datos: dict) -> str:
    controles = datos["controles"]
    filas = []
    for fila in controles.get("controles_fisicos") or []:
        if not isinstance(fila, dict):
            continue
        filas.append(
            f"<tr><td>{_html(fila.get('gen'))}</td><td>{_html(fila.get('estado'))}</td><td>{_entero(fila.get('n_obs'))}</td><td>{_num(fila.get('log2fc'))}</td></tr>"
        )
    veredicto = controles.get("veredicto_v3_4")
    texto = veredicto.get("estado") if isinstance(veredicto, dict) else veredicto
    return f"""
    <article id="r-controles">
      <header><h3>Controles del tejido</h3>{_estado(texto)}</header>
      <p>El bloque pasa si hay al menos un control evaluable, ningún evaluable falla y el control de calidad de la mucosa pasa en el contraste principal. El veredicto de esta regla es {_html(texto)}. La mucosa pasa el control de calidad en todos los participantes del contraste: {_html("sí" if controles.get("qc_mucosa_todos_pasan") is True else controles.get("qc_mucosa_todos_pasan"))}.</p>
      <p>CA1 baja en el tumor y también es marcador de eritrocito, así que su descenso no basta solo. AQP8 y GUCA2A no tienen fila. SLC26A3 tiene menos pares de los que la regla exige.</p>
      <table><thead><tr><th>Gen</th><th>Estado</th><th>Pares</th><th>log2FC</th></tr></thead><tbody>{"".join(filas)}</tbody></table>
    </article>"""


def _h6(datos: dict) -> str:
    h6 = datos["h6"]
    fisher = h6.get("fisher") or {}
    tabla = fisher.get("tabla_mismo_signo_luego_distinto") or []
    if isinstance(tabla, list) and len(tabla) == 2 and all(isinstance(fila, list) and len(fila) == 2 for fila in tabla):
        texto_tabla = (
            f"{_entero(tabla[0][0])} epiteliales y {_entero(tabla[1][0])} no epiteliales con el mismo signo; "
            f"{_entero(tabla[0][1])} y {_entero(tabla[1][1])} con signo distinto"
        )
    else:
        texto_tabla = _html(tabla)
    conteo = h6.get("conteo") or {}
    return f"""
    <article id="r-h6">
      <header><h3>H6. ¿La concordancia depende del tipo celular?</h3>{_estado(h6.get("estado"))}</header>
      <p>{_html(h6.get("motivo"))}. La clase sale de {_html(h6.get("fuente"))}. Pelka queda {_html(h6.get("pelka"))}: la matriz de esa referencia no está en este disco.</p>
      <p>Se clasifican {_entero(h6.get("n_genes"))} genes: {_entero(conteo.get("epitelial"))} epiteliales, {_entero(conteo.get("no_epitelial"))} no epiteliales y {_entero(conteo.get("mixta"))} mixtos. Los mixtos no entran en la prueba. La tabla queda en {texto_tabla}. La razón de odds no es finita cuando una de esas casillas es cero. El estado es {_html(h6.get("estado"))}.</p>
    </article>"""


def _h7(datos: dict) -> str:
    h7 = datos["h7"]
    por = h7.get("por_recurso") or {}
    barras = "".join(
        _auc(nombre, bloque) for nombre, bloque in por.items() if isinstance(bloque, dict) and bloque.get("auc") is not None
    )
    otros = []
    for nombre, bloque in por.items():
        if isinstance(bloque, dict) and bloque.get("auc") is None:
            otros.append(f"<li>{_html(nombre)}: {_html(bloque.get('estado'))}, {_entero(bloque.get('n_pos'))} con etiqueta y {_entero(bloque.get('n_neg'))} sin ella.</li>")
    kappas = []
    for fila in h7.get("concordancia_de_etiquetas") or []:
        if isinstance(fila, dict):
            kappas.append(f"<li>{_html(fila.get('a'))} con {_html(fila.get('b'))}: kappa {_num(fila.get('kappa'))} en {_entero(fila.get('n'))} proteínas.</li>")
    constructos = h7.get("constructo_por_recurso") or {}
    lista_c = "".join(f"<li>{_html(nombre)}: {_html(texto)}</li>" for nombre, texto in constructos.items())
    variables = _variables(h7.get("variables_usadas"))
    sec = datos["secuencia"]
    indep = h7.get("independientes") or {}
    return f"""
    <article id="r-h7">
      <header><h3>H7. ¿Lo que se sabe del tejido predice el cambio en sangre?</h3>{_estado(h7.get("estado"))}</header>
      <p>Esta es la pregunta central. El modelo principal usa solo variables anteriores a la sangre: {variables}. El modelo está completo: {_si(h7.get("modelo_completo"))}. Faltantes: {_html(h7.get("variables_faltantes") or "ninguna")}.</p>
      <p>La regla de refutación es que el intervalo del AUC incluya 0,50 en todos los recursos. La de sostén pide el límite inferior por encima de 0,55 en al menos dos. La línea vertical de cada figura es el azar.</p>
      {barras}
      <ul>{"".join(otros)}</ul>
      <p>Los tres recursos no miden lo mismo.</p>
      <ul>{lista_c}</ul>
      <p>{_html(h7.get("nota_alvez_hpa"))} {"Los casos están comprobados como independientes." if indep.get("comprobado") is True else "Los casos no están comprobados como independientes."} {_oracion(indep.get("motivo"))}</p>
      <ul>{"".join(kappas)}</ul>
      <p>El modelo secundario, en {_html(sec.get("recurso"))}, pasa de AUC {_num(sec.get("auc_principal"))} a {_num(sec.get("auc_secundario"), 4)}. El incremento es {_num(sec.get("incremento"), 4)}. {_html(sec.get("nota"))} Hay {_entero(sec.get("n_con_secuencia"))} genes con secuencia y {_entero(sec.get("n_sin_secuencia"))} sin ella. {_html(h7.get("nota_tmbed"))}</p>
      <p>La fase de Pepe de este resultado es {_html(h7.get("fase_pepe"))}.</p>
    </article>"""


def _resto(datos: dict) -> str:
    h8, h9, h10 = datos["h8"], datos["h9"], datos["h10"]
    multi = datos["multiverso"]
    celda = multi.get("celda_preregistrada") or {}
    r1 = multi.get("r1") or {}
    pre = datos["preanalitica"]
    tri = datos["triangulacion"]
    proteinas = []
    for fila in h9.get("proteinas") or []:
        if not isinstance(fila, dict):
            continue
        proteinas.append(
            f"<tr><td>{_html(fila.get('proteina'))}</td><td>{_html(fila.get('estado'))}</td><td>{_html(fila.get('rsid') or '—')}</td><td>{_num(fila.get('beta'), 4)}</td><td>{_num(fila.get('p'))}</td></tr>"
        )
    conteo = h10.get("conteo") or {}
    if tri.get("sin_pruebas_nuevas") is True:
        frase_tri = "No añade una prueba nueva."
    elif tri.get("sin_pruebas_nuevas") is False:
        frase_tri = "El JSON marca pruebas nuevas en esta pieza."
    else:
        frase_tri = ""
    return f"""
    <article id="r-multiverso">
      <header><h3>Sangre. La rejilla y la celda escrita de antemano</h3>{_estado("calculado")}</header>
      <p>El suero se reanaliza en {_entero(multi.get("n_celdas_suero"))} celdas. El aprendizaje, en {_entero(multi.get("n_celdas_ml"))}. Hay {_entero(multi.get("n_r1"))} normalizaciones del tejido; la concordancia de la sensibilidad pasa de {_num(r1.get("antes"), 1)} a {_num(r1.get("despues"), 1)}.</p>
      <p>La celda preregistrada usa {_html(celda.get("cuantificacion"))}, contaminantes {_html(celda.get("contaminantes"))}, paneles {_html(celda.get("paneles"))} y presencia {_html(celda.get("presencia"))}. Sobre {_entero(celda.get("n"))} proteínas hay {_entero(celda.get("diferente"))} diferentes, {_entero(celda.get("equivalente"))} equivalentes y {_entero(celda.get("indeterminada"))} indeterminadas entre cáncer y colonoscopia negativa.</p>
    </article>
    <article id="r-h8">
      <header><h3>H8. ¿El estadio IV concentra proteínas del tumor?</h3>{_estado(h8.get("estado"))}</header>
      <p>{_html(h8.get("motivo"))} La corrida unió {_entero(h8.get("n_muestras_unidas"))} muestras al metadato. Mapa uno a uno: {_si(h8.get("mapa_uno_a_uno"))}.</p>
    </article>
    <article id="r-h9">
      <header><h3>H9. Dirección con un instrumento genético</h3>{_estado(h9.get("estado"))}</header>
      <p>{_html(h9.get("motivo"))} El estudio es {_html(h9.get("gwas"))}. Proteínas exploradas: {_entero(h9.get("n_exploratorios"))}.</p>
      <table><thead><tr><th>Proteína</th><th>Estado</th><th>Instrumento</th><th>Beta</th><th>p</th></tr></thead><tbody>{"".join(proteinas)}</tbody></table>
    </article>
    <article id="r-h10">
      <header><h3>H10. Especificidad en el atlas de sangre</h3>{_estado(h10.get("estado"))}</header>
      <p>{_html(h10.get("motivo_estado"))} Proteínas específicas de cáncer colorrectal: {_entero(h10.get("n_especificas_crc"))}. El resto del conteo: {_conteo(conteo)}. Fase de Pepe: {_html(h10.get("fase_pepe"))}.</p>
    </article>
    <article id="r-pre">
      <header><h3>Orden de adquisición y triangulación</h3>{_estado(pre.get("estado"))}</header>
      <p>La preanalítica es exploratoria y no cambia el titular. {_oracion(pre.get("motivo"))} Crudos contados por la corrida: {_entero(pre.get("n_crudos"))}.</p>
      <p>La triangulación es {_html(str(tri.get("tipo") or "").replace("_", " "))}, con {_entero(tri.get("n_genes"))} genes. {_html(tri.get("nota"))} {frase_tri}</p>
    </article>"""


def _galeria(datos: dict) -> str:
    tarjetas = []
    for figura in datos.get("figuras") or []:
        ident = _html(figura["id"])
        tarjetas.append(
            f'<button class="figura" type="button" data-src="figuras/{ident}.png" data-titulo="{_html(figura["titulo"])}">'
            f'<img src="figuras/{ident}.png" alt="{_html(figura["titulo"])}">'
            f"<span>{ident}. {_html(figura['titulo'])}</span></button>"
        )
    if not tarjetas:
        return "<p>Esta carpeta no trae las figuras de la corrida v3.3.</p>"
    return "".join(tarjetas)


def _certificado(datos: dict) -> str:
    filas = []
    for casilla in datos["certificado"].get("casillas") or []:
        filas.append(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
                _html(casilla.get("id")), _estado(casilla.get("estado")), _html(casilla.get("motivo")), _html(casilla.get("fuente"))
            )
        )
    return "".join(filas)


def pagina(datos: dict) -> str:
    certificado = datos["certificado"]
    sello = datos["sello"]
    commit = str(sello.get("commit") or "")
    corto = commit[:7]
    h7 = datos["h7"]
    celda = (datos["multiverso"].get("celda_preregistrada") or {})
    return f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>OncoGO-CRC · la tesis y la corrida v3.4</title>
<style>
:root {{
  color-scheme: light;
  --ink:#1c1915; --muted:#5c564c; --paper:#f6f1e8; --card:#fffdf8; --line:#e2d8c8;
  --teal:#0f3d38; --teal-2:#1f6f64; --gold:#a15c07; --good:#1e6b45; --bad:#8c3a2f; --soft:#24577a;
}}
* {{ box-sizing:border-box; }}
body {{ margin:0; font-family:"Avenir Next","Segoe UI",sans-serif; background:var(--paper); color:var(--ink); line-height:1.55; }}
a {{ color:var(--teal-2); }}
.aviso {{ background:var(--teal); color:#f6f1e8; padding:.75rem 1.25rem; }}
.aviso p {{ max-width:74rem; margin:0 auto; }}
header.porta {{ max-width:74rem; margin:0 auto; padding:2.2rem 1.25rem 1rem; }}
header.porta h1 {{ font-family:Palatino,"Iowan Old Style",Georgia,serif; font-weight:600; font-size:clamp(2rem,4vw,3.3rem); line-height:1.05; margin:0 0 .6rem; }}
header.porta .lede {{ font-size:1.15rem; max-width:42rem; }}
nav {{ position:sticky; top:0; z-index:2; background:rgba(246,241,232,.94); border-bottom:1px solid var(--line); }}
nav ul {{ max-width:74rem; margin:0 auto; padding:.4rem 1.25rem; display:flex; gap:.9rem; list-style:none; overflow:auto; }}
nav a {{ color:var(--ink); text-decoration:none; font-size:.92rem; }}
main {{ max-width:74rem; margin:0 auto; padding:1rem 1.25rem 3rem; }}
section {{ margin:2.2rem 0; }}
h2 {{ font-family:Palatino,"Iowan Old Style",Georgia,serif; font-size:1.8rem; margin-bottom:.4rem; }}
.resumen {{ display:grid; grid-template-columns:repeat(3,1fr); gap:.8rem; }}
.resumen article, article.bloque, .paso {{ background:var(--card); border:1px solid var(--line); border-radius:14px; padding:1rem 1.05rem; }}
.resumen h3 {{ margin:.2rem 0; font-size:1rem; }}
.recorrido {{ display:grid; grid-template-columns:repeat(4,1fr); gap:.7rem; }}
.paso {{ color:inherit; text-decoration:none; display:block; }}
.paso .n {{ display:inline-grid; place-items:center; width:1.6rem; height:1.6rem; border-radius:99px; background:var(--teal); color:#fff; font-size:.85rem; }}
.paso p {{ color:var(--muted); font-size:.92rem; margin:.45rem 0 0; }}
article.bloque {{ margin:.8rem 0; }}
article.bloque header, .resumen header {{ display:flex; justify-content:space-between; gap:.6rem; align-items:flex-start; }}
article.bloque h3 {{ margin:0; font-size:1.15rem; }}
.estado {{ border-radius:99px; padding:.1rem .55rem; font-size:.78rem; letter-spacing:.02em; background:#eee; white-space:nowrap; }}
.e-sostiene, .e-verde, .e-pasa, .e-calculado {{ background:#e5f4ea; color:var(--good); }}
.e-refuta, .e-rojo {{ background:#f8e6e2; color:var(--bad); }}
.e-no_sostiene, .e-limitacion_declarada, .e-no_evaluable, .e-no_disponible {{ background:#f8efd8; color:var(--gold); }}
.e-exploratorio, .e-descriptiva {{ background:#e6f0f6; color:var(--soft); }}
table {{ width:100%; border-collapse:collapse; font-size:.92rem; margin:.6rem 0; }}
th, td {{ text-align:left; padding:.35rem .3rem; border-bottom:1px solid var(--line); vertical-align:top; }}
.auc svg {{ width:100%; height:2.4rem; }}
.auc .ref {{ stroke:#8c3a2f; stroke-width:.4; }}
.auc .ic {{ stroke:var(--teal); stroke-width:1.1; }}
.auc .punto {{ fill:var(--teal); }}
.eje {{ display:flex; justify-content:space-between; color:var(--muted); font-size:.78rem; margin:0; }}
.galeria {{ display:grid; grid-template-columns:repeat(3,1fr); gap:.7rem; }}
.figura {{ background:var(--card); border:1px solid var(--line); border-radius:12px; padding:.4rem; text-align:left; cursor:zoom-in; color:inherit; }}
.figura img {{ width:100%; height:9rem; object-fit:contain; background:#fff; }}
.figura span {{ display:block; font-size:.82rem; padding:.35rem; }}
dialog {{ width:min(920px, 94vw); border:0; border-radius:16px; padding:1rem; }}
dialog img {{ width:100%; }}
dialog::backdrop {{ background:rgba(28,25,21,.55); }}
footer {{ color:var(--muted); font-size:.9rem; }}
@media (max-width: 900px) {{
  .resumen, .recorrido, .galeria {{ grid-template-columns:1fr; }}
}}
</style>
</head>
<body>
<div class="aviso"><p>{_html(certificado.get("advertencia"))} Fase de Pepe {_html(h7.get("fase_pepe"))}. El certificado de esta corrida es {_html(certificado.get("estado"))}.</p></div>
<header class="porta">
  <p>OncoGO-CRC · corrida v3.4 · commit {_html(corto)} · {_html(sello.get("fecha_utc"))}</p>
  <h1>Qué pregunta la tesis, en qué orden se responde y qué salió.</h1>
  <p class="lede">Se busca saber si las proteínas que cambian en el tumor de colon dejan una señal reconocible en la sangre, y si esa señal se puede anticipar con lo que ya se sabe del tejido. Cada número de esta página sale del JSON de la corrida firmada.</p>
</header>
<nav><ul>
  <li><a href="#pregunta">Pregunta</a></li>
  <li><a href="#recorrido">Recorrido</a></li>
  <li><a href="#fuentes">Fuentes</a></li>
  <li><a href="#resultados">Resultados</a></li>
  <li><a href="#figuras">Figuras</a></li>
  <li><a href="#conclusion">Conclusión</a></li>
  <li><a href="#certificado">Certificado</a></li>
</ul></nav>
<main>
<section id="pregunta">
  <h2>La pregunta</h2>
  <div class="resumen">
    <article><h3>Lo que se busca</h3><p>Una regla, escrita antes de mirar estos resultados, que diga si el conocimiento del tejido predice qué proteínas tumorales cambian en sangre.</p></article>
    <article><h3>Lo que hizo esta corrida</h3><p>H7 queda {_html(h7.get("estado"))}. En la celda de suero escrita de antemano hay {_entero(celda.get("diferente"))} proteínas diferentes entre cáncer y colonoscopia negativa.</p></article>
    <article><h3>Cómo se lee un estado</h3><p>Sostiene y refuta son los dos finales de una regla fija. Exploratorio y descriptivo no deciden. No evaluable significa que el dato pedido no se pudo calcular.</p></article>
  </div>
</section>
<section id="recorrido">
  <h2>El orden del trabajo</h2>
  <p>Primero el tejido, después la calidad de ese contraste, después la sangre y, al final, la predicción. La sangre no se usa para elegir las proteínas del tejido.</p>
  <div class="recorrido">{_recorrido(datos)}</div>
</section>
<section id="fuentes">
  <h2>De dónde sale cada dato</h2>
  <p>Cada fuente entra por una razón escrita en el preregistro. Una cohorte no se fusiona con otra: no hay personas compartidas que autoricen esa fusión.</p>
  {_fuentes()}
</section>
<section id="resultados">
  <h2>Qué se hizo y qué se concluye en cada paso</h2>
  {_h5(datos)}
  {_controles(datos)}
  {_h6(datos)}
  {_h7(datos)}
  {_resto(datos)}
</section>
<section id="figuras">
  <h2>Figuras de la corrida v3.3</h2>
  <p>Estas catorce figuras pertenecen a la corrida anterior, con su sello de v3.3. Sirven para ver el camino. La lectura vigente es la de arriba: en v3.4, H5 queda sostiene porque cambió la regla del universo, y H7 sigue en refuta. El título de una figura de v3.3 no sustituye el JSON de v3.4.</p>
  <div class="galeria">{_galeria(datos)}</div>
</section>
<section id="conclusion">
  <h2>La conclusión de esta corrida</h2>
  <article class="bloque">
    <p>Con las variables del tejido que el preregistro fijó, el modelo no predice qué proteínas tumorales cambian en sangre: los intervalos de Álvez 2023 y del atlas HPA incluyen el azar, y Papier 2024 no se pudo evaluar. En el plasma de PXD055218, la celda escrita de antemano no encuentra proteínas diferentes entre cáncer y colonoscopia negativa.</p>
    <p>H5 queda sostiene porque la refutación de especificidad de órgano no se cumple entera. H6 no muestra una diferencia entre proteínas epiteliales y no epiteliales en la prueba bilateral. H8 y el orden de adquisición quedan sin el cálculo que la enmienda pidió, y la página lo dice con el motivo de cada JSON.</p>
    <p>{_html(certificado.get("advertencia"))}</p>
  </article>
</section>
<section id="certificado">
  <h2>Certificado para publicar</h2>
  <p>El estado publicable exige que ninguna casilla esté en rojo. Commit completo {_html(commit)}.</p>
  <table><thead><tr><th>Casilla</th><th>Estado</th><th>Qué se comprobó</th><th>Archivo</th></tr></thead><tbody>{_certificado(datos)}</tbody></table>
</section>
<footer>
  <p><a href="v3_3.html">Página de la corrida v3.3.</a> <a href="datos/v3_4.json">JSON que alimenta esta página.</a> Lo escribe <code>scripts/v3/exportar_sitio_v3_4.py</code> desde <code>results/v3_4/</code>.</p>
</footer>
</main>
<dialog id="visor"><img alt=""><p></p><form method="dialog"><button>Cerrar</button></form></dialog>
<script>
const visor = document.getElementById("visor");
document.querySelectorAll(".figura").forEach((boton) => {{
  boton.addEventListener("click", () => {{
    visor.querySelector("img").src = boton.dataset.src;
    visor.querySelector("img").alt = boton.dataset.titulo;
    visor.querySelector("p").textContent = boton.dataset.titulo;
    visor.showModal();
  }});
}});
</script>
</body>
</html>
"""


def escribir(base: Path, sitio: Path = SITIO) -> None:
    datos = datos_publicos(base, sitio)
    sitio.mkdir(parents=True, exist_ok=True)
    indice = sitio / "index.html"
    historico = sitio / "v3_3.html"
    if indice.is_file() and not historico.is_file() and "corrida v3.3" in indice.read_text(encoding="utf-8"):
        historico.write_text(indice.read_text(encoding="utf-8"), encoding="utf-8")
    (sitio / "datos").mkdir(parents=True, exist_ok=True)
    (sitio / "datos" / "v3_4.json").write_text(json.dumps(datos, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    indice.write_text(pagina(datos), encoding="utf-8")


def main() -> None:
    base = Path(sys.argv[1]) if len(sys.argv) > 1 else RESULTADOS
    try:
        escribir(base)
    except (RuntimeError, FileNotFoundError) as error:
        print(error)
        sys.exit(1)
    print(f"[v3.4] sitio escrito desde {base}")


if __name__ == "__main__":
    main()
