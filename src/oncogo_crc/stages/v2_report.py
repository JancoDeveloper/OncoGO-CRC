"""Informe de v2: figuras, Markdown, HTML y PDF, todo desde las tablas.

Lo que esta etapa no hace. No recalcula nada. Lee los JSON de auditoría y los
TSV que escribieron las etapas anteriores, y de ahí salen todas las cifras. Si
una etapa no ha corrido, su sección dice que no ha corrido en lugar de repetir el
número de v1.

Por qué importa. En v1 los números del PDF se escribían a mano. Cuando el
pipeline cambiaba, el texto se quedaba viejo y nadie lo notaba. Aquí el texto se
construye con los valores que acaban de leerse, y
``scripts/v2/verificar_numeros.py`` falla si alguno no aparece en una tabla.

El idioma es español, con frases cortas y voz activa, y está escrito para alguien
que no conoce el tema. El cierre obligatorio es el mismo en los tres formatos:
ningún resultado de esta tesis sirve para decidir sobre un paciente.
"""

from __future__ import annotations

import html
import json
import shutil
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.v2 import figuras as F

CIERRE = "Ningún resultado de esta tesis sirve para decidir sobre un paciente."


def _leer_json(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError):
        return None


def _leer_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except (pd.errors.ParserError, UnicodeDecodeError):
        return pd.DataFrame()


def _coma(valor, digitos: int = 3) -> str:
    """Número en formato español, con signo menos de verdad."""
    if valor is None or (isinstance(valor, float) and valor != valor):
        return "no disponible"
    return f"{float(valor):.{digitos}f}".replace(".", ",").replace("-", "\u2212")


def _pct(valor, digitos: int = 1) -> str:
    if valor is None or (isinstance(valor, float) and valor != valor):
        return "no disponible"
    return f"{float(valor) * 100:.{digitos}f}".replace(".", ",") + " %"


def _ic(d: dict, digitos: int = 3) -> str:
    low, high = d.get("ci_low"), d.get("ci_high")
    if low is None or high is None:
        return "sin intervalo"
    return f"IC 95 % {_coma(low, digitos)}–{_coma(high, digitos)}"


def _frac(d: dict) -> str:
    return f"{int(d['numerador'])} de {int(d['denominador'])}"


@dataclass
class Seccion:
    """Una sección del informe, con su figura opcional."""

    numero: int
    titulo: str
    cuerpo: list[str]
    figura: F.Figura | None = None


# --------------------------------------------------------------------- secciones
def _sec_que_es_esto(cert: dict | None) -> Seccion:
    cuerpo = [
        "Esta tesis busca proteínas que distingan un tumor de colon de tejido "
        "normal del mismo paciente, y después comprueba si esas mismas proteínas "
        "siguen distinguiéndolos en otras cohortes y si alguna se puede medir en "
        "sangre.",
        "La versión 1 del análisis llegó a unas conclusiones. Esta versión 2 las "
        "vuelve a comprobar con el método corregido. No sustituye a v1: v1 queda "
        "congelada, con su manifiesto de hashes, y lo que cambia se documenta "
        "número por número.",
        "Lo importante de leer un informe así es distinguir tres cosas. Un "
        "resultado que se sostiene. Un resultado que se cae al corregir el método. "
        "Y un resultado que no se puede evaluar con los datos que hay. Las tres "
        "aparecen aquí, y la tercera no es un fracaso: es la respuesta honesta "
        "cuando el dato no alcanza.",
    ]
    if cert:
        cuerpo.append(
            "El certificado de esta corrida usa alfa "
            f"{_coma(cert['alpha'], 2)}, un umbral de magnitud de "
            f"{_coma(cert['treat_lfc'], 3)} en log₂ y la semilla base "
            f"{int(cert['semilla_base'])}. Todo análisis con azar deriva su semilla "
            "de esa y queda escrito en la tabla correspondiente."
        )
    return Seccion(1, "Qué se busca y cómo leer esto", cuerpo)


def _sec_gradiente(axes: dict | None, fig: F.Figura | None) -> Seccion:
    if not axes:
        return Seccion(2, "El gradiente normal–adenoma–carcinoma",
                       ["La etapa de ejes no ha corrido, así que esta sección no "
                        "tiene resultado."])
    c5 = axes["c5_adenoma_carcinoma"]
    comp = c5["veredicto"]["comparacion"]
    indep = c5["seleccion_independiente_cptac"]["intermedio_dado_mismo_signo"]
    nulo = c5["nulo_que_reselecciona"]
    direccion = c5["seleccion_independiente_cptac"]["direccion"]
    cuerpo = [
        "La idea es intuitiva. Un pólipo es un paso intermedio entre el tejido "
        "normal y el tumor, así que sus proteínas deberían quedar a medio camino. "
        "v1 dijo que eso pasaba en el "
        f"{_pct(comp['mismo_brazo_v1'])} de los casos.",
        "El problema es cómo se eligieron las proteínas. v1 las eligió por ser las "
        "más alteradas en el brazo del tumor, y después midió el gradiente en ese "
        "mismo brazo. Elegir y medir en el mismo sitio garantiza encontrar lo que "
        "se busca, porque una proteína entra en la lista en parte por ruido y ese "
        "ruido vuelve a aparecer cuando se la mide.",
        "v2 lo hace de cuatro formas distintas para separar el efecto real del "
        "artefacto de selección. Con selección independiente, que es el análisis "
        "principal preregistrado, las proteínas se eligen en CPTAC y se miden en la "
        f"réplica: el gradiente baja a {_coma(indep['valor'])} "
        f"({_frac(indep)} proteínas, {_ic(indep)}), que no se distingue del azar.",
        "El nulo lo confirma por el otro lado. Si se sortea al azar qué brazo hace "
        "de tumor y se vuelve a seleccionar, el gradiente que sale solo por "
        f"regresión a la media ya es {_coma(nulo['nulo_mediana'])}. Lo observado "
        f"en el mismo brazo, {_coma(nulo['observado_intermedio'])}, está por encima "
        f"de eso (p = {_coma(nulo['p_empirico'], 6)}), pero la mayor parte del "
        "número de v1 era el artefacto, no el gradiente.",
        "Lo que sí sobrevive es la **dirección**: de las proteínas que la réplica "
        f"alcanza a medir, {_frac(direccion)} suben o bajan en el pólipo en el "
        f"mismo sentido que en el tumor ({_coma(direccion['valor'])}, "
        f"{_ic(direccion)}). Lo que no sobrevive es la magnitud ordenada. El pólipo "
        "se parece al tumor en el sentido del cambio, no en quedar a medio camino.",
        "Esto obligó a una enmienda. El nulo preregistrado estaba mal especificado: "
        "intercambiaba las etiquetas por paciente, lo que destruye justo la ventaja "
        "de selección que se quería medir, y daba un resultado significativo con "
        "ruido puro. Queda escrito en `docs/enmiendas_v2.md` como enmienda E7, con "
        "la fecha y el motivo. El preregistro no se edita.",
    ]
    return Seccion(2, "El gradiente normal–adenoma–carcinoma no se sostiene", cuerpo, fig)


def _sec_modelos(axes: dict | None, fig: F.Figura | None) -> Seccion:
    if not axes:
        return Seccion(3, "Comparar modelos", ["La etapa de ejes no ha corrido."])
    c6 = axes["c6_seleccion_de_modelo"]
    pac = c6["bootstrap_por_paciente"]
    ver = c6["veredicto"]
    disc = c6["discrepancia_entre_criterios"]
    mejor = pac["por_modelo"][ver["mejor"]]
    cuerpo = [
        "Aquí se entrena un clasificador en una cohorte y se prueba en otra. La "
        f"validación externa tiene {int(c6['n_validacion'])} muestras de "
        f"{int(c6['n_pacientes_validacion'])} pacientes, que es el detalle que lo "
        "cambia todo: cada paciente aporta su tumor y su normal, así que hay "
        f"{int(c6['n_pacientes_validacion'])} personas independientes y no "
        f"{int(c6['n_validacion'])}.",
        "v1 remuestreaba especímenes. Eso trata el tumor y el normal del mismo "
        "paciente como dos observaciones independientes, y no lo son. v2 remuestrea "
        "pacientes: cada uno entra con sus dos muestras o con ninguna.",
        "Conviene decir algo que va contra la intuición. Agrupar casi siempre "
        "ensancha el intervalo, pero aquí lo estrecha, porque el diseño es "
        "emparejado y el efecto del paciente se cancela dentro del par. La razón "
        "para agrupar no es obtener un intervalo más ancho: es remuestrear la "
        "unidad que de verdad se repetiría si se repitiera el estudio.",
        f"El mejor modelo es {ver['mejor'].replace('_', ' ')}, con AUC "
        f"{_coma(mejor['mediana'])} ({_ic(mejor)}). Pero la regla de decisión "
        "preregistrada dice que solo hay ganador si su intervalo no se solapa con "
        "el del segundo, y se solapan. Así que no hay ganador declarado.",
        f"Hay una discrepancia que merece quedar escrita. {disc['lectura']}",
        "Son dos preguntas distintas, y la diferencia por pares es la más potente "
        "porque aprovecha que los modelos se evalúan en las mismas muestras. Pero "
        "cambiar el criterio después de ver el resultado sería elegir el que da la "
        "respuesta que gusta. Se reporta la discrepancia y se mantiene la regla.",
    ]
    empatados = [
        f"{nombre.replace('_', ' ')} ({int(d['n_hiperparametros_empatados'])} valores)"
        for nombre, d in c6["seleccion_interna_agrupada"].items()
        if d.get("eleccion_arbitraria")
    ]
    if empatados:
        cuerpo.append(
            "Y una limitación del propio ajuste. Con tan pocos pacientes el bucle "
            "interno de validación no distingue entre hiperparámetros: empatan en "
            f"los tres modelos ({', '.join(empatados)}). Se toma el primero de la "
            "rejilla, así que esa elección no está respaldada por el dato y queda "
            "declarada como arbitraria en la tabla."
        )
    return Seccion(3, "Comparar modelos sin repetir pacientes", cuerpo, fig)


def _sec_vias(vias: dict | None, fig: F.Figura | None) -> Seccion:
    if not vias:
        return Seccion(4, "La cláusula de vías", ["La etapa de vías no ha corrido."])
    cl = vias["clausula_via"]
    rob = cl["robustez"]
    v1 = vias["comparacion_con_v1"]
    reac = vias["por_coleccion"]["reactome"]["concordancia_entre_cohortes"]
    uni = vias["universo"]
    cuerpo = [
        "Una vía es un conjunto de proteínas que trabajan juntas. La pregunta es si "
        "las vías alteradas en una cohorte también lo están en otra, lo que sería "
        "más robusto que mirar proteínas una a una.",
        f"v1 dio esta cláusula por satisfecha con una concordancia de "
        f"{_coma(v1['v1_concordancia_via'], 1)}, que suena perfecta. Detrás había "
        f"{int(v1['v1_n_terminos'])} términos y un intervalo que bajaba hasta "
        f"{_coma(v1['v1_ci_low'])}. Con cuatro términos un intervalo así no "
        "distingue una vía real de una moneda al aire.",
        f"v2 exige {int(vias['n_min_terminos'])} términos como mínimo, prueba tres "
        "colecciones independientes y usa una prueba competitiva, que compara el "
        "conjunto contra el resto de proteínas en lugar de contra nada. La prueba "
        "competitiva además infla la varianza según la correlación entre genes del "
        "conjunto: dejarla en cero haría la prueba anticonservadora.",
        f"De las tres colecciones solo Reactome llega al mínimo, con "
        f"{int(reac['n_terminos_evaluados'])} términos evaluables de "
        f"{int(reac['n_terminos_significativos_descubrimiento'])} significativos en "
        f"el descubrimiento. Su concordancia es {_coma(reac['concordancia']['valor'])} "
        f"({_frac(reac['concordancia'])}, {_ic(reac['concordancia'])}) y supera el "
        f"umbral de {_coma(vias['umbral_concordancia_via'], 2)} por "
        f"{_coma(reac['margen_del_ic_inferior_sobre_el_umbral'])}.",
        f"Eso es poco margen y hay que decirlo: {rob['advertencia']}",
        f"El cuello de botella no es el resultado, es la cobertura. La réplica mapea "
        f"a símbolo de gen solo {_pct(uni['fraccion_de_la_replica_con_simbolo_de_gen'])} "
        f"de sus proteínas, {int(uni['n_genes_replica'])} de "
        f"{int(uni['n_proteinas_en_la_tabla_de_la_replica'])}. Con tan pocos genes "
        "casi ningún conjunto llega al mínimo de términos.",
    ]
    return Seccion(4, "La cláusula de vías, con un mínimo de términos", cuerpo, fig)


def _sec_deteccion(det: dict | None, fig: F.Figura | None) -> Seccion:
    if not det:
        return Seccion(5, "La brecha tejido–plasma", ["La etapa de detección no ha corrido."])
    curva = det["curva_de_deteccion"]
    conc = det["concordancia"]
    ver = conc["veredicto"]
    fuente = det["fuente_concentraciones"]
    cuerpo = [
        "Una proteína alterada en el tumor sería útil si se pudiera medir en un "
        "análisis de sangre. En los datos no es así: la dirección del cambio en "
        "tejido y en plasma concuerda en "
        f"{_coma(conc['todas']['valor'])} de los casos "
        f"({_frac(conc['todas'])}), que es azar.",
        "Antes de concluir algo hay que descartar una explicación aburrida. Puede "
        "que la concordancia falle no porque la biología sea distinta, sino porque "
        "el experimento de suero no llega a medir las proteínas poco abundantes. Eso "
        "sería un artefacto de detección, no un resultado.",
        f"Para comprobarlo se usan las concentraciones en sangre del Human Protein "
        f"Atlas ({int(fuente['n_con_concentracion_y_uniprot'])} proteínas con "
        "concentración y accesión), prefiriendo las medidas por espectrometría de "
        "masas porque es la tecnología del propio experimento de suero.",
        "La ventana de detección existe y es fuerte: por cada orden de magnitud de "
        f"concentración la probabilidad de cuantificarse sube con pendiente "
        f"{_coma(curva['pendiente_por_orden_de_magnitud'])} en escala logit "
        f"(z = {_coma(curva['z_pendiente'], 1)}). El umbral donde la probabilidad "
        f"pasa de la mitad está en {_coma(curva['umbral_log10_pg_l'], 2)} en log₁₀ "
        "de pg/L.",
        f"Y sin embargo la ventana no explica la brecha. Dentro de la ventana la "
        f"concordancia es {_coma(conc['dentro_de_la_ventana']['valor'])} "
        f"({_frac(conc['dentro_de_la_ventana'])}, "
        f"{_ic(conc['dentro_de_la_ventana'])}), que sigue siendo compatible con el "
        f"azar. {ver['lectura']}",
        "Esto es una afirmación más fuerte que la de v1 y mejor acotada. v1 decía "
        "que no había traslación. v2 dice que no la hay incluso restringiendo a las "
        "proteínas que el experimento sí podía ver, lo que descarta la explicación "
        "técnica y deja la biológica: lo que cambia en el tumor no es lo que llega "
        "medible a la sangre.",
    ]
    return Seccion(5, "La brecha tejido–plasma no es un límite de detección", cuerpo, fig)


def _sec_meta(meta: dict | None, fig: F.Figura | None) -> Seccion:
    if not meta:
        return Seccion(6, "Heterogeneidad entre cohortes", ["La etapa de meta no ha corrido."])
    diag = meta["diagnostico_de_tecnologia"]
    cuerpo = [
        "R3 mide si las cohortes concuerdan en la **magnitud** del efecto, no solo "
        "en su dirección. En v1 era no calculable: el preregistro pide al menos "
        f"{int(meta['min_cohortes_para_r3'])} cohortes de tejido y solo había dos.",
        f"v2 incorpora una tercera y ya se puede calcular. Sale "
        f"{_coma(meta['mediana_1_menos_i2'])} sobre "
        f"{int(meta['n_proteinas_evaluables'])} proteínas presentes en las tres, de "
        f"{int(meta['n_proteinas_con_al_menos_dos_cohortes'])} que están en al menos "
        "dos.",
        f"Ese número conviene leerlo del revés. R3 es uno menos I², así que un R3 de "
        f"{_coma(meta['mediana_1_menos_i2'])} significa un I² mediano de "
        f"{_coma(meta['mediana_i2'])}. Es heterogeneidad altísima: las tres cohortes "
        "casi no concuerdan en cuánto cambia cada proteína, aunque sí concuerden en "
        "si sube o baja.",
        "La explicación fácil sería la tecnología. CPTAC usa marcaje isobárico y las "
        "otras dos usan cuantificación sin marcaje, y comparar razones logarítmicas "
        "con intensidades podría generar heterogeneidad por sí solo. Se comprobó "
        "repitiendo el cálculo en las mismas proteínas sin CPTAC: el I² mediano baja "
        f"de {_coma(diag['i2_mediano_con_las_tres'])} a "
        f"{_coma(diag['i2_mediano_solo_lfq'])}, es decir "
        f"{_coma(diag['caida_al_quitar_cptac'])}. Apenas cambia.",
        f"{diag['advertencia']} {diag['nota']}",
        f"De las {int(meta['n_proteinas_evaluables'])} proteínas evaluables, "
        f"{int(meta['n_con_efecto_combinado_distinto_de_cero'])} tienen un efecto "
        "combinado cuyo intervalo excluye el cero. Hay señal compartida; lo que no "
        "hay es acuerdo en su tamaño.",
    ]
    return Seccion(6, "Las cohortes concuerdan en dirección, no en magnitud", cuerpo, fig)


def _sec_independencia(meta: dict | None, fig: F.Figura | None) -> Seccion:
    if not meta:
        return Seccion(7, "Independencia de las cohortes", ["La etapa de meta no ha corrido."])
    ind = meta["independencia_de_cohortes"]
    cuerpo = [
        "Añadir una cohorte solo vale si de verdad es otra cohorte. Si dos depósitos "
        "públicos contienen los mismos pacientes, sumarlos infla la confianza sin "
        "añadir información, y eso es peor que no sumar nada.",
        "La comprobación tiene dos niveles. El primero compara los checksums de los "
        "ficheros en PRIDE. El segundo correlaciona cada muestra de un depósito con "
        "cada muestra del otro sobre las proteínas comunes. El segundo es el que "
        "decide, porque dos depósitos pueden ser ficheros distintos del mismo "
        "paciente y entonces los hashes no coinciden pero los datos sí.",
        f"Resultado: ningún par de muestras llega al umbral de sospecha de "
        f"{_coma(ind['umbral_de_sospecha'], 2)}. El máximo es "
        f"{_coma(ind['r_maximo_entre_muestras'])} y la mediana "
        f"{_coma(ind['r_mediano_entre_muestras'])}, una diferencia pequeña que es "
        "justo lo que se espera de cohortes distintas del mismo tejido.",
        "El primer nivel sí sirvió para descartar algo. PXD031705 "
        f"{ind['pxd031705_excluida'][0].lower()}{ind['pxd031705_excluida'][1:]}",
        f"Veredicto: {ind['veredicto']}.",
    ]
    return Seccion(7, "Comprobar que una cohorte nueva es de verdad nueva", cuerpo, fig)


def _sec_composicion(certified: dict | None, fig: F.Figura | None) -> Seccion:
    if not certified:
        return Seccion(8, "Composición del tejido",
                       ["La etapa certificada de v2 no ha corrido."])
    c7 = certified["c7_composicion"]
    principal = next((a for a in c7["ajustes"] if a.get("es_principal")), None)
    tras = certified["c15_r1_direccion"]["tras_el_ajuste_de_composicion"]
    base = certified["c15_r1_direccion"]["global"]
    cuerpo = [
        "Un tumor no es solo células tumorales. Lleva tejido de soporte, vasos y "
        "células inmunes, y la proporción cambia de una muestra a otra. Una proteína "
        "puede parecer alterada simplemente porque el trozo de tumor tenía más "
        "estroma que el trozo de tejido normal.",
        "v1 ajustaba por un índice estromal construido con las proteínas de matriz "
        "extracelular, y después contrastaba esas mismas proteínas contra ese "
        "índice. Cada una contribuía a su propia covariable, así que parte del "
        "ajuste era la proteína explicándose a sí misma.",
        "v2 recalcula el índice sin la proteína que se está ajustando. Como el "
        "índice es la media de valores estandarizados, quitar una es aritmética "
        "simple y no hay que recomputar nada.",
    ]
    if principal:
        f = principal["fraccion"]
        cuerpo.append(
            f"Con ese ajuste sobreviven {_frac(f)} de los hallazgos certificados, "
            f"{_coma(f['valor'])} ({_ic(f)}). Menos de la mitad. "
            f"{c7['lenguaje_obligatorio']}"
        )
    if tras.get("estado") == "calculado":
        ft = tras["fraccion"]
        cuerpo.append(
            "Queda una pregunta distinta: ¿la concordancia entre cohortes venía del "
            "estroma? Se recalculó R1 usando el efecto ya ajustado tanto para elegir "
            f"las proteínas como para su signo. Pasa de {_coma(base['valor'])} a "
            f"{_coma(ft['valor'])} ({_frac(ft)}, {_ic(ft)}), con el denominador "
            f"reducido en {int(tras['n_proteinas_que_pierde_el_denominador'])} "
            f"proteínas. {tras['lectura']} La concordancia no venía del estroma."
        )
        cuerpo.append(
            "Con una limitación que hay que declarar: la réplica entra sin ajustar, "
            "porque no tiene suficientes proteínas de matriz extracelular medidas en "
            "los dos brazos para construirle su propio índice. Así que la "
            "comparación es efecto ajustado contra efecto crudo, y es una cota de lo "
            "que el ajuste puede romper, no un contraste simétrico."
        )
    return Seccion(8, "Composición del tejido, sin circularidad", cuerpo, fig)


def _sec_indice(cert: dict | None, fig: F.Figura | None) -> Seccion:
    if not cert:
        return Seccion(9, "El índice R", ["El certificado de v2 no se ha generado."])
    cuerpo = [f"El índice R no es una nota. {cert['regla_de_lectura']}"]
    for s in cert["subindices"]:
        if s["estado"] == "calculado":
            detalle = _coma(s["valor"])
            if s.get("numerador") is not None and s.get("denominador"):
                detalle += f" ({int(s['numerador'])} de {int(s['denominador'])})"
            cuerpo.append(f"**{s['nombre']}**: {detalle}. {s['definicion']}")
        elif s["estado"] == "calculado_como_tope":
            cuerpo.append(
                f"**{s['nombre']}**: {_coma(s['valor'], 6)}, y es un tope, no una "
                f"estimación. {s.get('motivo') or s['definicion']}"
            )
        else:
            cuerpo.append(f"**{s['nombre']}**: {s['estado']}. {s.get('motivo', '')}")
    h1 = cert["h1"]
    cuerpo.append(
        f"**H1** queda en `{h1['veredicto']}`. Es la conjunción de dos cláusulas: la "
        f"de proteína {h1['clausula_proteina']} y la de vía {h1['clausula_via']}. "
        f"{h1['motivo_clausula_proteina']}"
    )
    h2 = cert["h2"]
    cuerpo.append(
        f"**H2** queda en `{h2['estado']}`. {h2['motivo']} Para evaluarla haría "
        f"falta {h2['requisito_para_evaluar']}. {h2['que_hacia_v1']}"
    )
    return Seccion(9, "El índice R y las hipótesis", cuerpo, fig)


def _sec_limitaciones(auditorias: dict) -> Seccion:
    cuerpo = [
        "Lo que sigue no es un trámite. Cada punto es una razón concreta por la que "
        "algún número de este informe podría estar diciendo menos de lo que parece.",
    ]
    vistas: set[str] = set()
    for nombre, aud in auditorias.items():
        if not aud:
            continue
        for lim in aud.get("limitaciones", []) or []:
            texto = lim if isinstance(lim, str) else str(lim)
            if texto not in vistas:
                vistas.add(texto)
                cuerpo.append(f"- {texto}")
    cuerpo.append(
        "Y la limitación que engloba a todas las demás: esto es un trabajo "
        "metodológico sobre datos públicos, con tamaños de muestra de decenas de "
        f"pacientes. {CIERRE}"
    )
    return Seccion(10, "Limitaciones", cuerpo)


# --------------------------------------------------------------------- formatos
def _comparacion(changelog: pd.DataFrame, cert: dict | None) -> str:
    """Tabla de v1 contra v2, número por número, desde `changelog_numeros.tsv`."""
    lineas = [
        "# Qué cambia de v1 a v2, número por número",
        "",
        "Esta tabla sale de `certification/changelog_numeros.tsv`. Cada fila pone el "
        "valor de v1 al lado del de v2, con el motivo y el código de la corrección "
        "que lo movió. Si una fila dice que no cambia, el número es el mismo y lo "
        "que cambia es lo que se publica junto a él.",
        "",
        f"> **{CIERRE}**",
        "",
    ]
    if changelog.empty:
        lineas.append("El changelog no está disponible: `v2_certificado` no ha corrido.")
        return "\n".join(lineas)

    lineas += [
        "| Número | v1 | v2 | ¿Cambia? | Corrección | Motivo |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for _, r in changelog.iterrows():
        def fmt(v):
            num = pd.to_numeric(v, errors="coerce")
            if num != num:
                return "no disponible" if str(v) in ("nan", "") else str(v)
            return _coma(num, 6 if abs(float(num)) < 0.005 and num != 0 else 4)

        cambia = str(r.get("cambia", ""))
        marca = {"True": "sí", "False": "no"}.get(cambia, "no aplica")
        lineas.append(
            f"| `{r['numero']}` | {fmt(r['v1'])} | {fmt(r['v2'])} | {marca} | "
            f"{r.get('correccion', '')} | {str(r.get('motivo', '')).strip()} |"
        )
    lineas.append("")
    if cert:
        n_cambian = int((changelog["cambia"].astype(str) == "True").sum())
        lineas += [
            f"De los {len(changelog)} números del changelog, {n_cambian} cambian.",
            "",
            "## Lo que no cambia",
            "",
            "v1 no se toca. Su manifiesto de hashes está en "
            "`results/v1_manifest.sha256` y `make manifest-verify` falla si algún "
            "fichero de v1 o alguna matriz de entrada se ha movido. Las matrices "
            "crudas se leen en modo solo lectura y su manifiesto está en "
            "`results/v2/manifest_inputs.sha256`.",
            "",
            "## Enmiendas",
            "",
            "El preregistro no se edita. Cuando un umbral o un método cambió después "
            "de ver los datos de v1, queda escrito como enmienda en "
            "`docs/enmiendas_v2.md`, con la fecha y la justificación.",
        ]
    return "\n".join(lineas)


def _readme(figs: dict, auditorias: dict, *, hay_pdf: bool) -> str:
    """Cómo reproducir v2 desde cero, en el orden en el que las etapas dependen."""
    etapas = [
        ("make manifest-write", "Congela los hashes de v1 y de las matrices de entrada."),
        ("make session-info", "Guarda versiones de Python y de cada dependencia."),
        ("make regresion-v1", "Re-ejecuta el contraste de v1 en un workspace temporal."),
        ("make v2-external", "Descargas idempotentes. Cada fuente queda marcada."),
        ("make v2-certified", "Contraste primario, nulos, C7 composición, R1 y R2."),
        ("make v2-axes", "C5 pólipos y C6 selección de modelo."),
        ("make v2-vias", "C4 cláusula de vías y C13 nulo de etiquetas de gen."),
        ("make v2-deteccion", "C8 ventana de detección con el Human Protein Atlas."),
        ("make v2-meta", "C10 meta-análisis REML y C11 tercera cohorte."),
        ("make v2-certificado", "Certificado, índice R, changelog e informe de brecha."),
        ("make v2-report", "Figuras, Markdown, HTML y PDF."),
        ("make verificar-numeros", "Falla si un número del informe no sale de una tabla."),
        ("make manifest-verify", "Comprueba que v1 y las entradas siguen intactas."),
    ]
    lineas = [
        "# results/v2 — qué hay aquí y cómo se rehace",
        "",
        f"> **{CIERRE}**",
        "",
        "Esta carpeta es la versión 2 del análisis. v1 queda congelada en "
        "`results/certified/` y no se sobrescribe nunca.",
        "",
        "## Reproducir desde cero",
        "",
        "`make all-v2` ejecuta todo en orden. Por separado:",
        "",
        "| Orden | Comando | Qué hace |",
        "| --- | --- | --- |",
    ]
    for i, (cmd, que) in enumerate(etapas, start=1):
        lineas.append(f"| {i} | `{cmd}` | {que} |")
    lineas += [
        "",
        "Las etapas de v2 están fuera de `--stage all` a propósito: una corrida del "
        "pipeline antiguo no puede sobrescribir nada de v2.",
        "",
        "## Qué hay en cada carpeta",
        "",
        "| Carpeta | Contenido |",
        "| --- | --- |",
        "| `certified/` | Contraste primario de v2, nulos, ajustes por composición. |",
        "| `axes/` | C5 pólipos y C6 modelos. |",
        "| `vias/` | C4 por colección, prueba competitiva y nulo de etiquetas. |",
        "| `deteccion/` | C8 curva de detección y concordancia por ventana. |",
        "| `meta/` | C10 meta-análisis por proteína e I² sin CPTAC. |",
        "| `certification/` | Certificado, índice R, changelog, informe de brecha. |",
        "| `report/` | Informe en Markdown, HTML y PDF, con sus figuras. |",
        "",
        "## Figuras",
        "",
        f"Las {len(figs)} figuras están en `report/figuras/`, cada una en PNG a "
        f"{F.DPI} dpi y en SVG, y cada una con el TSV con el que se dibujó.",
        "",
        "| Figura | Tabla de origen |",
        "| --- | --- |",
    ]
    for fig in figs.values():
        lineas.append(f"| `{fig.png.name}` | `{fig.tsv.name}` |")
    lineas += [
        "",
        "## Reglas que esta carpeta cumple",
        "",
        "- Ninguna matriz cruda se modifica. Se leen en modo solo lectura y su "
        "manifiesto de hashes se verifica al final.",
        "- No se imputan valores faltantes, no se fusionan cohortes en un solo "
        "modelo y el cero de MaxQuant se trata como ausencia.",
        "- Ningún número de los informes está escrito a mano. "
        "`make verificar-numeros` falla si alguno no sale de una tabla.",
        "- Los datos sintéticos solo existen dentro de `tests/`, con el prefijo "
        "`SINTETICO_`. Nunca entran en un resultado.",
        "- Si una fuente no se pudo descargar, queda marcada `no_disponible` con el "
        "motivo en `docs/fuentes_datos.md`. No se rellena con nada.",
        "- Ningún dato clínico personal entra en el repositorio ni en los informes.",
        "",
        "## Despliegue",
        "",
        "Es local. `make serve-v2` sirve `report/` en `http://localhost:8000`. "
        + ("El PDF está en `report/OncoGO-CRC_informe_v2.pdf`." if hay_pdf
           else "El PDF no se generó: falta `reportlab`."),
        "",
    ]
    faltan = [k for k, v in auditorias.items() if v is None]
    if faltan:
        lineas += [
            "## Etapas sin auditoría en esta corrida",
            "",
            f"{', '.join(faltan)}. Sus secciones del informe lo dicen en lugar de "
            "repetir el número de v1.",
            "",
        ]
    return "\n".join(lineas)


def _a_markdown(secciones: list[Seccion], *, figuras_rel: str) -> str:
    lineas = [
        "# OncoGO-CRC · Informe de la versión 2",
        "",
        "> Este informe se genera desde las tablas de `results/v2/`. "
        "Ningún número está escrito a mano.",
        "",
        f"> **{CIERRE}**",
        "",
        "## Índice",
        "",
    ]
    for s in secciones:
        lineas.append(f"{s.numero}. [{s.titulo}](#sec-{s.numero})")
    lineas.append("")
    for s in secciones:
        lineas += [f'<a id="sec-{s.numero}"></a>', "", f"## {s.numero}. {s.titulo}", ""]
        for p in s.cuerpo:
            lineas += [p, ""]
        if s.figura:
            lineas += [
                f"![{s.figura.titulo}]({figuras_rel}/{s.figura.png.name})",
                "",
                f"**Figura {s.numero}.** {s.figura.titulo}. {s.figura.pie}",
                "",
                f"*Tabla de origen: `{s.figura.tsv.name}`.*",
                "",
            ]
    return "\n".join(lineas)


def _a_html(secciones: list[Seccion], *, figuras_rel: str) -> str:
    def p(texto: str) -> str:
        # El cuerpo solo usa **negrita** y `código`; se convierten a mano para no
        # arrastrar una dependencia de Markdown solo por dos marcas.
        t = html.escape(texto)
        partes = t.split("**")
        t = "".join(x if i % 2 == 0 else f"<strong>{x}</strong>" for i, x in enumerate(partes))
        partes = t.split("`")
        return "".join(x if i % 2 == 0 else f"<code>{x}</code>" for i, x in enumerate(partes))

    cuerpo = []
    for s in secciones:
        cuerpo.append(f'<section id="sec-{s.numero}">')
        cuerpo.append(f"<h2>{s.numero}. {html.escape(s.titulo)}</h2>")
        for texto in s.cuerpo:
            if texto.startswith("- "):
                cuerpo.append(f"<ul><li>{p(texto[2:])}</li></ul>")
            else:
                cuerpo.append(f"<p>{p(texto)}</p>")
        if s.figura:
            cuerpo.append(
                f'<figure><img src="{figuras_rel}/{s.figura.svg.name}" '
                f'alt="{html.escape(s.figura.titulo)}">'
                f"<figcaption><strong>Figura {s.numero}.</strong> "
                f"{html.escape(s.figura.titulo)}. {html.escape(s.figura.pie)} "
                f"<em>Tabla de origen: {s.figura.tsv.name}.</em></figcaption></figure>"
            )
        cuerpo.append("</section>")

    nav = "".join(
        f'<li><a href="#sec-{s.numero}">{html.escape(s.titulo)}</a></li>' for s in secciones
    )
    return f"""<!doctype html>
<html lang="es"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>OncoGO-CRC · Informe v2</title>
<style>
 :root {{ --tinta:#1d1d1f; --suave:#5a5a5a; --linea:#e4e4e4; --azul:#2a5d8f; }}
 * {{ box-sizing:border-box; }}
 body {{ margin:0; font:16px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
        color:var(--tinta); background:#fff; }}
 .envoltorio {{ display:grid; grid-template-columns:260px minmax(0,1fr); gap:3rem;
        max-width:1180px; margin:0 auto; padding:2.5rem 1.5rem 5rem; }}
 nav {{ position:sticky; top:2.5rem; align-self:start; font-size:13.5px; }}
 nav ol {{ margin:0; padding-left:1.2rem; color:var(--suave); }}
 nav li {{ margin:.45rem 0; }}
 nav a {{ color:var(--suave); text-decoration:none; }}
 nav a:hover {{ color:var(--azul); }}
 h1 {{ font-size:1.9rem; line-height:1.25; margin:0 0 .5rem; letter-spacing:-.02em; }}
 h2 {{ font-size:1.25rem; margin:3rem 0 .9rem; letter-spacing:-.01em;
       border-top:1px solid var(--linea); padding-top:1.6rem; }}
 section:first-of-type h2 {{ border-top:0; padding-top:0; margin-top:1.5rem; }}
 p {{ margin:0 0 1rem; }}
 code {{ font:13.5px/1.4 ui-monospace,SFMono-Regular,Menlo,monospace;
         background:#f4f4f6; padding:.1rem .32rem; border-radius:4px; }}
 figure {{ margin:1.8rem 0; }}
 figure img {{ width:100%; height:auto; border:1px solid var(--linea); border-radius:6px; }}
 figcaption {{ font-size:13.5px; color:var(--suave); margin-top:.7rem; }}
 .aviso {{ border-left:3px solid var(--azul); background:#f5f8fb; padding:.9rem 1.1rem;
           margin:1.2rem 0; font-size:14.5px; }}
 .clinico {{ border-left:3px solid #a4243b; background:#fdf5f6; padding:.9rem 1.1rem;
             margin:1.2rem 0 2rem; font-weight:600; }}
 ul {{ margin:0 0 1rem; padding-left:1.3rem; }}
 @media (max-width:860px) {{ .envoltorio {{ grid-template-columns:1fr; gap:1rem; }}
   nav {{ position:static; }} }}
</style></head><body><div class="envoltorio">
<nav><ol>{nav}</ol></nav>
<main>
<h1>OncoGO-CRC · Informe de la versión 2</h1>
<div class="aviso">Este informe se genera desde las tablas de <code>results/v2/</code>.
Ningún número está escrito a mano.</div>
<div class="clinico">{CIERRE}</div>
{"".join(cuerpo)}
</main></div></body></html>
"""


def _a_pdf(secciones: list[Seccion], destino: Path) -> bool:
    """PDF con reportlab. Devuelve False si la dependencia no está."""
    try:
        from reportlab.lib.enums import TA_JUSTIFY
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import ParagraphStyle
        from reportlab.lib.units import cm
        from reportlab.platypus import (
            Image,
            PageBreak,
            Paragraph,
            SimpleDocTemplate,
            Spacer,
        )
    except ImportError:
        return False

    from reportlab.lib import colors

    cuerpo = ParagraphStyle(
        "cuerpo", fontName="Helvetica", fontSize=9.6, leading=14.2,
        alignment=TA_JUSTIFY, spaceAfter=7,
    )
    h1 = ParagraphStyle(
        "h1", fontName="Helvetica-Bold", fontSize=19, leading=23, spaceAfter=10,
    )
    h2 = ParagraphStyle(
        "h2", fontName="Helvetica-Bold", fontSize=12.6, leading=16,
        spaceBefore=16, spaceAfter=7, textColor=colors.HexColor("#1d1d1f"),
    )
    pie = ParagraphStyle(
        "pie", fontName="Helvetica-Oblique", fontSize=8.2, leading=11.4,
        alignment=TA_JUSTIFY, textColor=colors.HexColor("#5a5a5a"), spaceAfter=10,
    )
    aviso = ParagraphStyle(
        "aviso", fontName="Helvetica-Bold", fontSize=9.8, leading=14,
        textColor=colors.HexColor("#a4243b"), spaceAfter=14,
    )

    def marcas(texto: str) -> str:
        t = texto.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        partes = t.split("**")
        t = "".join(x if i % 2 == 0 else f"<b>{x}</b>" for i, x in enumerate(partes))
        partes = t.split("`")
        return "".join(
            x if i % 2 == 0 else f'<font face="Courier">{x}</font>'
            for i, x in enumerate(partes)
        )

    flujo = [
        Paragraph("OncoGO-CRC", h1),
        Paragraph("Informe de la versión 2", h2),
        Paragraph(
            "Todas las cifras de este documento salen de las tablas de "
            "<font face='Courier'>results/v2/</font>. Ninguna está escrita a mano.",
            cuerpo,
        ),
        Paragraph(CIERRE, aviso),
    ]
    for s in secciones:
        flujo.append(Paragraph(f"{s.numero}. {marcas(s.titulo)}", h2))
        for texto in s.cuerpo:
            flujo.append(Paragraph(marcas(texto.removeprefix("- ")), cuerpo))
        if s.figura and s.figura.png.exists():
            ancho = 16.0 * cm
            img = Image(str(s.figura.png))
            escala = ancho / img.imageWidth
            img.drawWidth = ancho
            img.drawHeight = img.imageHeight * escala
            flujo += [Spacer(1, 5), img, Spacer(1, 4)]
            flujo.append(
                Paragraph(
                    f"<b>Figura {s.numero}.</b> {marcas(s.figura.titulo)}. "
                    f"{marcas(s.figura.pie)} Tabla de origen: {s.figura.tsv.name}.",
                    pie,
                )
            )
    flujo += [PageBreak(), Paragraph("Aviso final", h2), Paragraph(CIERRE, aviso)]

    destino.parent.mkdir(parents=True, exist_ok=True)
    SimpleDocTemplate(
        str(destino), pagesize=A4,
        leftMargin=2.1 * cm, rightMargin=2.1 * cm,
        topMargin=2.0 * cm, bottomMargin=1.9 * cm,
        title="OncoGO-CRC · Informe v2", author="OncoGO-CRC",
    ).build(flujo)
    return True


# --------------------------------------------------------------------- run
def run(settings: Settings) -> None:
    base = settings.resolve("results") / "v2"
    salida = base / "report"
    dir_figs = salida / "figuras"
    salida.mkdir(parents=True, exist_ok=True)

    auditorias = {
        "certified": _leer_json(base / "certified" / "inference_audit_v2.json"),
        "axes": _leer_json(base / "axes" / "axes_audit_v2.json"),
        "vias": _leer_json(base / "vias" / "vias_audit_v2.json"),
        "deteccion": _leer_json(base / "deteccion" / "deteccion_audit_v2.json"),
        "meta": _leer_json(base / "meta" / "meta_audit_v2.json"),
        "certificado": _leer_json(base / "certification" / "certificado_v2.json"),
    }
    faltan = [k for k, v in auditorias.items() if v is None]
    if faltan:
        print(f"[v2-report] etapas sin auditoría: {', '.join(faltan)}. Sus secciones lo dirán.")

    figs: dict[str, F.Figura] = {}
    if auditorias["axes"]:
        figs["c5"] = F.fig_c5_gradiente(auditorias["axes"], dir_figs)
        figs["c6"] = F.fig_c6_modelos(auditorias["axes"], dir_figs)
    if auditorias["vias"]:
        figs["c4"] = F.fig_c4_vias(auditorias["vias"], dir_figs)
    if auditorias["deteccion"]:
        figs["c8"] = F.fig_c8_deteccion(
            auditorias["deteccion"],
            _leer_tsv(base / "deteccion" / "ventana_deteccion_por_proteina.tsv"),
            dir_figs,
        )
    if auditorias["meta"]:
        figs["c10"] = F.fig_c10_heterogeneidad(
            auditorias["meta"],
            _leer_tsv(base / "meta" / "meta_por_proteina.tsv"),
            _leer_tsv(base / "meta" / "i2_solo_cohortes_lfq.tsv"),
            dir_figs,
        )
        dir_ind = settings.resolve("external") / "independencia"
        correlaciones = {
            p.stem.replace("correlaciones_", ""): _leer_tsv(p)
            for p in sorted(dir_ind.glob("correlaciones_*.tsv"))
            if not p.name.startswith("._")
        }
        figs["c11"] = F.fig_c11_independencia(
            auditorias["meta"], correlaciones, dir_figs
        )
    if auditorias["certified"]:
        figs["c7"] = F.fig_c7_c9_composicion(auditorias["certified"], dir_figs)
    if auditorias["certificado"]:
        figs["r"] = F.fig_indice_r(
            _leer_tsv(base / "certification" / "changelog_numeros.tsv"),
            _leer_tsv(base / "certification" / "indice_R_v2.tsv"),
            dir_figs,
        )
    print(f"[v2-report] {len(figs)} figuras en PNG a {F.DPI} dpi y SVG, cada una con su TSV.")

    secciones = [
        _sec_que_es_esto(auditorias["certificado"]),
        _sec_gradiente(auditorias["axes"], figs.get("c5")),
        _sec_modelos(auditorias["axes"], figs.get("c6")),
        _sec_vias(auditorias["vias"], figs.get("c4")),
        _sec_deteccion(auditorias["deteccion"], figs.get("c8")),
        _sec_meta(auditorias["meta"], figs.get("c10")),
        _sec_independencia(auditorias["meta"], figs.get("c11")),
        _sec_composicion(auditorias["certified"], figs.get("c7")),
        _sec_indice(auditorias["certificado"], figs.get("r")),
        _sec_limitaciones(auditorias),
    ]

    (salida / "informe_v2.md").write_text(
        _a_markdown(secciones, figuras_rel="figuras"), encoding="utf-8"
    )
    (salida / "index.html").write_text(
        _a_html(secciones, figuras_rel="figuras"), encoding="utf-8"
    )
    hecho = _a_pdf(secciones, salida / "OncoGO-CRC_informe_v2.pdf")
    if not hecho:
        print("[v2-report] reportlab no está instalado: no hay PDF, sí HTML y Markdown.")

    changelog = _leer_tsv(base / "certification" / "changelog_numeros.tsv")
    (salida / "comparacion_v1_v2.md").write_text(
        _comparacion(changelog, auditorias["certificado"]), encoding="utf-8"
    )
    (base / "README.md").write_text(
        _readme(figs, auditorias, hay_pdf=hecho), encoding="utf-8"
    )

    # El informe de brecha ya lo genera v2_certificado; se copia al directorio del
    # despliegue para que todo lo publicable esté en un sitio.
    brecha = base / "certification" / "informe_brecha_v2.md"
    if brecha.exists():
        shutil.copy2(brecha, salida / "informe_brecha_v2.md")

    print(f"[v2-report] escrito en {salida}")
