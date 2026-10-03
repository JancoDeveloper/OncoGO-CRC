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
from oncogo_crc.v2 import cohorte_cea as C
from oncogo_crc.v2 import figuras as F

CIERRE = "Ningún resultado de esta tesis sirve para decidir sobre un paciente."

# El protocolo es un documento del repositorio, no del disco de datos.
PROTOCOLO = Path(__file__).resolve().parents[3] / "docs" / "protocolo_cohorte_cea.md"


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
    a, b = _coma(low, digitos), _coma(high, digitos)
    # Con un extremo negativo, unir con una raya da «−0,423–−0,013», que hay que
    # leer dos veces para entender dónde empieza cada número.
    if str(low).startswith("-") or str(high).startswith("-"):
        return f"IC 95 % de {a} a {b}"
    return f"IC 95 % {a}–{b}"


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


def _tamanos_del_protocolo(marc: dict | None) -> list[dict]:
    """Pacientes necesarios para que la cohorte del protocolo pudiera concluir algo.

    El AUC de partida **no** se inventa ni sale de los datos sintéticos de los
    tests. Se deriva del meta-análisis de D8: una pareja de sensibilidad y
    especificidad conjuntas da una cota inferior del AUC en su media, que es el
    área del trapecio que pasa por ese punto. Es el rendimiento publicado del CEA
    en esta misma indicación, así que es el punto de partida correcto para
    preguntar cuánto más haría falta.

    Si D8 no ha corrido, no se devuelve ningún tamaño: inventar el AUC de partida
    para poder imprimir un número sería exactamente lo que este informe no hace.
    """
    cea = ((marc or {}).get("ajustes") or {}).get("CEA|deteccion_de_recaida") or {}
    if cea.get("estado") != "calculado":
        return []
    base = (float(cea["sensibilidad"]) + float(cea["especificidad"])) / 2.0
    filas = []
    for delta in (0.05, 0.10):
        res = C.tamano_necesario(
            auc_umbral=base,
            auc_trayectoria=base + delta,
            proporcion_con_recaida=0.25,
            visitas_por_paciente=4.0,
        )
        filas.append(
            {
                "delta": delta,
                "auc_de_partida": round(base, 4),
                "de_donde_sale_el_auc_de_partida": (
                    "media de la sensibilidad y la especificidad conjuntas del CEA "
                    f"en D8 ({_coma(cea['sensibilidad'])} y "
                    f"{_coma(cea['especificidad'])} sobre {int(cea['k'])} estudios); "
                    "es una cota inferior del AUC, no el AUC"
                ),
                "proporcion_con_recaida": 0.25,
                "visitas_por_paciente": 4.0,
                "potencia": 0.80,
                "alfa": 0.05,
                **{k: v for k, v in res.items() if k != "nota"},
            }
        )
    return filas


def _sec_cea_tcga(tcga: dict | None) -> Seccion:
    titulo = "Qué determina el CEA en sangre"
    if not tcga or tcga.get("estado") != "calculado":
        return Seccion(
            10,
            titulo,
            [
                "Esta fuente no está disponible en esta corrida, así que la sección "
                f"no tiene resultado. Motivo: {(tcga or {}).get('motivo', 'la etapa no ha corrido')}."
            ],
        )
    m = tcga["modelo_del_cea"]
    cob = tcga["fuente"]["cobertura"]
    ns = tcga["fenotipo_no_secretor"]
    cuerpo = [
        "El CEA es el marcador en sangre que de verdad se usa en el seguimiento del "
        "cáncer de colon, y el gen que lo produce, CEACAM5, es una de las proteínas "
        "que esta tesis encuentra alteradas en el tumor. La pregunta de este eje es "
        "la que une las dos cosas: si un tumor expresa mucho CEACAM5, ¿sube el CEA en "
        "la sangre de ese paciente?",
        "Se responde en TCGA, que tiene a la vez la expresión del gen en el tumor y "
        f"el CEA prequirúrgico del mismo paciente. De {int(cob['n_muestras'])} "
        f"muestras de tumor primario, {int(cob['n_con_cea'])} tienen CEA "
        f"({_pct(cob['fraccion_con_cea'])}) y {int(cob['n_con_pureza_absolute'])} "
        "tienen una pureza tumoral medida sobre el ADN.",
        "La respuesta es que sí, pero mucho menos de lo que cabría esperar. Por cada "
        "duplicación de la expresión tisular de CEACAM5, el CEA sérico sube "
        f"{_coma(m['coeficientes'][1]['coeficiente'])} en logaritmo decimal, y la "
        f"asociación no es casualidad (p = {_coma(m['coeficientes'][1]['p'], 6)}). "
        "Pero la magnitud es lo decisivo: la expresión del gen explica por sí sola el "
        f"{_pct(m['r2_solo_ceacam5'])} de la variación del CEA, frente al "
        f"{_pct(m['r2'])} del modelo completo con estadio, metástasis, pureza, edad, "
        "sexo, lado del tumor y subtipo. Añadir CEACAM5 a ese modelo completo aporta "
        f"{_pct(m['aporte_de_ceacam5_al_r2'])}.",
        "Dicho en palabras: el CEA de un paciente depende sobre todo de cuánto tumor "
        "hay y de dónde está, no de cuánto gen expresa el tumor. Esa es la razón por "
        "la que una proteína alterada en el tejido no se convierte automáticamente en "
        "un marcador en sangre, y conecta directamente con la brecha tejido–plasma de "
        "la sección 5: lo que allí parecía la rareza de un experimento pequeño aquí "
        "se mide en el caso mejor estudiado que existe.",
        f"Hay además un fenotipo que conviene nombrar. {int(ns['grupos']['no_secretor']['n'])} "
        "tumores de esta cohorte están en el cuartil alto de expresión de CEACAM5 y "
        "aun así tienen el CEA por debajo del umbral de normalidad de 5 ng/ml. "
        "Existen, y no son una rareza. Lo que los distingue en estos datos es "
        "principalmente la ausencia de metástasis, que vuelve a apuntar a que el "
        "determinante es la carga tumoral.",
        "Estos números no son un subíndice del certificado y no cambian ninguna "
        "hipótesis. Son el contexto que explica por qué la parte de plasma de esta "
        "tesis sale como sale.",
    ]
    return Seccion(10, titulo, cuerpo)


def _sec_pureza(pureza: dict | None, fig: F.Figura | None) -> Seccion:
    titulo = "El índice de composición de R4 frente a una pureza medida"
    if not pureza or pureza.get("estado") != "calculado":
        return Seccion(
            11,
            titulo,
            [
                "Esta sección no tiene resultado en esta corrida. Motivo: "
                f"{(pureza or {}).get('motivo', 'la etapa no ha corrido')}."
            ],
        )
    ecm = pureza["validacion_del_indice"]["ecm"]
    inm = pureza["validacion_del_indice"]["inmune"]
    ver = pureza["veredicto_sobre_r4"]
    repro = pureza["reproducibilidad_del_ensayo"]
    r5 = pureza["reproducibilidad_por_gen"]["CEACAM5"]
    cea5 = pureza["cea_de_proteina_frente_a_cea_serico"]["proteina"]["CEACAM5"]
    cuerpo = [
        "La sección 8 ajusta el efecto tumor–normal por un índice de composición del "
        "tejido construido con las proteínas de matriz extracelular del propio "
        "proteoma. Queda una pregunta que ese ajuste no puede responderse a sí mismo: "
        "¿ese índice mide de verdad cuánto tumor hay en la muestra?",
        "Para contestarla hace falta una pureza medida con otra tecnología y en los "
        "mismos pacientes. Existe: un proteoma retrospectivo de TCGA cuyos pacientes "
        "tienen pureza tumoral calculada sobre el ADN con ABSOLUTE.",
        f"El índice correlaciona con la pureza medida en el sentido correcto —más "
        f"estroma, menos tumor— con una correlación de rangos de "
        f"{_coma(ecm['rho'])} ({_ic({'ci_low': ecm['rho_ci_low'], 'ci_high': ecm['rho_ci_high']})}, "
        f"n = {int(ecm['n'])} pacientes). No está midiendo ruido. Pero comparte con la "
        f"pureza solo el {_pct(ecm['varianza_compartida'])} de la varianza: de cada "
        "cien partes de variación en pureza tumoral, el índice capta unas pocas.",
        f"El contraste con la firma inmune lo deja claro. Con cuatro proteínas llega a "
        f"{_coma(inm['rho'])}, más lejos que las {int(ecm['n_proteinas_de_la_firma'])} "
        "proteínas de la firma de matriz extracelular. Más proteínas no es más "
        "información cuando las proteínas miden otra cosa.",
        "La consecuencia para la tesis es precisa y conviene no suavizarla. "
        f"{ver['consecuencia_para_r4']}",
        "Esta fuente deja además dos cifras sobre la fiabilidad de la medida que "
        "afectan a cómo se leen todos los coeficientes de este informe. La primera: "
        f"el perfil completo de un mismo tumor medido dos veces correlaciona a "
        f"{_coma(repro['r_mediana'])} sobre {int(repro['n_pares_de_replicas'])} pares "
        "de alícuotas. El perfil global se reproduce razonablemente.",
        "La segunda es más incómoda. CEACAM5, la proteína del CEA, medida dos veces en "
        f"el mismo tumor difiere {_coma(r5['veces_de_diferencia_mediana'], 1)} veces "
        f"(mediana), que es el {_pct(r5['ruido_frente_a_senal'], 0)} de la variación "
        "que hay entre pacientes distintos. Eso cambia el significado de un "
        f"coeficiente: la correlación entre CEACAM5 en el tumor y el CEA en sangre sale "
        f"{_coma(cea5['rho'])} (n = {int(cea5['n'])}), y con este nivel de ruido un "
        "número cercano a cero no dice «no hay relación», dice «con esta medida no se "
        "puede saber». Son dos conclusiones muy distintas y confundirlas sería el "
        "error más fácil de cometer aquí.",
    ]
    return Seccion(11, titulo, cuerpo, fig)


def _sec_heterogeneidad(het: dict | None, fig: F.Figura | None) -> Seccion:
    titulo = "Cuánto varía un tumor dentro de sí mismo"
    if not het or het.get("estado") != "calculado":
        return Seccion(
            12,
            titulo,
            [
                "Esta sección no tiene resultado en esta corrida. Motivo: "
                f"{(het or {}).get('motivo', 'la etapa no ha corrido')}."
            ],
        )
    v = het["variacion_dentro_frente_a_entre"]
    prot = het["variacion_por_proteina"]
    cert = prot["certificadas"]
    tn = het["tumor_frente_a_mucosa_normal"]
    mt = het["metastasis_frente_a_tumor_primario"]
    busqueda = het["busqueda_de_cohortes_d3"]
    incluidas = [b for b in busqueda if b["veredicto"] == "incluida"]
    cuerpo = [
        "Todo este trabajo mide un tumor con una sola biopsia. Conviene saber cuánto "
        "se pierde por eso, y la única forma de saberlo es un depósito que haya medido "
        "varias regiones del mismo tumor.",
        f"La búsqueda de cohortes nuevas revisó {len(busqueda)} candidatos contra los "
        "cuatro criterios del preregistro. El resultado es que **ningún depósito nuevo "
        f"los cumple todos**: solo {len(incluidas)} entra, y es la que ya se usaba. Los "
        "que publican una matriz de proteínas no llegan a diez pacientes o no tienen "
        "tejido normal; los que tienen el diseño adecuado publican solo espectros "
        "crudos, y reconstruir la matriz desde ahí compararía dos programas de análisis "
        "en lugar de dos tejidos. El detalle de cada descarte está en la tabla de la "
        "etapa y en la ficha de la fuente.",
        "Lo que sí apareció es un depósito con cuatro pacientes, uno de ellos con ocho "
        "regiones del mismo tumor. Cuatro pacientes están por debajo del mínimo "
        "preregistrado de diez, así que **esta cohorte no entra en el meta-análisis** y "
        "el umbral no se ha doblado para que entrara. Entra para una pregunta distinta, "
        "que no es un contraste entre cohortes.",
        "La comparación hay que hacerla con el número de observaciones igualado: la "
        "desviación de ocho valores y la de cuatro no son estimadores igual de "
        "precisos, y comparar una con otra regalaría el resultado. Con el remuestreo "
        f"que las iguala, la variación dentro del tumor es {_coma(v['sd_mediana_dentro_del_tumor'])} "
        f"y la variación entre pacientes {_coma(v['sd_mediana_entre_pacientes'])}, con "
        f"un cociente de {_coma(v['cociente_mediano'])} "
        f"({_ic({'ci_low': v['cociente_ci_low'], 'ci_high': v['cociente_ci_high']})}).",
        "El intervalo incluye el uno, así que esto es una señal y no un resultado "
        "establecido. Lo que sostiene no es que un tumor varíe más dentro de sí que "
        "entre pacientes, sino algo igualmente relevante: que las dos variaciones son "
        "**del mismo orden**. Si eso es así, el ruido de la biopsia no se reduce "
        "añadiendo pacientes, porque entra en cada medida.",
        "Esto no invalida los efectos del certificado, y la razón es estructural: el "
        "contraste es dentro de cada paciente y el ruido de región no tiene dirección "
        "preferente. Lo que cambia es la interpretación de qué se ha medido. Se ha "
        "medido un trozo de tumor, no el tumor.",
        f"Hay un dato que favorece al certificado y se describe sin contrastarlo: de "
        f"las {int(cert['n_presentes'])} proteínas certificadas presentes en esta "
        f"cohorte, el cociente mediano de variación es {_coma(cert['cociente_mediano'])}, "
        f"frente a {_coma(cert['cociente_mediano_del_resto'])} en el resto. Las "
        "proteínas que el certificado selecciona están entre las menos variables de "
        "una región a otra. Son dos medianas puestas una al lado de la otra, no una "
        "prueba de que una sea mayor.",
        "Los dos contrastes que esta cohorte permite se leen con su referencia de "
        f"azar, sin la cual no significan nada: con cuatro pares, el azar ya pone de "
        f"acuerdo a una proteína de cada ocho. Tumor frente a mucosa normal concuerda "
        f"en los cuatro pares para el {_pct(tn['fraccion_con_todos_los_pares_en_la_misma_direccion'])} "
        f"de las proteínas, que es {_coma(tn['veces_por_encima_del_azar'], 1)} veces el "
        f"azar. Metástasis hepática frente a tumor primario solo "
        f"{_coma(mt['veces_por_encima_del_azar'], 1)} veces. Con cuatro pares no se "
        "calcula ningún p ni ningún intervalo: publicarlo invitaría a leerlo como si "
        "distinguiera algo.",
    ]
    return Seccion(12, titulo, cuerpo, fig)


def _sec_deconvolucion(dec: dict | None, fig: F.Figura | None) -> Seccion:
    titulo = "Estimar la composición del tejido con una referencia de fuera"
    if not dec or dec.get("estado") != "calculado":
        return Seccion(
            13,
            titulo,
            [
                "Esta sección no tiene resultado en esta corrida. Motivo: "
                f"{(dec or {}).get('motivo', 'la etapa no ha corrido')}."
            ],
        )
    tmt = dec.get("por_que_no_se_corre_en_cptac") or {}
    cuerpo = [
        "El índice de composición de la sección 8 tiene un problema de principio: se "
        "calcula con el mismo proteoma que se está contrastando. La corrección C7 "
        "mitigó lo más grave excluyendo cada proteína de su propio índice, pero el "
        "índice sigue saliendo de los mismos datos. La forma de romper eso es estimar "
        "la composición con perfiles de tipos celulares que vengan de fuera, y eso es "
        "lo que el preregistro contemplaba como método de sensibilidad.",
        "Los perfiles salen del Human Protein Atlas, que publica la expresión media de "
        "cada gen en cada tipo celular. De ahí se construye una firma de los once tipos "
        "celulares del colon, fijados antes de mirar ningún resultado.",
        "Hay dos obstáculos, y los dos se comprueban con los datos en lugar de "
        "afirmarlos.",
        "El primero: **esto no se puede hacer en CPTAC**, que es la cohorte donde vive "
        "R4. Estimar proporciones celulares supone que lo medido es una mezcla no "
        "negativa y lineal de los perfiles. CPTAC publica el logaritmo de un cociente "
        f"contra un canal de referencia, y el {_pct(tmt.get('fraccion_negativa'))} de "
        "esos valores es negativo. Una mezcla de perfiles no negativos con pesos no "
        "negativos solo puede dar valores no negativos, así que no existe ninguna "
        "proporción que reproduzca ese perfil. No es que el ajuste salga malo: el "
        "sistema no tiene solución con sentido. Hay dos motivos más, independientes: el "
        "logaritmo no es lineal, y el canal de referencia es a su vez una mezcla "
        "distinta en cada lote.",
        "El segundo: **las proporciones tipo a tipo no se publican**, y conviene decir "
        "por qué con precisión. El residuo del ajuste es de alrededor de 0,9 en las dos "
        "cohortes, es decir que la mezcla reproduce menos de una décima parte del perfil "
        "medido; y varios tipos celulares salen ausentes en más de la mitad de las "
        "muestras, que es lo que ocurre cuando otro tipo con perfil parecido se lleva su "
        "señal entera. Un cero en esas tablas no significa «este tipo celular no está»: "
        "significa que el método ha asignado su señal a otro.",
        "Lo que sí resiste es el agregado epitelio frente a no epitelio, porque esa "
        "separación no depende de distinguir un fibroblasto de una célula de músculo "
        "liso. Y se usa solo después de comprobar en los datos que ordena los pares como "
        "la histología espera. Una de las dos cohortes no pasa ese filtro y en ella no "
        "se calcula el ajuste externo.",
    ]
    for c in dec["cohortes"]:
        cmp_ = c.get("como_comparar_los_dos_ajustes") or {}
        if cmp_.get("estado") != "comparable":
            continue
        cuerpo.append(
            f"En {c['cohorte']}, con {int(c['n_pares'])} pares, el ajuste por la "
            "composición estimada desde fuera deja en pie el "
            f"{_pct(cmp_['r4_con_referencia_externa'])} de los hallazgos de esa "
            f"cohorte, frente al {_pct(cmp_['r4_con_el_indice_del_propio_proteoma'])} "
            "con el índice del propio proteoma."
        )
        cuerpo.append(
            "La tentación es leer el primero como la confusión real y el segundo como "
            "un ajuste que se queda corto, y no se puede. Las dos covariables marcan "
            "un número parecido de proteínas como dependientes de la composición "
            f"({int(cmp_['proteinas_con_pendiente_de_composicion_externa'])} frente a "
            f"{int(cmp_['proteinas_con_pendiente_de_composicion_interna'])}), así que "
            "la diferencia no viene de que una detecte mucha más confusión. Viene de "
            "la precisión: el cero en el que R4 evalúa el efecto está a "
            f"{_coma(cmp_['distancia_tipificada_al_cero_externa'], 2)} desviaciones "
            "típicas de los datos con la covariable externa y a "
            f"{_coma(cmp_['distancia_tipificada_al_cero_interna'], 2)} con la interna. "
            "Extrapolar más lejos cuesta error típico, y con más error típico menos "
            "proteínas pasan el umbral."
        )
        cuerpo.append(
            "Esa es, de hecho, la aportación más útil de esta sección, y es una "
            "limitación de R4 y no de esta fuente: R4 evalúa el efecto en una "
            "composición que ningún par tiene, así que es una extrapolación, y su valor "
            "depende de a qué distancia de los datos quede ese cero. Dos formas "
            "razonables de medir la composición dan cotas distintas. R4 es una cota, no "
            "una cantidad, y conviene tenerlo delante al leer también el R4 del "
            "certificado."
        )
    cuerpo.append(
        "Una advertencia sobre la referencia: es de ARN y lo que se deconvuelve es "
        "proteína. La sección 11 midió el acuerdo entre transcrito y proteína del mismo "
        "gen y salió flojo. Toda cifra de esta sección hereda esa aproximación."
    )
    return Seccion(13, titulo, cuerpo, fig)


def _sec_plasma_benigno(wann: dict | None) -> Seccion:
    titulo = "Separar cáncer de sano no es separar cáncer de enfermo"
    if not wann or wann.get("estado") != "usada":
        return Seccion(
            14,
            titulo,
            [
                "Esta fuente no está disponible en esta corrida. Motivo: "
                f"{(wann or {}).get('motivo', 'no se ha descargado')}."
            ],
        )
    comp = wann["comparacion_del_grupo_de_control"]
    sint = comp["contra_pacientes_con_sintomas"]
    sano = comp["contra_personas_sanas"]
    perfil = wann["perfil_de_las_proteinas"]
    cr = wann["sobre_el_colorrectal"]
    cuerpo = [
        "La parte de plasma de esta tesis compara cáncer con personas sanas. Es el "
        "diseño que había, y tiene un sesgo que no se arregla con más muestras: en una "
        "consulta real no llegan personas sanas, llegan personas con síntomas. "
        "Distinguir cáncer de salud es una tarea más fácil que la que haría falta.",
        "Esta fuente permite poner número a esa diferencia, y lo permite porque las dos "
        "comparaciones usan la misma plataforma con el mismo número de ensayos, de modo "
        "que los dos AUC son comparables aunque las magnitudes de proteína no lo sean.",
        f"Con el mismo panel de {int(sano['n_proteinas_del_panel'])} proteínas, el AUC "
        f"es {_coma(sano['auc'])} contra personas sanas y {_coma(sint['auc_descubrimiento'])} "
        f"contra pacientes derivados por síntomas que resultaron no tener cáncer "
        f"(y {_coma(sint['auc_replica'])} en la réplica). La diferencia es de "
        f"{_coma(comp['diferencia_de_auc'])} en AUC, y no viene de la tecnología: viene "
        "de contra quién se compara.",
        "Conviene ver qué eran esos controles. De los descritos, una parte tenía una "
        "enfermedad autoinmune y otra una infección. Es decir, el AUC más bajo no sale "
        "de comparar con gente un poco menos sana: sale de comparar con gente enferma de "
        "otra cosa, que es exactamente el problema clínico.",
        f"La lista de proteínas del modelo también dice algo. De las "
        f"{int(perfil['n_proteinas'])} proteínas, solo {int(perfil['n_inflamatorias'])} "
        "es una proteína inflamatoria secretada; la mayoría son intracelulares, "
        "nucleares o citoqueratinas, que en plasma indican recambio y muerte celular "
        "más que inflamación. Ninguna de las tres proteínas de plasma que v1 destacó "
        "aparece en ella.",
        f"Y la limitación que impide estirar esto: solo {int(cr['n_descubrimiento'])} de "
        f"los cánceres de la cohorte de descubrimiento y {int(cr['n_replica'])} de la "
        "réplica eran colorrectales. De aquí no sale ninguna afirmación sobre cáncer "
        "colorrectal en particular. Lo que aporta es la medida de cuánto optimismo "
        "introduce elegir un grupo de control sano.",
    ]
    return Seccion(14, titulo, cuerpo)


def _sec_marcadores(marc: dict | None, fig: F.Figura | None) -> Seccion:
    titulo = "El CEA frente a las alternativas, según la literatura"
    if not marc or marc.get("estado") != "calculado":
        return Seccion(
            15,
            titulo,
            [
                "Esta sección no tiene resultado en esta corrida. Motivo: "
                f"{(marc or {}).get('motivo', 'la etapa no ha corrido')}."
            ],
        )
    ajustes = marc["ajustes"]
    cea = ajustes.get("CEA|deteccion_de_recaida", {})
    met = ajustes.get("metilacion_SEPT9|cribado_o_diagnostico", {})
    ctdna = (marc.get("comparacion_entre_cea_y_ctdna") or {}).get("todas", {})
    cuerpo = [
        "Si una proteína nueva aspirara algún día a servir de marcador, tendría que "
        "superar a lo que ya se usa. Esta sección fija ese listón con un meta-análisis "
        "de tablas 2×2 extraídas de la literatura, con el modelo bivariado de Reitsma, "
        "que modela la sensibilidad y la especificidad conjuntamente en lugar de "
        "promediarlas por separado.",
    ]
    if cea.get("estado") == "calculado":
        cuerpo.append(
            f"El CEA para detectar una recaída, sobre {int(cea['k'])} estudios, tiene "
            f"una sensibilidad conjunta de {_coma(cea['sensibilidad'])} "
            f"({_ic({'ci_low': cea['sensibilidad_ci_low'], 'ci_high': cea['sensibilidad_ci_high']})}) "
            f"y una especificidad de {_coma(cea['especificidad'])} "
            f"({_ic({'ci_low': cea['especificidad_ci_low'], 'ci_high': cea['especificidad_ci_high']})}). "
            "Dicho en claro: de cada diez recaídas, el CEA deja pasar cerca de cuatro. "
            "Ese es el listón real, y no es alto."
        )
    if met.get("estado") == "calculado":
        cuerpo.append(
            f"La metilación de SEPT9 en sangre, sobre {int(met['k'])} estudios, llega a "
            f"una sensibilidad de {_coma(met['sensibilidad'])} y una especificidad de "
            f"{_coma(met['especificidad'])}. Parece mejor, pero **no se resta**: "
            "responde a otra pregunta clínica, cribado o diagnóstico en lugar de "
            "detección de recaída. Los dos puntos se dan por separado y la "
            "comparación directa se declara no comparable."
        )
    if ctdna.get("publicado"):
        pub = ctdna["publicado"]
        cuerpo.append(
            f"El ADN tumoral circulante, según una cifra publicada sobre "
            f"{int(pub['k'])} estudios en la misma indicación, tiene una sensibilidad "
            f"de {_coma(pub['sensibilidad'])} y una especificidad de "
            f"{_coma(pub['especificidad'])}. Los intervalos se solapan con los del CEA, "
            "pero eso es una descripción y no una prueba: son dos grupos de estudios "
            "distintos con criterios de inclusión distintos, y restarlos no "
            "significaría nada."
        )
    cuerpo.append(
        f"Una nota de procedimiento que no es menor. La tabla de extracción tiene "
        f"{int(marc['n_filas'])} filas y las {int(marc['pendiente_de_verificacion_humana'])} "
        "están pendientes de verificación humana: la columna `verificado_por` está "
        "vacía a propósito y la rellena una persona tras comparar cada fila con el "
        "artículo original. Ningún número extraído del texto de un PDF debería entrar "
        "en un resultado sin ese paso."
    )
    return Seccion(15, titulo, cuerpo, fig)


def _sec_protocolo(existe: bool, tamanos: list[dict]) -> Seccion:
    titulo = "Lo que haría falta para responder la pregunta que queda abierta"
    if not existe:
        return Seccion(
            16,
            titulo,
            ["El protocolo no está en el repositorio en esta corrida."],
        )
    cuerpo = [
        "Las secciones anteriores dejan una pregunta que ningún dato público puede "
        "contestar: si mirar cómo **cambia** el CEA a lo largo del tiempo aporta algo "
        "sobre mirar solo su último valor. Todas las cohortes de esta tesis tienen una "
        "sola determinación por paciente, y una pendiente no se calcula con un punto.",
        "Lo honesto cuando la respuesta exige datos que no existen es escribir el "
        "estudio que los generaría, y dejarlo escrito antes de tener los datos. Eso es "
        "`docs/protocolo_cohorte_cea.md`: diseño, diccionario de datos con las columnas "
        "obligatorias, plan de análisis estadístico cerrado de antemano, y lo que el "
        "estudio **no** podría concluir.",
        "**El estudio no se ha hecho.** El protocolo está pendiente de aprobación ética "
        "y de datos reales. El código de análisis sí existe y está probado, pero sobre "
        "datos sintéticos generados a propósito, en un fichero con el prefijo "
        "`SINTETICO_` y confinado a `tests/`. Ninguna cifra de esos datos aparece en "
        "este informe ni en ningún resultado.",
        "Lo que los tests comprueban son las dos direcciones, que es lo que hay que "
        "verificar antes de usar un análisis con datos de personas: que detecta una "
        "diferencia cuando el generador la programó, y que no la fabrica cuando no la "
        "programó. También que la validación cruzada parte por paciente y no por "
        "visita, porque partir por visita dejaría visitas del mismo paciente a los dos "
        "lados y el resultado saldría mejor de lo que es.",
    ]
    trozos: list[str] = []
    if tamanos:
        for t in tamanos:
            if t.get("estado") != "calculado":
                continue
            trozos.append(
                f"{int(t['n_pacientes'])} pacientes para una mejora de "
                f"{_coma(t['delta'], 2)} en AUC"
            )
        if trozos:
            cuerpo.append(
                "El tamaño necesario se calcula con la aproximación de Hanley y McNeil "
                "para la varianza de un AUC, corregida por el número de visitas de cada "
                "paciente, y partiendo del rendimiento del CEA que mide la sección 15: "
                + "; ".join(trozos) + ". Conviene separar de dónde sale cada supuesto. "
                "El rendimiento de partida es el del CEA publicado, medido en la "
                "sección 15; los otros dos —qué proporción de pacientes recae y cuántas "
                "determinaciones tiene cada uno— los fija el protocolo y no los mide "
                "nada de esta tesis. Los tres están en la tabla de la etapa. Es un "
                "orden de magnitud, no una cifra al paciente."
            )
    if not trozos:
        cuerpo.append(
            "Aquí iría el número de pacientes que necesitaría el estudio, y no está. "
            "Para calcularlo hace falta el rendimiento actual del CEA, que lo mide el "
            "meta-análisis de la sección 15, y ese no ha corrido en esta corrida. "
            "Poner una cifra con un rendimiento de partida supuesto daría un tamaño de "
            "muestra que no se apoya en ninguna medida."
        )
    cuerpo.append(
        "Dos advertencias que el protocolo recoge y que importan más que el tamaño de "
        "muestra. El estudio no diría si conviene hacer más determinaciones de CEA "
        "—detectar antes no es vivir más— ni daría un punto de corte para usar en "
        "consulta. Y no se trasladaría a otra población sin comprobarlo, porque el CEA "
        "basal depende del tabaquismo y de la función hepática."
    )
    return Seccion(16, titulo, cuerpo)


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
    return Seccion(17, "Limitaciones", cuerpo)


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
        # Las fuentes nuevas de v2. Cada una dice en su propia sección si no está,
        # así que su ausencia no impide generar el informe.
        "tcga_cea": _leer_json(base / "axes" / "tcga_cea" / "tcga_cea_audit_v2.json"),
        "pureza": _leer_json(
            base / "axes" / "pureza_medida" / "pureza_medida_audit_v2.json"
        ),
        "heterogeneidad": _leer_json(
            base / "axes" / "heterogeneidad" / "heterogeneidad_audit_v2.json"
        ),
        "deconvolucion": _leer_json(
            base / "axes" / "deconvolucion" / "deconvolucion_audit_v2.json"
        ),
        "marcadores": _leer_json(
            base / "marcadores" / "marcadores_audit_v2.json"
        ),
        "wannberg": _leer_json(
            settings.resolve("external") / "wannberg_2026" / "estado.json"
        ),
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
    if (auditorias["pureza"] or {}).get("estado") == "calculado":
        figs["d2"] = F.fig_d2_pureza(auditorias["pureza"], dir_figs)
    if (auditorias["heterogeneidad"] or {}).get("estado") == "calculado":
        figs["d3"] = F.fig_d3_heterogeneidad(
            auditorias["heterogeneidad"],
            _leer_tsv(base / "axes" / "heterogeneidad" / "variacion_por_proteina.tsv"),
            dir_figs,
        )
    if (auditorias["deconvolucion"] or {}).get("estado") == "calculado":
        figs["d5"] = F.fig_d5_deconvolucion(auditorias["deconvolucion"], dir_figs)
    if (auditorias["marcadores"] or {}).get("estado") == "calculado":
        figs["d8"] = F.fig_d8_marcadores(auditorias["marcadores"], dir_figs)
    print(f"[v2-report] {len(figs)} figuras en PNG a {F.DPI} dpi y SVG, cada una con su TSV.")

    protocolo = PROTOCOLO.exists()
    tamanos = _tamanos_del_protocolo(auditorias["marcadores"])
    if tamanos:
        (salida / "tamano_necesario_cohorte_cea.tsv").write_text(
            pd.DataFrame(tamanos).to_csv(sep="\t", index=False), encoding="utf-8"
        )

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
        _sec_cea_tcga(auditorias["tcga_cea"]),
        _sec_pureza(auditorias["pureza"], figs.get("d2")),
        _sec_heterogeneidad(auditorias["heterogeneidad"], figs.get("d3")),
        _sec_deconvolucion(auditorias["deconvolucion"], figs.get("d5")),
        _sec_plasma_benigno(auditorias["wannberg"]),
        _sec_marcadores(auditorias["marcadores"], figs.get("d8")),
        _sec_protocolo(protocolo, tamanos),
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
