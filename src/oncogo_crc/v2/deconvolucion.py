"""D5. Deconvolución de la composición celular con una referencia externa.

Qué problema resuelve. R4 mide cuánto del efecto tumor-normal sobrevive al
ajustar por la composición del tejido. En v2 esa composición se estima con un
índice de matriz extracelular construido del propio proteoma. C7 corrigió la
circularidad más grave (cada proteína se evalúa contra un índice que la excluye),
pero el índice sigue saliendo de los mismos datos. El preregistro lista
``deconvolucion_referencia_externa`` como método de sensibilidad precisamente
para romper eso: los perfiles de cada tipo celular vienen de fuera.

Cómo se estima. Mínimos cuadrados no negativos. Se busca el vector de
proporciones ``x >= 0`` que mejor explica el perfil medido ``b`` como mezcla
``A x`` de los perfiles de referencia ``A``.

Y aquí está la condición que casi nunca se comprueba: eso solo tiene sentido si
``b`` es una mezcla no negativa y **lineal** de los perfiles. Una intensidad LFQ
lo es, aproximadamente: si una muestra tiene el doble de fibroblastos, aporta el
doble de sus proteínas. Un **log** de un **cociente** contra un canal de
referencia, que es lo que publica CPTAC, no lo es por dos motivos
independientes: el logaritmo no es lineal, y el cociente ya ha dividido por otra
mezcla distinta en cada canal. Por eso la deconvolución **no** se corre sobre
CPTAC. La función ``por_que_no_en_tmt`` lo comprueba con los datos delante en vez
de afirmarlo.

Lo que esta deconvolución sostiene y lo que no. Las proporciones tipo a tipo no
son publicables: varios tipos comparten marcadores, y cuando dos columnas de
``A`` se parecen, los mínimos cuadrados no negativos reparten la señal entre
ellas de forma inestable y suelen dejar una en cero exacto. Lo que sí resulta
estable es el agregado epitelio frente a no epitelio, porque esa separación no
depende de distinguir fibroblasto de músculo liso. ``diagnosticar`` comprueba las
dos cosas en los datos y devuelve el veredicto; no se da por supuesto.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import nnls

# Mínimo de genes marcadores medidos para que la mezcla tenga con qué resolverse.
MIN_MARCADORES_MEDIDOS = 100

# Un tipo celular que queda por debajo de esto se declara ausente en esa muestra.
FRACCION_DESPRECIABLE = 0.005

# Si un tipo sale ausente en más de esta proporción de muestras, el método no lo
# está estimando en cero: lo está descartando. Un linfocito T no falta en la
# mitad de los tumores de colon, así que un cero repetido describe el método.
PROPORCION_DE_MUESTRAS_EN_CERO = 0.5

# Correlación entre columnas de la referencia por encima de la cual dos tipos
# celulares se declaran no separables entre sí.
COLINEALIDAD_ALTA = 0.8

# Residuo relativo por encima del cual la mezcla no explica el perfil medido.
# Con 0,9 la mezcla reproduce una décima parte de la norma del perfil: las
# proporciones no son una descripción de la muestra, son el mejor ajuste de un
# modelo que no encaja.
RESIDUO_ALTO = 0.5

# Estos tres umbrales no estaban en el preregistro y se declaran como parte de la
# enmienda E9.
ENMIENDA_DE_LOS_UMBRALES = "E9"

# Mínimo de muestras para correlacionar la fracción con otra cosa.
MIN_PARA_CORRELAR = 20


class FaltaColumna(KeyError):
    """Una columna que el análisis necesita no está en la tabla."""


@dataclass(frozen=True)
class Fracciones:
    """Proporciones estimadas y lo que costó estimarlas.

    ``residuo_relativo`` es la norma del residuo dividida por la norma del perfil
    medido. Cerca de cero significa que la mezcla explica el perfil; cerca de uno
    significa que no lo explica y que las proporciones no describen nada.
    """

    tabla: pd.DataFrame
    residuo_relativo: pd.Series
    genes_usados: list[str]
    tipos: list[str]


def _normalizar_columnas(a: np.ndarray) -> np.ndarray:
    """Cada tipo celular suma uno.

    Sin esto, un tipo con genes muy expresados se lleva una proporción mayor solo
    por la escala de su perfil, no por estar más presente.
    """
    s = a.sum(axis=0, keepdims=True)
    return a / np.where(s > 0, s, 1.0)


def fracciones(
    matriz: pd.DataFrame,
    columnas: list[str],
    firma: pd.DataFrame,
) -> Fracciones:
    """Resuelve la mezcla en cada muestra.

    ``matriz`` tiene que estar en escala **lineal** de abundancia, indexada por
    símbolo de gen en mayúsculas. Pasarle log2 daría proporciones de un problema
    que no es este.
    """
    faltan = [c for c in columnas if c not in matriz.columns]
    if faltan:
        raise FaltaColumna(
            f"no están en la matriz las columnas {faltan[:5]}; sin ellas no hay "
            "muestras que deconvolver"
        )
    genes = sorted(set(matriz.index.astype(str)) & set(firma.index.astype(str)))
    if len(genes) < MIN_MARCADORES_MEDIDOS:
        raise ValueError(
            f"solo {len(genes)} genes marcadores están medidos (mínimo "
            f"{MIN_MARCADORES_MEDIDOS}); con menos la mezcla no se resuelve"
        )
    tipos = [str(t) for t in firma.columns]
    a = _normalizar_columnas(firma.loc[genes, firma.columns].to_numpy(dtype=float))

    filas: dict[str, np.ndarray] = {}
    residuos: dict[str, float] = {}
    for c in columnas:
        b = matriz.loc[genes, c].to_numpy(dtype=float)
        # El cero de MaxQuant es ausencia, no cantidad cero medida, pero para la
        # mezcla una proteína ausente aporta cero: es el mismo número por una
        # razón distinta, y conviene que quede dicho.
        b = np.where(np.isfinite(b) & (b > 0.0), b, 0.0)
        total = float(b.sum())
        if total <= 0.0:
            continue
        b = b / total
        x, resid = nnls(a, b)
        suma = float(x.sum())
        filas[c] = x / suma if suma > 0 else x
        nb = float(np.linalg.norm(b))
        residuos[c] = float(resid / nb) if nb > 0 else float("nan")

    if not filas:
        raise ValueError(
            "ninguna muestra tenía señal en los genes marcadores; no hay nada "
            "que deconvolver"
        )
    return Fracciones(
        tabla=pd.DataFrame(filas, index=tipos).T,
        residuo_relativo=pd.Series(residuos, name="residuo_relativo"),
        genes_usados=genes,
        tipos=tipos,
    )


def eje_epitelial(fr: Fracciones, tipos_epiteliales: list[str]) -> pd.Series:
    """Suma de las proporciones epiteliales: el agregado que sí es estable.

    Agrupar es lo que lo salva. Que el método no sepa si una proteína de matriz
    viene de un fibroblasto o de una célula de músculo liso da igual para esta
    suma, porque ninguno de los dos es epitelio.
    """
    presentes = [t for t in tipos_epiteliales if t in fr.tabla.columns]
    if not presentes:
        raise FaltaColumna(
            f"ninguno de los tipos epiteliales {tipos_epiteliales} está en la "
            "deconvolución; sin ellos no hay eje epitelial"
        )
    return fr.tabla[presentes].sum(axis=1).rename("fraccion_epitelial")


def eje_epitelial_sin_gen(
    matriz: pd.DataFrame,
    columnas: list[str],
    firma: pd.DataFrame,
    gen: str,
    tipos_epiteliales: list[str],
) -> pd.Series:
    """El eje epitelial recalculado dejando fuera un gen.

    Es la misma cautela que C7 aplicó al índice estromal. Si se contrasta una
    proteína y se ajusta por una composición que esa proteína ayudó a estimar,
    parte del ajuste es la proteína ajustándose a sí misma.
    """
    g = str(gen).strip().upper()
    if g in firma.index:
        firma = firma.drop(index=g)
    return eje_epitelial(fracciones(matriz, columnas, firma), tipos_epiteliales)


def diagnosticar(fr: Fracciones, firma: pd.DataFrame) -> dict:
    """Comprueba si las proporciones tipo a tipo se pueden publicar.

    Tres comprobaciones, las tres sobre los datos reales.

    La primera es el residuo. Si la mezcla ajustada reproduce una parte pequeña
    del perfil medido, las proporciones no describen la muestra: son el mejor
    ajuste de un modelo que no encaja. Esta es la comprobación que primero hay
    que mirar y la que más veces se omite.

    La segunda es de colinealidad. Si dos columnas de la referencia están muy
    correlacionadas, los dos tipos no son separables y la proporción de cada uno
    no es un número estimable, aunque el método devuelva uno.

    La tercera es de colapso. Si un tipo sale en cero en buena parte de las
    muestras, no está siendo estimado en cero: está siendo descartado por el
    reparto entre tipos parecidos. Mirar solo el máximo no lo detecta, porque
    basta una muestra con señal para que el máximo no sea cero.
    """
    a = firma.loc[fr.genes_usados, firma.columns].to_numpy(dtype=float)
    a = _normalizar_columnas(a)
    tipos = list(fr.tabla.columns)
    pares = []
    for i in range(len(tipos)):
        for j in range(i + 1, len(tipos)):
            r = float(np.corrcoef(a[:, i], a[:, j])[0, 1])
            if np.isfinite(r) and r >= COLINEALIDAD_ALTA:
                pares.append(
                    {"tipo_a": tipos[i], "tipo_b": tipos[j], "correlacion": round(r, 3)}
                )
    pares.sort(key=lambda d: -d["correlacion"])

    n = int(fr.tabla.shape[0])
    ausente = (fr.tabla < FRACCION_DESPRECIABLE).sum(axis=0) / max(n, 1)
    colapsados = {
        str(t): round(float(ausente[t]), 3)
        for t in ausente.index[ausente > PROPORCION_DE_MUESTRAS_EN_CERO]
    }
    residuo = float(np.nanmedian(fr.residuo_relativo.to_numpy(dtype=float)))
    no_encaja = bool(np.isfinite(residuo) and residuo > RESIDUO_ALTO)

    usable = not pares and not colapsados and not no_encaja
    return {
        "residuo_relativo_mediano": round(residuo, 4),
        "umbral_de_residuo": RESIDUO_ALTO,
        "la_mezcla_no_explica_el_perfil": no_encaja,
        "pares_no_separables": pares,
        "umbral_de_colinealidad": COLINEALIDAD_ALTA,
        "tipos_ausentes_en_demasiadas_muestras": colapsados,
        "umbral_de_proporcion_en_cero": PROPORCION_DE_MUESTRAS_EN_CERO,
        "umbral_de_fraccion_despreciable": FRACCION_DESPRECIABLE,
        "los_umbrales_no_estaban_preregistrados": (
            f"Los tres umbrales de este diagnóstico se fijan aquí y quedan "
            f"declarados como enmienda {ENMIENDA_DE_LOS_UMBRALES}."
        ),
        "proporciones_por_tipo_publicables": bool(usable),
        "veredicto": (
            "Las proporciones tipo a tipo se pueden leer: la mezcla explica el "
            "perfil, ningún par de tipos está confundido y ninguno desaparece de "
            "la mayoría de las muestras."
            if usable
            else (
                "Las proporciones tipo a tipo **no** se publican. "
                + (
                    f"El residuo relativo mediano es {residuo:.2f}: la mezcla "
                    "ajustada reproduce una parte pequeña del perfil medido, así "
                    "que las proporciones no describen la muestra sino el mejor "
                    "ajuste de un modelo que no encaja. "
                    if no_encaja
                    else ""
                )
                + (
                    f"Hay {len(pares)} pares de tipos celulares cuyos perfiles de "
                    "referencia están tan correlacionados que no son separables, y "
                    "los mínimos cuadrados no negativos reparten la señal entre "
                    "ellos de forma inestable. "
                    if pares
                    else ""
                )
                + (
                    f"Y {len(colapsados)} tipos salen ausentes en más de la mitad "
                    "de las muestras, que es lo que ocurre cuando otro tipo con "
                    "perfil parecido se lleva su señal entera. "
                    if colapsados
                    else ""
                )
            )
        ),
        "como_se_lee": (
            "Un cero en una proporción no significa «ese tipo celular no está». "
            "Significa que el método ha asignado su señal a otro tipo con un "
            "perfil parecido. Leerlo como ausencia sería leer un artefacto del "
            "método como biología."
        ),
    }


def veredicto_del_eje_epitelial(
    epi: pd.Series,
    pares: list[tuple[str, str]],
    comparacion: dict,
    diagnostico: dict,
) -> dict:
    """Decide si el eje epitelial se puede usar como covariable, y como qué.

    Con un residuo alto no se puede sostener que el número sea una proporción de
    células. La pregunta que queda es más modesta y sí tiene respuesta: ¿ordena
    las muestras? Dos indicios lo dicen, y los dos salen de los datos. El primero
    es que el tumor quede por debajo de su propio normal en la mayoría de los
    pares, que es lo que la histología dice del colon. El segundo es que el eje
    coincida con el índice estromal, que es otra lectura del mismo eje.

    Ninguno de los dos demuestra que las proporciones sean correctas. Lo que
    sostienen es que el orden lleva señal, y para ajustar por una covariable el
    orden y la escala relativa es lo que se usa.
    """
    bajan = [
        1
        for a, b in pares
        if np.isfinite(epi.get(a, np.nan))
        and np.isfinite(epi.get(b, np.nan))
        and float(epi[a]) < float(epi[b])
    ]
    n = len([1 for a, b in pares if a in epi.index and b in epi.index])
    coincide = bool(comparacion.get("coincide_con_lo_esperado"))
    direccion_ok = bool(n > 0 and len(bajan) / n >= 0.75)
    return {
        "n_pares": n,
        "pares_con_menos_epitelio_en_el_tumor": len(bajan),
        "proporcion": round(len(bajan) / n, 3) if n else None,
        "coincide_con_el_indice_estromal": coincide,
        "usable_como_covariable": bool(direccion_ok and coincide),
        "no_es_una_proporcion_de_celulas": (
            "El número no se puede leer como «este tanto por ciento de las "
            f"células es epitelio». El residuo del ajuste es "
            f"{diagnostico['residuo_relativo_mediano']:.2f}, así que la mezcla no "
            "reproduce el perfil medido y el nivel absoluto no significa nada. Lo "
            "que se usa es el orden entre muestras."
        ),
        "por_que_se_puede_usar_igual": (
            f"El tumor tiene menos epitelio que su propio normal en "
            f"{len(bajan)} de {n} pares, que es lo que la histología dice del "
            "colon, y el eje coincide en sentido con el índice estromal. Dos "
            "indicios de que el orden lleva señal. Para ajustar por una "
            "covariable eso es lo que hace falta."
            if direccion_ok and coincide
            else (
                "No se usa como covariable: el eje no ordena los pares en el "
                "sentido que la histología espera, o no coincide con el índice "
                "estromal. Sin ninguno de los dos indicios, no hay motivo para "
                "creer que el orden lleve señal."
            )
        ),
    }


def por_que_no_en_tmt(
    matriz_log_ratio: pd.DataFrame,
    columnas: list[str],
    firma: pd.DataFrame,
) -> dict:
    """Comprueba sobre los datos que la deconvolución no vale en un log-cociente.

    No es una advertencia de manual. Es una comprobación: si el perfil medido
    tiene valores negativos, ninguna combinación con pesos no negativos de
    perfiles no negativos puede reproducirlos, porque todo lo que ``A x`` puede
    producir es no negativo. El método devolvería un número de todas formas, y
    ese número no sería una proporción de nada.
    """
    genes = sorted(set(matriz_log_ratio.index.astype(str)) & set(firma.index.astype(str)))
    if len(genes) < MIN_MARCADORES_MEDIDOS:
        raise ValueError(
            f"solo {len(genes)} genes marcadores están medidos en la matriz TMT; "
            "con menos la comprobación no es informativa"
        )
    cols = [c for c in columnas if c in matriz_log_ratio.columns]
    if not cols:
        raise FaltaColumna("ninguna de las columnas pedidas está en la matriz TMT")
    sub = matriz_log_ratio.loc[genes, cols].apply(pd.to_numeric, errors="coerce")
    v = sub.to_numpy(dtype=float)
    finitos = v[np.isfinite(v)]
    n_neg = int((finitos < 0).sum())
    frac_neg = float(n_neg / finitos.size) if finitos.size else float("nan")

    a = _normalizar_columnas(firma.loc[genes, firma.columns].to_numpy(dtype=float))
    forzados = []
    for c in cols[:10]:
        b = sub[c].to_numpy(dtype=float)
        b = np.where(np.isfinite(b), b, 0.0)
        nb = float(np.linalg.norm(b))
        if nb <= 0:
            continue
        _, resid = nnls(a, np.clip(b, 0.0, None))
        forzados.append(float(resid / nb))

    return {
        "matriz": "log-cociente contra el canal de referencia (CPTAC, TMT)",
        "n_valores_finitos": int(finitos.size),
        "n_valores_negativos": n_neg,
        "fraccion_negativa": round(frac_neg, 4),
        "primer_motivo": (
            f"El {frac_neg:.1%} de los valores medidos es negativo. La mezcla "
            "A·x con pesos no negativos y perfiles no negativos solo puede "
            "producir valores no negativos, así que no existe ninguna "
            "proporción que reproduzca este perfil. No es que el ajuste salga "
            "malo: es que el sistema no tiene solución con sentido."
        ),
        "segundo_motivo": (
            "El logaritmo no es lineal. Aunque se recortaran los negativos, el "
            "doble de fibroblastos no da el doble de log-intensidad, y la "
            "deconvolución supone exactamente eso."
        ),
        "tercer_motivo": (
            "Es un cociente contra un canal de referencia, y ese canal es a su "
            "vez una mezcla distinta en cada lote. Dividir por él ya ha mezclado "
            "la composición que se quiere estimar con la del canal."
        ),
        "residuo_relativo_si_se_fuerza": (
            [round(x, 3) for x in forzados] if forzados else None
        ),
        "que_se_hace_en_su_lugar": (
            "La deconvolución se corre solo en las cohortes LFQ, cuyas "
            "intensidades sí están en escala de abundancia. Para CPTAC, R4 se "
            "queda con el índice estromal y con el lenguaje obligatorio del "
            "preregistro: es una cota del efecto de composición, no una pureza "
            "medida."
        ),
        "por_que_se_comprueba_y_no_se_afirma": (
            "Afirmarlo sin mirar dejaría la duda de si se ha evitado el análisis "
            "por comodidad. El dato de arriba lo cierra."
        ),
    }


def comparar_con_indice(
    fraccion: pd.Series, indice: pd.Series, nombre_indice: str = "indice_estromal"
) -> dict:
    """Confronta el eje epitelial con el índice estromal del propio proteoma.

    Se espera que vayan en sentidos opuestos: más epitelio, menos matriz
    extracelular. Que coincidan apoya que ambos miden el mismo eje. Lo que **no**
    hace es validar ninguno de los dos, y el motivo hay que decirlo: los dos se
    calculan del mismo proteoma. Dos lecturas del mismo termómetro coinciden
    aunque el termómetro esté descalibrado. La validación de verdad necesita una
    medida de otra tecnología, y eso es lo que hizo D2 con la pureza ABSOLUTE.
    """
    comun = sorted(set(fraccion.dropna().index) & set(indice.dropna().index))
    if len(comun) < MIN_PARA_CORRELAR:
        return {
            "estado": "no_calculable",
            "motivo": (
                f"solo {len(comun)} muestras tienen las dos medidas (mínimo "
                f"{MIN_PARA_CORRELAR})"
            ),
            "n": len(comun),
        }
    x = fraccion.loc[comun].to_numpy(dtype=float)
    y = indice.loc[comun].to_numpy(dtype=float)
    rho, p = stats.spearmanr(x, y)
    return {
        "estado": "calculada",
        "n": len(comun),
        "indice": nombre_indice,
        "rho_spearman": round(float(rho), 4),
        "p": float(p),
        "signo_esperado": "negativo",
        "signo_observado": "negativo" if rho < 0 else "positivo",
        "coincide_con_lo_esperado": bool(rho < 0),
        "no_es_una_validacion": (
            "Las dos medidas salen del mismo proteoma, así que coincidir no "
            "demuestra que ninguna sea correcta: demuestra que son dos lecturas "
            "del mismo dato. Validar de verdad exige otra tecnología, como la "
            "pureza ABSOLUTE del ADN que se usó en D2."
        ),
    }


def recuperacion_en_mezclas_conocidas(
    firma: pd.DataFrame,
    proporciones: np.ndarray,
    semilla: int,
    ruido_log: float = 0.0,
) -> dict:
    """Mide si el método recupera proporciones que se conocen de antemano.

    Es la calibración del método, no un resultado sobre ningún tejido. Las
    mezclas se construyen a partir de la propia referencia, así que este es el
    caso más favorable posible: si aquí no recupera las proporciones, en un
    proteoma real tampoco. Los datos de esta función son sintéticos y por eso
    solo se usa desde ``tests/``; ninguna cifra suya entra en los resultados.
    """
    rng = np.random.default_rng(semilla)
    a = _normalizar_columnas(firma.to_numpy(dtype=float))
    errores = []
    for fila in np.atleast_2d(proporciones):
        x = np.asarray(fila, dtype=float)
        x = x / x.sum()
        b = a @ x
        if ruido_log > 0:
            b = b * np.exp(rng.normal(0.0, ruido_log, size=b.shape))
        b = b / b.sum()
        est, _ = nnls(a, b)
        est = est / (est.sum() or 1.0)
        errores.append(float(np.abs(est - x).sum() / 2.0))
    return {
        "n_mezclas": int(np.atleast_2d(proporciones).shape[0]),
        "ruido_log": ruido_log,
        "error_total_mediano": round(float(np.median(errores)), 4),
        "error_total_maximo": round(float(np.max(errores)), 4),
        "como_se_lee": (
            "El error total es la mitad de la suma de diferencias absolutas, así "
            "que va de 0 a 1: 0 es recuperación exacta y 1 es haber asignado "
            "toda la proporción a los tipos equivocados."
        ),
        "datos_sinteticos": True,
    }
