"""Meta-análisis bivariado de sensibilidad y especificidad (corrección D8).

Sustituye a ``mada::reitsma`` y a ``metafor`` porque no hay R instalado (enmienda
E2). Implementa el modelo de Reitsma, que es el estándar para combinar pruebas
diagnósticas:

    y_i = (logit sens_i, logit espec_i)'
    y_i ~ N(mu, Sigma + S_i)

donde ``S_i`` es la varianza dentro del estudio ``i``, que sale directamente del
binomio, y ``Sigma`` es la varianza entre estudios, con su correlación.

Por qué bivariado y no dos meta-análisis separados. La sensibilidad y la
especificidad de una prueba no se mueven independientemente: un umbral más bajo
sube una y baja la otra. Si se combinan por separado, esa correlación se pierde y
el resultado puede ser un punto que ningún umbral real produce. El modelo
bivariado la estima, y el signo de esa correlación es un resultado por sí mismo:
negativa significa que los estudios se diferencian sobre todo por dónde pusieron
el corte.

Dos decisiones que condicionan los números y que conviene tener delante:

1. La varianza dentro del estudio se trata como conocida. Es la aproximación
   normal de Reitsma 2005, no el binomio exacto de Chu y Cole 2006. Con celdas
   pequeñas la aproximación empeora, y por eso se avisa cuando alguna celda es
   menor que cinco.
2. El intervalo del punto combinado sale de la covarianza de mínimos cuadrados
   generalizados, que trata ``Sigma`` como si fuera conocida. Eso lo estrecha un
   poco. Es lo que hace ``mada`` por defecto y se declara en lugar de disimularse.

Referencias: Reitsma JB et al. 2005 (Journal of Clinical Epidemiology 58:982);
Chu H y Cole SR 2006 (Journal of Clinical Epidemiology 59:1331) para la versión
binomial exacta; Harbord RM et al. 2007 para la equivalencia con HSROC.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logit
from scipy.stats import norm

# Corrección de continuidad. Se aplica a las cuatro celdas y solo cuando alguna
# es cero, porque con un cero el logit no existe. Añadirla siempre sesgaría hacia
# el centro todos los estudios, también los que no la necesitan.
CORRECCION = 0.5

# Por debajo de esta cuenta en alguna celda, la aproximación normal del logit
# deja de ser buena y el resultado se marca.
CELDA_PEQUENA = 5

MAX_ITER = 2000


@dataclass(frozen=True)
class ResultadoReitsma:
    k: int
    sensibilidad: float
    sensibilidad_ci_low: float
    sensibilidad_ci_high: float
    especificidad: float
    especificidad_ci_low: float
    especificidad_ci_high: float
    logit_sensibilidad: float
    logit_especificidad: float
    se_logit_sensibilidad: float
    se_logit_especificidad: float
    tau_sensibilidad: float
    tau_especificidad: float
    correlacion: float
    dor: float
    dor_ci_low: float
    dor_ci_high: float
    convergio: bool
    n_estudios_con_correccion: int
    n_estudios_con_celda_pequena: int
    metodo: str

    def to_dict(self) -> dict:
        return asdict(self)


def _a_logit(tp, fp, fn, tn) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int]:
    """Pasa las tablas 2x2 a logits con sus varianzas.

    La varianza del logit de una proporción p = a/(a+b) es 1/a + 1/b. Para la
    sensibilidad las dos celdas son TP y FN; para la especificidad, TN y FP. No
    hay covarianza dentro del estudio porque los enfermos y los sanos son
    personas distintas: una no informa de la otra.
    """
    tp = np.asarray(tp, dtype=float)
    fp = np.asarray(fp, dtype=float)
    fn = np.asarray(fn, dtype=float)
    tn = np.asarray(tn, dtype=float)
    celdas = np.vstack([tp, fp, fn, tn])

    pequena = int((celdas.min(axis=0) < CELDA_PEQUENA).sum())
    necesita = celdas.min(axis=0) == 0
    n_corregidos = int(necesita.sum())
    ajuste = np.where(necesita, CORRECCION, 0.0)
    tp, fp, fn, tn = (c + ajuste for c in (tp, fp, fn, tn))

    y_sens = logit(tp / (tp + fn))
    y_espec = logit(tn / (tn + fp))
    v_sens = 1.0 / tp + 1.0 / fn
    v_espec = 1.0 / tn + 1.0 / fp
    return y_sens, y_espec, v_sens, v_espec, n_corregidos, pequena


def _gls(y: np.ndarray, sigmas: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Media y covarianza por mínimos cuadrados generalizados.

    ``y`` es (k, 2) y ``sigmas`` es (k, 2, 2) con la covarianza total de cada
    estudio. Devuelve mu (2,) y su covarianza (2, 2).
    """
    acum = np.zeros((2, 2))
    peso_y = np.zeros(2)
    for yi, si in zip(y, sigmas, strict=True):
        inv = np.linalg.inv(si)
        acum += inv
        peso_y += inv @ yi
    cov = np.linalg.inv(acum)
    return cov @ peso_y, cov


def _menos_log_verosimilitud(
    parametros: np.ndarray, y: np.ndarray, v: np.ndarray
) -> float:
    """Verosimilitud con mu perfilada fuera.

    Los parámetros llegan sin restricciones: los dos primeros son logaritmos de
    las desviaciones típicas entre estudios y el tercero es la arcotangente
    hiperbólica de la correlación. Así el optimizador no tiene que respetar que
    una varianza sea positiva ni que una correlación esté entre −1 y 1, y nunca
    se sale del espacio válido.
    """
    log_t1, log_t2, z_rho = parametros
    t1, t2 = np.exp(log_t1), np.exp(log_t2)
    rho = np.tanh(z_rho)
    entre = np.array([[t1**2, rho * t1 * t2], [rho * t1 * t2, t2**2]])

    sigmas = np.array([entre + np.diag(vi) for vi in v])
    mu, _ = _gls(y, sigmas)

    total = 0.0
    for yi, si in zip(y, sigmas, strict=True):
        signo, logdet = np.linalg.slogdet(si)
        if signo <= 0:
            return np.inf
        d = yi - mu
        total += logdet + float(d @ np.linalg.solve(si, d))
    return 0.5 * total


def reitsma(tp, fp, fn, tn) -> ResultadoReitsma:
    """Ajusta el modelo bivariado a un conjunto de tablas 2x2.

    Hacen falta al menos cuatro estudios. Con menos, los cinco parámetros del
    modelo no están identificados y lo que saldría sería ruido con intervalos de
    aspecto respetable, que es peor que no dar nada.
    """
    y_sens, y_espec, v_sens, v_espec, n_corr, n_peq = _a_logit(tp, fp, fn, tn)
    k = len(y_sens)
    if k < 4:
        raise ValueError(
            f"el modelo bivariado tiene cinco parámetros y hay {k} estudios; "
            "hacen falta al menos 4"
        )

    y = np.column_stack([y_sens, y_espec])
    v = np.column_stack([v_sens, v_espec])

    # El arranque: las desviaciones observadas de los logits, y correlación cero.
    inicio = np.array(
        [
            np.log(max(y_sens.std(ddof=1), 1e-3)),
            np.log(max(y_espec.std(ddof=1), 1e-3)),
            0.0,
        ]
    )
    opt = minimize(
        _menos_log_verosimilitud,
        inicio,
        args=(y, v),
        method="Nelder-Mead",
        options={"maxiter": MAX_ITER, "xatol": 1e-8, "fatol": 1e-10},
    )
    log_t1, log_t2, z_rho = opt.x
    t1, t2 = float(np.exp(log_t1)), float(np.exp(log_t2))
    rho = float(np.tanh(z_rho))
    entre = np.array([[t1**2, rho * t1 * t2], [rho * t1 * t2, t2**2]])
    sigmas = np.array([entre + np.diag(vi) for vi in v])
    mu, cov = _gls(y, sigmas)

    z = norm.ppf(0.975)
    se1, se2 = float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))

    # El cociente de odds diagnóstico es el producto de los odds de los dos
    # aciertos, así que en escala logarítmica es la suma de los dos logits y su
    # varianza incluye la covarianza entre ellos. Olvidarla daría un intervalo
    # equivocado en la dirección que depende del signo de la correlación.
    log_dor = float(mu[0] + mu[1])
    var_log_dor = float(cov[0, 0] + cov[1, 1] + 2 * cov[0, 1])
    se_log_dor = float(np.sqrt(max(var_log_dor, 0.0)))

    return ResultadoReitsma(
        k=k,
        sensibilidad=float(expit(mu[0])),
        sensibilidad_ci_low=float(expit(mu[0] - z * se1)),
        sensibilidad_ci_high=float(expit(mu[0] + z * se1)),
        especificidad=float(expit(mu[1])),
        especificidad_ci_low=float(expit(mu[1] - z * se2)),
        especificidad_ci_high=float(expit(mu[1] + z * se2)),
        logit_sensibilidad=float(mu[0]),
        logit_especificidad=float(mu[1]),
        se_logit_sensibilidad=se1,
        se_logit_especificidad=se2,
        tau_sensibilidad=t1,
        tau_especificidad=t2,
        correlacion=rho,
        dor=float(np.exp(log_dor)),
        dor_ci_low=float(np.exp(log_dor - z * se_log_dor)),
        dor_ci_high=float(np.exp(log_dor + z * se_log_dor)),
        convergio=bool(opt.success),
        n_estudios_con_correccion=n_corr,
        n_estudios_con_celda_pequena=n_peq,
        metodo=(
            "reitsma_bivariado_aproximacion_normal; mu por GLS; Sigma por "
            "maxima verosimilitud con mu perfilada"
        ),
    )


def diferencia(a: ResultadoReitsma, b: ResultadoReitsma, *, comparable: bool) -> dict:
    """Diferencia entre dos pruebas, en escala logit y con su intervalo.

    El argumento ``comparable`` no tiene valor por defecto a propósito. Dos
    pruebas evaluadas en indicaciones distintas —una para cribado y otra para
    detectar recaída— no se pueden comparar aunque los números se resten sin
    protestar, y quien llame a esta función tiene que decidirlo y dejarlo
    escrito. Si no son comparables, se devuelven los dos puntos y ninguna
    diferencia.
    """
    base = {
        "comparable": comparable,
        "a": {
            "sensibilidad": a.sensibilidad,
            "especificidad": a.especificidad,
            "k": a.k,
        },
        "b": {
            "sensibilidad": b.sensibilidad,
            "especificidad": b.especificidad,
            "k": b.k,
        },
    }
    if not comparable:
        base["motivo"] = (
            "las dos estimaciones no responden a la misma pregunta clínica, así "
            "que restarlas no significa nada. Se dan los dos puntos por separado."
        )
        return base

    z = norm.ppf(0.975)
    for nombre, ma, sa, mb, sb in (
        (
            "sensibilidad",
            a.logit_sensibilidad,
            a.se_logit_sensibilidad,
            b.logit_sensibilidad,
            b.se_logit_sensibilidad,
        ),
        (
            "especificidad",
            a.logit_especificidad,
            a.se_logit_especificidad,
            b.logit_especificidad,
            b.se_logit_especificidad,
        ),
    ):
        d = ma - mb
        se = float(np.sqrt(sa**2 + sb**2))
        base[f"diferencia_logit_{nombre}"] = {
            "diferencia": float(d),
            "se": se,
            "ci_low": float(d - z * se),
            "ci_high": float(d + z * se),
            # El veredicto sale del intervalo de la diferencia, no de si los dos
            # intervalos originales se solapan, que no es lo mismo.
            "difieren": bool((d - z * se) * (d + z * se) > 0),
        }
    return base


def comparar_con_publicado(
    propio: ResultadoReitsma,
    *,
    sensibilidad: float,
    sensibilidad_ci: tuple[float, float],
    especificidad: float,
    especificidad_ci: tuple[float, float],
    fuente: str,
    k_publicado: int,
) -> dict:
    """Pone el ajuste propio al lado de una estimación combinada ya publicada.

    No es un contraste y no se presenta como tal. Cuando de una prueba no hay
    tablas 2x2 por estudio y sí un resumen publicado, lo único honesto es poner
    los dos puntos con sus intervalos y decir que se solapan o no. Comparar un
    ajuste propio con el ajuste de otro grupo arrastra todas las diferencias de
    criterio de inclusión, y eso no entra en ningún intervalo.
    """
    solapan_sens = (
        propio.sensibilidad_ci_low <= sensibilidad_ci[1]
        and sensibilidad_ci[0] <= propio.sensibilidad_ci_high
    )
    solapan_espec = (
        propio.especificidad_ci_low <= especificidad_ci[1]
        and especificidad_ci[0] <= propio.especificidad_ci_high
    )
    return {
        "propio": {
            "k": propio.k,
            "sensibilidad": propio.sensibilidad,
            "sensibilidad_ci": [
                propio.sensibilidad_ci_low,
                propio.sensibilidad_ci_high,
            ],
            "especificidad": propio.especificidad,
            "especificidad_ci": [
                propio.especificidad_ci_low,
                propio.especificidad_ci_high,
            ],
        },
        "publicado": {
            "fuente": fuente,
            "k": k_publicado,
            "sensibilidad": sensibilidad,
            "sensibilidad_ci": list(sensibilidad_ci),
            "especificidad": especificidad,
            "especificidad_ci": list(especificidad_ci),
        },
        "intervalos_de_sensibilidad_se_solapan": bool(solapan_sens),
        "intervalos_de_especificidad_se_solapan": bool(solapan_espec),
        "no_es_un_contraste": (
            "Son dos ajustes de grupos distintos, con criterios de inclusión "
            "distintos. Que los intervalos se solapen o no es una descripción, "
            "no una prueba estadística."
        ),
    }
