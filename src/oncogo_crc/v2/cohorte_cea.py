"""D9. Análisis de la cohorte de seguimiento con CEA seriado.

Este módulo es el código del protocolo de ``docs/protocolo_cohorte_cea.md``. El
estudio **no se ha hecho**: está pendiente de aprobación ética y de datos reales.
Lo que hay aquí es el análisis escrito y probado sobre una cohorte sintética, que
vive en ``tests/SINTETICO_cohorte_cea.csv`` y en ningún otro sitio.

La pregunta viene de D8. En ctDNA, medir una sola vez después de la cirugía da
una sensibilidad del 51 % y medir repetidamente da el 84 %. Si lo que importa es
la trayectoria y no el valor, al CEA debería pasarle lo mismo, y eso explicaría
por qué el CEA tiene mala fama midiéndolo como se mide: comparando cada valor
contra un umbral fijo de 5 ng/ml.

Así que se comparan tres modelos sobre los mismos datos:

1. **Umbral.** El CEA de esta visita pasa de 5 ng/ml o no. Es lo que se hace hoy.
2. **Valor.** El CEA de esta visita como número, sin umbral.
3. **Trayectoria.** El valor más cómo ha cambiado respecto a las visitas
   anteriores.

Dos decisiones de diseño que no son adornos:

- **La unidad de remuestreo es el paciente, no la visita.** Un paciente aporta
  varias visitas y están correlacionadas entre sí. Remuestrear visitas daría
  intervalos más estrechos de lo que corresponde, porque trataría como
  independientes medidas de la misma persona.
- **El valor de referencia de cada paciente es el suyo.** El cambio se mide
  contra las propias visitas anteriores de esa persona, no contra la media de la
  cohorte. Esa es toda la diferencia entre una trayectoria y un umbral.

Nada de este módulo se aplica a ninguna persona concreta. Ningún resultado de
esta tesis sirve para decidir sobre un paciente.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import expit

# Umbral clínico del CEA en ng/ml. Es el que usan los 8 estudios de D8.
UMBRAL_CLINICO = 5.0

# Ventana de la variable de resultado: ¿hay recaída en los meses siguientes a
# esta visita? Doce meses es el horizonte en el que una detección cambia el
# manejo; con una ventana mucho más corta casi ninguna visita tendría el suceso.
VENTANA_MESES = 12.0

# Cuántas visitas previas entran en la pendiente. Tres puntos es el mínimo para
# que una pendiente no sea simplemente la diferencia entre dos valores.
VISITAS_PARA_LA_PENDIENTE = 3

# Visitas mínimas por paciente para que entre en el análisis de trayectoria. Con
# una sola visita no hay trayectoria que medir.
VISITAS_MINIMAS = 2

COLUMNAS_OBLIGATORIAS = (
    "paciente",
    "visita_meses",
    "cea",
    "recaida",
    "mes_de_la_recaida",
    "ultimo_seguimiento_meses",
)


class FaltaColumna(KeyError):
    """La tabla de entrada no tiene el diccionario de datos del protocolo."""


def comprobar_entrada(d: pd.DataFrame) -> dict:
    """Comprueba que la tabla cumple el diccionario de datos del protocolo.

    Se comprueba antes de calcular nada. Un análisis que arranca con una columna
    ausente y la rellena con el valor por defecto produce resultados que parecen
    normales.
    """
    faltan = [c for c in COLUMNAS_OBLIGATORIAS if c not in d.columns]
    if faltan:
        raise FaltaColumna(
            f"faltan las columnas {faltan}; el diccionario de datos está en "
            "docs/protocolo_cohorte_cea.md"
        )
    problemas = []
    if (d["cea"] < 0).any():
        problemas.append("hay valores de CEA negativos")
    if (d["visita_meses"] < 0).any():
        problemas.append("hay visitas con tiempo negativo")

    # Una recaída declarada sin mes, o un mes de recaída en alguien sin recaída,
    # son los dos errores de captura más frecuentes en una base de seguimiento.
    por_paciente = d.groupby("paciente").agg(
        recaida=("recaida", "max"),
        mes=("mes_de_la_recaida", "max"),
    )
    sin_mes = por_paciente[(por_paciente["recaida"] == 1) & por_paciente["mes"].isna()]
    if not sin_mes.empty:
        problemas.append(
            f"{len(sin_mes)} pacientes con recaída declarada y sin mes de recaída"
        )
    mes_sobrante = por_paciente[
        (por_paciente["recaida"] == 0) & por_paciente["mes"].notna()
    ]
    if not mes_sobrante.empty:
        problemas.append(
            f"{len(mes_sobrante)} pacientes sin recaída y con mes de recaída"
        )
    pocas = por_paciente.index[
        d.groupby("paciente").size() < VISITAS_MINIMAS
    ].tolist()
    return {
        "n_filas": int(len(d)),
        "n_pacientes": int(d["paciente"].nunique()),
        "n_con_recaida": int(por_paciente["recaida"].sum()),
        "problemas": problemas,
        "coherente": not problemas,
        "pacientes_con_menos_de_dos_visitas": sorted(map(str, pocas)),
    }


def preparar_visitas(d: pd.DataFrame) -> pd.DataFrame:
    """Una fila por visita, con las variables de los tres modelos y el resultado.

    La variable de resultado de una visita es si el paciente recae en los
    ``VENTANA_MESES`` siguientes. Una visita que no llega a tener esos meses de
    seguimiento y en la que no ha pasado nada **no entra**: no se sabe si recayó
    después, y meterla como «no recayó» inventaría información. Es censura, y la
    forma de tratarla es excluir la visita, no adivinarla.
    """
    comprobar_entrada(d)
    t = d.sort_values(["paciente", "visita_meses"]).copy()
    t["log2_cea"] = np.log2(t["cea"].clip(lower=0.1))

    g = t.groupby("paciente", sort=False)
    t["visita_numero"] = g.cumcount() + 1
    t["log2_cea_anterior"] = g["log2_cea"].shift(1)
    t["meses_desde_la_anterior"] = g["visita_meses"].diff()
    # El cambio respecto a la visita anterior, en duplicaciones por año. En
    # duplicaciones porque el CEA sube de forma multiplicativa, y por año para
    # que dos visitas separadas seis meses y dos separadas tres no den el mismo
    # número por el mero hecho de estar más juntas.
    t["cambio_anual"] = (
        (t["log2_cea"] - t["log2_cea_anterior"])
        / t["meses_desde_la_anterior"].replace(0, np.nan)
        * 12.0
    )
    t["pendiente"] = g["log2_cea"].transform(
        lambda s: s.rolling(VISITAS_PARA_LA_PENDIENTE, min_periods=2).apply(
            _pendiente_local, raw=True
        )
    )
    # El valor de referencia de cada paciente es el suyo: la primera visita.
    t["log2_cea_basal"] = g["log2_cea"].transform("first")
    t["cambio_sobre_el_basal"] = t["log2_cea"] - t["log2_cea_basal"]

    t["supera_el_umbral"] = (t["cea"] > UMBRAL_CLINICO).astype(int)

    recae = t["recaida"] == 1
    hasta_la_recaida = t["mes_de_la_recaida"] - t["visita_meses"]
    t["recae_en_la_ventana"] = (
        recae & (hasta_la_recaida > 0) & (hasta_la_recaida <= VENTANA_MESES)
    ).astype(int)

    # Visitas evaluables: o el suceso ocurre dentro de la ventana, o hay
    # seguimiento suficiente para afirmar que no ocurrió.
    seguimiento_restante = t["ultimo_seguimiento_meses"] - t["visita_meses"]
    suficiente = seguimiento_restante >= VENTANA_MESES
    ya_recayo = recae & (hasta_la_recaida <= 0)
    t["evaluable"] = (
        (t["recae_en_la_ventana"] == 1) | (suficiente & ~ya_recayo)
    ) & ~ya_recayo
    return t


def _pendiente_local(valores: np.ndarray) -> float:
    """Pendiente de una recta por mínimos cuadrados sobre los últimos valores.

    Se calcula sobre el índice de la visita y no sobre el tiempo porque la
    ventana rodante de pandas no pasa las dos series a la vez. La diferencia la
    recoge ``cambio_anual``, que sí usa el tiempo; esta variable mide la forma de
    la trayectoria y no su velocidad.
    """
    y = np.asarray(valores, dtype=float)
    y = y[np.isfinite(y)]
    if len(y) < 2:
        return np.nan
    x = np.arange(len(y), dtype=float)
    return float(np.polyfit(x, y, 1)[0])


MODELOS = {
    "umbral": ["supera_el_umbral"],
    "valor": ["log2_cea"],
    "trayectoria": ["log2_cea", "cambio_anual", "pendiente", "cambio_sobre_el_basal"],
}


def _ajustar_logistica(
    X: np.ndarray, y: np.ndarray, *, max_iter: int = 200, tol: float = 1e-8
) -> np.ndarray:
    """Regresión logística por Newton-Raphson con una cresta mínima.

    La cresta es de 1e-6 y está solo para que la matriz sea invertible cuando una
    variable es casi constante. No es regularización con intención: con ese valor
    no mueve los coeficientes de forma apreciable, y evita que el ajuste falle
    por una columna degenerada en un remuestreo.
    """
    n, p = X.shape
    beta = np.zeros(p)
    for _ in range(max_iter):
        eta = X @ beta
        mu = expit(eta)
        w = np.clip(mu * (1 - mu), 1e-10, None)
        gradiente = X.T @ (y - mu) - 1e-6 * beta
        hessiana = (X * w[:, None]).T @ X + 1e-6 * np.eye(p)
        paso = np.linalg.solve(hessiana, gradiente)
        beta = beta + paso
        if np.max(np.abs(paso)) < tol:
            break
    return beta


def _auc(y: np.ndarray, puntuacion: np.ndarray) -> float:
    """AUC por la U de Mann-Whitney, con los empates valiendo medio.

    Los empates importan aquí más que de costumbre: el modelo del umbral solo
    produce dos puntuaciones distintas, así que casi todo son empates. Contarlos
    como aciertos infla su AUC y no contarlos lo hunde.
    """
    y = np.asarray(y)
    pos = puntuacion[y == 1]
    neg = puntuacion[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    orden = np.argsort(np.concatenate([pos, neg]), kind="mergesort")
    rangos = np.empty(len(orden), dtype=float)
    valores = np.concatenate([pos, neg])[orden]
    i = 0
    while i < len(valores):
        j = i
        while j + 1 < len(valores) and valores[j + 1] == valores[i]:
            j += 1
        rangos[i : j + 1] = (i + j) / 2.0 + 1.0
        i = j + 1
    rangos_originales = np.empty(len(orden), dtype=float)
    rangos_originales[orden] = rangos
    suma_pos = rangos_originales[: len(pos)].sum()
    return float(
        (suma_pos - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))
    )


def _diseno(t: pd.DataFrame, variables: list[str]) -> tuple[np.ndarray, np.ndarray]:
    sub = t[variables + ["recae_en_la_ventana"]].dropna()
    X = np.column_stack([np.ones(len(sub))] + [sub[v].to_numpy() for v in variables])
    return X, sub["recae_en_la_ventana"].to_numpy(dtype=float)


def _auc_por_validacion_cruzada(
    t: pd.DataFrame, variables: list[str], *, semilla: int, pliegues: int = 5
) -> float:
    """AUC fuera de muestra, con los pliegues partidos por paciente.

    Partir por visita dejaría visitas del mismo paciente a los dos lados, y el
    modelo vería en entrenamiento a la misma persona que luego evalúa. El AUC
    saldría más alto y la causa sería la fuga, no el modelo.
    """
    sub = t[variables + ["recae_en_la_ventana", "paciente"]].dropna()
    if sub.empty or sub["recae_en_la_ventana"].nunique() < 2:
        return np.nan
    pacientes = sub["paciente"].unique()
    rng = np.random.default_rng(semilla)
    asignacion = dict(
        zip(pacientes, rng.integers(0, pliegues, len(pacientes)), strict=True)
    )
    pliegue_de = sub["paciente"].map(asignacion).to_numpy()

    predicciones = np.full(len(sub), np.nan)
    for p in range(pliegues):
        fuera = pliegue_de == p
        dentro = ~fuera
        if not fuera.any() or sub["recae_en_la_ventana"][dentro].nunique() < 2:
            continue
        X = np.column_stack(
            [np.ones(len(sub))] + [sub[v].to_numpy() for v in variables]
        )
        beta = _ajustar_logistica(
            X[dentro], sub["recae_en_la_ventana"].to_numpy(dtype=float)[dentro]
        )
        predicciones[fuera] = X[fuera] @ beta
    ok = np.isfinite(predicciones)
    if ok.sum() == 0:
        return np.nan
    return _auc(
        sub["recae_en_la_ventana"].to_numpy()[ok], predicciones[ok]
    )


def comparar_modelos(
    t: pd.DataFrame, *, semilla: int, n_bootstrap: int = 2000
) -> dict:
    """Compara los tres modelos por AUC, con un bootstrap agrupado por paciente.

    La comparación es de diferencias pareadas y no de intervalos que se solapan.
    En cada réplica se remuestrean pacientes, se recalculan los tres AUC sobre
    los mismos pacientes y se guarda la diferencia. El intervalo de la diferencia
    es lo que decide, porque dos intervalos que se solapan pueden corresponder a
    una diferencia que no incluye el cero.
    """
    evaluables = t[t["evaluable"]].copy()
    if evaluables.empty or evaluables["recae_en_la_ventana"].nunique() < 2:
        return {
            "estado": "no_calculable",
            "motivo": (
                "no hay visitas evaluables con las dos clases del resultado"
            ),
        }

    auc_observado = {
        nombre: _auc_por_validacion_cruzada(
            evaluables, variables, semilla=semilla + i
        )
        for i, (nombre, variables) in enumerate(MODELOS.items())
    }

    pacientes = evaluables["paciente"].unique()
    rng = np.random.default_rng(semilla)
    diferencias: dict[str, list[float]] = {
        "trayectoria_menos_umbral": [],
        "trayectoria_menos_valor": [],
        "valor_menos_umbral": [],
    }
    fallos = 0
    for b in range(n_bootstrap):
        elegidos = rng.choice(pacientes, size=len(pacientes), replace=True)
        # Las visitas de un paciente elegido dos veces entran dos veces, que es
        # lo que hace que el bootstrap agrupado respete la correlación interna.
        trozos = [evaluables[evaluables["paciente"] == p] for p in elegidos]
        replica = pd.concat(trozos, ignore_index=True)
        if replica["recae_en_la_ventana"].nunique() < 2:
            fallos += 1
            continue
        aucs = {
            nombre: _auc_por_validacion_cruzada(
                replica, variables, semilla=semilla + 10_000 + b
            )
            for nombre, variables in MODELOS.items()
        }
        if any(not np.isfinite(v) for v in aucs.values()):
            fallos += 1
            continue
        diferencias["trayectoria_menos_umbral"].append(
            aucs["trayectoria"] - aucs["umbral"]
        )
        diferencias["trayectoria_menos_valor"].append(
            aucs["trayectoria"] - aucs["valor"]
        )
        diferencias["valor_menos_umbral"].append(aucs["valor"] - aucs["umbral"])

    resumen = {}
    for nombre, valores in diferencias.items():
        if len(valores) < 100:
            resumen[nombre] = {
                "estado": "no_calculable",
                "n_replicas_utiles": len(valores),
                "motivo": "menos de 100 réplicas útiles",
            }
            continue
        v = np.asarray(valores)
        bajo, alto = np.percentile(v, [2.5, 97.5])
        resumen[nombre] = {
            "diferencia_mediana": float(np.median(v)),
            "ci_low": float(bajo),
            "ci_high": float(alto),
            "mejora": bool(bajo > 0),
            "n_replicas_utiles": int(len(v)),
        }
    return {
        "estado": "calculado",
        "auc_por_modelo": {k: (None if not np.isfinite(v) else float(v)) for k, v in auc_observado.items()},
        "diferencias_pareadas": resumen,
        "n_visitas_evaluables": int(len(evaluables)),
        "n_visitas_excluidas_por_censura": int((~t["evaluable"]).sum()),
        "n_pacientes": int(len(pacientes)),
        "n_replicas_pedidas": n_bootstrap,
        "n_replicas_fallidas": fallos,
        "unidad_de_remuestreo": "paciente",
        "por_que_el_paciente": (
            "Un paciente aporta varias visitas y están correlacionadas. "
            "Remuestrear visitas daría intervalos más estrechos de lo que "
            "corresponde, porque trataría como independientes medidas de la "
            "misma persona."
        ),
        "variables_de_cada_modelo": MODELOS,
        "ventana_del_resultado_meses": VENTANA_MESES,
        "umbral_clinico_ng_ml": UMBRAL_CLINICO,
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
    }


def tamano_necesario(
    *,
    auc_umbral: float,
    auc_trayectoria: float,
    proporcion_con_recaida: float,
    visitas_por_paciente: float,
    potencia: float = 0.80,
    alfa: float = 0.05,
) -> dict:
    """Pacientes necesarios para detectar una diferencia de AUC dada.

    Usa la aproximación de Hanley y McNeil para la varianza de un AUC y la
    corrige por el número de visitas por paciente mediante el efecto de diseño
    de un conglomerado. La correlación dentro del paciente se fija en 0,5, que es
    un valor conservador: cuanto más se parecen las visitas de una persona, menos
    información aporta cada una.

    Es una aproximación y se declara como tal: el cálculo exacto exigiría
    simular, y lo que hace falta aquí es el orden de magnitud para pedir el
    estudio, no la cifra al paciente.
    """
    from scipy.stats import norm

    rho_intra = 0.5
    efecto_diseno = 1 + (visitas_por_paciente - 1) * rho_intra

    def var_auc(a, n_pos, n_neg):
        q1 = a / (2 - a)
        q2 = 2 * a**2 / (1 + a)
        return (
            a * (1 - a)
            + (n_pos - 1) * (q1 - a**2)
            + (n_neg - 1) * (q2 - a**2)
        ) / (n_pos * n_neg)

    z_alfa = norm.ppf(1 - alfa / 2)
    z_beta = norm.ppf(potencia)
    delta = auc_trayectoria - auc_umbral
    if delta <= 0:
        return {
            "estado": "no_calculable",
            "motivo": "la diferencia de AUC que se quiere detectar no es positiva",
        }

    # Se busca el n de pacientes que hace que la diferencia sea detectable. Las
    # dos AUC se calculan sobre las mismas visitas, así que la varianza de la
    # diferencia se reduce por su correlación; se toma 0,6, que es la habitual
    # entre dos modelos anidados.
    rho_entre_auc = 0.6
    for n_pacientes in range(20, 20001, 10):
        visitas = n_pacientes * visitas_por_paciente / efecto_diseno
        n_pos = visitas * proporcion_con_recaida
        n_neg = visitas - n_pos
        if n_pos < 5 or n_neg < 5:
            continue
        var = (
            var_auc(auc_umbral, n_pos, n_neg)
            + var_auc(auc_trayectoria, n_pos, n_neg)
            - 2
            * rho_entre_auc
            * np.sqrt(
                var_auc(auc_umbral, n_pos, n_neg)
                * var_auc(auc_trayectoria, n_pos, n_neg)
            )
        )
        if delta >= (z_alfa + z_beta) * np.sqrt(max(var, 1e-12)):
            return {
                "estado": "calculado",
                "n_pacientes": int(n_pacientes),
                "n_visitas_efectivas": float(visitas),
                "diferencia_de_auc_objetivo": float(delta),
                "potencia": potencia,
                "alfa": alfa,
                "efecto_diseno": float(efecto_diseno),
                "correlacion_intrapaciente_supuesta": rho_intra,
                "correlacion_entre_las_dos_auc_supuesta": rho_entre_auc,
                "es_una_aproximacion": (
                    "Varianza de Hanley y McNeil corregida por efecto de "
                    "diseño. El cálculo exacto exigiría simular; esto da el "
                    "orden de magnitud para pedir el estudio."
                ),
            }
    return {
        "estado": "no_calculable",
        "motivo": "harían falta más de 20 000 pacientes con esos supuestos",
    }
