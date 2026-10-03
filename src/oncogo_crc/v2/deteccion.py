"""Ventana de detección tejido→plasma (corrección C8).

Qué estaba mal en v1. Se midió la concordancia de dirección entre tejido y suero
sobre todas las proteínas que aparecían en las dos matrices y se concluyó que la
dirección no se transfiere. Pero el suero tiene un rango dinámico de diez
órdenes de magnitud y la espectrometría sin enriquecimiento solo alcanza la
parte alta. Una proteína de tejido que circula a picogramos por litro no es que
vaya en dirección contraria: es que no se mide. Mezclar «no se detecta» con «no
está» convierte un límite del método en una conclusión biológica.

Qué se hace. Con las concentraciones en sangre del Human Protein Atlas se
estima la probabilidad de detección en función de la concentración, se localiza
el punto donde esa probabilidad cruza la mitad, y la concordancia de signo se
mide **solo por encima** de ese punto. El resultado deja de ser «la dirección no
se transfiere» y pasa a ser «la dirección no se transfiere entre las proteínas
que el método alcanza», que es una afirmación mucho más estrecha y además
comprobable.

El modelo es una regresión logística de una sola variable, el logaritmo decimal
de la concentración. No se ajusta nada más: con unos pocos cientos de proteínas
con concentración conocida, cualquier modelo mayor sobreajusta.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from oncogo_crc.v2.metricas import Fraccion


def _logistica(x: np.ndarray, y: np.ndarray, *, max_iter: int = 200, tol: float = 1e-9):
    """Newton-Raphson para una logística con intercepto y una pendiente.

    Se implementa a mano para no arrastrar un ajustador con regularización por
    defecto: aquí la pendiente y el punto de corte son el resultado, así que una
    penalización silenciosa los sesgaría.
    """
    X = np.column_stack([np.ones_like(x), x])
    beta = np.zeros(2)
    for _ in range(max_iter):
        eta = X @ beta
        p = 1.0 / (1.0 + np.exp(-np.clip(eta, -500, 500)))
        w = np.clip(p * (1.0 - p), 1e-10, None)
        gradiente = X.T @ (y - p)
        hessiana = X.T @ (X * w[:, None])
        try:
            paso = np.linalg.solve(hessiana, gradiente)
        except np.linalg.LinAlgError:
            return None
        beta = beta + paso
        if np.max(np.abs(paso)) < tol:
            break
    else:
        return None
    eta = X @ beta
    p = 1.0 / (1.0 + np.exp(-np.clip(eta, -500, 500)))
    w = np.clip(p * (1.0 - p), 1e-10, None)
    try:
        cov = np.linalg.inv(X.T @ (X * w[:, None]))
    except np.linalg.LinAlgError:
        return None
    return beta, np.sqrt(np.diag(cov))


def curva_de_deteccion(
    concentracion_pg_l: pd.Series,
    detectada: pd.Series,
    *,
    n_min: int,
) -> dict:
    """Probabilidad de detección frente a log10 de la concentración en sangre.

    ``concentracion_pg_l`` y ``detectada`` comparten índice. Solo entran las
    proteínas con concentración conocida: las demás no informan de la curva.
    """
    comun = concentracion_pg_l.dropna().index.intersection(detectada.dropna().index)
    conc = pd.to_numeric(concentracion_pg_l.loc[comun], errors="coerce")
    det = detectada.loc[comun].astype(bool)
    ok = conc.notna() & (conc > 0)
    conc, det = conc[ok], det[ok]
    n = int(conc.size)
    if n < n_min or det.nunique() < 2:
        return {
            "estado": "no_evaluable",
            "motivo": (
                f"{n} proteínas con concentración conocida y detección informativa "
                f"no llegan al mínimo de {n_min}, o todas caen del mismo lado."
            ),
            "n": n,
        }
    x = np.log10(conc.to_numpy(dtype=float))
    y = det.to_numpy(dtype=float)
    ajuste = _logistica(x, y)
    if ajuste is None:
        return {"estado": "no_evaluable", "motivo": "la logística no converge", "n": n}
    beta, se = ajuste
    # El umbral es la concentración donde la probabilidad cruza la mitad.
    umbral_log10 = None if beta[1] == 0 else float(-beta[0] / beta[1])
    z = float(beta[1] / se[1]) if se[1] > 0 else None
    return {
        "estado": "calculado",
        "n": n,
        "n_detectadas": int(det.sum()),
        "intercepto": float(beta[0]),
        "pendiente_por_orden_de_magnitud": float(beta[1]),
        "se_pendiente": float(se[1]),
        "z_pendiente": z,
        "umbral_log10_pg_l": umbral_log10,
        "umbral_pg_l": None if umbral_log10 is None else float(10.0**umbral_log10),
        "rango_log10_observado": [float(x.min()), float(x.max())],
        "interpretacion": (
            "Una pendiente positiva y grande significa que la detección depende de la "
            "concentración, es decir que el límite es del método. El umbral es la "
            "concentración a la que la proteína tiene la mitad de probabilidad de "
            "medirse."
        ),
    }


def concordancia_dentro_de_la_ventana(
    lfc_tejido: pd.Series,
    lfc_plasma: pd.Series,
    concentracion_pg_l: pd.Series,
    *,
    umbral_log10: float | None,
    n_min: int,
) -> dict:
    """Concordancia de signo por encima y por debajo del umbral de detección."""
    idx = lfc_tejido.dropna().index.intersection(lfc_plasma.dropna().index)
    tejido = pd.to_numeric(lfc_tejido.loc[idx], errors="coerce")
    plasma = pd.to_numeric(lfc_plasma.loc[idx], errors="coerce")
    conc = pd.to_numeric(concentracion_pg_l.reindex(idx), errors="coerce")
    valido = tejido.notna() & plasma.notna() & (tejido != 0) & (plasma != 0)
    tejido, plasma, conc = tejido[valido], plasma[valido], conc[valido]

    def _fraccion(mascara: pd.Series, nombre: str, nota: str) -> dict:
        t, p = tejido[mascara], plasma[mascara]
        f = Fraccion(
            nombre=nombre,
            numerador=int((np.sign(t) == np.sign(p)).sum()),
            denominador=int(t.size),
            unidad_denominador="proteinas_unicas",
            nota=nota,
        )
        d = f.to_dict()
        d["suficiente"] = bool(t.size >= n_min)
        return d

    todas = _fraccion(
        pd.Series(True, index=tejido.index),
        "concordancia_tejido_plasma_todas",
        "Sin filtrar por concentración. Es el cálculo de v1.",
    )
    salida = {
        "n_con_las_dos_medidas": int(tejido.size),
        "n_con_concentracion_conocida": int(conc.notna().sum()),
        "todas": todas,
        "umbral_log10_pg_l": umbral_log10,
    }
    if umbral_log10 is None:
        salida["estado"] = "no_evaluable"
        salida["motivo"] = "sin umbral de detección no se puede definir la ventana"
        return salida
    arriba = conc.notna() & (np.log10(conc.where(conc > 0)) >= umbral_log10)
    abajo = conc.notna() & (np.log10(conc.where(conc > 0)) < umbral_log10)
    salida["dentro_de_la_ventana"] = _fraccion(
        arriba,
        "concordancia_tejido_plasma_dentro_de_la_ventana",
        "Solo proteínas cuya concentración en sangre está por encima del umbral.",
    )
    salida["fuera_de_la_ventana"] = _fraccion(
        abajo,
        "concordancia_tejido_plasma_fuera_de_la_ventana",
        "Por debajo del umbral. Aquí «no concuerda» puede ser solo «no se mide».",
    )
    dentro = salida["dentro_de_la_ventana"]
    salida["estado"] = "calculado" if dentro["suficiente"] else "no_evaluable"
    if not dentro["suficiente"]:
        salida["motivo"] = (
            f"{dentro['denominador']} proteínas dentro de la ventana no llegan al "
            f"mínimo de {n_min}."
        )
    # El veredicto sale de los dos intervalos, no de comparar los puntos.
    compatible_con_azar = (
        dentro["ci_low"] is not None
        and dentro["ci_low"] <= 0.5 <= dentro["ci_high"]
    )
    se_solapan = (
        dentro["ci_low"] is not None
        and todas["ci_low"] is not None
        and dentro["ci_low"] <= todas["ci_high"]
        and todas["ci_low"] <= dentro["ci_high"]
    )
    salida["veredicto"] = {
        "dentro_de_la_ventana_es_compatible_con_azar": bool(compatible_con_azar),
        "los_dos_intervalos_se_solapan": bool(se_solapan),
        "la_ventana_explica_la_brecha": bool(not compatible_con_azar and not se_solapan),
        "lectura": (
            "Restringir a la ventana de detección no cambia el resultado: la "
            "concordancia sigue siendo compatible con el azar. Así que la brecha "
            "tejido→plasma no es un artefacto del límite de detección, y eso es una "
            "conclusión más fuerte que la de v1, no más débil: ahora se sabe que "
            "tampoco se transfiere entre las proteínas que el método sí alcanza."
            if compatible_con_azar
            else "Dentro de la ventana la concordancia deja de ser compatible con el "
            "azar, así que parte de la brecha de v1 era el método y no la biología."
        ),
    }
    return salida


def tabla_por_proteina(
    lfc_tejido: pd.Series,
    lfc_plasma: pd.Series,
    concentracion_pg_l: pd.Series,
    detectada: pd.Series,
    *,
    umbral_log10: float | None,
) -> pd.DataFrame:
    """Una fila por proteína, para que cada número del informe salga de aquí."""
    idx = sorted(set(lfc_tejido.index) | set(lfc_plasma.index))
    conc = pd.to_numeric(concentracion_pg_l.reindex(idx), errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        log10 = np.log10(conc.where(conc > 0))
    tabla = pd.DataFrame(
        {
            "identificador": idx,
            "log2fc_tejido": pd.to_numeric(lfc_tejido.reindex(idx), errors="coerce"),
            "log2fc_plasma": pd.to_numeric(lfc_plasma.reindex(idx), errors="coerce"),
            "concentracion_sangre_pg_l": conc.to_numpy(),
            "log10_concentracion": log10.to_numpy(),
            "detectada_en_plasma": detectada.reindex(idx).to_numpy(),
        }
    )
    if umbral_log10 is not None:
        tabla["dentro_de_la_ventana"] = tabla["log10_concentracion"] >= umbral_log10
    tabla["mismo_signo"] = np.where(
        tabla["log2fc_tejido"].notna() & tabla["log2fc_plasma"].notna(),
        np.sign(tabla["log2fc_tejido"]) == np.sign(tabla["log2fc_plasma"]),
        None,
    )
    return tabla
