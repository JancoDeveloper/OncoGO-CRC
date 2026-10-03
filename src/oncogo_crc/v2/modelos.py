"""Remuestreo y validación agrupados por paciente (corrección C6).

El problema de v1. Los 32 especímenes de prueba de PXD034575 son 16 pacientes
con dos muestras cada uno: su tumor y su tejido normal. El bootstrap de v1
remuestreaba especímenes, lo que trata las dos muestras de un paciente como si
fueran de dos personas distintas. La información real son 16 unidades
independientes, no 32, y el intervalo tiene que reflejarlo.

Una advertencia sobre la dirección del cambio. Es tentador decir que agrupar
siempre ensancha el intervalo. Aquí no: el diseño es emparejado y la
discriminación ocurre **dentro** de cada paciente, así que el efecto de paciente
se cancela al comparar su tumor con su normal. Mantener la pareja junta conserva
esa cancelación; romperla mezcla tumores de un paciente con normales de otro y
añade una varianza que el diseño no tiene. Por eso el bootstrap por
conglomerado puede salir más estrecho. Lo que lo justifica no es que dé un
intervalo mayor, sino que remuestrea la unidad que de verdad se repetiría si se
volviera a hacer el estudio: el paciente.

Dos remedios.

- **Bootstrap por conglomerado**: se remuestrean pacientes; cada paciente entra
  con sus dos muestras o con ninguna.
- **Validación cruzada agrupada**: ningún paciente aparece a la vez en
  entrenamiento y en validación, y los hiperparámetros se eligen en un bucle
  interno que también agrupa.

Se añade la diferencia de AUC por pares con su intervalo, porque comparar dos
intervalos que se solapan no es lo mismo que mirar el intervalo de la
diferencia: el segundo usa que los dos modelos se evalúan en las mismas
muestras.
"""

from __future__ import annotations

import numpy as np
from sklearn.metrics import roc_auc_score


def bootstrap_por_conglomerado(
    y: np.ndarray,
    grupos: np.ndarray,
    puntajes: dict[str, np.ndarray],
    *,
    n_bootstrap: int,
    semilla: int,
) -> dict:
    """Intervalos de AUC remuestreando pacientes, no especímenes.

    ``grupos`` dice a qué paciente pertenece cada muestra. En cada repetición se
    sortean pacientes con reemplazo y se juntan todas sus muestras.
    """
    y = np.asarray(y)
    grupos = np.asarray(grupos)
    if y.size != grupos.size:
        raise ValueError("y y grupos deben tener el mismo tamaño")
    unicos = np.unique(grupos)
    indices_por_grupo = {g: np.flatnonzero(grupos == g) for g in unicos}
    rng = np.random.default_rng(semilla)

    nombres = list(puntajes)
    muestras: dict[str, list[float]] = {n: [] for n in nombres}
    # Las diferencias se calculan en la misma remuestra para los dos modelos.
    diferencias: dict[tuple[str, str], list[float]] = {
        (a, b): [] for i, a in enumerate(nombres) for b in nombres[i + 1 :]
    }
    usadas = 0
    for _ in range(n_bootstrap):
        elegidos = rng.choice(unicos, size=unicos.size, replace=True)
        idx = np.concatenate([indices_por_grupo[g] for g in elegidos])
        yb = y[idx]
        if np.unique(yb).size < 2:
            continue
        usadas += 1
        aucs = {}
        for nombre in nombres:
            aucs[nombre] = float(roc_auc_score(yb, np.asarray(puntajes[nombre])[idx]))
            muestras[nombre].append(aucs[nombre])
        for (a, b) in diferencias:
            diferencias[(a, b)].append(aucs[a] - aucs[b])

    def _resumen(arr: list[float]) -> dict:
        a = np.asarray(arr)
        if a.size == 0:
            return {"estado": "no_calculable"}
        lo, hi = np.percentile(a, [2.5, 97.5])
        return {
            "mediana": float(np.median(a)),
            "ci_low": float(lo),
            "ci_high": float(hi),
            "n_bootstrap_utiles": int(a.size),
        }

    return {
        "unidad_de_remuestreo": "paciente",
        "n_pacientes": int(unicos.size),
        "n_muestras": int(y.size),
        "n_bootstrap_pedidos": int(n_bootstrap),
        "n_bootstrap_utiles": usadas,
        "semilla": int(semilla),
        "por_modelo": {n: _resumen(muestras[n]) for n in nombres},
        "diferencias_por_pares": {
            f"{a}_menos_{b}": {
                **_resumen(v),
                "incluye_cero": (
                    None
                    if not v
                    else bool(
                        np.percentile(v, 2.5) <= 0.0 <= np.percentile(v, 97.5)
                    )
                ),
            }
            for (a, b), v in diferencias.items()
        },
        "nota": (
            "Cada paciente entra con sus dos muestras o con ninguna. Remuestrear "
            "especímenes trataría el tumor y el normal del mismo paciente como dos "
            "personas. En un diseño emparejado eso no estrecha el intervalo: lo "
            "ensancha, porque rompe la cancelación del efecto de paciente."
        ),
    }


def veredicto_sin_solape(por_modelo: dict[str, dict], aucs: dict[str, float]) -> dict:
    """Regla preregistrada: solo hay ganador si los intervalos no se solapan."""
    orden = sorted(aucs, key=lambda k: aucs[k], reverse=True)
    if len(orden) < 2:
        return {"veredicto": "no_evaluable", "ranking": orden}
    mejor, segundo = orden[0], orden[1]
    a = por_modelo.get(mejor, {})
    b = por_modelo.get(segundo, {})
    if "ci_low" not in a or "ci_high" not in b:
        return {"veredicto": "no_evaluable", "ranking": orden}
    disjunto = a["ci_low"] > b["ci_high"]
    return {
        "veredicto": f"gana {mejor}" if disjunto else "sin ganador: los intervalos se solapan",
        "ranking": orden,
        "mejor": mejor,
        "segundo": segundo,
        "intervalos_disjuntos": bool(disjunto),
    }


def grupos_de_cptac(columnas: list[str], meta) -> np.ndarray:
    """Participante de cada columna de CPTAC, para agrupar la validación cruzada."""
    mapa = meta.set_index("sample_id")["participant_id"].astype(str)
    return np.array([str(mapa.get(c, c)) for c in columnas])
