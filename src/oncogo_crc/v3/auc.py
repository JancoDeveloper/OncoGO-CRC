"""AUC con signo, intervalo de DeLong, bootstrap y nulo del máximo.

El AUC no se pliega para reportar una proteína. El plegado solo entra en el
nulo del máximo del catálogo, que es el estadístico de haber buscado entre
muchas candidatas.
"""

from __future__ import annotations

import numpy as np

Z95 = 1.959963984540054


def auc_firmado(positivo: np.ndarray, negativo: np.ndarray) -> float:
    """AUC con la primera clase como caso. Los empates valen 0,5."""
    pos = np.asarray(positivo, dtype=float)
    neg = np.asarray(negativo, dtype=float)
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    ordenado = np.sort(neg)
    izquierda = np.searchsorted(ordenado, pos, side="left")
    derecha = np.searchsorted(ordenado, pos, side="right")
    aciertos = float(izquierda.sum() + 0.5 * (derecha - izquierda).sum())
    return aciertos / float(pos.size * neg.size)


def delong_ic(positivo: np.ndarray, negativo: np.ndarray) -> tuple[float, float, float]:
    """AUC y su intervalo al 95 % con la varianza de DeLong para una curva."""
    pos = np.asarray(positivo, dtype=float)
    neg = np.asarray(negativo, dtype=float)
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    if pos.size < 2 or neg.size < 2:
        auc = auc_firmado(pos, neg)
        return auc, float("nan"), float("nan")
    mayor = pos[:, None] > neg[None, :]
    igual = pos[:, None] == neg[None, :]
    psi = mayor.astype(float) + 0.5 * igual.astype(float)
    v10 = psi.mean(axis=1)
    v01 = psi.mean(axis=0)
    auc = float(v10.mean())
    varianza = float(v10.var(ddof=1) / pos.size + v01.var(ddof=1) / neg.size)
    if not np.isfinite(varianza) or varianza < 0:
        return auc, float("nan"), float("nan")
    margen = Z95 * float(np.sqrt(varianza))
    return auc, float(np.clip(auc - margen, 0.0, 1.0)), float(np.clip(auc + margen, 0.0, 1.0))


def bootstrap_auc(
    positivo: np.ndarray,
    negativo: np.ndarray,
    *,
    n_boot: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    """Intervalo percentil estratificado. Remuestrea cada grupo por su lado."""
    pos = np.asarray(positivo, dtype=float)
    neg = np.asarray(negativo, dtype=float)
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    if pos.size < 2 or neg.size < 2:
        return float("nan"), float("nan")
    aucs = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        aucs[i] = auc_firmado(
            pos[rng.integers(0, pos.size, pos.size)],
            neg[rng.integers(0, neg.size, neg.size)],
        )
    bajo, alto = np.nanquantile(aucs, [0.025, 0.975])
    return float(bajo), float(alto)


def _auc_de_fila(valores: np.ndarray, etiquetas: np.ndarray, *, min_n: int) -> float:
    pos = valores[etiquetas == 1]
    neg = valores[etiquetas == 0]
    if np.isfinite(pos).sum() < min_n or np.isfinite(neg).sum() < min_n:
        return float("nan")
    return auc_firmado(pos, neg)


def nulo_del_maximo_plegado(
    matriz: np.ndarray,
    etiquetas: np.ndarray,
    *,
    n_perm: int,
    min_n: int,
    rng: np.random.Generator,
) -> dict[str, float | int]:
    """Máximo del AUC plegado del catálogo, y su nulo al permutar etiquetas.

    Una permutación mueve las etiquetas de todas las proteínas a la vez.
    """
    matriz = np.asarray(matriz, dtype=float)
    etiquetas = np.asarray(etiquetas, dtype=int)
    observados = np.array([_auc_de_fila(matriz[i], etiquetas, min_n=min_n) for i in range(matriz.shape[0])])
    plegado = np.maximum(observados, 1.0 - observados)
    plegado = plegado[np.isfinite(plegado)]
    if plegado.size == 0:
        return {
            "maximo_observado": float("nan"),
            "n_proteinas": 0,
            "n_permutaciones": n_perm,
            "b": n_perm,
            "p": 1.0,
            "media_nula": float("nan"),
            "percentil95_nulo": float("nan"),
        }
    maximo = float(plegado.max())
    nulos = np.empty(n_perm, dtype=float)
    n = etiquetas.size
    for b in range(n_perm):
        movidas = etiquetas[rng.permutation(n)]
        aucs = np.array([_auc_de_fila(matriz[i], movidas, min_n=min_n) for i in range(matriz.shape[0])])
        doblados = np.maximum(aucs, 1.0 - aucs)
        doblados = doblados[np.isfinite(doblados)]
        nulos[b] = float(doblados.max()) if doblados.size else float("nan")
    validos = nulos[np.isfinite(nulos)]
    b = int(np.sum(validos >= maximo))
    return {
        "maximo_observado": maximo,
        "n_proteinas": int(plegado.size),
        "n_permutaciones": int(n_perm),
        "b": b,
        "p": float((b + 1) / (n_perm + 1)),
        "media_nula": float(np.mean(validos)) if validos.size else float("nan"),
        "percentil95_nulo": float(np.quantile(validos, 0.95)) if validos.size else float("nan"),
    }


def bootstrap_diferencia(
    crc: np.ndarray,
    enfermedad: np.ndarray,
    sanos: np.ndarray,
    *,
    n_boot: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    """Diferencia de AUC entre dos contrastes que comparten un grupo.

    DeLong no sirve: no son las mismas dos muestras. Se remuestrean los tres
    grupos y se recalculan los dos AUC.
    """
    a = np.asarray(crc, dtype=float)
    b = np.asarray(enfermedad, dtype=float)
    c = np.asarray(sanos, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    c = c[np.isfinite(c)]
    auc_crc = auc_firmado(a, b)
    auc_enf = auc_firmado(b, c)
    if not (np.isfinite(auc_crc) and np.isfinite(auc_enf)):
        return {
            "auc_crc_vs_enfermedad": auc_crc,
            "auc_enfermedad_vs_sano": auc_enf,
            "diferencia": float("nan"),
            "ci_bajo": float("nan"),
            "ci_alto": float("nan"),
        }
    diferencias = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        aa = a[rng.integers(0, a.size, a.size)]
        bb = b[rng.integers(0, b.size, b.size)]
        cc = c[rng.integers(0, c.size, c.size)]
        diferencias[i] = auc_firmado(aa, bb) - auc_firmado(bb, cc)
    bajo, alto = np.nanquantile(diferencias, [0.025, 0.975])
    return {
        "auc_crc_vs_enfermedad": float(auc_crc),
        "auc_enfermedad_vs_sano": float(auc_enf),
        "diferencia": float(auc_crc - auc_enf),
        "ci_bajo": float(bajo),
        "ci_alto": float(alto),
    }
