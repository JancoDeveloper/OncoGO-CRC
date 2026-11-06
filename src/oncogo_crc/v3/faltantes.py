"""Presencia y ausencia, separadas de la abundancia.

Una proteína con menos del umbral de observaciones no se rescata imputando.
Se pregunta si el hueco depende del grupo, y el contraste de abundancia queda
solo en los valores observados.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def patron_por_grupo(valores: np.ndarray, etiquetas: np.ndarray) -> dict[str, int]:
    """Cuántos valores finitos y cuántos huecos hay en cada etiqueta."""
    salida: dict[str, int] = {}
    for etiqueta in pd.unique(etiquetas):
        bloque = valores[etiquetas == etiqueta]
        presentes = int(np.isfinite(bloque).sum())
        salida[f"presentes_{etiqueta}"] = presentes
        salida[f"ausentes_{etiqueta}"] = int(bloque.size - presentes)
    return salida


def fisher_presencia(
    matriz: np.ndarray,
    etiquetas: np.ndarray,
    *,
    grupo_a: str,
    grupo_b: str,
    alpha: float,
) -> pd.DataFrame:
    """Fisher exacto de presente frente a ausente, y BH sobre esas proteínas."""
    etiquetas = np.asarray(etiquetas)
    en_a = etiquetas == grupo_a
    en_b = etiquetas == grupo_b
    p_valores = np.full(matriz.shape[0], np.nan)
    odds = np.full(matriz.shape[0], np.nan)
    tablas = []
    for i in range(matriz.shape[0]):
        fila = matriz[i]
        a_si = int(np.isfinite(fila[en_a]).sum())
        a_no = int(en_a.sum() - a_si)
        b_si = int(np.isfinite(fila[en_b]).sum())
        b_no = int(en_b.sum() - b_si)
        tablas.append((a_si, a_no, b_si, b_no))
        if min(a_si + a_no, b_si + b_no) == 0:
            continue
        razon, p = fisher_exact([[a_si, a_no], [b_si, b_no]])
        odds[i] = float(razon)
        p_valores[i] = float(p)
    padj = np.full(p_valores.shape, np.nan)
    ok = np.isfinite(p_valores)
    if ok.any():
        _, q, _, _ = multipletests(p_valores[ok], alpha=alpha, method="fdr_bh")
        padj[ok] = q
    return pd.DataFrame(
        {
            f"presentes_{grupo_a}": [t[0] for t in tablas],
            f"ausentes_{grupo_a}": [t[1] for t in tablas],
            f"presentes_{grupo_b}": [t[2] for t in tablas],
            f"ausentes_{grupo_b}": [t[3] for t in tablas],
            "odds_ratio_presencia": odds,
            "p_fisher": p_valores,
            "padj_fisher": padj,
            "presencia_significativa": padj < alpha,
        }
    )
