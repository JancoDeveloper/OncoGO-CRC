"""Composición de tejido como covariable, sin circularidad (corrección C7).

El problema de v1. El índice estromal se construyó con las 282 proteínas
anotadas en matriz extracelular (GO:0031012), y esas mismas 282 proteínas se
contrastaron contra ese índice. Una proteína de matriz extracelular contribuye a
su propia covariable, así que parte del ajuste es la proteína explicándose a sí
misma.

La corrección. Al ajustar la proteína ``i``, su índice se recalcula sin ella.
Como el índice es la media de los valores estandarizados, quitar una proteína es
aritmética simple:

    indice_sin_i = (suma_z - z_i) / (m - 1)

Eso vuelve el cálculo exacto y barato: no hay que recomputar nada por proteína.
Para las proteínas que no están en la firma, el índice no cambia.

Las firmas alternativas (músculo liso, inmune, epitelial) están preregistradas
en ``config/preregistro_v2.yaml`` y sirven para ver si la conclusión depende de
qué compartimento se mide o no.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class IndiceComposicion:
    """Índice de composición por espécimen, con el detalle para la exclusión."""

    nombre: str
    score: pd.Series            # un valor por espécimen
    z: pd.DataFrame             # proteínas de la firma x especímenes, estandarizadas
    genes_presentes: list[str]

    @property
    def m(self) -> int:
        return len(self.genes_presentes)

    def score_sin(self, gen: str) -> pd.Series | None:
        """Índice recalculado sin ``gen``. ``None`` si el gen no está en la firma."""
        clave = str(gen).strip().upper()
        if clave not in self._indice_z:
            return None
        if self.m <= 2:
            # Con dos proteínas, quitar una deja un índice de una sola: no es un índice.
            return None
        fila = self.z.loc[self._indice_z[clave]]
        suma = self.z.sum(axis=0, skipna=True)
        n_ok = self.z.notna().sum(axis=0)
        presente = fila.notna()
        suma_sin = suma - fila.fillna(0.0)
        n_sin = n_ok - presente.astype(int)
        return (suma_sin / n_sin.replace(0, np.nan)).rename(f"{self.nombre}_sin_{clave}")

    def __post_init__(self) -> None:
        self._indice_z = {str(g).strip().upper(): g for g in self.z.index}


def construir_indice(
    matriz: pd.DataFrame,
    columnas: list[str],
    genes: set[str],
    *,
    nombre: str,
) -> IndiceComposicion | None:
    """Media de los valores estandarizados de la firma, por espécimen.

    ``matriz`` va indexada por identificador de proteína y sin filas de resumen.
    Estandarizar por proteína evita que las proteínas más abundantes dominen la
    media. Devuelve ``None`` si la firma no llega a tres proteínas presentes.
    """
    presentes = [g for g in matriz.index.astype(str) if str(g).strip().upper() in genes]
    if len(presentes) < 3:
        return None
    bloque = matriz.loc[presentes, columnas].apply(pd.to_numeric, errors="coerce")
    mu = bloque.mean(axis=1)
    sd = bloque.std(axis=1, ddof=1).replace(0.0, np.nan)
    z = bloque.sub(mu, axis=0).div(sd, axis=0)
    return IndiceComposicion(
        nombre=nombre,
        score=z.mean(axis=0, skipna=True),
        z=z,
        genes_presentes=presentes,
    )


def delta_por_par(score: pd.Series, pares: list[tuple[str, str]]) -> np.ndarray:
    """Contraste de composición dentro de cada par: tumor menos su normal."""
    return np.array(
        [float(score.get(a, np.nan)) - float(score.get(b, np.nan)) for a, b in pares],
        dtype=float,
    )


def regresion_por_proteina(
    diffs: np.ndarray,
    covariables: np.ndarray,
    *,
    min_pares: int,
    evaluar_en: float = 0.0,
) -> tuple[np.ndarray, ...]:
    """Ajusta ``diff = a + b·x`` por proteína, con una covariable distinta por proteína.

    ``covariables`` tiene la misma forma que ``diffs``: una fila por proteína y
    una columna por par. Así cada proteína puede llevar su índice sin ella misma.
    El efecto que se devuelve es la predicción en ``x = evaluar_en``, que con cero
    es el efecto que quedaría si el par no tuviera diferencia de composición.

    Evaluar en la media de ``x`` devolvería la media de ``y`` y no ajustaría nada:
    ese fue el primer error de v1 y por eso el parámetro es explícito.
    """
    diffs = np.asarray(diffs, dtype=float)
    covariables = np.asarray(covariables, dtype=float)
    if covariables.shape != diffs.shape:
        raise ValueError(
            f"covariables {covariables.shape} no coincide con diffs {diffs.shape}: "
            "se necesita una covariable por proteína y par"
        )
    n_prot = diffs.shape[0]
    coef = np.full(n_prot, np.nan)
    s2 = np.full(n_prot, np.nan)
    df = np.zeros(n_prot)
    stdev = np.full(n_prot, np.nan)
    n_used = np.zeros(n_prot, dtype=int)
    pendiente = np.full(n_prot, np.nan)
    pendiente_t = np.full(n_prot, np.nan)
    minimo = max(int(min_pares), 3)

    for i in range(n_prot):
        y = diffs[i]
        x = covariables[i]
        ok = np.isfinite(y) & np.isfinite(x)
        n = int(ok.sum())
        n_used[i] = n
        if n < minimo:
            continue
        yi = y[ok]
        xi = x[ok]
        xbar = float(xi.mean())
        xc = xi - xbar
        sxx = float(np.dot(xc, xc))
        if sxx <= 0:
            # Covariable constante: no hay nada que ajustar, el efecto es la media.
            coef[i] = float(yi.mean())
            resid = yi - yi.mean()
            df[i] = n - 1
            if df[i] > 0:
                s2[i] = float(np.dot(resid, resid) / df[i])
            stdev[i] = float(np.sqrt(1.0 / n))
            pendiente[i] = 0.0
            continue
        b = float(np.dot(xc, yi) / sxx)
        ybar = float(yi.mean())
        resid = yi - (ybar + b * xc)
        dof = n - 2
        if dof <= 0:
            continue
        var = float(np.dot(resid, resid) / dof)
        coef[i] = ybar + b * (float(evaluar_en) - xbar)
        s2[i] = var
        df[i] = dof
        stdev[i] = float(np.sqrt(1.0 / n + (float(evaluar_en) - xbar) ** 2 / sxx))
        pendiente[i] = b
        if var > 0:
            pendiente_t[i] = b / float(np.sqrt(var / sxx))
    return coef, s2, df, stdev, n_used, pendiente, pendiente_t


def regresion_multiple_por_proteina(
    diffs: np.ndarray,
    covariables: list[np.ndarray],
    *,
    min_pares: int,
) -> tuple[np.ndarray, ...]:
    """Ajusta ``diff = a + b1·x1 + ... + bk·xk`` por proteína.

    El efecto reportado es el intercepto, que es la predicción con todas las
    covariables en cero: el efecto que quedaría si el par no tuviera ninguna
    diferencia de composición en ninguno de los compartimentos medidos.

    Sirve para separar compartimentos que se solapan. El índice de matriz
    extracelular y el de músculo liso comparten proteínas de pared intestinal, así
    que ajustar por uno solo no dice cuál manda.
    """
    diffs = np.asarray(diffs, dtype=float)
    k = len(covariables)
    if k == 0:
        raise ValueError("se necesita al menos una covariable")
    covs = [np.asarray(c, dtype=float) for c in covariables]
    for c in covs:
        if c.shape != diffs.shape:
            raise ValueError(f"covariable {c.shape} no coincide con diffs {diffs.shape}")

    n_prot = diffs.shape[0]
    coef = np.full(n_prot, np.nan)
    s2 = np.full(n_prot, np.nan)
    df = np.zeros(n_prot)
    stdev = np.full(n_prot, np.nan)
    n_used = np.zeros(n_prot, dtype=int)
    pendientes = np.full((n_prot, k), np.nan)
    # Hacen falta k pendientes más el intercepto, y al menos un grado para el error.
    minimo = max(int(min_pares), k + 2)

    for i in range(n_prot):
        y = diffs[i]
        xs = [c[i] for c in covs]
        ok = np.isfinite(y)
        for x in xs:
            ok = ok & np.isfinite(x)
        n = int(ok.sum())
        n_used[i] = n
        if n < minimo:
            continue
        design = np.column_stack([np.ones(n)] + [x[ok] for x in xs])
        yi = y[ok]
        try:
            xtx_inv = np.linalg.inv(design.T @ design)
        except np.linalg.LinAlgError:
            # Covariables colineales: no se puede separar su contribución.
            continue
        beta = xtx_inv @ (design.T @ yi)
        resid = yi - design @ beta
        dof = n - design.shape[1]
        if dof <= 0:
            continue
        var = float(resid @ resid / dof)
        coef[i] = float(beta[0])
        s2[i] = var
        df[i] = dof
        stdev[i] = float(np.sqrt(max(xtx_inv[0, 0], 0.0)))
        pendientes[i] = beta[1:]
    return coef, s2, df, stdev, n_used, pendientes


def covariables_leave_one_out(
    indice: IndiceComposicion,
    pares: list[tuple[str, str]],
    protein_ids: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Matriz de covariables con la proteína excluida de su propio índice.

    Devuelve la matriz y cuántas proteínas recibieron un índice distinto del
    general, que son las que pertenecen a la firma.
    """
    base = delta_por_par(indice.score, pares)
    n_prot = len(protein_ids)
    matriz = np.tile(base, (n_prot, 1))
    en_firma = {str(g).strip().upper() for g in indice.genes_presentes}
    n_excluidas = 0
    cache: dict[str, np.ndarray] = {}
    for i, pid in enumerate(protein_ids):
        clave = str(pid).strip().upper()
        if clave not in en_firma:
            continue
        if clave not in cache:
            score_sin = indice.score_sin(clave)
            if score_sin is None:
                continue
            cache[clave] = delta_por_par(score_sin, pares)
        matriz[i] = cache[clave]
        n_excluidas += 1
    return matriz, n_excluidas
