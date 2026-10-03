"""Tests de C7: el índice de composición no puede explicarse a sí mismo.

Matrices sintéticas, solo aquí. El preregistro prohíbe datos simulados en
results/.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.v2.composicion import (
    construir_indice,
    covariables_leave_one_out,
    delta_por_par,
    regresion_por_proteina,
)


def _matriz(genes: list[str], columnas: list[str], valores: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(valores, index=genes, columns=columnas)


def test_el_indice_necesita_al_menos_tres_proteinas():
    cols = ["t1", "n1", "t2", "n2"]
    m = _matriz(["A", "B"], cols, np.random.default_rng(1).normal(size=(2, 4)))
    assert construir_indice(m, cols, {"A", "B"}, nombre="firma") is None


def test_el_indice_estandariza_por_proteina():
    """Una proteína mil veces más abundante no debe dominar la media."""
    cols = ["t1", "n1", "t2", "n2"]
    base = np.array([[1.0, -1.0, 1.0, -1.0], [2.0, -2.0, 2.0, -2.0], [0.5, -0.5, 0.5, -0.5]])
    chica = construir_indice(_matriz(["A", "B", "C"], cols, base), cols, {"A", "B", "C"}, nombre="f")
    grande = construir_indice(
        _matriz(["A", "B", "C"], cols, base * np.array([[1.0], [1000.0], [1.0]])),
        cols,
        {"A", "B", "C"},
        nombre="f",
    )
    assert chica is not None and grande is not None
    # Al estandarizar, multiplicar una fila por mil no cambia el índice.
    assert np.allclose(chica.score.to_numpy(), grande.score.to_numpy())


def test_quitar_una_proteina_da_la_media_de_las_demas():
    cols = ["t1", "n1", "t2", "n2", "t3", "n3"]
    rng = np.random.default_rng(20260914)
    vals = rng.normal(size=(5, 6))
    genes = ["A", "B", "C", "D", "E"]
    idx = construir_indice(_matriz(genes, cols, vals), cols, set(genes), nombre="f")
    assert idx is not None
    sin_c = idx.score_sin("C")
    assert sin_c is not None
    esperado = idx.z.drop(index="C").mean(axis=0, skipna=True)
    assert np.allclose(sin_c.to_numpy(), esperado.to_numpy())


def test_una_proteina_fuera_de_la_firma_no_tiene_indice_propio():
    cols = ["t1", "n1", "t2", "n2"]
    genes = ["A", "B", "C"]
    idx = construir_indice(
        _matriz(genes, cols, np.random.default_rng(2).normal(size=(3, 4))),
        cols,
        set(genes),
        nombre="f",
    )
    assert idx is not None
    assert idx.score_sin("ZZZ") is None


def test_solo_las_proteinas_de_la_firma_cambian_de_covariable():
    cols = ["t1", "n1", "t2", "n2", "t3", "n3"]
    pares = [("t1", "n1"), ("t2", "n2"), ("t3", "n3")]
    genes_firma = ["A", "B", "C", "D"]
    rng = np.random.default_rng(5)
    idx = construir_indice(
        _matriz(genes_firma, cols, rng.normal(size=(4, 6))), cols, set(genes_firma), nombre="f"
    )
    assert idx is not None
    protein_ids = np.array(["A", "B", "C", "D", "FUERA1", "FUERA2"])
    matriz, n_excluidas = covariables_leave_one_out(idx, pares, protein_ids)
    assert n_excluidas == 4
    base = delta_por_par(idx.score, pares)
    # Las dos de fuera llevan el índice general.
    assert np.allclose(matriz[4], base)
    assert np.allclose(matriz[5], base)
    # Las cuatro de la firma llevan un índice distinto del general.
    for i in range(4):
        assert not np.allclose(matriz[i], base)


def test_un_efecto_explicado_por_la_composicion_se_cae_en_cero():
    """Si la diferencia es proporcional al contraste de composición, al igualar desaparece."""
    rng = np.random.default_rng(20260914)
    n = 40
    x = rng.normal(0, 1, size=n)
    y = 1.5 * x  # el efecto es solo composición
    coef, _s2, _df, _stdev, _n, pend, _t = regresion_por_proteina(
        y.reshape(1, -1), x.reshape(1, -1), min_pares=3, evaluar_en=0.0
    )
    assert abs(coef[0]) < 1e-9
    assert pend[0] == pytest.approx(1.5, abs=1e-9)


def test_un_efecto_independiente_de_la_composicion_se_conserva():
    rng = np.random.default_rng(7)
    n = 40
    x = rng.normal(0, 1, size=n)
    y = np.full(n, 2.0) + rng.normal(0, 0.01, size=n)
    coef, *_rest = regresion_por_proteina(
        y.reshape(1, -1), x.reshape(1, -1), min_pares=3, evaluar_en=0.0
    )
    assert coef[0] == pytest.approx(2.0, abs=0.02)


def test_evaluar_en_la_media_no_ajusta_nada():
    """El bug original: con la covariable centrada, el intercepto es la media cruda."""
    rng = np.random.default_rng(11)
    n = 30
    x = rng.normal(0, 1, size=n)
    y = 1.0 + 2.0 * x
    xc = (x - x.mean()).reshape(1, -1)
    coef_media, *_ = regresion_por_proteina(y.reshape(1, -1), xc, min_pares=3, evaluar_en=0.0)
    # Con x centrada, evaluar en 0 es evaluar en la media: devuelve la media de y.
    assert coef_media[0] == pytest.approx(float(y.mean()), abs=1e-9)
    # Con x sin centrar, evaluar en 0 sí ajusta.
    coef_cero, *_ = regresion_por_proteina(
        y.reshape(1, -1), x.reshape(1, -1), min_pares=3, evaluar_en=0.0
    )
    assert coef_cero[0] == pytest.approx(1.0, abs=1e-9)
    assert abs(coef_cero[0] - float(y.mean())) > 0.1


def test_la_forma_de_la_covariable_tiene_que_coincidir():
    with pytest.raises(ValueError, match="no coincide"):
        regresion_por_proteina(np.zeros((3, 10)), np.zeros((3, 9)), min_pares=3)


def test_pocos_pares_observados_no_se_imputan():
    diffs = np.array([[1.0, np.nan, np.nan, np.nan]])
    cov = np.array([[0.1, 0.2, 0.3, 0.4]])
    coef, _s2, _df, _stdev, n_used, *_ = regresion_por_proteina(diffs, cov, min_pares=3)
    assert n_used[0] == 1
    assert np.isnan(coef[0])


def test_covariable_constante_devuelve_la_media_sin_pendiente():
    diffs = np.array([[1.0, 1.2, 0.8, 1.1]])
    cov = np.zeros((1, 4))
    coef, _s2, df, _stdev, _n, pend, _t = regresion_por_proteina(diffs, cov, min_pares=3)
    assert coef[0] == pytest.approx(float(diffs.mean()))
    assert pend[0] == 0.0
    assert df[0] == 3  # n - 1, no n - 2: no se gastó un grado en una pendiente
