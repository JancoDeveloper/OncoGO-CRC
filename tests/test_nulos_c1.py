"""Tests de C1: el p de permutación no puede quedarse en su piso sin decirlo."""

from __future__ import annotations

import numpy as np
import pytest

from oncogo_crc.utils.moderated import (
    fit_f_prior,
    n_significant_after_signs,
    paired_observed_stats,
    posterior_scale,
)
from oncogo_crc.v2.metricas import p_empirico
from oncogo_crc.v2.nulos import nulo_aleatorio, nulo_exacto, preparar


def _datos(n_prot=120, n_pares=10, efecto=0.0, semilla=20260914):
    rng = np.random.default_rng(semilla)
    d = rng.normal(0.0, 0.5, size=(n_prot, n_pares))
    if efecto:
        d[: n_prot // 4] += efecto
    return d


def test_la_forma_cerrada_coincide_con_recalcular_la_tabla():
    """El atajo vectorizado tiene que dar el mismo conteo que el camino largo de v1."""
    diffs = _datos(efecto=1.5)
    prep = preparar(diffs, min_pares=3)
    _coef, s2, df, stdev, _n = paired_observed_stats(diffs, min_pairs=3)
    d0, s02 = fit_f_prior(s2, df)
    s2_post, df_total = posterior_scale(s2, df, d0, s02)

    rng = np.random.default_rng(1)
    for _ in range(8):
        signos = rng.choice(np.array([-1.0, 1.0]), size=diffs.shape[1])
        esperado = n_significant_after_signs(
            diffs,
            signos,
            s2_post=s2_post,
            df_total=df_total,
            stdev_unscaled=stdev,
            lfc=0.585,
            alpha=0.05,
        )
        from oncogo_crc.v2.nulos import _conteos_del_bloque

        obtenido = int(
            _conteos_del_bloque(prep, signos.reshape(1, -1), lfc=0.585, alpha=0.05)[0]
        )
        assert obtenido == esperado


def test_el_observado_es_el_vector_de_signos_todo_positivo():
    diffs = _datos(efecto=2.0)
    res = nulo_aleatorio(
        diffs, lfc=0.585, alpha=0.05, min_pares=3, n_permutaciones=50, semilla=20260914
    )
    _coef, s2, df, stdev, _n = paired_observed_stats(diffs, min_pairs=3)
    d0, s02 = fit_f_prior(s2, df)
    s2_post, df_total = posterior_scale(s2, df, d0, s02)
    esperado = n_significant_after_signs(
        diffs,
        np.ones(diffs.shape[1]),
        s2_post=s2_post,
        df_total=df_total,
        stdev_unscaled=stdev,
        lfc=0.585,
        alpha=0.05,
    )
    assert res["observado"] == esperado


def test_la_enumeracion_exacta_da_dos_elevado_a_k():
    diffs = _datos(n_prot=60, n_pares=8, efecto=2.0)
    res = nulo_exacto(diffs, lfc=0.585, alpha=0.05, min_pares=3)
    assert res["tipo"] == "exacto"
    assert res["n_permutaciones"] == 2**8 == 256
    assert res["nulo"].size == 256
    assert res["semilla"] is None


def test_el_p_exacto_no_puede_bajar_de_tres_sobre_b_mas_uno():
    """La identidad y su complemento siempre empatan con lo observado."""
    diffs = _datos(n_prot=60, n_pares=8, efecto=4.0)
    res = nulo_exacto(diffs, lfc=0.585, alpha=0.05, min_pares=3)
    assert res["b_al_menos_observado"] >= res["b_triviales"] == 2
    assert res["p_minimo_alcanzable_exacto"] == pytest.approx(3 / (res["n_permutaciones"] + 1))
    p = p_empirico(res["b_al_menos_observado"], res["n_permutaciones"])
    assert p["p"] >= res["p_minimo_alcanzable_exacto"]
    assert res["b_no_triviales"] == res["b_al_menos_observado"] - 2
    assert res["es_minimo_alcanzable"] == (res["b_no_triviales"] == 0)


def test_la_distribucion_exacta_es_simetrica_en_el_signo_global():
    """Invertir todos los signos no cambia |t|, así que cada conteo aparece dos veces."""
    diffs = _datos(n_prot=40, n_pares=6, efecto=1.5)
    res = nulo_exacto(diffs, lfc=0.585, alpha=0.05, min_pares=3)
    nulo = res["nulo"]
    # El orden de itertools.product empareja el índice i con 2^k - 1 - i por complemento.
    assert np.array_equal(nulo, nulo[::-1])


def test_un_efecto_real_queda_en_el_extremo_del_nulo():
    diffs = _datos(n_prot=200, n_pares=12, efecto=3.0)
    res = nulo_aleatorio(
        diffs, lfc=0.585, alpha=0.05, min_pares=3, n_permutaciones=200, semilla=20260914
    )
    assert res["observado"] > res["nulo_p95"]
    p = p_empirico(res["b_al_menos_observado"], res["n_permutaciones"])
    assert p["p"] < 0.05
    # Con un efecto constante en todos los pares, invertir signos al azar todavía
    # deja medias grandes, así que el nulo no se vacía. Por eso el p no baja al piso.
    assert p["es_tope"] is False


def test_sin_efecto_el_observado_no_destaca():
    diffs = _datos(n_prot=200, n_pares=12, efecto=0.0)
    res = nulo_aleatorio(
        diffs, lfc=0.585, alpha=0.05, min_pares=3, n_permutaciones=200, semilla=20260914
    )
    p = p_empirico(res["b_al_menos_observado"], res["n_permutaciones"])
    assert p["p"] > 0.05


def test_mas_permutaciones_bajan_el_piso_del_p():
    """Es el motivo de C1: con 100 el p no puede bajar de 0,0099."""
    assert p_empirico(0, 100)["p_minimo_alcanzable"] == pytest.approx(1 / 101)
    assert p_empirico(0, 10000)["p_minimo_alcanzable"] == pytest.approx(1 / 10001)
    assert p_empirico(0, 65536)["p_minimo_alcanzable"] == pytest.approx(1 / 65537)
    assert p_empirico(0, 10000)["p"] < p_empirico(0, 100)["p"]


def test_el_nulo_aleatorio_es_reproducible():
    diffs = _datos(efecto=1.0)
    kw = dict(lfc=0.585, alpha=0.05, min_pares=3, n_permutaciones=60, semilla=20260914)
    a = nulo_aleatorio(diffs, **kw)
    b = nulo_aleatorio(diffs, **kw)
    assert np.array_equal(a["nulo"], b["nulo"])
    assert a["b_al_menos_observado"] == b["b_al_menos_observado"]


def test_demasiados_pares_no_se_enumeran():
    diffs = _datos(n_prot=10, n_pares=25)
    with pytest.raises(ValueError, match="demasiadas"):
        nulo_exacto(diffs, lfc=0.585, alpha=0.05, min_pares=3)


def test_el_prior_de_la_varianza_se_estima_una_sola_vez():
    """Invertir signos no cambia la varianza, así que el prior es el de los datos."""
    diffs = _datos(efecto=1.0)
    prep = preparar(diffs, min_pares=3)
    _coef, s2, df, _stdev, _n = paired_observed_stats(diffs, min_pairs=3)
    d0, s02 = fit_f_prior(s2, df)
    assert prep.prior_df == pytest.approx(d0)
    assert prep.prior_s02 == pytest.approx(s02)


def test_las_ausencias_no_se_imputan_en_el_nulo():
    diffs = _datos(n_prot=30, n_pares=10, efecto=2.0)
    diffs[0, :8] = np.nan  # solo 2 pares observados: por debajo del mínimo
    prep = preparar(diffs, min_pares=3)
    assert prep.n_obs[0] == 2
    assert not prep.evaluable[0]
    res = nulo_aleatorio(
        diffs, lfc=0.585, alpha=0.05, min_pares=3, n_permutaciones=30, semilla=20260914
    )
    # La proteína sin pares suficientes no puede contarse como significativa.
    assert res["nulo"].max() <= int(prep.evaluable.sum())
