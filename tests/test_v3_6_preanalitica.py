from math import comb

from oncogo_crc.stages.v3_5_preanalitica import cota_con_60, p_exacta_mann_whitney, presentar_ic


def test_el_ic_degenerado_queda_nulo_y_la_p_es_la_exacta():
    ic, motivo = presentar_ic(0.0, 0.0, 0.0)
    assert ic is None
    assert "degenera" in motivo
    p = p_exacta_mann_whitney(list(range(1, 20)), list(range(20, 37)))
    assert abs(p - 2 / comb(36, 17)) < 1e-18


def test_la_cota_con_los_extremos_es_0_y_77_sobre_400():
    cota = cota_con_60()
    assert cota["auc_min"] == 0.0
    assert cota["auc_max"] == 77 / 400
    assert cota["n_parejas_conocidas"] == 323
    assert cota["n_parejas_desconocidas"] == 77
