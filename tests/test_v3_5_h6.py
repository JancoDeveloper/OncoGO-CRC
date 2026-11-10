"""H6 marca una tabla sin variación en el desenlace y no cambia la regla del estado."""

from __future__ import annotations

from scipy.stats import fisher_exact

from oncogo_crc.v3.h6 import tabla_degenerada


def test_una_columna_en_cero_marca_la_tabla_degenerada():
    tabla = [[20, 0], [98, 0]]
    _, p = fisher_exact(tabla, alternative="two-sided")
    assert tabla_degenerada(tabla) is True
    assert (p < 0.05) is False
    assert tabla_degenerada([[20, 3], [10, 8]]) is False
