"""El veredicto de v3.3 sale de n_fallos, no del texto estado."""

from __future__ import annotations

from oncogo_crc.stages.v3_4 import veredicto_desde_n_fallos


def test_los_dos_json_de_controles():
    assert veredicto_desde_n_fallos({"estado": "calculado", "n_fallos": 3}) == "falla"
    assert veredicto_desde_n_fallos({"estado": "calculado", "n_fallos": 0}) == "pasa"
