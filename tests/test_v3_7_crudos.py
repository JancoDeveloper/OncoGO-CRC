"""La integridad de un crudo se decide con el SHA-1 de PRIDE, no con el tamaño."""

from __future__ import annotations

import importlib.util
from pathlib import Path

_RUTA = Path(__file__).resolve().parents[1] / "scripts" / "download" / "pxd055218_crudos.py"
_SPEC = importlib.util.spec_from_file_location("pxd055218_crudos", _RUTA)
_MOD = importlib.util.module_from_spec(_SPEC)
assert _SPEC.loader is not None
_SPEC.loader.exec_module(_MOD)


def test_un_sha1_distinto_pide_redescarga(tmp_path):
    crudo = tmp_path / "CN_01_DIA.raw"
    crudo.write_bytes(b"abc")
    fila = _MOD.evaluar_crudo(crudo, sha1_pride="0" * 40, bytes_pride=99)
    assert fila["igual"] is False
    assert fila["pide_redescarga"] is True
    assert fila["bytes_local"] == 3
