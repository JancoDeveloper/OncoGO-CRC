"""La integridad de un crudo se decide con el SHA-1 de PRIDE, no con el tamaño."""

from __future__ import annotations

import importlib.util
from pathlib import Path

_RUTA = Path(__file__).resolve().parents[1] / "scripts" / "download" / "pxd055218_crudos.py"
_SPEC = importlib.util.spec_from_file_location("pxd055218_crudos", _RUTA)
_MOD = importlib.util.module_from_spec(_SPEC)
assert _SPEC.loader is not None
_SPEC.loader.exec_module(_MOD)


def test_el_checksum_de_pride_trae_la_ruta_delante_del_sha1():
    texto = "D:\\Akhila\\PhD\\DIA\\PRIDE\\CN_01_DIA.raw\t991006918db00e872d13525fa23eb4bae459b53a\n"
    tabla = _MOD.parsear_checksum(texto)
    assert tabla["CN_01_DIA.raw"] == "991006918db00e872d13525fa23eb4bae459b53a"


def test_colocar_el_crudo_nuevo_deja_el_viejo_en_cuarentena(tmp_path):
    crudo = tmp_path / "CN_01_DIA.raw"
    crudo.write_bytes(b"viejo")
    temporal = tmp_path / "cache" / "CN_01_DIA.raw"
    temporal.parent.mkdir()
    temporal.write_bytes(b"nuevo")
    _MOD.colocar(temporal, crudo, "2026-11-10")
    assert crudo.read_bytes() == b"nuevo"
    assert (tmp_path / "cuarentena" / "CN_01_DIA.raw.danado-2026-11-10").read_bytes() == b"viejo"
    assert not temporal.exists()


def test_un_sha1_distinto_pide_redescarga(tmp_path):
    crudo = tmp_path / "CN_01_DIA.raw"
    crudo.write_bytes(b"abc")
    fila = _MOD.evaluar_crudo(crudo, sha1_pride="0" * 40, bytes_pride=99)
    assert fila["igual"] is False
    assert fila["pide_redescarga"] is True
    assert fila["bytes_local"] == 3
