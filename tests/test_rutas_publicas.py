from pathlib import Path

from oncogo_crc.settings import load_settings
from oncogo_crc.v3.rutas_publicas import ruta_publica, texto_con_ruta_local

RAIZ = Path(__file__).resolve().parents[1]


def test_una_ruta_del_repositorio_no_lleva_el_disco():
    settings = load_settings(RAIZ)
    info = ruta_publica(RAIZ / "config/entradas_v3_6.yaml", settings)
    assert not info["ruta"].startswith("/")
    assert not texto_con_ruta_local(info["ruta"])
    assert info["sha256"]


def test_detecta_rutas_locales():
    assert texto_con_ruta_local("/Volumes/JANCODATA/OncoGO-CRC/results/v3_3/x.json")
    assert texto_con_ruta_local("C:\\Sonia\\SH_01.mzML")
    assert not texto_con_ruta_local("external/modelos/h7_tmbed.pred")
