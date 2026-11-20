from pathlib import Path

from oncogo_crc.settings import load_settings
from oncogo_crc.v3.rutas_publicas import ruta_publica, texto_con_ruta_local

RAIZ = Path(__file__).resolve().parents[1]


def test_una_ruta_del_repositorio_no_lleva_el_disco():
    settings = load_settings(RAIZ)
    info = ruta_publica(RAIZ / "config/entradas_v3_7.yaml", settings)
    assert not info["ruta"].startswith("/")
    assert not texto_con_ruta_local(info["ruta"])
    assert info["sha256"]


def test_detecta_rutas_locales():
    assert texto_con_ruta_local("/Volumes/JANCODATA/OncoGO-CRC/results/v3_3/x.json")
    assert texto_con_ruta_local("C:\\Sonia\\SH_01.mzML")
    assert not texto_con_ruta_local("external/modelos/h7_tmbed.pred")


def test_un_enlace_conserva_la_ruta_del_repositorio(tmp_path, monkeypatch):
    datos = tmp_path / "datos"
    hpa = datos / "external/hpa/hpa_concentracion_sangre.tsv"
    pred = datos / "external/modelos/h7_tmbed.pred"
    hpa.parent.mkdir(parents=True)
    pred.parent.mkdir(parents=True)
    hpa.write_text("gene\tlocalizacion_secretoma\nTP53\tx\n", encoding="utf-8")
    pred.write_text("pred\n", encoding="utf-8")
    salida = tmp_path / "salida"
    (salida / "external/hpa").mkdir(parents=True)
    (salida / "external/modelos").mkdir(parents=True)
    (salida / "external/hpa/hpa_concentracion_sangre.tsv").symlink_to(hpa)
    (salida / "external/modelos/h7_tmbed.pred").symlink_to(pred)
    monkeypatch.delenv("ONCOGO_DATA_DIR", raising=False)
    monkeypatch.setenv("ONCOGO_CRC_DATA_WORKSPACE", str(salida))
    settings = load_settings()
    assert ruta_publica(settings.resolve("external/hpa/hpa_concentracion_sangre.tsv"), settings)["ruta"] == (
        "external/hpa/hpa_concentracion_sangre.tsv"
    )
    assert ruta_publica(hpa.resolve(), settings)["ruta"] == "external/hpa/hpa_concentracion_sangre.tsv"
    assert ruta_publica(settings.resolve("external/modelos/h7_tmbed.pred"), settings)["ruta"] == (
        "external/modelos/h7_tmbed.pred"
    )


def test_data_dir_conserva_el_prefijo_de_la_entrada(tmp_path, monkeypatch):
    datos = tmp_path / "datos"
    pred = datos / "external/modelos/h7_tmbed.pred"
    pred.parent.mkdir(parents=True)
    pred.write_text("pred\n", encoding="utf-8")
    salida = tmp_path / "salida"
    monkeypatch.setenv("ONCOGO_DATA_DIR", str(datos))
    monkeypatch.setenv("ONCOGO_CRC_DATA_WORKSPACE", str(salida))
    settings = load_settings()
    info = ruta_publica(settings.resolve("external/modelos/h7_tmbed.pred"), settings)
    assert info["ruta"] == "external/modelos/h7_tmbed.pred"
