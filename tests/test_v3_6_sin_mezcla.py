from pathlib import Path

from oncogo_crc.stages.v3_5_preanalitica import _sin_ruta_local
from oncogo_crc.v3.vigilante import clasificar, excepciones_de


def test_una_lectura_de_v3_4_esta_prohibida_y_la_corrida_propia_no():
    assert clasificar("/tmp/repo/results/v3_4/x.json") == "prohibida"
    assert clasificar("/Volumes/JANCODATA/OncoGO-CRC/results/v3_6/h8.json") == "permitida"
    assert clasificar("/Volumes/JANCODATA/OncoGO-CRC/results/certified/inference_audit.json") == "permitida"
    assert clasificar(
        "/Volumes/JANCODATA/OncoGO-CRC/results/certification/certificado_resumen.json",
        ("results/certification/certificado_resumen.json",),
    ) == "excepcion"


def test_la_lista_de_entradas_cubre_el_resumen_de_certificacion():
    raiz = Path(__file__).resolve().parents[1]
    assert "results/certification/certificado_resumen.json" in excepciones_de(raiz / "config" / "entradas_v3_6.yaml")


def test_el_error_del_analizador_no_guarda_una_ruta_absoluta():
    texto = _sin_ruta_local("Started parsing /Volumes/JANCODATA/OncoGO-CRC/data/raw/a.raw y falló")
    assert "/Volumes/" not in texto
    assert "JANCODATA" not in texto
