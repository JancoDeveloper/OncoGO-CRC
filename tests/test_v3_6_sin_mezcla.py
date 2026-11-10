from oncogo_crc.v3.vigilante import clasificar


def test_una_lectura_de_v3_4_esta_prohibida_y_la_corrida_propia_no():
    assert clasificar("/tmp/repo/results/v3_4/x.json") == "prohibida"
    assert clasificar("/Volumes/JANCODATA/OncoGO-CRC/results/v3_6/h8.json") == "permitida"
    assert clasificar("/Volumes/JANCODATA/OncoGO-CRC/results/certified/inference_audit.json") == "permitida"
    assert clasificar(
        "/Volumes/JANCODATA/OncoGO-CRC/results/certification/certificado_resumen.json",
        ("results/certification/certificado_resumen.json",),
    ) == "excepcion"
