from pathlib import Path

from oncogo_crc.v3.salida import activar_corrida_limpia, cerrar_corrida_limpia, entrada_v3


class _Ajustes:
    def resolve(self, rel: str) -> Path:
        return Path("/no/abrir") / rel


def test_la_corrida_limpia_no_cae_en_v3_3(tmp_path):
    salida = tmp_path / "v3_6"
    salida.mkdir()
    tokens = activar_corrida_limpia(salida)
    try:
        ruta = entrada_v3(_Ajustes(), salida, "k1_primario_frente_a_mucosa.tsv")
    finally:
        cerrar_corrida_limpia(tokens)
    assert ruta == salida / "k1_primario_frente_a_mucosa.tsv"
    assert "v3_3" not in str(ruta)
