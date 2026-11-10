from pathlib import Path

from oncogo_crc.v3.vigilante import clasificar


def test_durante_v3_7_leer_v3_6_esta_prohibido():
    assert clasificar("/tmp/repo/results/v3_6/x.json", version="v3_7") == "prohibida"
    assert clasificar("/tmp/repo/results/v3_7/x.json", version="v3_7") == "permitida"


def test_h8_no_asigna_el_gen_de_tanaka_solo_con_uniprot_revisado():
    texto = Path("src/oncogo_crc/v3/h8.py").read_text(encoding="utf-8")
    assert "uni.get(acceso_base" not in texto
