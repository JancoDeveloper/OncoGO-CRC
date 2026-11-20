"""La descarga de Zenodo comprueba la huella. Estas pruebas no abren la red."""

from __future__ import annotations

import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path

_RUTA = Path(__file__).resolve().parents[1] / "scripts" / "descargar_datos.py"
_SPEC = importlib.util.spec_from_file_location("descargar_datos", _RUTA)
_MOD = importlib.util.module_from_spec(_SPEC)
assert _SPEC.loader is not None
_SPEC.loader.exec_module(_MOD)

RAIZ = Path(__file__).resolve().parents[1]
PAYLOAD = b"oncogo-zenodo-prueba\n"
DIGEST = hashlib.sha256(PAYLOAD).hexdigest()
ESPERADAS = {
    "data/raw/cptac/cptac_crc_matrix.tsv": (
        "a1272baf3bb84218b60d48e971f576f3b5aee6f3fdce4c25a46ad85a352117d1",
        "cptac_crc_matrix.tsv",
    ),
    "results/certified/cptac_unshared_treat.tsv": (
        "e5e73c2bdf543ec5493846e71b0a92143bde59fd891550856e632d451eed7cf6",
        "cptac_unshared_treat.tsv",
    ),
    "results/certified/inference_audit.json": (
        "1fd5e8cee6ef044b12d0b0b0a05f522bd788ad6f9330099fb23dc9e57053f9ce",
        "inference_audit.json",
    ),
    "results/certification/certificado_resumen.json": (
        "ece8a9079cb1e7db29403ce984fc6322fac7b0387b4e7e87259b8127db4d9590",
        "certificado_resumen.json",
    ),
    "external/modelos/h7_esm2_embeddings.tsv": (
        "2b1e48e3d90476bd1cd9f2c8a2091b18be9b0dd119c61c87ac912729bd64e6b7",
        "h7_esm2_embeddings.tsv",
    ),
}


def _item(rel: str) -> dict:
    return {"rel": rel, "url": "https://zenodo.test/fichero", "sha256": DIGEST, "bytes": len(PAYLOAD)}


def _escribir(destino: Path, contenido: bytes) -> None:
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_bytes(contenido)


def _argv(tmp_path: Path, rel: str) -> list[str]:
    return ["descargar_datos.py", "--destino", str(tmp_path), "--solo", rel]


def test_el_deposito_publicado_tiene_cinco_rutas_y_el_doi():
    assert _MOD.ZENODO_DOI == "10.5281/zenodo.24398807"
    assert _MOD.ZENODO == "https://zenodo.org/records/24398807"
    assert _MOD.ZENODO_LICENCIA == "CC-BY-4.0"
    assert "PENDIENTE" not in _MOD.ZENODO
    assert len(_MOD.ZENODO_DESCARGAS) == 5
    for item in _MOD.ZENODO_DESCARGAS:
        digest, nombre = ESPERADAS[item["rel"]]
        assert item["sha256"] == digest
        assert item["url"] == f"https://zenodo.org/records/24398807/files/{nombre}?download=1"
        assert item in _MOD.FUENTES
    assert "PDC000116" in _MOD.MENSAJE_ZENODO
    assert "Vasaikar et al. 2019" in _MOD.MENSAJE_ZENODO
    assert _MOD.ZENODO_DOI in _MOD.MENSAJE_ZENODO


def test_las_huellas_coinciden_con_el_manifiesto():
    texto = (RAIZ / "data" / "checksums.sha256").read_text(encoding="utf-8")
    for item in _MOD.ZENODO_DESCARGAS:
        assert f"{item['sha256']}  {item['rel']}" in texto


def test_no_vuelve_a_bajar_si_la_huella_coincide(tmp_path):
    rel = "data/raw/cptac/cptac_crc_matrix.tsv"
    destino = tmp_path / rel
    _escribir(destino, PAYLOAD)
    llamadas: list[str] = []
    assert _MOD.asegurar(_item(rel), tmp_path, bajar=lambda url, dest: llamadas.append(url))
    assert llamadas == []
    assert destino.read_bytes() == PAYLOAD


def test_baja_y_acepta_la_huella(tmp_path):
    rel = "results/certified/cptac_unshared_treat.tsv"

    def bajar(url, destino):
        assert url == "https://zenodo.test/fichero"
        _escribir(destino, PAYLOAD)

    assert _MOD.asegurar(_item(rel), tmp_path, bajar=bajar)
    assert (tmp_path / rel).read_bytes() == PAYLOAD


def test_una_huella_distinta_falla_y_deja_el_fichero(tmp_path, capsys):
    rel = "results/certified/inference_audit.json"

    def bajar(url, destino):
        _escribir(destino, b"otro")

    assert _MOD.asegurar(_item(rel), tmp_path, bajar=bajar) is False
    err = capsys.readouterr().err
    assert "[falla]" in err
    assert DIGEST in err
    assert rel in err
    assert "inspección" in err
    assert (tmp_path / rel).read_bytes() == b"otro"


def test_una_huella_previa_incorrecta_vuelve_a_bajar(tmp_path):
    rel = "results/certification/certificado_resumen.json"
    destino = tmp_path / rel
    _escribir(destino, b"viejo")

    def bajar(url, destino_nuevo):
        destino_nuevo.write_bytes(PAYLOAD)

    assert _MOD.asegurar(_item(rel), tmp_path, bajar=bajar)
    assert destino.read_bytes() == PAYLOAD


def test_un_fallo_de_curl_se_informa(tmp_path, capsys):
    rel = "external/modelos/h7_esm2_embeddings.tsv"

    def bajar(url, destino):
        raise subprocess.CalledProcessError(22, ["curl", url])

    assert _MOD.asegurar(_item(rel), tmp_path, bajar=bajar) is False
    assert "[falla]" in capsys.readouterr().err
    assert not (tmp_path / rel).exists()


def test_sin_huella_sellada_no_rebaja_si_el_tamano_coincide(tmp_path):
    item = {
        "rel": "external/pxd064655/Sample_Metadata.txt",
        "url": "https://example.invalid/m",
        "sha256": "",
        "bytes": 3,
    }
    destino = tmp_path / item["rel"]
    _escribir(destino, b"abc")
    llamadas: list[str] = []
    assert _MOD.asegurar(item, tmp_path, bajar=lambda url, dest: llamadas.append(url))
    assert llamadas == []


def test_solo_baja_la_ruta_pedida_y_usa_la_url_publicada(tmp_path, monkeypatch):
    rel = "results/certified/inference_audit.json"
    item = next(fila for fila in _MOD.ZENODO_DESCARGAS if fila["rel"] == rel)
    vistas: list[str] = []

    def bajar(url, destino):
        vistas.append(url)
        _escribir(destino, b"x")

    monkeypatch.setattr(_MOD, "descargar", bajar)
    monkeypatch.setattr(_MOD, "sha256", lambda path: item["sha256"])
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, rel))
    assert _MOD.main() == 0
    assert vistas == [item["url"]]
    assert not (tmp_path / "data/raw/cptac/cptac_crc_matrix.tsv").exists()
    assert not (tmp_path / "external/modelos/h7_esm2_embeddings.tsv").exists()


def test_main_rechaza_una_huella_que_no_coincide(tmp_path, monkeypatch, capsys):
    rel = "data/raw/cptac/cptac_crc_matrix.tsv"
    item = next(fila for fila in _MOD.ZENODO_DESCARGAS if fila["rel"] == rel)

    def bajar(url, destino):
        _escribir(destino, b"matriz-falsa")

    monkeypatch.setattr(_MOD, "descargar", bajar)
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, rel))
    assert _MOD.main() == 1
    err = capsys.readouterr().err
    assert item["sha256"] in err
    assert rel in err
    assert (tmp_path / rel).read_bytes() == b"matriz-falsa"


def test_main_no_rebaja_cuando_el_archivo_ya_tiene_la_huella(tmp_path, monkeypatch):
    rel = "external/modelos/h7_esm2_embeddings.tsv"
    item = next(fila for fila in _MOD.ZENODO_DESCARGAS if fila["rel"] == rel)
    destino = tmp_path / rel
    _escribir(destino, b"ya-esta")
    llamadas: list[str] = []
    monkeypatch.setattr(_MOD, "descargar", lambda url, dest: llamadas.append(url))
    monkeypatch.setattr(_MOD, "sha256", lambda path: item["sha256"])
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, rel))
    assert _MOD.main() == 0
    assert llamadas == []
    assert destino.read_bytes() == b"ya-esta"


def test_al_bajar_cptac_se_recuerda_la_atribucion(tmp_path, monkeypatch, capsys):
    rel = "data/raw/cptac/cptac_crc_matrix.tsv"
    item = next(fila for fila in _MOD.ZENODO_DESCARGAS if fila["rel"] == rel)

    def bajar(url, destino):
        _escribir(destino, b"x")

    monkeypatch.setattr(_MOD, "descargar", bajar)
    monkeypatch.setattr(_MOD, "sha256", lambda path: item["sha256"])
    monkeypatch.setattr(sys, "argv", _argv(tmp_path, rel))
    assert _MOD.main() == 0
    salida = capsys.readouterr().out
    assert "PDC000116" in salida
    assert "Vasaikar et al. 2019" in salida
    assert "10.5281/zenodo.24398807" in salida
    assert "atribución" in salida
