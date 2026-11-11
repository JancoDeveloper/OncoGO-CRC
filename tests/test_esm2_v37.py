"""ESM-2 se calcula fuera de la corrida y la corrida solo lee el TSV."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.v3 import secuencia_features as sec

RUTA = Path(__file__).resolve().parents[1] / "scripts" / "v3" / "esm2_embeddings.py"


def script():
    spec = importlib.util.spec_from_file_location("esm2_embeddings", RUTA)
    modulo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(modulo)
    return modulo


def _paquete(tmp_path: Path, *, hash_modelo: str | None = None, hash_fasta: str | None = None,
             longitud: int = 4, gen_extra: str | None = None) -> dict:
    fasta = tmp_path / "humano.fasta"
    fasta.write_text(">sp|P00001|X_HUMAN GN=ABC\nACDE\n", encoding="utf-8")
    modelo = tmp_path / "modelo"
    modelo.mkdir()
    (modelo / "peso.txt").write_text("abc", encoding="utf-8")
    columnas = [f"esm2_{i:03d}" for i in range(480)]
    fila = {"gen": "ABC", "longitud": longitud, "n_ventanas": 1, **{columna: 0.1 for columna in columnas}}
    tsv = tmp_path / "h7_esm2_embeddings.tsv"
    pd.DataFrame([fila]).to_csv(tsv, sep="\t", index=False, float_format="%.9f")
    omitidas = tmp_path / "h7_esm2_omitidas.tsv"
    if gen_extra:
        omitidas.write_text(f"gen\tlongitud\tmotivo\n{gen_extra}\t\tsin secuencia en el FASTA revisado\n", encoding="utf-8")
    else:
        omitidas.write_text("gen\tlongitud\tmotivo\n", encoding="utf-8")
    sello = {
        "sha256_fasta": hash_fasta if hash_fasta is not None else sec.sha256_fichero(fasta),
        "sha256_modelo": hash_modelo if hash_modelo is not None else sec.sha256_directorio(modelo),
        "sha256_tsv": sec.sha256_fichero(tsv),
    }
    sello_ruta = tmp_path / "h7_esm2_embeddings.tsv.sello.json"
    sello_ruta.write_text(json.dumps(sello), encoding="utf-8")
    return {"tsv": tsv, "sello": sello_ruta, "omitidas": omitidas, "fasta": fasta, "modelo": modelo}


def test_la_corrida_no_importa_torch():
    texto = Path(sec.__file__).read_text(encoding="utf-8")
    assert "import torch" not in texto
    assert "embeddings_esm2" not in texto


def test_un_sello_que_no_coincide_deja_el_secundario_no_disponible(tmp_path):
    for campo, valor in (("modelo", "0" * 64), ("fasta", "1" * 64)):
        carpeta = tmp_path / campo
        carpeta.mkdir()
        hecho = _paquete(
            carpeta,
            hash_modelo=valor if campo == "modelo" else None,
            hash_fasta=valor if campo == "fasta" else None,
        )
        frame, info = sec.leer_esm2(
            hecho["tsv"], hecho["sello"], hecho["omitidas"], hecho["fasta"], hecho["modelo"],
            {"ABC"}, {"ABC": "ACDE"},
        )
        assert frame is None
        assert info["esm2"] == "no_disponible"
        assert "no coincide" in info["esm2_motivo"]


def test_un_tsv_coherente_si_se_lee(tmp_path):
    hecho = _paquete(tmp_path)
    frame, info = sec.leer_esm2(
        hecho["tsv"], hecho["sello"], hecho["omitidas"], hecho["fasta"], hecho["modelo"],
        {"ABC"}, {"ABC": "ACDE"},
    )
    assert info["esm2"] == "calculado"
    assert info["n_con_esm2"] == 1
    assert info["n_omitidas_esm2"] == 0
    assert frame is not None and frame.shape == (1, 480)


def test_una_proteina_sin_vector_ni_omision_no_se_lee(tmp_path):
    hecho = _paquete(tmp_path)
    frame, info = sec.leer_esm2(
        hecho["tsv"], hecho["sello"], hecho["omitidas"], hecho["fasta"], hecho["modelo"],
        {"ABC", "DEF"}, {"ABC": "ACDE"},
    )
    assert frame is None
    assert info["esm2"] == "no_disponible"


def test_un_error_queda_en_las_omitidas(tmp_path):
    mod = script()

    class Tokens:
        def __init__(self):
            import torch
            self.input_ids = torch.zeros(1, 3, dtype=torch.long)

        def __getitem__(self, clave):
            return self.input_ids

    class Tokenizador:
        def __call__(self, secuencia, return_tensors="pt"):
            return Tokens()

    _linea, omision = mod.procesar_gen("GEN", "ACDE", None, Tokenizador(), lambda: None)
    assert omision is not None and omision.startswith("GEN\t4\tValueError:")
    destino = tmp_path / "h7_esm2_omitidas.tsv"
    mod.anexar(destino, omision, "gen\tlongitud\tmotivo")
    texto = destino.read_text(encoding="utf-8")
    assert "GEN" in texto and "ValueError" in texto


def test_la_reanudacion_no_duplica_ni_salta(tmp_path):
    mod = script()
    vector = [0.25] * 480
    parcial = tmp_path / "parcial.tsv"
    encabezado = "\t".join(mod.COLUMNAS)
    mod.anexar(parcial, mod.texto_vector("AAA", 4, 1, vector), encabezado)
    hechos = mod.leer_parcial(parcial, mod.ANCHO)
    assert mod.genes_pendientes(["AAA", "BBB", "AAA"], set(hechos)) == ["BBB"]
    mod.anexar(parcial, mod.texto_vector("BBB", 2044, 2, vector), encabezado)
    hechos = mod.leer_parcial(parcial, mod.ANCHO)
    assert set(hechos) == {"AAA", "BBB"}
    final = tmp_path / "final.tsv"
    mod.escribir_final(final, encabezado, hechos)
    genes = [linea.split("\t")[0] for linea in final.read_text(encoding="utf-8").splitlines()[1:]]
    assert genes == ["AAA", "BBB"]
    mod.anexar(parcial, mod.texto_vector("AAA", 4, 1, vector), encabezado)
    with pytest.raises(RuntimeError, match="repetido"):
        mod.leer_parcial(parcial, mod.ANCHO)


def test_usa_cpu_aunque_mps_este_disponible(monkeypatch):
    import torch
    monkeypatch.setattr(torch.backends.mps, "is_available", lambda: True)
    mod = script()
    assert mod.dispositivo_de_calculo().type == "cpu"
    with pytest.raises(RuntimeError, match="CPU"):
        mod.comprobar_cpu(torch.device("mps"))


@pytest.fixture(scope="module")
def modelo_real():
    import torch
    from oncogo_crc.settings import load_settings
    anterior = torch.are_deterministic_algorithms_enabled()
    hilos = torch.get_num_threads()
    mod = script()
    mod.preparar(4)
    settings = load_settings()
    ruta = settings.resolve(sec.MODELO_ESM)
    modelo, tokenizador = mod.cargar_modelo(ruta)
    fabrica = mod.fabrica_atencion()
    yield mod, modelo, tokenizador, fabrica
    torch.use_deterministic_algorithms(anterior)
    torch.set_num_threads(hilos)


def test_una_secuencia_larga_da_480_y_al_menos_dos_ventanas(modelo_real):
    mod, modelo, tokenizador, fabrica = modelo_real
    vector, n_ventanas = mod.vector_de_secuencia(modelo, tokenizador, "A" * 1023, fabrica)
    assert n_ventanas >= 2
    assert vector.shape == (480,)
    assert len(sec.ventanas("A" * 1023)) >= 2


def test_dos_llamadas_dan_el_mismo_vector(modelo_real):
    mod, modelo, tokenizador, fabrica = modelo_real
    secuencia = "ACDEFGHIKLMNPQRSTVWY"
    uno, _ = mod.vector_de_secuencia(modelo, tokenizador, secuencia, fabrica)
    dos, _ = mod.vector_de_secuencia(modelo, tokenizador, secuencia, fabrica)
    assert np.array_equal(uno, dos)
