"""ONCOGO_DATA_DIR se lee; las salidas van al workspace."""

from __future__ import annotations

import pytest

from oncogo_crc.settings import load_settings


def test_data_dir_separa_lectura_y_escritura(tmp_path, monkeypatch):
    datos = tmp_path / "datos"
    salida = tmp_path / "salida"
    matriz = datos / "data/raw/cptac/cptac_crc_matrix.tsv"
    matriz.parent.mkdir(parents=True)
    matriz.write_text("protein_id\n", encoding="utf-8")
    (datos / "results/certified").mkdir(parents=True)
    monkeypatch.setenv("ONCOGO_DATA_DIR", str(datos))
    monkeypatch.setenv("ONCOGO_CRC_DATA_WORKSPACE", str(salida))
    settings = load_settings()
    assert settings.resolve("data/raw/cptac/cptac_crc_matrix.tsv") == matriz.resolve()
    assert settings.resolve("results/certified/cptac_unshared_treat.tsv").is_relative_to(datos.resolve())
    procesada = settings.resolve("data/processed/plasma/pxd055218_pg_maxlfq_filtrado.tsv")
    assert procesada.is_relative_to(salida.resolve())
    derivados = settings.resolve("external/pxd064655/derivados/pr_matrix_SH.aq_reformat.tsv")
    assert derivados.is_relative_to(salida.resolve())
    assert settings.resolve("results/v3_7/certificado.json").is_relative_to(salida.resolve())


def test_data_dir_y_workspace_no_pueden_ser_el_mismo(tmp_path, monkeypatch):
    monkeypatch.setenv("ONCOGO_DATA_DIR", str(tmp_path))
    monkeypatch.setenv("ONCOGO_CRC_DATA_WORKSPACE", str(tmp_path))
    with pytest.raises(RuntimeError):
        load_settings()
