"""Extremo a extremo de la etapa v3_7 con datos que no son la cohorte.

La etapa de suero rechaza cualquier informe que no sea el DIA-NN sellado
(213 313 241 bytes y su SHA-256). Un report sintético tiene que detener la
corrida ahí. El cierre, por su lado, emite 22 casillas cuando recibe 22
fuentes que existen.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.e2e


def test_un_report_sintetico_detiene_la_etapa(tmp_path, monkeypatch):
    from oncogo_crc.settings import load_settings
    from oncogo_crc.stages import v3_7
    from oncogo_crc.stages import v3_3_suero, v3_3_tejido

    report = tmp_path / "data" / "raw" / "plasma" / "DIA-NN_report.tsv"
    report.parent.mkdir(parents=True)
    report.write_text("Run\tProtein.Group\nCP_01_DIA\tP00000\n", encoding="utf-8")

    def _no_exige(_sello):
        return None

    monkeypatch.setattr(v3_7, "exigir_arbol_limpio", _no_exige)
    monkeypatch.setattr(v3_3_suero, "exigir_arbol_limpio", _no_exige)
    monkeypatch.setattr(v3_3_tejido, "exigir_arbol_limpio", _no_exige)
    monkeypatch.setattr(v3_7, "correr_tests", lambda _root, _salida: 0)

    settings = load_settings(data_workspace=tmp_path)
    with pytest.raises(RuntimeError, match="no coincide con la enmienda"):
        v3_7.run(settings)


def test_el_certificado_tiene_22_casillas(tmp_path):
    from oncogo_crc.v3.cierre_v37 import emitir

    casillas = []
    for indice in range(22):
        ruta = tmp_path / f"casilla_{indice}.json"
        ruta.write_text("{}\n", encoding="utf-8")
        casillas.append(
            {
                "nombre": f"c{indice}",
                "fuente": "sintetica",
                "criterio": "la fuente existe",
                "valor": indice,
                "estado": "pasa",
                "ruta": ruta,
            }
        )
    certificado = emitir(casillas, [])
    assert certificado["estado"] == "publicable"
    assert len(certificado["casillas"]) == 22
    assert all(casilla["estado"] == "pasa" for casilla in certificado["casillas"])
