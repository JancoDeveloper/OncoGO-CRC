"""Parser del suero. El report sintético no es una cohorte. El real se salta si no está."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import pytest

from oncogo_crc.utils.dia_nn_plasma import (
    BYTES_REPORT,
    SHA256_REPORT,
    controles_del_report,
    especies_de_nombres,
    matriz_filtrada,
    postprocess_dia_nn_plasma_download,
)

REPORT = Path("/Volumes/JANCODATA/OncoGO-CRC/data/raw/plasma/DIA-NN_report.tsv")


def _report_sintetico() -> pd.DataFrame:
    filas = []
    for corrida, q_pg in (("CP_1", 0.001), ("CN_1", 0.2)):
        filas.append({
            "Run": corrida,
            "Protein.Group": "HUMANO",
            "Protein.Ids": "P1",
            "Protein.Names": "ALBU_HUMAN",
            "Genes": "ALB",
            "PG.MaxLFQ": 10.0,
            "Genes.MaxLFQ": 10.0,
            "Q.Value": 0.001,
            "PG.Q.Value": q_pg,
            "Global.PG.Q.Value": 0.001,
            "GG.Q.Value": 0.001,
            "Proteotypic": 1,
        })
    filas.append({
        "Run": "CP_1",
        "Protein.Group": "RATON",
        "Protein.Ids": "P2",
        "Protein.Names": "ALBU_MOUSE",
        "Genes": "Alb",
        "PG.MaxLFQ": 8.0,
        "Genes.MaxLFQ": 8.0,
        "Q.Value": 0.001,
        "PG.Q.Value": 0.001,
        "Global.PG.Q.Value": 0.001,
        "GG.Q.Value": 0.001,
        "Proteotypic": 1,
    })
    return pd.DataFrame(filas)


def test_el_filtro_deja_fuera_al_raton_y_la_corrida_que_no_pasa():
    assert especies_de_nombres("ALBU_HUMAN;TRYP_PIG") == ["HUMAN", "PIG"]
    ancha = matriz_filtrada(_report_sintetico())
    assert list(ancha["protein_id"]) == ["HUMANO"]
    assert ancha.loc[0, "CP_1"] == 10.0
    assert pd.isna(ancha.loc[0, "CN_1"])


def test_el_postproceso_no_toca_el_fichero(tmp_path: Path):
    origen = tmp_path / "SINTETICO_report.tsv"
    origen.write_text("Run\tProtein.Group\nCP_1\tHUMANO\n", encoding="utf-8")
    antes = origen.read_bytes()
    postprocess_dia_nn_plasma_download(origen, tmp_path / "meta.tsv")
    assert origen.read_bytes() == antes
    meta = pd.read_csv(tmp_path / "meta.tsv", sep="\t")
    assert meta.loc[0, "condition"] == "crc"


@pytest.mark.datos_reales
def test_el_report_real_reproduce_los_controles_de_la_revision():
    if not REPORT.is_file():
        pytest.skip("DIA-NN_report.tsv no está en el disco")
    if REPORT.stat().st_size != BYTES_REPORT:
        pytest.skip(f"report incompleto: {REPORT.stat().st_size} de {BYTES_REPORT} bytes")
    resumen = hashlib.sha256()
    with REPORT.open("rb") as entrada:
        while True:
            bloque = entrada.read(1 << 20)
            if not bloque:
                break
            resumen.update(bloque)
    if resumen.hexdigest().upper() != SHA256_REPORT:
        pytest.fail(f"SHA-256 no coincide: {resumen.hexdigest()}")
    from oncogo_crc.utils.dia_nn_plasma import leer_report

    cifras = controles_del_report(leer_report(REPORT))
    assert cifras["filas"] == 283_488
    assert cifras["corridas"] == 60
    assert cifras["q_value_max"] <= 0.01
    assert cifras["grupos_global_pg_q_0_01"] == 329
    assert cifras["grupos_sin_pg_q_0_01_en_ninguna_corrida"] == 432
    assert cifras["humanos_pg_y_global"] == 312
    assert cifras["humanos_pg_y_global_con_gen"] == 307
