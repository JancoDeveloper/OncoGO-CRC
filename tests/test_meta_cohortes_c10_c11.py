"""Tests de C10 y C11: parseo de la tercera cohorte y filtro de independencia.

Las matrices de este fichero son sintéticas y viven solo aquí, con el prefijo
que exige el preregistro. Ninguna entra en results/.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from oncogo_crc.stages.v2_meta import (
    _contraste_019504,
    _cohortes_permitidas,
    _diagnostico_tecnologia,
    _matriz_maxquant,
    _por_uniprot,
)


def _escribir_proteingroups(path, *, sufijos, n_prot=40, semilla=7, con_banderas=True):
    rng = np.random.default_rng(semilla)
    columnas = {
        "Majority protein IDs": [f"sp|P{i:05d}|GEN{i}_HUMAN" for i in range(n_prot)],
    }
    for s in sufijos:
        base = rng.lognormal(20, 1, size=n_prot)
        columnas[f"LFQ intensity {s}"] = base
    if con_banderas:
        columnas["Reverse"] = [None] * n_prot
        columnas["Potential contaminant"] = [None] * n_prot
        # Las dos últimas son basura de MaxQuant y tienen que desaparecer.
        columnas["Reverse"][-1] = "+"
        columnas["Potential contaminant"][-2] = "+"
    pd.DataFrame(columnas).to_csv(path, sep="\t", index=False)


def test_se_descartan_las_banderas_de_maxquant(tmp_path):
    path = tmp_path / "SINTETICO_proteinGroups.txt"
    _escribir_proteingroups(path, sufijos=["CRC1N", "CRC1T"], n_prot=40)
    matriz, pacientes = _matriz_maxquant(path)
    assert matriz.shape[0] == 38
    assert set(pacientes) == {"CRC1"}
    assert set(pacientes["CRC1"]) == {"N", "T"}


def test_el_tumor_de_pxd019504_es_t_y_no_m(tmp_path):
    """El detalle que invertiría el signo de toda la cohorte si se confundiera."""
    path = tmp_path / "SINTETICO_proteinGroups.txt"
    _escribir_proteingroups(path, sufijos=["CRC1M", "CRC1N"], n_prot=40)
    _matriz, pacientes = _matriz_maxquant(path)
    # M no es una etiqueta de esta cohorte, así que el paciente queda incompleto.
    assert pacientes.get("CRC1", {}).get("T") is None
    with pytest.raises(RuntimeError, match="el tumor es T, no M"):
        _contraste_019504(path, lfc=0.585, alpha=0.05)


def test_el_cero_de_maxquant_no_se_convierte_en_log2(tmp_path):
    path = tmp_path / "SINTETICO_proteinGroups.txt"
    tabla = pd.DataFrame(
        {
            "Majority protein IDs": ["sp|P00001|A_HUMAN", "sp|P00002|B_HUMAN"],
            "LFQ intensity CRC1N": [0.0, 1024.0],
            "LFQ intensity CRC1T": [2048.0, 0.0],
        }
    )
    tabla.to_csv(path, sep="\t", index=False)
    matriz, _ = _matriz_maxquant(path)
    assert np.isnan(matriz.loc["P00001", "LFQ intensity CRC1N"])
    assert matriz.loc["P00002", "LFQ intensity CRC1N"] == pytest.approx(10.0)
    assert np.isnan(matriz.loc["P00002", "LFQ intensity CRC1T"])


def test_el_contraste_recupera_un_efecto_conocido(tmp_path):
    """Un desplazamiento de +2 en log2 tiene que salir como log2fc cercano a 2."""
    path = tmp_path / "SINTETICO_proteinGroups.txt"
    rng = np.random.default_rng(20260914)
    n_prot, n_pac = 60, 8
    filas = {"Majority protein IDs": [f"sp|P{i:05d}|G{i}_HUMAN" for i in range(n_prot)]}
    base = rng.lognormal(20, 0.3, size=(n_prot, n_pac))
    for j in range(n_pac):
        filas[f"LFQ intensity CRC{j}N"] = base[:, j]
        # Las 10 primeras suben 4 veces (2 en log2) en el tumor.
        factor = np.where(np.arange(n_prot) < 10, 4.0, 1.0)
        filas[f"LFQ intensity CRC{j}T"] = base[:, j] * factor
    pd.DataFrame(filas).to_csv(path, sep="\t", index=False)
    tabla, n_pares = _contraste_019504(path, lfc=0.585, alpha=0.05)
    assert n_pares == 8
    arriba = tabla.set_index("protein_id").loc[[f"P{i:05d}" for i in range(10)]]
    assert arriba["log2fc"].median() == pytest.approx(2.0, abs=0.01)
    assert bool(arriba["significant"].all())


def test_por_uniprot_se_queda_con_el_error_estandar_menor():
    tabla = pd.DataFrame(
        {
            "protein_id": ["P00001", "P00001", "P00002"],
            "log2fc": [1.0, 2.0, 3.0],
            "se": [0.5, 0.1, 0.2],
        }
    )
    res = _por_uniprot(tabla)
    assert len(res) == 2
    assert res.set_index("uniprot").loc["P00001", "log2fc"] == 2.0


def test_por_uniprot_descarta_errores_estandar_no_validos():
    tabla = pd.DataFrame(
        {
            "protein_id": ["P00001", "P00002", "P00003", "P00004"],
            "log2fc": [1.0, 2.0, 3.0, np.nan],
            "se": [0.0, np.inf, 0.2, 0.3],
        }
    )
    res = _por_uniprot(tabla)
    assert set(res["uniprot"]) == {"P00003"}


def test_sin_auditoria_de_independencia_no_entra_ninguna_cohorte(tmp_path):
    class _S:
        def resolve(self, _):
            return tmp_path

    permitidas, resumen = _cohortes_permitidas(_S())
    assert permitidas == []
    assert resumen["estado"] == "no_disponible"


def test_una_cohorte_con_ficheros_compartidos_queda_fuera(tmp_path):
    base = tmp_path / "independencia"
    base.mkdir()
    (base / "estado.json").write_text(
        json.dumps(
            {
                "nivel_a_hashes_de_pride": [
                    {"par": "PXD019504_vs_PXD034575", "n_hashes_compartidos": 3}
                ],
                "nivel_b_correlacion_entre_muestras": {
                    "PXD019504_vs_PXD034575": {
                        "estado": "evaluado",
                        "n_pares_sospechosos": 0,
                        "r_maximo": 0.8,
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    class _S:
        def resolve(self, _):
            return tmp_path

    permitidas, resumen = _cohortes_permitidas(_S())
    assert permitidas == []
    assert resumen["nivel_a_sin_ficheros_compartidos"] is False
    assert "no entra" in resumen["veredicto"]


def test_una_muestra_repetida_tambien_deja_fuera_la_cohorte(tmp_path):
    base = tmp_path / "independencia"
    base.mkdir()
    (base / "estado.json").write_text(
        json.dumps(
            {
                "nivel_a_hashes_de_pride": [
                    {"par": "PXD019504_vs_PXD034575", "n_hashes_compartidos": 0}
                ],
                "nivel_b_correlacion_entre_muestras": {
                    "PXD019504_vs_PXD034575": {
                        "estado": "evaluado",
                        "n_pares_sospechosos": 2,
                        "r_maximo": 0.999,
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    class _S:
        def resolve(self, _):
            return tmp_path

    permitidas, resumen = _cohortes_permitidas(_S())
    assert permitidas == []
    assert resumen["nivel_b_sin_muestras_repetidas"] is False


def test_el_diagnostico_de_tecnologia_avisa_cuando_la_caida_es_grande():
    con_tres = pd.DataFrame({"uniprot": ["A", "B", "C"], "i2": [0.9, 0.92, 0.88]})
    solo_lfq = pd.DataFrame(
        {"uniprot": ["A", "B", "C"], "i2_solo_lfq": [0.2, 0.25, 0.1], "tau2_solo_lfq": [0.1] * 3}
    )
    res = _diagnostico_tecnologia(con_tres, solo_lfq, ["a", "b"])
    assert res["caida_al_quitar_cptac"] > 0.2
    assert "diferencia de tecnología" in res["advertencia"]

    sin_caida = pd.DataFrame(
        {"uniprot": ["A", "B", "C"], "i2_solo_lfq": [0.85, 0.9, 0.86], "tau2_solo_lfq": [0.1] * 3}
    )
    res2 = _diagnostico_tecnologia(con_tres, sin_caida, ["a", "b"])
    assert res2["caida_al_quitar_cptac"] < 0.2
    assert "apenas cambia" in res2["advertencia"]
