"""El contaminante se cruza por acceso y especie, no por el símbolo GN."""

from pathlib import Path

import pandas as pd
import pytest

from oncogo_crc.utils.dia_nn_plasma import acceso_base, leer_libreria, trazabilidad

FASTA = """>sp|Cont_P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4
MKWVTFISLL
>sp|Cont_P04264|K2C1_HUMAN Keratin, type II cytoskeletal 1 OS=Homo sapiens OX=9606 GN=KRT1 PE=1 SV=6
MSRQFSSRSG
>sp|Cont_P05787|K2C8_HUMAN Keratin, type II cytoskeletal 8 OS=Homo sapiens OX=9606 GN=KRT8 PE=1 SV=7
MSIRVTQKSY
>sp|Cont_P00722|BGAL_ECOLI Beta-galactosidase OS=Escherichia coli (strain K12) OX=83333 GN=lacZ PE=1 SV=2
MTMITDSLAV
"""
PLASMA_CLASICO = {"ALB", "APOA1", "C3", "TF", "SERPINA1", "A2M", "APOE", "FGA", "FGB", "FGG", "LRG1"}


def _fila(grupo, ids, nombres, genes):
    return {
        "Run": "CP_1",
        "Protein.Group": grupo,
        "Protein.Ids": ids,
        "Protein.Names": nombres,
        "Genes": genes,
        "PG.Q.Value": 0.001,
        "Global.PG.Q.Value": 0.001,
    }


def test_acceso_base():
    assert acceso_base("Cont_P02768") == "P02768"
    assert acceso_base("P02768-2") == "P02768"
    assert acceso_base("sp|Cont_P02769|ALBU_BOVIN") == "P02769"


def test_albumina_humana_no_es_contaminante(tmp_path):
    fasta = tmp_path / "c.fasta"
    fasta.write_text(FASTA, encoding="utf-8")
    lib = leer_libreria(fasta)
    assert lib.n_entradas == 4
    assert lib.genes_humanos == {"KRT1", "KRT8"}
    tabla = pd.DataFrame([
        _fila("P02768", "P02768", "ALBU_HUMAN", "ALB"),
        _fila("P02769", "Cont_P02769", "ALBU_BOVIN", "ALB"),
        _fila("P04264", "P04264", "K2C1_HUMAN", "KRT1"),
        _fila("P05787", "P05787", "K2C8_HUMAN", "KRT8"),
        _fila("P02768;P02769", "P02768;Cont_P02769", "ALBU_HUMAN;ALBU_BOVIN", "ALB"),
    ])
    traza = trazabilidad(tabla, libreria=lib).set_index("protein_group")
    assert not traza.at["P02768", "contaminante"]
    assert traza.at["P02769", "contaminante"]
    assert traza.at["P04264", "contaminante"]
    assert traza.at["P05787", "contaminante"] and traza.at["P05787", "queratina_epitelial"]
    assert not traza.at["P02768;P02769", "contaminante"]
    assert traza.at["P02768;P02769", "comparte_con_contaminante"]


def test_libreria_real_no_trae_plasma_humano():
    fasta = Path("external/frankenfield_2022/0602_Universal_Contaminants.fasta")
    if not fasta.is_file():
        fasta = Path("/Volumes/JANCODATA/OncoGO-CRC/external/frankenfield_2022/0602_Universal_Contaminants.fasta")
    if not fasta.is_file():
        pytest.skip("la librería de Frankenfield no está en este equipo")
    lib = leer_libreria(fasta)
    assert lib.n_entradas == 381
    assert lib.por_taxon.get("9606") == 151
    assert not (lib.genes_humanos & PLASMA_CLASICO)
