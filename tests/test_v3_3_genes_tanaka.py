import pandas as pd

from oncogo_crc.stages.v3_3_controles import _genes_bajos
from oncogo_crc.v3.genes_tanaka import base, gen_unico, genes_de, mapa_genes


class _Ajustes:
    def __init__(self, root):
        self.root = root

    def resolve(self, rel):
        return self.root / rel


def test_el_mapa_encuentra_el_gen_por_las_tres_vias(tmp_path):
    grupos = tmp_path / "grupos.tsv"
    pd.DataFrame([
        {"Protein IDs": "P00915", "Gene names": "CA1", "Fasta headers": ""},
        {"Protein IDs": "Q9UNE7", "Gene names": "", "Fasta headers": "sp|Q9UNE7|AQP8_HUMAN GN=AQP8 PE=1"},
        {"Protein IDs": "A0A000", "Gene names": "", "Fasta headers": ""},
        {"Protein IDs": "P99999", "Gene names": "FOO;BAR", "Fasta headers": ""},
    ]).to_csv(grupos, sep="\t", index=False)
    uni = tmp_path / "uni.tsv"
    pd.DataFrame({"gene_symbol": ["SLC26A3"], "uniprot_accession": ["A0A000"]}).to_csv(uni, sep="\t", index=False)
    mapa = mapa_genes(grupos, uni)
    assert "CA1" in genes_de("P00915", mapa)
    assert gen_unico("Q9UNE7", mapa) == "AQP8"
    assert gen_unico("A0A000", mapa) == "SLC26A3"
    assert gen_unico("P99999", mapa) is None
    assert base("CON__P02768") == "P02768"


def test_control_positivo_no_evaluable_es_falla(tmp_path):
    # Enmienda v3.3: un control positivo que no se encuentra ya no pasa en silencio.
    grupos = tmp_path / "data/raw/pride"
    grupos.mkdir(parents=True)
    pd.DataFrame([{"Protein IDs": "P00915", "Gene names": "OTRO", "Fasta headers": ""}]).to_csv(
        grupos / "Cohort1_proteinGroups.txt", sep="\t", index=False,
    )
    uni = tmp_path / "metadata/protein_id_maps"
    uni.mkdir(parents=True)
    pd.DataFrame({"gene_symbol": ["ZZZ"], "uniprot_accession": ["P00000"]}).to_csv(
        uni / "gene_uniprot_human_reviewed_bulk.tsv", sep="\t", index=False,
    )
    salida = tmp_path / "salida"
    salida.mkdir()
    pd.DataFrame({"protein_id": ["P00915"], "log2fc": [-1.0], "n_obs": [8]}).to_csv(
        salida / "k1_primario_frente_a_mucosa.tsv", sep="\t", index=False,
    )
    filas = _genes_bajos(_Ajustes(tmp_path), salida)
    assert filas
    assert all(fila["pasa"] == "falla" for fila in filas)
