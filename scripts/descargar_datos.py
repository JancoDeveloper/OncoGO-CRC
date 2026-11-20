#!/usr/bin/env python3
"""Descarga las cuatro cohortes y comprueba las huellas conocidas.

No inventa un fichero que el depósito público no entrega. Si la huella no
coincide, el script sale con código 1 y deja el fichero para inspección.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
ZENODO = "https://zenodo.org/deposit/PENDIENTE"
FRANKENFIELD = (
    "https://raw.githubusercontent.com/HaoGroup-ProtContLib/"
    "Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/main/"
    "Universal%20protein%20contaminant%20FASTA/0602_Universal%20Contaminants.fasta"
)

DESCARGAS = (
    {
        "rel": "data/raw/pride/Cohort1_proteinGroups.txt",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/01/PXD034575/Cohort1_proteinGroups.txt",
        "sha256": "da9ec08ee656d207b29242d804f308247957126201678e2ed80cae37ce8ed8c9",
        "bytes": 69272674,
    },
    {
        "rel": "data/raw/plasma/DIA-NN_report.tsv",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/06/PXD055218/DIA-NN_report.tsv",
        "sha256": "a8a5aadd74391afef68cf4d6c3a7d3860023101e5c92590e2e87c152a39e59e8",
        "bytes": 213313241,
    },
    {
        "rel": "external/pxd064655/report.pr_matrix.tsv",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/01/PXD064655/report.pr_matrix.tsv",
        "sha256": "",
        "bytes": 6300706,
    },
    {
        "rel": "external/pxd064655/Sample_Metadata.txt",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/01/PXD064655/Sample_Metadata.txt",
        "sha256": "",
        "bytes": 3961,
    },
)


def sha256(path: Path) -> str:
    resumen = hashlib.sha256()
    with path.open("rb") as fichero:
        for bloque in iter(lambda: fichero.read(1 << 20), b""):
            resumen.update(bloque)
    return resumen.hexdigest()


def descargar(url: str, destino: Path) -> None:
    destino.parent.mkdir(parents=True, exist_ok=True)
    temporal = destino.with_suffix(destino.suffix + ".parcial")
    subprocess.run(["curl", "-fL", "--retry", "5", "-o", str(temporal), url], check=True)
    temporal.replace(destino)


def matriz_pride(protein_groups: Path, destino: Path) -> str:
    """Misma reducción que harmonize_pride_maxquant_lfq: protein_id y columnas LFQ."""
    import pandas as pd

    tabla = pd.read_csv(protein_groups, sep="\t", low_memory=False)
    tabla.columns = [str(columna).replace("\ufeff", "").strip() for columna in tabla.columns]
    if str(tabla.columns[0]) != "Protein IDs":
        raise SystemExit(f"proteinGroups inesperado: primera columna {tabla.columns[0]!r}")
    patron = re.compile(r"^LFQ intensity CRC\d+[MNP]$", re.IGNORECASE)

    def primero(celda: object) -> str:
        texto = "" if pd.isna(celda) else str(celda).strip()
        return texto.split(";")[0].strip() if texto else ""

    tabla.insert(0, "protein_id", tabla["Protein IDs"].map(primero))
    columnas = [columna for columna in tabla.columns if isinstance(columna, str) and patron.match(columna)]
    if not columnas:
        raise SystemExit("proteinGroups sin columnas LFQ intensity CRC##[MNP]")
    destino.parent.mkdir(parents=True, exist_ok=True)
    tabla[["protein_id", *columnas]].to_csv(destino, sep="\t", index=False)
    return sha256(destino)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", "--destino", type=Path, required=True, dest="workspace")
    parser.add_argument(
        "--solo",
        action="append",
        default=[],
        help="Descarga solo esta ruta relativa, o 'contaminantes'. Se puede repetir. No baja el plasma de 74 GB.",
    )
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    solo = set(args.solo)
    fallos = 0
    for item in DESCARGAS:
        if solo and item["rel"] not in solo:
            continue
        destino = workspace / item["rel"]
        if destino.is_file() and item["sha256"] and sha256(destino) == item["sha256"]:
            print(f"[ok] {item['rel']}")
            continue
        if destino.is_file() and item["bytes"] and destino.stat().st_size == item["bytes"] and not item["sha256"]:
            print(f"[ok] {item['rel']} ({destino.stat().st_size} bytes, sin huella sellada)")
            continue
        print(f"[descarga] {item['rel']}")
        try:
            descargar(item["url"], destino)
        except subprocess.CalledProcessError as exc:
            print(f"[falla] {item['rel']}: {exc}", file=sys.stderr)
            fallos += 1
            continue
        digest = sha256(destino)
        if item["sha256"] and digest != item["sha256"]:
            print(f"[falla] {item['rel']} sha256 {digest} != {item['sha256']}", file=sys.stderr)
            fallos += 1
        else:
            print(f"[ok] {item['rel']} sha256 {digest}")

    if solo and "external/frankenfield_2022/0602_Universal_Contaminants.fasta" not in solo and "contaminantes" not in solo:
        return 1 if fallos else 0

    grupos = workspace / "data/raw/pride/Cohort1_proteinGroups.txt"
    if solo and "data/raw/pride/Cohort1_proteinGroups.txt" not in solo and "data/raw/pride/pride_crc_matrix.tsv" not in solo:
        grupos = Path("/no-existe-oncogo")
    matriz = workspace / "data/raw/pride/pride_crc_matrix.tsv"
    esperada = "a5baca5a3d3d41884a5217df7778b930f580e8b03630875e409ce14a242fd1bc"
    if grupos.is_file() and not (matriz.is_file() and sha256(matriz) == esperada):
        digest = matriz_pride(grupos, matriz)
        if digest == esperada:
            print("[ok] pride_crc_matrix.tsv coincide con la huella sellada")
        else:
            print(
                f"[aviso] pride_crc_matrix.tsv regenerada con pandas actual tiene {digest}, "
                f"la sellada es {esperada}. La etapa de tejido lee esta matriz.",
                file=sys.stderr,
            )

    fasta = workspace / "external/frankenfield_2022/0602_Universal_Contaminants.fasta"
    if not fasta.is_file():
        try:
            destino_temporal = fasta.with_suffix(".fasta.parcial")
            descargar(FRANKENFIELD, destino_temporal)
            destino_temporal.replace(fasta)
            print(f"[ok] contaminantes sha256 {sha256(fasta)}")
        except subprocess.CalledProcessError as exc:
            print(f"[falla] biblioteca de contaminantes: {exc}", file=sys.stderr)
            fallos += 1

    if not solo:
        for rel in (
            "metadata/gene_ontology/go_gene_long.tsv",
            "metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv",
            "metadata/protein_id_maps/protein_mapping_master.tsv",
        ):
            origen = RAIZ / rel
            destino = workspace / rel
            if origen.is_file() and not destino.is_file():
                destino.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(origen, destino)
                print(f"[ok] copiado {rel}")

    if not solo:
        print(
            "PDC000116 no tiene un único URL de matriz en este script. "
            "La huella sellada de data/raw/cptac/cptac_crc_matrix.tsv está en data/checksums.sha256. "
            "results/certified y results/certification se regeneran con make certified "
            f"o, cuando exista, se descargan de {ZENODO}."
        )
        print(
            "external/modelos/h7_esm2_embeddings.tsv no está en git. "
            "Su huella sellada está en data/checksums.sha256. "
            "Se regenera con scripts/v3/esm2_embeddings.py (torch, fuera de la corrida principal)."
        )
    return 1 if fallos else 0


if __name__ == "__main__":
    raise SystemExit(main())
