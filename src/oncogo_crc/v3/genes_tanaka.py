"""Símbolo de gen para cada fila de PXD034575, desde el proteinGroups de MaxQuant."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

_GN = re.compile(r"\bGN=(\S+)")
_ISOFORMA = re.compile(r"-\d+$")


def base(identificador: str) -> str:
    texto = str(identificador).strip()
    for prefijo in ("CON__", "REV__"):
        if texto.startswith(prefijo):
            texto = texto[len(prefijo):]
    if "|" in texto:
        partes = texto.split("|")
        texto = partes[1] if len(partes) > 1 else partes[0]
    return _ISOFORMA.sub("", texto)


def mapa_genes(grupos_path: Path, mapa_uniprot_path: Path) -> dict[str, frozenset[str]]:
    cabecera = set(pd.read_csv(grupos_path, sep="\t", nrows=0).columns)
    if "Protein IDs" not in cabecera:
        raise RuntimeError("proteinGroups sin la columna 'Protein IDs'")
    usar = [c for c in ("Protein IDs", "Gene names", "Fasta headers") if c in cabecera]
    grupos = pd.read_csv(grupos_path, sep="\t", usecols=usar, dtype=str, low_memory=False).fillna("")
    uni = pd.read_csv(mapa_uniprot_path, sep="\t", usecols=["gene_symbol", "uniprot_accession"], dtype=str).dropna()
    uni_a_gen: dict[str, set[str]] = {}
    for gen, acceso in zip(uni["gene_symbol"].str.upper(), uni["uniprot_accession"], strict=True):
        uni_a_gen.setdefault(base(acceso), set()).add(gen)
    mapa: dict[str, frozenset[str]] = {}
    for fila in grupos.to_dict("records"):
        ids = [x.strip() for x in fila["Protein IDs"].split(";") if x.strip()]
        if not ids:
            continue
        genes = {g.strip().upper() for g in fila.get("Gene names", "").split(";") if g.strip()}
        genes |= {g.upper() for g in _GN.findall(fila.get("Fasta headers", ""))}
        for x in ids:
            genes |= uni_a_gen.get(base(x), set())
        if not genes:
            continue
        for clave in (fila["Protein IDs"].strip(), ids[0], base(ids[0])):
            mapa.setdefault(clave, frozenset(genes))
    return mapa


def genes_de(pid: object, mapa: dict[str, frozenset[str]]) -> frozenset[str]:
    texto = str(pid).strip()
    primero = texto.split(";")[0].strip()
    for clave in (texto, primero, base(primero)):
        if clave in mapa:
            return mapa[clave]
    return frozenset()


def gen_unico(pid: object, mapa: dict[str, frozenset[str]]) -> str | None:
    genes = genes_de(pid, mapa)
    return next(iter(genes)) if len(genes) == 1 else None


def filas_de_gen(indice, mapa: dict[str, frozenset[str]], gen: str) -> list[str]:
    buscado = gen.upper()
    return [str(pid) for pid in indice if buscado in genes_de(pid, mapa)]
