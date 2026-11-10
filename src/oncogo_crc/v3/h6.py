"""H6. Origen celular con el HPA si Pelka no tiene la matriz de expresión."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from oncogo_crc.settings import Settings
from oncogo_crc.v3.genes_tanaka import gen_unico, mapa_genes
from oncogo_crc.v3.salida import entrada_v3, volcar_json, volcar_tsv


def clase_de_expresion(media_epitelial: float, maxima_no_epitelial: float) -> str:
    if not np.isfinite(media_epitelial) or not np.isfinite(maxima_no_epitelial):
        return "no_evaluable"
    if media_epitelial >= 2.0 * maxima_no_epitelial:
        return "epitelial"
    if maxima_no_epitelial >= 2.0 * media_epitelial:
        return "no_epitelial"
    return "mixta"


def ejecutar(settings: Settings, salida: Path, sello: dict) -> dict:
    estado_ref = settings.resolve("external/referencia_unicelular/estado.json")
    matriz = settings.resolve("external/referencia_unicelular/crudos/rna_single_cell_type.tsv")
    if not estado_ref.is_file() or not matriz.is_file():
        payload = {"estado": "no_evaluable", "motivo": "no está la leyenda o la matriz de tipo celular del HPA"}
        volcar_json(payload, salida / "h6_origen_celular.json", sello)
        return payload
    ficha = json.loads(estado_ref.read_text(encoding="utf-8"))
    pelka = (ficha.get("referencia_especifica_de_colon") or {}).get("estado")
    tipos = ficha.get("tipos_celulares") or {}
    epiteliales = [t.lower() for t in tipos.get("tipos_epiteliales") or []]
    presentes = [t.lower() for t in tipos.get("tipos_presentes") or []]
    otros = [t for t in presentes if t not in epiteliales]
    if not epiteliales or not otros or pelka == "usada":
        motivo = "Pelka está marcada como usada, pero esta etapa no reabre una matriz que no está" if pelka == "usada" else "la ficha del HPA no separa tipos epiteliales y no epiteliales"
        payload = {"estado": "no_evaluable", "motivo": motivo}
        volcar_json(payload, salida / "h6_origen_celular.json", sello)
        return payload
    tabla = pd.read_csv(matriz, sep="\t", usecols=["Gene name", "Cell type", "nCPM"])
    tabla["tipo"] = tabla["Cell type"].astype(str).str.lower()
    tabla["gen"] = tabla["Gene name"].astype(str).str.upper()
    tabla["ncpm"] = pd.to_numeric(tabla["nCPM"], errors="coerce")
    filas = []
    for gen, bloque in tabla.groupby("gen"):
        epi = bloque.loc[bloque["tipo"].isin(epiteliales), "ncpm"]
        no = bloque.loc[bloque["tipo"].isin(otros), "ncpm"]
        if epi.empty or no.empty:
            continue
        filas.append({
            "gen": gen,
            "media_epitelial": float(epi.mean()),
            "maxima_no_epitelial": float(no.max()),
            "clase": clase_de_expresion(float(epi.mean()), float(no.max())),
        })
    clases = pd.DataFrame(filas)
    if clases.empty:
        payload = {"estado": "no_evaluable", "motivo": "el HPA no produjo una clase para ningún gen"}
        volcar_json(payload, salida / "h6_origen_celular.json", sello)
        return payload
    volcar_tsv(clases, salida / "h6_origen_celular.tsv", sello)
    prueba = _fisher(settings, salida, clases)
    payload = {
        "estado": prueba["estado"],
        "motivo": prueba.get("motivo", ""),
        "fuente": "HPA RNA single cell type",
        "pelka": pelka,
        "tipos_epiteliales": epiteliales,
        "tipos_no_epiteliales": otros,
        "n_genes": int(len(clases)),
        "conteo": {str(k): int(v) for k, v in clases["clase"].value_counts().items()},
        "fisher": {k: v for k, v in prueba.items() if k != "estado"},
        "nota": "Pelka 2021 no tiene la matriz h5 en este disco. La clase sale del HPA y las células tumorales no están en ese atlas.",
    }
    volcar_json(payload, salida / "h6_origen_celular.json", sello)
    return payload


def _fisher(settings: Settings, salida: Path, clases: pd.DataFrame) -> dict:
    cptac = settings.resolve("results/certified/cptac_unshared_treat.tsv")
    tanaka_ruta = entrada_v3(settings, salida, "k1_primario_frente_a_mucosa.tsv")
    grupos = settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt")
    uni = settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv")
    if not cptac.is_file() or not tanaka_ruta.is_file() or not grupos.is_file() or not uni.is_file():
        return {"estado": "no_evaluable", "motivo": "faltan las tablas certificadas o el mapa de genes"}
    cert = pd.read_csv(cptac, sep="\t")
    cert = cert.loc[cert["significant"].astype(bool)].copy()
    cert["gen"] = cert["protein_id"].astype(str).str.upper()
    cert = cert.drop_duplicates("gen").set_index("gen")
    tanaka = pd.read_csv(tanaka_ruta, sep="\t")
    mapa = mapa_genes(grupos, uni)
    tanaka["gen"] = tanaka["protein_id"].map(lambda pid: gen_unico(str(pid), mapa))
    tanaka["log2fc"] = pd.to_numeric(tanaka["log2fc"], errors="coerce")
    tanaka = tanaka.dropna(subset=["gen"]).drop_duplicates("gen").set_index("gen")
    clase = clases.drop_duplicates("gen").set_index("gen")["clase"]
    mismo = {"epitelial": 0, "no_epitelial": 0}
    distinto = {"epitelial": 0, "no_epitelial": 0}
    for gen in cert.index.intersection(tanaka.index).intersection(clase.index):
        grupo = clase.at[gen]
        if grupo not in mismo:
            continue
        a = float(tanaka.at[gen, "log2fc"])
        b = float(cert.at[gen, "log2fc"])
        if not np.isfinite(a) or not np.isfinite(b) or a == 0 or b == 0:
            continue
        if np.sign(a) == np.sign(b):
            mismo[grupo] += 1
        else:
            distinto[grupo] += 1
    tabla = np.array([
        [mismo["epitelial"], distinto["epitelial"]],
        [mismo["no_epitelial"], distinto["no_epitelial"]],
    ])
    if tabla.sum() == 0 or (tabla.sum(axis=1) == 0).any():
        return {"estado": "no_evaluable", "motivo": "no hay epiteliales y no epiteliales con signo en las dos cohortes", "tabla": tabla.tolist()}
    razon, p = fisher_exact(tabla, alternative="two-sided")
    return {
        "estado": "sostiene" if p < 0.05 else "no_sostiene",
        "motivo": "Fisher bilateral, alfa 0,05; el sentido no se predice",
        "odds_ratio": float(razon),
        "p": float(p),
        "tabla_mismo_signo_luego_distinto": tabla.tolist(),
    }
