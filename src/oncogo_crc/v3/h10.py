"""H10. Especificidad clínica en el atlas de sangre del HPA.

Una proteína es específica si el colorrectal frente a sanos pasa el umbral y
ninguna otra enfermedad medida lo pasa. Si no está medida en colorrectal, o
no está medida en ningún otro cáncer y en ninguna enfermedad no cancerosa,
queda no_medida. El umbral es el alfa del preregistro.
"""

from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.v3.prereg import PreregistroV3
from oncogo_crc.v3.salida_v32 import volcar_json, volcar_tsv

CANCER = "Colorectal cancer"


def clasificar_gen(bloque: pd.DataFrame, alpha: float) -> str:
    propio = bloque.loc[bloque["Disease"].eq(CANCER)]
    if propio.empty or not np.isfinite(propio["padj"].iloc[0]):
        return "no_medida"
    otros = bloque.loc[~bloque["Disease"].eq(CANCER)]
    hay_otro_cancer = bool(otros["Class"].eq("Cancer").any())
    hay_no_cancer = bool((~otros["Class"].eq("Cancer")).any())
    if otros.empty or not (hay_otro_cancer and hay_no_cancer):
        return "no_medida"
    if not bool(propio["padj"].iloc[0] < alpha):
        return "no_especifica"
    if bool((otros["padj"] < alpha).any()):
        return "no_especifica"
    return "especifica"


def ejecutar(settings: Settings, salida: Path, sello: dict, prereg: PreregistroV3) -> dict:
    base = settings.resolve("external/hpa")
    de = base / "blood_pea_disease_de.tsv.zip"
    enfermedades = base / "blood_pea_diseases.tsv.zip"
    if not de.is_file() or not enfermedades.is_file():
        payload = {"estado": "no_disponible", "motivo": f"faltan {de.name} o {enfermedades.name}"}
        volcar_json(payload, salida / "h10_especificidad.json", sello)
        return payload
    with zipfile.ZipFile(de) as zf:
        nombre = zf.namelist()[0]
        tabla = pd.read_csv(zf.open(nombre), sep="\t")
    with zipfile.ZipFile(enfermedades) as zf:
        catalogo = pd.read_csv(zf.open(zf.namelist()[0]), sep="\t")
    tabla = tabla.rename(columns={"p-value adjusted": "padj", "Gene": "gen"})
    tabla["padj"] = pd.to_numeric(tabla["padj"], errors="coerce")
    alpha = float(prereg.get("inferencia.alpha"))
    filas = []
    for gen, bloque in tabla.groupby("gen"):
        filas.append({"gen": gen, "clase": clasificar_gen(bloque, alpha)})
    clases = pd.DataFrame(filas)
    piel = set(prereg.get("enmienda_v3_2.contaminante_piel.genes"))
    paneles = set()
    for nombre_panel in ("eritrocitos", "plaquetas", "coagulacion"):
        paneles.update(prereg.get(f"enmienda_v3_2.paneles_geyer_2019.{nombre_panel}"))
    candidatas = settings.resolve("metadata/tesis/candidatas_literatura.tsv")
    if candidatas.is_file():
        lista = pd.read_csv(candidatas, sep="\t")
        auditoria = lista.merge(clases, left_on="gene", right_on="gen", how="left")
        auditoria["clase"] = auditoria["clase"].fillna("no_medida")
        auditoria["contaminante_piel"] = auditoria["gene"].isin(piel)
        auditoria["panel_geyer"] = auditoria["gene"].isin(paneles)
        volcar_tsv(auditoria, salida / "h10_candidatas.tsv", sello)
    conteo = clases["clase"].value_counts().to_dict()
    payload = {
        "estado": "calculado",
        "n_enfermedades_en_el_catalogo": int(len(catalogo)),
        "nota_catalogo": "El prompt contaba 59 enfermedades. El fichero descargado trae las filas de este catalogo.",
        "alfa": alpha,
        "n_especificas": int(conteo.get("especifica", 0)),
        "n_no_especificas": int(conteo.get("no_especifica", 0)),
        "n_no_medidas": int(conteo.get("no_medida", 0)),
        "fase_pepe": 1,
    }
    volcar_tsv(clases, salida / "h10_especificidad.tsv", sello)
    volcar_json(payload, salida / "h10_especificidad.json", sello)
    return payload
