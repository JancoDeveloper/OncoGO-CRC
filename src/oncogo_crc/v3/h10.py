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
from oncogo_crc.v3.salida import volcar_json, volcar_tsv

CANCER = "Colorectal cancer"
# "Class" contiene la subcadena "all". El mapeo es por el texto entero, no por subcadena.
CONTROLES_EXPLICITOS = {
    "healthy": "sanos",
    "class": "misma_clase",
    "same class": "misma_clase",
    "all other diseases": "todas_las_demas",
    "all other": "todas_las_demas",
}
CONTROLES = ("sanos", "misma_clase", "todas_las_demas")
PARTIDA = ("SAA1", "LBP", "CHMP4A", "BAIAP3", "FN1", "TPM4")


def tipo_de_control(valor: object) -> str:
    texto = str(valor).strip().lower()
    if texto in CONTROLES_EXPLICITOS:
        return CONTROLES_EXPLICITOS[texto]
    raise ValueError(f"control del HPA sin mapear: {valor!r}")


def clasificar_gen(bloque: pd.DataFrame, alpha: float, lfc_min: float) -> dict:
    propio = bloque.loc[bloque["Disease"].eq(CANCER)]
    filas = {t: propio.loc[propio["tipo_control"].eq(t)] for t in CONTROLES}
    if any(len(f) != 1 for f in filas.values()):
        return {"clase": "no_medida"}
    fila = {t: f.iloc[0] for t, f in filas.items()}

    def signif(f):
        return bool(np.isfinite(f["padj"]) and f["padj"] < alpha and abs(f["logFC"]) >= lfc_min)

    signo = np.sign(fila["sanos"]["logFC"])
    detalle = {f"padj_{t}": float(fila[t]["padj"]) for t in CONTROLES} | {f"logfc_{t}": float(fila[t]["logFC"]) for t in CONTROLES}
    if not signif(fila["sanos"]):
        return {"clase": "no_elevada_frente_a_sanos", **detalle}
    todas = signif(fila["todas_las_demas"]) and np.sign(fila["todas_las_demas"]["logFC"]) == signo
    clase = signif(fila["misma_clase"]) and np.sign(fila["misma_clase"]["logFC"]) == signo
    if todas and clase:
        nombre = "especifica_crc"
    elif todas:
        nombre = "especifica_frente_a_enfermedades"
    elif clase:
        nombre = "especifica_frente_a_canceres"
    else:
        nombre = "no_especifica"
    return {"clase": nombre, **detalle}


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
    tabla["logFC"] = pd.to_numeric(tabla["logFC"], errors="coerce")
    if "Control" not in tabla.columns:
        raise RuntimeError("el atlas del HPA no trae la columna Control")
    tabla["tipo_control"] = tabla["Control"].map(tipo_de_control)
    alpha = float(prereg.get("inferencia.alpha"))
    lfc_min = float(prereg.get("inferencia.treat_lfc"))
    filas = []
    for gen, bloque in tabla.groupby("gen"):
        filas.append({"gen": gen, **clasificar_gen(bloque, alpha, lfc_min)})
    clases = pd.DataFrame(filas)
    orden = {gen: i for i, gen in enumerate(PARTIDA)}
    partida = clases.loc[clases["gen"].isin(PARTIDA)].copy()
    faltan = [gen for gen in PARTIDA if gen not in set(partida["gen"])]
    if faltan:
        partida = pd.concat([partida, pd.DataFrame({"gen": faltan, "clase": "no_medida"})], ignore_index=True)
    partida = partida.sort_values("gen", key=lambda s: s.map(orden), kind="stable")
    volcar_tsv(partida, salida / "h10_proteinas_de_partida.tsv", sello)
    go = settings.resolve("metadata/gene_ontology/go_gene_long.tsv")
    if go.is_file():
        anotacion = pd.read_csv(go, sep="\t", usecols=["gene_symbol", "go_id"], dtype=str)
        genes_aguda = set(anotacion.loc[anotacion["go_id"].eq("GO:0006953"), "gene_symbol"].str.upper())
        aguda = clases.loc[clases["gen"].astype(str).str.upper().isin(genes_aguda)].copy()
    else:
        aguda = pd.DataFrame(columns=clases.columns)
    volcar_tsv(aguda, salida / "h10_fase_aguda.tsv", sello)
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
    columna_enfermedad = next((c for c in catalogo.columns if str(c).lower() in {"disease", "disease name", "name"}), "Disease")
    nombres = catalogo[columna_enfermedad].astype(str) if columna_enfermedad in catalogo.columns else tabla["Disease"].astype(str)
    n_enfermedades = int(nombres.str.lower().ne("healthy").sum())
    payload = {
        "estado": "descriptiva",
        "fase_pepe": 1,
        "motivo_estado": "El preregistro no fijó umbral. Fijarlo después de ver 0/1165 sería elegir por resultado.",
        "valores_de_control": sorted(tabla["Control"].astype(str).unique()),
        "n_enfermedades": n_enfermedades,
        "n_enfermedades_en_el_catalogo": int(len(catalogo)),
        "nota_catalogo": "El catálogo incluye Healthy. n_enfermedades lo deja fuera.",
        "alfa": alpha,
        "treat_lfc": lfc_min,
        "n_especificas_crc": int(conteo.get("especifica_crc", 0)),
        "conteo": {k: int(v) for k, v in conteo.items()},
    }
    volcar_tsv(clases, salida / "h10_especificidad.tsv", sello)
    volcar_json(payload, salida / "h10_especificidad.json", sello)
    return payload
