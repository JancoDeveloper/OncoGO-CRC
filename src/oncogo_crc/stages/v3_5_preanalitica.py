"""Orden de adquisición de PXD055218. Con menos de 60 crudos el estado es no_evaluable."""

from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from oncogo_crc.settings import Settings
from oncogo_crc.v3.ml_suero import ic_ledell
from oncogo_crc.v3.salida import volcar_json, volcar_tsv

_CRUDO = re.compile(r"(CP|CN|HC)_\d{2}_DIA\.raw")
_PRIDE_BYTES = {
    "CN_01_DIA.raw": 1049764670,
    "CN_08_DIA.raw": 1000697970,
    "CN_16_DIA.raw": 1205043346,
    "CP_18_DIA.raw": 1226017242,
    "HC_17_DIA.raw": 1065175806,
    "HC_20_DIA.raw": 1088171988,
}


def contar_crudos(nombres: list[str]) -> tuple[list[str], int]:
    validos = sorted(n for n in nombres if _CRUDO.fullmatch(n))
    apple = sum(1 for n in nombres if n.startswith("._"))
    return validos, apple


def fecha_de_metadato(meta: dict, campo: str) -> str | None:
    for bloque in meta.get("FileProperties", []):
        if bloque.get("name") == campo and bloque.get("value"):
            momento = datetime.strptime(str(bloque["value"]), "%m/%d/%Y %H:%M:%S")
            return momento.strftime("%Y-%m-%dT%H:%M:%S")
    return None


def _valor(meta: dict, fragmento: str) -> str | None:
    for bloque in meta.get("InstrumentProperties", []):
        nombre = str(bloque.get("name", ""))
        if fragmento.lower() in nombre.lower() and bloque.get("value"):
            return str(bloque["value"])
    return None


def resumir_legibles(*, n_validos: int, n_con_fecha: int, auc, spearman: dict) -> dict:
    return {
        "estado": "no_evaluable" if n_con_fecha < 60 else "calculado",
        "motivo": f"{n_con_fecha} de {n_validos} crudos legibles; la prueba preregistrada pide 60" if n_con_fecha < 60 else "60 crudos con fecha",
        "calculado_con_legibles": {"auc": auc, "n_con_fecha": n_con_fecha, "spearman": spearman},
    }


def ejecutar(settings: Settings, salida: Path, sello: dict, prereg) -> dict:
    campo = str(prereg.get("enmienda_v3_5.preanalitica.campo_fecha"))
    carpeta = settings.resolve("data/raw/plasma/pxd055218_raw")
    metadatos = settings.resolve("external/pxd055218_metadata")
    nombres = [p.name for p in carpeta.iterdir()] if carpeta.is_dir() else []
    validos, apple = contar_crudos(nombres)
    filas = []
    sin_fecha = []
    for nombre in validos:
        crudo = carpeta / nombre
        json_path = metadatos / f"{crudo.stem}.json"
        tamaño = crudo.stat().st_size if crudo.is_file() else None
        if not json_path.is_file():
            pride = _PRIDE_BYTES.get(nombre)
            sin_fecha.append({
                "crudo": nombre,
                "bytes_local": tamaño,
                "bytes_pride": pride,
                "tamaño_igual_a_pride": tamaño == pride if pride else None,
                "motivo": "ThermoRawFileParser no escribió JSON. El tamaño coincide con PRIDE, así que la descarga no quedó corta." if tamaño == pride else "sin JSON",
            })
            continue
        meta = json.loads(json_path.read_text(encoding="utf-8"))
        filas.append({
            "crudo": nombre,
            "grupo": nombre[:2],
            "fecha": fecha_de_metadato(meta, campo),
            "instrumento": _valor(meta, "instrument model"),
            "serie": _valor(meta, "serial number"),
        })
    tabla = pd.DataFrame(filas)
    con_fecha = tabla.dropna(subset=["fecha"]).sort_values(["fecha", "crudo"]) if len(tabla) else tabla
    if len(con_fecha):
        con_fecha = con_fecha.copy()
        con_fecha["orden"] = range(1, len(con_fecha) + 1)
    else:
        con_fecha["orden"] = []
    volcar_tsv(con_fecha, salida / "preanalitica_orden.tsv", sello)
    auc = None
    ic = None
    n_cp = n_cn = 0
    if len(con_fecha):
        pares = con_fecha.loc[con_fecha["grupo"].isin(["CP", "CN"])]
        n_cp = int((pares["grupo"] == "CP").sum())
        n_cn = int((pares["grupo"] == "CN").sum())
        if n_cp >= 3 and n_cn >= 3:
            y = (pares["grupo"] == "CP").to_numpy().astype(int)
            auc, bajo, alto = ic_ledell(y, pares["orden"].to_numpy(dtype=float))
            ic = [bajo, alto]
    spearman = {}
    calidad = settings.resolve("results/v3_3/suero_calidad_muestra.tsv")
    if calidad.is_file() and len(con_fecha):
        qc = pd.read_csv(calidad, sep="\t")
        qc["crudo"] = qc["sample_id"].astype(str) + ".raw"
        unido = con_fecha.merge(qc, on="crudo", how="inner")
        if "n_precursores" in unido.columns:
            unido["log2_n_precursores"] = np.log2(pd.to_numeric(unido["n_precursores"], errors="coerce"))
        for columna in ("eritrocitos", "plaquetas", "coagulacion", "piel", "log2_n_precursores"):
            if columna not in unido.columns:
                spearman[columna] = {"estado": "no_disponible"}
                continue
            x = pd.to_numeric(unido["orden"], errors="coerce")
            y = pd.to_numeric(unido[columna], errors="coerce")
            ok = x.notna() & y.notna()
            if int(ok.sum()) < 5:
                spearman[columna] = {"estado": "no_evaluable", "n": int(ok.sum())}
                continue
            rho, p = spearmanr(x[ok], y[ok])
            spearman[columna] = {"rho": float(rho), "p": float(p), "n": int(ok.sum())}
    resumen = resumir_legibles(
        n_validos=len(validos),
        n_con_fecha=int(len(con_fecha)),
        auc={"valor": auc, "ic95_ledell": ic, "n_cp": n_cp, "n_cn": n_cn},
        spearman=spearman,
    )
    fechas = con_fecha["fecha"].tolist() if len(con_fecha) else []
    por_dia: dict[str, int] = {}
    for fecha in fechas:
        por_dia[fecha[:10]] = por_dia.get(fecha[:10], 0) + 1
    hc = con_fecha.loc[con_fecha["grupo"].eq("HC"), "orden"] if len(con_fecha) else pd.Series(dtype=float)
    otros = con_fecha.loc[con_fecha["grupo"].isin(["CP", "CN"]), "orden"] if len(con_fecha) else pd.Series(dtype=float)
    if len(hc) and len(otros):
        dentro = int(((hc >= otros.min()) & (hc <= otros.max())).sum())
        posicion_hc = "entre_CP_y_CN" if dentro else "aparte"
    else:
        dentro = 0
        posicion_hc = "no_evaluable"
    payload = {
        **resumen,
        "n_crudos_validos": len(validos),
        "n_appledouble": apple,
        "n_con_fecha": int(len(con_fecha)),
        "sin_fecha": sin_fecha,
        "campo_fecha": campo,
        "formato_fecha": "%m/%d/%Y %H:%M:%S",
        "inicio": fechas[0] if fechas else None,
        "fin": fechas[-1] if fechas else None,
        "corridas_por_dia": por_dia,
        "corridas_por_grupo": con_fecha["grupo"].value_counts().to_dict() if len(con_fecha) else {},
        "hc_dentro_del_rango_cp_cn": dentro,
        "posicion_hc": posicion_hc,
        "instrumento": con_fecha["instrumento"].dropna().mode().iloc[0] if len(con_fecha) and con_fecha["instrumento"].notna().any() else None,
        "no_cambia_el_titular": True,
    }
    volcar_json(payload, salida / "preanalitica_v3_5.json", sello)
    return payload
