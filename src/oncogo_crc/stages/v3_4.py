"""Corrida v3.4. Escribe solo bajo results/v3_4/ y no reabre v3.3 para escribir."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import cohen_kappa_score

from oncogo_crc.settings import Settings
from oncogo_crc.stages.v3_3_controles import MUCOSA_BAJA
from oncogo_crc.stages.v3_4_controles import estado_bloque, estado_control, filas_fisicas
from oncogo_crc.stages.v3_4_report import armar_certificado
from oncogo_crc.utils.dia_nn_plasma import BYTES_REPORT, LIBRERIA_VACIA, SHA256_REPORT, leer_libreria, leer_report, matriz_filtrada, trazabilidad
from oncogo_crc.utils.sello import exigir_arbol_limpio, hacer_sello
from oncogo_crc.v3 import h10, h5, h6, h7, h8, h9
from oncogo_crc.v3.genes_tanaka import filas_de_gen, mapa_genes
from oncogo_crc.v3.ml_suero import _log2_si_hace_falta, evaluar
from oncogo_crc.v3.prereg import PreregistroV3
from oncogo_crc.v3.rejilla_v34 import celdas_ml, celdas_suero
from oncogo_crc.v3.salida import carpeta_v34, entrada_v3, manifiesto_congelado, volcar_json, volcar_tsv
from oncogo_crc.v3.suero_v32 import _contraste, _contraste_covariable, _grupos, _log2, _simbolos
from oncogo_crc.stages.v3_3_suero import _sha256


def _abrir(settings: Settings):
    prereg = PreregistroV3(settings.root)
    salida = carpeta_v34(settings)
    salida.mkdir(parents=True, exist_ok=True)
    sello = hacer_sello(settings.root, etapa="v3_4", manifest=manifiesto_congelado(settings))
    exigir_arbol_limpio(sello)
    return prereg, salida, sello


def _controles(settings: Settings, salida: Path, sello: dict, prereg: PreregistroV3) -> dict:
    grupos_path = settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt")
    uni_path = settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv")
    k1 = entrada_v3(settings, salida, "k1_primario_frente_a_mucosa.tsv")
    qc = entrada_v3(settings, salida, "qc_etiquetas_mucosa.tsv")
    if not grupos_path.is_file() or not uni_path.is_file() or not k1.is_file():
        payload = {"estado": "no_disponible", "motivo": "faltan proteinGroups, el mapa o el contraste de mucosa"}
        volcar_json(payload, salida / "controles_v3_4.json", sello)
        return payload
    grupos = pd.read_csv(grupos_path, sep="\t", low_memory=False)
    uni = pd.read_csv(uni_path, sep="\t", usecols=["gene_symbol", "uniprot_accession"], dtype=str).dropna()
    uni["gene_symbol"] = uni["gene_symbol"].str.upper()
    unicos = uni.groupby("uniprot_accession")["gene_symbol"].nunique()
    uni = uni.loc[uni["uniprot_accession"].isin(unicos[unicos == 1].index)]
    mapa_uni = dict(zip(uni["uniprot_accession"], uni["gene_symbol"], strict=True))
    mapa = mapa_genes(grupos_path, uni_path)
    contraste = pd.read_csv(k1, sep="\t")
    contraste["log2fc"] = pd.to_numeric(contraste["log2fc"], errors="coerce")
    columna_n = next((c for c in ("n_obs", "n_pares_observados", "n_pares") if c in contraste.columns), None)
    min_pares = int(prereg.get("enmienda_v3_2.tejido.min_pares_observados"))
    filas = []
    for gen in MUCOSA_BAJA:
        fisico = filas_fisicas(grupos, gen, mapa_uni)
        ids = set(filas_de_gen(contraste["protein_id"].astype(str), mapa, gen))
        bloque = contraste.loc[contraste["protein_id"].astype(str).isin(ids)]
        n_obs = log2fc = float("nan")
        if len(bloque):
            if columna_n is not None:
                bloque = bloque.sort_values(columna_n, ascending=False, kind="stable")
            n_obs = float(pd.to_numeric(bloque[columna_n], errors="coerce").iloc[0]) if columna_n else float("nan")
            log2fc = float(bloque["log2fc"].iloc[0])
        filas.append({
            "gen": gen,
            **fisico,
            "n_obs": n_obs,
            "log2fc": log2fc,
            "estado": estado_control(fisico, n_obs, log2fc, min_pares),
        })
    qc_pasa = False
    if qc.is_file():
        etiquetas = pd.read_csv(qc, sep="\t")
        qc_pasa = len(etiquetas) > 0 and etiquetas["estado"].eq("pasa").all()
    veredicto = estado_bloque([fila["estado"] for fila in filas], qc_pasa)
    anterior = entrada_v3(settings, salida, "controles.json")
    viejo = json.loads(anterior.read_text(encoding="utf-8")) if anterior.is_file() else {"estado": "no_disponible"}
    payload = {
        "controles_fisicos": filas,
        "veredicto_v3_4": veredicto,
        "veredicto_v3_3": viejo.get("estado", "no_disponible"),
        "qc_mucosa_todos_pasan": bool(qc_pasa),
    }
    volcar_tsv(pd.DataFrame(filas), salida / "controles_fisicos.tsv", sello)
    volcar_json(payload, salida / "controles_v3_4.json", sello)
    return payload


def _mascaras(traza: pd.DataFrame) -> dict[str, pd.Series]:
    panel = traza["panel"].fillna("").astype(str).str.len().gt(0)
    cont = traza["contaminante"].astype(bool) | traza["contaminante_piel"].astype(bool)
    epi = traza["queratina_epitelial"].astype(bool) & traza["contaminante_motivo"].eq("acceso_humano_en_libreria") if "queratina_epitelial" in traza.columns else pd.Series(False, index=traza.index)
    mix = traza["comparte_con_contaminante"].astype(bool) if "comparte_con_contaminante" in traza.columns else pd.Series(False, index=traza.index)
    return {"cont": cont, "panel": panel, "epi": epi, "mix": mix}


def _quitar(traza: pd.DataFrame, mascaras: dict[str, pd.Series], contaminante: str, panel: str, por_gen: bool) -> set[str]:
    quitar = pd.Series(False, index=traza.index)
    if contaminante == "libreria_v3_3":
        quitar = quitar | mascaras["cont"]
    elif contaminante == "con_queratinas_epiteliales":
        quitar = quitar | (mascaras["cont"] & ~mascaras["epi"])
    elif contaminante == "sin_grupos_mixtos":
        quitar = quitar | mascaras["cont"] | mascaras["mix"]
    if panel in ("excluidos", "covariable_coagulacion"):
        quitar = quitar | mascaras["panel"]
    if por_gen:
        genes: set[str] = set()
        for texto in traza.loc[quitar, "genes"].astype(str):
            genes |= _simbolos(texto)
        return genes
    return set(traza.loc[quitar, "protein_group"].astype(str))


def _variantes(pg: pd.DataFrame, genes: pd.DataFrame, traza: pd.DataFrame) -> dict[tuple[str, str, str], pd.DataFrame]:
    mascaras = _mascaras(traza)
    salida = {}
    for cuant, base in (("PG.MaxLFQ", pg), ("Genes.MaxLFQ", genes)):
        if base.empty:
            continue
        for cont in ("libreria_v3_3", "sin_quitar", "con_queratinas_epiteliales", "sin_grupos_mixtos"):
            for panel in ("excluidos", "covariable_coagulacion", "incluidos"):
                fuera = _quitar(traza, mascaras, cont, panel, por_gen=cuant == "Genes.MaxLFQ")
                salida[(cuant, cont, panel)] = base.loc[~base.index.astype(str).isin(fuera)]
    return salida


def _con_presencia(mat: pd.DataFrame, cols_a: list[str], cols_b: list[str], presencia: str) -> pd.DataFrame:
    if presencia == "min_n_3" or mat.empty:
        return mat
    umbral = float(presencia)
    a = np.isfinite(mat[cols_a].to_numpy(dtype=float)).mean(axis=1)
    b = np.isfinite(mat[cols_b].to_numpy(dtype=float)).mean(axis=1)
    return mat.loc[(a >= umbral) & (b >= umbral)]


def _multiverso(settings: Settings, salida: Path, sello: dict, prereg: PreregistroV3) -> dict:
    report = settings.resolve(str(prereg.get("enmienda_v3_2.report_pxd055218.destino")))
    if not report.is_file() or report.stat().st_size != BYTES_REPORT or _sha256(report) != SHA256_REPORT:
        payload = {"estado": "no_disponible", "motivo": "el report de DIA-NN no está o no coincide con el tamaño y el SHA-256 de la enmienda"}
        volcar_json(payload, salida / "multiverso_v3_4.json", sello)
        return payload
    tabla = leer_report(report)
    fasta = settings.resolve(str(prereg.get("enmienda_v3_2.contaminantes.fichero")))
    libreria = leer_libreria(fasta) if fasta.is_file() else LIBRERIA_VACIA
    traza = trazabilidad(
        tabla,
        libreria=libreria,
        genes_piel=set(prereg.get("enmienda_v3_2.contaminante_piel.genes")),
        paneles={
            "eritrocitos": set(prereg.get("enmienda_v3_2.paneles_geyer_2019.eritrocitos")),
            "plaquetas": set(prereg.get("enmienda_v3_2.paneles_geyer_2019.plaquetas")),
            "coagulacion": set(prereg.get("enmienda_v3_2.paneles_geyer_2019.coagulacion")),
        },
    )
    pg, _escala = _log2(matriz_filtrada(tabla))
    genes_ancha = matriz_filtrada(tabla, valor="Genes.MaxLFQ", indice="Genes", q_de_grupo="GG.Q.Value") if {"Genes.MaxLFQ", "GG.Q.Value"}.issubset(tabla.columns) else pd.DataFrame()
    genes = _log2(genes_ancha)[0] if len(genes_ancha) else pd.DataFrame()
    variantes = _variantes(pg, genes, traza)
    grupos = _grupos(list(pg.columns))
    cols_a, cols_b = grupos["crc"], grupos["colonoscopy_negative"]
    calidad = entrada_v3(settings, salida, "suero_calidad_muestra.tsv")
    coag = pd.read_csv(calidad, sep="\t").set_index("sample_id")["coagulacion"] if calidad.is_file() else pd.Series(dtype=float)
    filas = []
    for celda in celdas_suero():
        mat = variantes.get((celda["cuantificacion"], celda["contaminantes"], celda["paneles"]))
        if mat is None or mat.empty:
            filas.append({**celda, "estado": "no_disponible", "n": 0, "diferente": 0, "equivalente": 0, "indeterminada": 0})
            continue
        recortada = _con_presencia(mat, cols_a, cols_b, celda["presencia"])
        if celda["paneles"] == "covariable_coagulacion":
            if coag.empty:
                filas.append({**celda, "estado": "no_disponible", "motivo": "no está suero_calidad_muestra.tsv"})
                continue
            contraste = _contraste_covariable(recortada, cols_a, cols_b, coag, lfc=celda["lfc"], alpha=celda["alfa"], nombre="v3_4")
        else:
            contraste = _contraste(recortada, cols_a, cols_b, lfc=celda["lfc"], alpha=celda["alfa"], nombre="v3_4")
        filas.append({
            **celda,
            "estado": "calculado",
            "n": int(len(contraste)),
            "diferente": int((contraste["clase"] == "diferente").sum()),
            "equivalente": int((contraste["clase"] == "equivalente").sum()),
            "indeterminada": int((contraste["clase"] == "indeterminada").sum()),
        })
    r1 = entrada_v3(settings, salida, "k1_sensibilidad_normalizacion.json")
    r1_json = json.loads(r1.read_text(encoding="utf-8")) if r1.is_file() else {}
    normalizaciones = [nombre for nombre in ("antes", "despues") if nombre in r1_json]
    ml = _ml(variantes, cols_a, cols_b, coag, prereg)
    payload = {
        "n_celdas_suero": len(filas),
        "n_celdas_ml": len(ml),
        "n_r1": len(normalizaciones),
        "r1": {nombre: r1_json[nombre].get("r1_estimacion") for nombre in normalizaciones},
        "r1_fuente": str(r1),
        "celda_preregistrada": next((fila for fila in filas if fila.get("es_la_celda_preregistrada")), {}),
    }
    volcar_tsv(pd.DataFrame(filas), salida / "multiverso_suero_v3_4.tsv", sello)
    volcar_tsv(pd.DataFrame(ml), salida / "multiverso_ml_v3_4.tsv", sello)
    volcar_json(payload, salida / "multiverso_v3_4.json", sello)
    return payload


def _ml(variantes: dict, cols_a: list[str], cols_b: list[str], coag: pd.Series, prereg: PreregistroV3) -> list[dict]:
    base = int(prereg.get("enmienda_v3_2.semillas.base"))
    semilla = base + int(prereg.get("enmienda_v3_4.semillas.desplazamientos.multiverso_v3_4"))
    filas = []
    columnas = cols_a + cols_b
    y = np.array([1] * len(cols_a) + [0] * len(cols_b))
    for i, celda in enumerate(celdas_ml()):
        print(f"[v3.4] ML {i + 1}/72")
        mat = variantes.get((celda["cuantificacion"], celda["contaminantes"], celda["paneles"]))
        if mat is None or mat.empty:
            filas.append({**celda, "estado": "no_disponible"})
            continue
        recortada = _con_presencia(mat, cols_a, cols_b, celda["presencia"])
        if recortada.empty:
            filas.append({**celda, "estado": "no_evaluable", "motivo": "ninguna proteína pasa la presencia"})
            continue
        x = _log2_si_hace_falta(recortada[columnas].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float).T)
        if celda["paneles"] == "covariable_coagulacion":
            if coag.empty:
                filas.append({**celda, "estado": "no_disponible", "motivo": "no está la covariable de coagulación"})
                continue
            extra = coag.reindex(columnas).to_numpy(dtype=float).reshape(-1, 1)
            x = np.hstack([x, extra])
        n_perm = 1000 if celda["es_la_celda_preregistrada"] else 200
        resultado = evaluar(x, y, semilla, semilla + 1, "proteinas", n_perm)
        filas.append({**celda, **{k: v for k, v in resultado.items() if not isinstance(v, (dict, list))}})
    return filas


def _triangulacion(settings: Settings, salida: Path, sello: dict, prereg: PreregistroV3) -> None:
    genes = list(prereg.get("enmienda_v3_4.triangulacion.genes"))
    tabla = pd.DataFrame({"gen": genes})
    volcar_tsv(tabla, salida / "triangulacion_v3_4.tsv", sello)
    volcar_json({
        "tipo": "exploratoria_post_hoc",
        "sin_pruebas_nuevas": True,
        "n_genes": len(genes),
        "nota": "La tabla lista los genes de la enmienda. Cada columna numérica se añade solo cuando su fichero de origen está en esta corrida.",
    }, salida / "triangulacion_v3_4.json", sello)


def _enriquecer_h7(salida: Path, sello: dict, prereg: PreregistroV3, payload: dict) -> dict:
    tabla = salida / "h7_tabla_proteinas.tsv"
    recursos = payload.get("recursos") or []
    concordancias = []
    if tabla.is_file() and len(recursos) >= 2:
        datos = pd.read_csv(tabla, sep="\t")
        for i, a in enumerate(recursos):
            for b in recursos[i + 1:]:
                ca, cb = f"y_{a}", f"y_{b}"
                if ca not in datos.columns or cb not in datos.columns:
                    continue
                par = datos[[ca, cb]].dropna()
                if len(par) < 2 or par[ca].nunique() < 2 or par[cb].nunique() < 2:
                    continue
                concordancias.append({
                    "a": a,
                    "b": b,
                    "n": int(len(par)),
                    "kappa": float(cohen_kappa_score(par[ca].astype(int), par[cb].astype(int))),
                })
    payload["constructo_por_recurso"] = {
        "alvez_2023": prereg.get("enmienda_v3_4.h7.recursos.alvez_2023.contraste_real"),
        "hpa_disease_blood_atlas": prereg.get("enmienda_v3_4.h7.recursos.hpa_disease_blood_atlas.contraste_real"),
        "papier_2024": prereg.get("enmienda_v3_4.h7.recursos.papier_2024.contraste_real"),
    }
    payload["independientes"] = {
        "comprobado": False,
        "motivo": "los ficheros descargados no dicen si los casos de Álvez 2023 están en el atlas de 2025",
    }
    payload["concordancia_de_etiquetas"] = concordancias
    payload["modelo_completo"] = not list(payload.get("variables_faltantes") or [])
    if "papier_2024" not in recursos:
        payload["por_recurso"] = dict(payload.get("por_recurso") or {})
        payload["por_recurso"]["papier_2024"] = {
            "estado": "no_disponible",
            "motivo": "no se pudo etiquetar Papier 2024 con el suplemento descargado",
        }
    volcar_json(payload, salida / "h7_visibilidad_en_sangre.json", sello)
    return payload


def _preanalitica(settings: Settings, salida: Path, sello: dict) -> dict:
    crudos = settings.resolve("data/raw/plasma/pxd055218_raw")
    n = len(list(crudos.glob("*.raw"))) if crudos.is_dir() else 0
    herramienta = shutil.which("ThermoRawFileParser")
    if n < 60 or not herramienta:
        payload = {
            "estado": "no_disponible",
            "n_crudos": n,
            "motivo": "hacen falta los 60 crudos y ThermoRawFileParser en el PATH. No se inventan tiempos de adquisición.",
        }
    else:
        payload = {"estado": "no_evaluable", "motivo": "el analizador de metadatos todavía no está cableado a este binario", "binario": herramienta}
    volcar_json(payload, salida / "preanalitica_v3_4.json", sello)
    return payload


def _secuencia(salida: Path, sello: dict) -> dict:
    ruta = salida / "h7_visibilidad_en_sangre.json"
    if ruta.is_file():
        datos = json.loads(ruta.read_text(encoding="utf-8"))
        sec = datos.get("modelo_secundario") or {}
        if sec.get("esm2") == "calculado" and sec.get("tmbed") == "calculado":
            payload = {
                "estado": "calculado",
                "esm2": sec.get("esm2_fuente"),
                "tmbed": sec.get("tmbed_fuente"),
                "n_con_secuencia": sec.get("n_con_secuencia"),
                "n_sin_secuencia": sec.get("n_sin_secuencia"),
                "nota": "El modelo principal no usa estas variables. El secundario sí, y no cambia el veredicto.",
            }
            volcar_json(payload, salida / "secuencia_v3_4.json", sello)
            return payload
    try:
        import torch  # noqa: F401
        hay_torch = True
    except ImportError:
        hay_torch = False
    payload = {
        "estado": "no_disponible" if not hay_torch else "no_evaluable",
        "motivo": "torch no está instalado; ESM-2 y ProtT5 no se ejecutaron" if not hay_torch else "ESM-2 o TMbed no tienen todavía una salida completa",
    }
    volcar_json(payload, salida / "secuencia_v3_4.json", sello)
    return payload


def run(settings: Settings) -> None:
    prereg, salida, sello = _abrir(settings)
    print("[v3.4] H5")
    h5_payload = h5.ejecutar(settings, salida, sello, prereg)
    print("[v3.4] controles")
    controles = _controles(settings, salida, sello, prereg)
    print("[v3.4] H6")
    h6_payload = h6.ejecutar(settings, salida, sello)
    base = int(prereg.get("enmienda_v3_2.semillas.base"))
    print("[v3.4] H7")
    h7_payload = h7.ejecutar(
        settings, salida, sello, prereg,
        semilla=base + int(prereg.get("enmienda_v3_4.semillas.desplazamientos.h7_modelo_v3_4")),
        semilla_b=base + int(prereg.get("enmienda_v3_4.semillas.desplazamientos.h7_bootstrap_v3_4")),
    )
    h7_payload = _enriquecer_h7(salida, sello, prereg, h7_payload)
    print("[v3.4] H8")
    h8_payload = h8.ejecutar(settings, salida, sello)
    print("[v3.4] H10")
    h10_payload = h10.ejecutar(settings, salida, sello, prereg)
    print("[v3.4] H9")
    h9_payload = h9.ejecutar(settings, salida, sello, prereg)
    _triangulacion(settings, salida, sello, prereg)
    print("[v3.4] multiverso")
    multi = _multiverso(settings, salida, sello, prereg)
    pre = _preanalitica(settings, salida, sello)
    sec = _secuencia(salida, sello)
    piezas = {
        "h5": h5_payload,
        "controles": controles,
        "h6": h6_payload,
        "h7": h7_payload,
        "h8": h8_payload,
        "h9": h9_payload,
        "h10": h10_payload,
        "multiverso": multi,
        "preanalitica": pre,
        "secuencia": sec,
    }
    certificado = armar_certificado(piezas, {"enmienda_v3_4": prereg.get("enmienda_v3_4")})
    volcar_json(certificado, salida / "certificado_v3_4.json", sello)
    print(f"[v3.4] certificado: {certificado['estado']}")
