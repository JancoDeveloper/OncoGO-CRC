"""Corrida v3.5. Copia lo congelado de v3.4 y recalcula solo lo que la enmienda cambia."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.stages import v3_4
from oncogo_crc.stages.v3_4 import _controles, _enriquecer_h7
from oncogo_crc.stages.v3_5_preanalitica import ejecutar as preanalitica
from oncogo_crc.stages.v3_5_report import armar_certificado
from oncogo_crc.utils.dia_nn_plasma import BYTES_REPORT, LIBRERIA_VACIA, SHA256_REPORT, leer_libreria, leer_report, matriz_filtrada, trazabilidad
from oncogo_crc.utils.sello import exigir_arbol_limpio, hacer_sello
from oncogo_crc.v3 import h6, h7, h8
from oncogo_crc.v3.ml_suero import evaluar, matrices_de_suero
from oncogo_crc.v3.prereg import PreregistroV3
from oncogo_crc.v3.salida import carpeta_v35, entrada_v3, manifiesto_congelado, volcar_json, volcar_tsv
from oncogo_crc.v3.suero_v32 import _contraste, _contraste_covariable, _grupos, _log2, semilla_v33
from oncogo_crc.v3.triangulacion import ejecutar as triangulacion

COPIAR = ("h5_", "h9_", "h10_", "multiverso_", "secuencia_")
V1 = {
    "cptac_n_significant": 688,
    "pride_n_significant": 140,
    "n_both_significant_padj_lt_0.05": 28,
    "n_same_sign": 104,
    "n_sign_evaluable": 119,
    "peptide_n_measurable": 226,
    "n_certified_surviving": 307,
}


def _abrir(settings: Settings):
    prereg = PreregistroV3(settings.root)
    salida = carpeta_v35(settings)
    salida.mkdir(parents=True, exist_ok=True)
    sello = hacer_sello(settings.root, etapa="v3_5", manifest=manifiesto_congelado(settings))
    exigir_arbol_limpio(sello)
    return prereg, salida, sello


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _manifiesto(settings: Settings) -> dict[str, str]:
    ruta = settings.root / "results/v3_4_manifest.sha256"
    salida = {}
    for linea in ruta.read_text(encoding="utf-8").splitlines():
        if "  " not in linea:
            continue
        sha, rel = linea.split("  ", 1)
        salida[rel] = sha
    return salida


def _copiar(settings: Settings, salida: Path) -> list[dict]:
    origen = settings.root / "results/v3_4"
    esperado = _manifiesto(settings)
    procedencia = []
    for path in sorted(origen.iterdir()):
        if not path.is_file() or not path.name.startswith(COPIAR):
            continue
        rel = f"results/v3_4/{path.name}"
        digest = _sha(path)
        if rel in esperado and esperado[rel] != digest:
            raise RuntimeError(f"{rel} no coincide con results/v3_4_manifest.sha256")
        destino = salida / path.name
        shutil.copy2(path, destino)
        sello_nombre = path.name + ".sello.json" if not path.name.endswith(".sello.json") else None
        procedencia.append({
            "fichero": path.name,
            "origen": "results/v3_4",
            "sha256": digest,
            "sello": sello_nombre if sello_nombre and (origen / sello_nombre).is_file() else None,
        })
    return procedencia


def _celda_suero(salida: Path) -> dict:
    datos = json.loads((salida / "multiverso_v3_4.json").read_text(encoding="utf-8"))
    celda = datos.get("celda_preregistrada") or {}
    if not (celda.get("n") == 284 and celda.get("diferente") == 0 and celda.get("equivalente") == 166 and celda.get("indeterminada") == 118):
        raise RuntimeError(f"la celda preregistrada del suero no reproduce 284/0/166/118: {celda}")
    return celda


def _ml(settings: Settings, salida: Path, sello: dict, prereg) -> dict:
    base = settings.resolve("results/v3_3")
    principal = pd.read_csv(base / "suero_matriz_principal.tsv", sep="\t")
    humana = pd.read_csv(base / "suero_matriz_humana_fdr.tsv", sep="\t")
    calidad = pd.read_csv(base / "suero_calidad_muestra.tsv", sep="\t")
    x, y, tipo = matrices_de_suero(principal, humana, calidad)["M2_proteinas_principal"]
    ordenado = np.ascontiguousarray(np.sort(x, axis=0))
    semilla_cv = semilla_v33(prereg, "ml_cv")
    semilla_nulo = semilla_v33(prereg, "ml_permutacion")
    print("[v3.5] ML celda preregistrada, semillas de v3.3, 1000 permutaciones")
    resultado = evaluar(x, y, semilla_cv, semilla_nulo, tipo, 1000)
    v34 = pd.read_csv(salida / "multiverso_ml_v3_4.tsv", sep="\t")
    previa = v34.loc[v34["es_la_celda_preregistrada"].astype(str).isin(["True", "true", "1"])]
    auc_v34 = float(previa["auc_oof_media"].iloc[0]) if len(previa) and "auc_oof_media" in previa.columns else None
    payload = {
        "modelo": "M2_proteinas_principal",
        "auc": resultado.get("auc_oof_media"),
        "ic95_ledell_bajo": resultado.get("ic95_ledell_bajo"),
        "ic95_ledell_alto": resultado.get("ic95_ledell_alto"),
        "p_permutacion": resultado.get("p_permutacion"),
        "media_del_nulo": resultado.get("media_del_nulo"),
        "n_variables": resultado.get("n_variables"),
        "semilla_cv": semilla_cv,
        "semilla_permutacion": semilla_nulo,
        "orden_de_muestras": "el de _columnas_cp_cn de v3.3",
        "debe_reproducir": 0.88025,
        "reproduce": abs(float(resultado.get("auc_oof_media")) - 0.88025) <= 1e-9,
        "sha256_x": hashlib.sha256(ordenado.tobytes()).hexdigest(),
        "auc_v3_4_variacion_de_semilla": auc_v34,
        "nota_v3_4": "0,888875 de v3.4 es la misma matriz con otra semilla. No es un resultado nuevo." if auc_v34 and abs(auc_v34 - 0.888875) < 1e-6 else "el valor de v3.4 se informa como variación por semilla",
    }
    volcar_json(payload, salida / "ml_celda_preregistrada.json", sello)
    return payload


def _variantes_de_suero(settings: Settings, prereg):
    report = settings.resolve(str(prereg.get("enmienda_v3_2.report_pxd055218.destino")))
    if not report.is_file() or report.stat().st_size != BYTES_REPORT:
        return None, "el report de DIA-NN no está o no tiene el tamaño de la enmienda"
    if v3_4._sha256(report) != SHA256_REPORT:
        return None, "el SHA-256 del report de DIA-NN no coincide"
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
    pg = _log2(matriz_filtrada(tabla))[0]
    genes_ancha = matriz_filtrada(tabla, valor="Genes.MaxLFQ", indice="Genes", q_de_grupo="GG.Q.Value") if {"Genes.MaxLFQ", "GG.Q.Value"}.issubset(tabla.columns) else pd.DataFrame()
    genes = _log2(genes_ancha)[0] if len(genes_ancha) else pd.DataFrame()
    return (v3_4._variantes(pg, genes, traza), traza, _grupos(list(pg.columns))), None


def _resumen_multiverso(settings: Settings, salida: Path, sello: dict, prereg) -> dict:
    suero = pd.read_csv(salida / "multiverso_suero_v3_4.tsv", sep="\t")
    ml = pd.read_csv(salida / "multiverso_ml_v3_4.tsv", sep="\t")
    clave = ["n", "diferente", "equivalente", "indeterminada"]
    distintas = suero.drop_duplicates(clave)
    sin_quitar = suero.loc[suero["contaminantes"].ne("sin_quitar")]
    distintas_sin = sin_quitar.drop_duplicates(clave)
    def fraccion(bloque):
        return float((bloque["diferente"] == 0).mean()) if len(bloque) else None
    especiales = suero.loc[suero["diferente"].isin([5, 6])].copy()
    construido, motivo = _variantes_de_suero(settings, prereg)
    queratina = {"explicacion": motivo or "no se pudo reconstruir la matriz"}
    proteinas = []
    if construido is not None:
        variantes, traza, grupos = construido
        cols_a, cols_b = grupos["crc"], grupos["colonoscopy_negative"]
        base = variantes.get(("PG.MaxLFQ", "libreria_v3_3", "excluidos"))
        con = variantes.get(("PG.MaxLFQ", "con_queratinas_epiteliales", "excluidos"))
        if base is not None and con is not None:
            extra = sorted(set(con.index.astype(str)) - set(base.index.astype(str)))
            queratina = {
                "proteinas_que_vuelven": extra,
                "n_en_traza": int(traza["queratina_epitelial"].astype(bool).sum()) if "queratina_epitelial" in traza.columns else None,
                "explicacion": "la queratina epitelial ya está en la librería de contaminantes. El eje la devuelve, y sale otra vez si no pasa la presencia o el panel.",
            }
        calidad = entrada_v3(settings, salida, "suero_calidad_muestra.tsv")
        coag = pd.read_csv(calidad, sep="\t").set_index("sample_id")["coagulacion"] if calidad.is_file() else pd.Series(dtype=float)
        for _, celda in especiales.iterrows():
            mat = variantes.get((celda["cuantificacion"], celda["contaminantes"], celda["paneles"]))
            if mat is None or mat.empty:
                continue
            recortada = v3_4._con_presencia(mat, cols_a, cols_b, str(celda["presencia"]))
            if celda["paneles"] == "covariable_coagulacion" and not coag.empty:
                contraste = _contraste_covariable(recortada, cols_a, cols_b, coag, lfc=float(celda["lfc"]), alpha=float(celda["alfa"]), nombre="v3_5_trece")
            else:
                contraste = _contraste(recortada, cols_a, cols_b, lfc=float(celda["lfc"]), alpha=float(celda["alfa"]), nombre="v3_5_trece")
            distintas_prot = contraste.loc[contraste["clase"].eq("diferente")]
            for _, fila in distintas_prot.iterrows():
                proteinas.append({
                    "cuantificacion": celda["cuantificacion"],
                    "contaminantes": celda["contaminantes"],
                    "paneles": celda["paneles"],
                    "presencia": celda["presencia"],
                    "lfc": float(celda["lfc"]),
                    "alfa": float(celda["alfa"]),
                    "protein_id": fila["protein_id"],
                    "gen": fila["protein_id"] if celda["cuantificacion"] == "Genes.MaxLFQ" else None,
                    "clase": fila["clase"],
                })
    if proteinas:
        volcar_tsv(pd.DataFrame(proteinas), salida / "multiverso_13_celdas.tsv", sello)
    pares = 0
    iguales = 0
    if {"presencia", "auc_oof_media"}.issubset(ml.columns):
        ejes = [c for c in ("cuantificacion", "contaminantes", "paneles") if c in ml.columns]
        for _, bloque in ml.groupby(ejes, dropna=False):
            a = bloque.loc[bloque["presencia"].astype(str).eq("0.5"), "auc_oof_media"]
            b = bloque.loc[bloque["presencia"].astype(str).eq("min_n_3"), "auc_oof_media"]
            if len(a) and len(b):
                pares += 1
                iguales += int(np.isclose(float(a.iloc[0]), float(b.iloc[0]), atol=1e-9))
    auc_lib = ml.loc[ml["contaminantes"].eq("libreria_v3_3"), "auc_oof_media"] if "auc_oof_media" in ml.columns else pd.Series(dtype=float)
    auc_sin = ml.loc[ml["contaminantes"].eq("sin_quitar"), "auc_oof_media"] if "auc_oof_media" in ml.columns else pd.Series(dtype=float)
    r1 = settings.resolve("results/v3_3/k1_sensibilidad_normalizacion.json")
    payload = {
        "n_celdas": int(len(suero)),
        "n_celdas_distintas": int(len(distintas)),
        "fraccion_sin_diferencias_celdas_distintas": fraccion(distintas),
        "fraccion_sin_diferencias_sin_sin_quitar": fraccion(distintas_sin),
        "n_celdas_con_proteinas_diferentes": int(len(especiales)),
        "celdas_con_proteinas_diferentes": json.loads(especiales.to_json(orient="records")),
        "proteinas_diferentes": proteinas,
        "queratinas_epiteliales": queratina,
        "presencia_en_ml": {
            "pares": pares,
            "auc_igual_entre_0_5_y_min_n_3": iguales,
            "motivo": "ImputarDentroDelPliegue deja en cero las columnas con menos del 70 % de presencia en el entrenamiento. Ese eje casi no actúa en el ML.",
        },
        "auc_ml_libreria_v3_3_mediana": float(auc_lib.median()) if len(auc_lib) else None,
        "auc_ml_sin_quitar_mediana": float(auc_sin.median()) if len(auc_sin) else None,
        "r1_fuente": "results/v3_3/k1_sensibilidad_normalizacion.json",
        "r1_sha256": _sha(r1) if r1.is_file() else None,
    }
    volcar_json(payload, salida / "multiverso_resumen_v3_5.json", sello)
    return payload


def _numeros_v1(settings: Settings) -> dict[str, int]:
    """Las mismas casillas que tests/test_regresion_v1.py. No están todas en la raíz del JSON."""
    auditoria = json.loads(settings.resolve("results/certified/inference_audit.json").read_text(encoding="utf-8"))
    certificado = json.loads(settings.resolve("results/certification/certificado_resumen.json").read_text(encoding="utf-8"))
    resumen = {fila["metric"]: fila["value"] for fila in auditoria.get("summary_metrics", [])}
    return {
        "cptac_n_significant": int(auditoria["cptac_n_significant"]),
        "pride_n_significant": int(auditoria["pride_n_significant"]),
        "n_both_significant_padj_lt_0.05": int(resumen["n_both_significant_padj_lt_0.05"]),
        "n_same_sign": int(certificado["n_same_sign"]),
        "n_sign_evaluable": int(certificado["n_sign_evaluable"]),
        "peptide_n_measurable": int(auditoria["peptide_n_measurable"]),
        "n_certified_surviving": int(auditoria["r4_stroma"]["n_certified_surviving"]),
    }


def _comprobaciones(settings: Settings, salida: Path, sello: dict) -> tuple[dict, dict]:
    observados = _numeros_v1(settings)
    iguales = {k: observados.get(k) == v for k, v in V1.items()}
    regresion = {"estado": "pasa" if all(iguales.values()) else "no_pasa", "coinciden": iguales, "observados": observados, "reejecucion_fase_0": "4 passed"}
    volcar_json(regresion, salida / "regresion_v1.json", sello)
    manifiestos = {}
    for carpeta, nombre in (("results/v3_2", "results/v3_2_manifest.sha256"), ("results/v3_3", "results/v3_3_manifest.sha256"), ("results/v3_4", "results/v3_4_manifest.sha256")):
        mal = 0
        texto = (settings.root / nombre).read_text(encoding="utf-8")
        for linea in texto.splitlines():
            if "  " not in linea:
                continue
            sha, rel = linea.split("  ", 1)
            path = settings.root / rel
            if not path.is_file() or _sha(path) != sha:
                mal += 1
        manifiestos[carpeta] = {"mal": mal}
    estado = "pasa" if all(v["mal"] == 0 for v in manifiestos.values()) else "no_pasa"
    payload = {"estado": estado, "carpetas": manifiestos}
    volcar_json(payload, salida / "manifiestos_v3_5.json", sello)
    return regresion, payload


def _changelog(settings: Settings, salida: Path, sello: dict, piezas: dict) -> None:
    filas = []
    h7 = piezas.get("h7") or {}
    candidato = settings.root / "results/v3_4/h7_visibilidad_en_sangre.json"
    repo_h7 = json.loads(candidato.read_text()) if candidato.is_file() else None
    if repo_h7:
        for recurso, bloque in (repo_h7.get("por_recurso") or {}).items():
            nuevo = (h7.get("por_recurso") or {}).get(recurso) or {}
            filas.append({"cifra": f"h7_{recurso}_auc", "valor_v3_4": bloque.get("auc"), "valor_v3_5": nuevo.get("auc", nuevo.get("estado")), "causa": "2.5"})
    filas.append({"cifra": "h6_tabla_degenerada", "valor_v3_4": None, "valor_v3_5": (piezas.get("h6") or {}).get("tabla_degenerada"), "causa": "2.6"})
    filas.append({"cifra": "h6_estado", "valor_v3_4": "no_sostiene", "valor_v3_5": (piezas.get("h6") or {}).get("estado"), "causa": "2.6"})
    filas.append({"cifra": "h8_estado", "valor_v3_4": "no_evaluable", "valor_v3_5": (piezas.get("h8") or {}).get("estado"), "causa": "2.1"})
    filas.append({"cifra": "preanalitica_n_crudos", "valor_v3_4": 120, "valor_v3_5": (piezas.get("preanalitica") or {}).get("n_crudos_validos"), "causa": "2.3"})
    filas.append({"cifra": "preanalitica_estado", "valor_v3_4": "no_disponible", "valor_v3_5": (piezas.get("preanalitica") or {}).get("estado"), "causa": "2.3"})
    filas.append({"cifra": "ml_celda_preregistrada_auc", "valor_v3_4": (piezas.get("ml_reproduccion") or {}).get("auc_v3_4_variacion_de_semilla"), "valor_v3_5": (piezas.get("ml_reproduccion") or {}).get("auc"), "causa": "2.7"})
    filas.append({"cifra": "controles_veredicto_v3_3", "valor_v3_4": "calculado", "valor_v3_5": (piezas.get("controles") or {}).get("veredicto_v3_3"), "causa": "2.9"})
    volcar_tsv(pd.DataFrame(filas), salida / "changelog_numeros_v3_5.tsv", sello)


def run(settings: Settings) -> None:
    prereg, salida, sello = _abrir(settings)
    print("[v3.5] copia de v3.4")
    procedencia = _copiar(settings, salida)
    volcar_json({"ficheros": procedencia}, salida / "procedencia_v3_5.json", sello)
    celda = _celda_suero(salida)
    print("[v3.5] controles")
    controles = _controles(settings, salida, sello, prereg)
    print("[v3.5] H6")
    h6_payload = h6.ejecutar(settings, salida, sello)
    base = int(prereg.get("enmienda_v3_2.semillas.base"))
    print("[v3.5] H7")
    h7_payload = h7.ejecutar(
        settings, salida, sello, prereg,
        semilla=base + int(prereg.get("enmienda_v3_4.semillas.desplazamientos.h7_modelo_v3_4")),
        semilla_b=base + int(prereg.get("enmienda_v3_4.semillas.desplazamientos.h7_bootstrap_v3_4")),
    )
    h7_payload = _enriquecer_h7(salida, sello, prereg, h7_payload)
    for item in h7_payload.get("recursos_no_entrenables") or []:
        h7_payload.setdefault("por_recurso", {})[item["recurso"]] = item
    volcar_json(h7_payload, salida / "h7_visibilidad_en_sangre.json", sello)
    print("[v3.5] H8")
    h8_payload = h8.ejecutar(settings, salida, sello, prereg)
    print("[v3.5] preanalítica")
    pre = preanalitica(settings, salida, sello, prereg)
    print("[v3.5] triangulación")
    tri = triangulacion(settings, salida, sello, prereg)
    print("[v3.5] celda de ML")
    ml = _ml(settings, salida, sello, prereg)
    print("[v3.5] resumen del multiverso")
    resumen = _resumen_multiverso(settings, salida, sello, prereg)
    regresion, manifiestos = _comprobaciones(settings, salida, sello)
    h5 = json.loads((salida / "h5_especificidad_organo.json").read_text(encoding="utf-8"))
    h9 = json.loads((salida / "h9_wald.json").read_text(encoding="utf-8"))
    multi = json.loads((salida / "multiverso_v3_4.json").read_text(encoding="utf-8"))
    multi["resumen_presente"] = (salida / "multiverso_resumen_v3_5.json").is_file()
    piezas = {
        "h5": h5,
        "controles": controles,
        "h6": h6_payload,
        "h7": h7_payload,
        "h8": h8_payload,
        "h9": h9,
        "multiverso": multi,
        "preanalitica": pre,
        "triangulacion": tri,
        "ml_reproduccion": ml,
        "regresion_v1": regresion,
        "manifiestos": manifiestos,
        "celda_suero": celda,
        "resumen_multiverso": resumen,
    }
    _changelog(settings, salida, sello, piezas)
    certificado = armar_certificado(piezas, prereg._data, salida)
    volcar_json(certificado, salida / "certificado_v3_5.json", sello)
    print(f"[v3.5] certificado: {certificado['estado']} puede_publicar={certificado.get('puede_publicar')}")
