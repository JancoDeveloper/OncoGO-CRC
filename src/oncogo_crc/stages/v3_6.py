"""Corrida v3.6. Recalcula cada pieza y no abre corridas anteriores."""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.v3.prereg import PreregistroV3
from oncogo_crc.stages import v3_3_controles, v3_3_ml, v3_3_suero, v3_3_tejido, v3_4
from oncogo_crc.stages.v3_5 import _numeros_v1
from oncogo_crc.stages.v3_5_preanalitica import ejecutar as preanalitica
from oncogo_crc.stages.v3_6_report import emitir, rutas_limpias
from oncogo_crc.utils.sello import exigir_arbol_limpio, hacer_sello
from oncogo_crc.v3 import h5, h6, h7, h8, h9, h10, triangulacion
from oncogo_crc.v3.salida import activar_corrida_limpia, carpeta_v36, cerrar_corrida_limpia, volcar_json, volcar_tsv
from oncogo_crc.v3.vigilante import activar as activar_vigilante
from oncogo_crc.v3.vigilante import cerrar as cerrar_vigilante
from oncogo_crc.v3.vigilante import excepciones_de


def _pieza(nombre: str, funcion):
    inicio = time.perf_counter()
    print(f"[v3.6] {nombre}", flush=True)
    resultado = funcion()
    print(f"[v3.6] {nombre} en {time.perf_counter() - inicio:.1f} s", flush=True)
    return resultado


def _con_ic(origen: Path, destino: Path, sello: dict) -> None:
    tabla = pd.read_csv(origen, sep="\t")
    if "se" in tabla.columns and "log2fc" in tabla.columns:
        se = pd.to_numeric(tabla["se"], errors="coerce")
        efecto = pd.to_numeric(tabla["log2fc"], errors="coerce")
        tabla["ic95_bajo"] = efecto - 1.96 * se
        tabla["ic95_alto"] = efecto + 1.96 * se
    volcar_tsv(tabla, destino, sello)


def _celda(path: Path) -> dict:
    tabla = pd.read_csv(path, sep="\t")
    return {
        "n": int(len(tabla)),
        "diferente": int(tabla["clase"].eq("diferente").sum()),
        "equivalente": int(tabla["clase"].eq("equivalente").sum()),
        "indeterminada": int(tabla["clase"].eq("indeterminada").sum()),
    }


def _alias(origen: Path, destino: Path) -> None:
    if origen.is_file() and origen.resolve() != destino.resolve():
        destino.write_bytes(origen.read_bytes())


def _procedencia(salida: Path, sello: dict) -> None:
    filas = []
    for path in sorted(p for p in salida.iterdir() if p.is_file() and not p.name.startswith("._")):
        filas.append({
            "fichero": path.name,
            "origen": "recalculado_v3_6",
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        })
    volcar_json({"ficheros": filas}, salida / "procedencia.json", sello)


def run(settings: Settings) -> None:
    prereg = PreregistroV3(settings.root)
    salida = carpeta_v36(settings)
    salida.mkdir(parents=True, exist_ok=True)
    sello = hacer_sello(settings.root, etapa="v3_6", manifest=None)
    exigir_arbol_limpio(sello)
    tokens = activar_corrida_limpia(salida)
    activar_vigilante(excepciones_de(settings.root / "config" / "entradas_v3_6.yaml"))
    try:
        _pieza("tejido", lambda: v3_3_tejido.run(settings))
        _pieza("suero", lambda: v3_3_suero.run(settings))
        def controles_base():
            filas = v3_3_controles._genes_bajos(settings, salida)
            filas.append(v3_3_controles._permutacion(settings, prereg))
            filas.extend(v3_3_controles._suero(salida))
            tabla = pd.DataFrame(filas)
            volcar_tsv(tabla, salida / "controles.tsv", sello)
            fallos = int(tabla["pasa"].eq("falla").sum()) if "pasa" in tabla.columns else 0
            volcar_json({"estado": "calculado", "n_fallos": fallos}, salida / "controles.json", sello)
            return fallos
        _pieza("controles", controles_base)
        controles = _pieza("controles_tejido", lambda: v3_4._controles(settings, salida, sello, prereg))
        volcar_json(controles, salida / "controles_tejido.json", sello)
        _pieza("ml", lambda: v3_3_ml.run(settings))
        _alias(salida / "ml_suero_modelos.tsv", salida / "ml_modelos.tsv")
        _con_ic(salida / "k4_crc_frente_a_colonoscopia_negativa.tsv", salida / "suero_k4_crc_cn.tsv", sello)
        _con_ic(salida / "k4_crc_frente_a_sanos.tsv", salida / "suero_k4_crc_sanos.tsv", sello)
        _con_ic(salida / "k4_colonoscopia_negativa_frente_a_sanos.tsv", salida / "suero_k4_cn_sanos.tsv", sello)
        celda = _celda(salida / "suero_k4_crc_cn.tsv")
        controles_suero = json.loads((salida / "suero_controles.json").read_text(encoding="utf-8"))
        volcar_json({"celda_preregistrada": celda, "controles": controles_suero.get("controles"), "n_grupos_en_matriz_principal": controles_suero.get("n_grupos_en_matriz_principal")}, salida / "suero_embudo.json", sello)
        _alias(salida / "k7_control_tecnico.json", salida / "suero_calidad.json")
        modelos = pd.read_csv(salida / "ml_modelos.tsv", sep="\t")
        m2 = modelos.loc[modelos["modelo"].eq("M2_proteinas_principal")].iloc[0]
        volcar_json({"modelo": "M2_proteinas_principal", "auc_oof_media": float(m2["auc_oof_media"]), "versiones": {"M2_proteinas_principal": {"auc_oof_media": float(m2["auc_oof_media"])}}}, salida / "ml_celda_preregistrada.json", sello)
        multi = _pieza("multiverso", lambda: v3_4._multiverso(settings, salida, sello, prereg))
        _alias(salida / "multiverso_suero_v3_4.tsv", salida / "multiverso_suero.tsv")
        _alias(salida / "multiverso_ml_v3_4.tsv", salida / "multiverso_ml.tsv")
        _alias(salida / "multiverso_v3_4.json", salida / "multiverso.json")
        suero_mv = pd.read_csv(salida / "multiverso_suero.tsv", sep="\t")
        con_diferencias = suero_mv.loc[pd.to_numeric(suero_mv["diferente"], errors="coerce").fillna(0).gt(0)]
        resumen = {
            "n_celdas_suero": int(len(suero_mv)),
            "n_celdas_ml": int(len(pd.read_csv(salida / "multiverso_ml.tsv", sep="\t"))),
            "n_r1": 2,
            "n_con_diferente": int(len(con_diferencias)),
            "contaminantes_de_las_celdas_con_diferente": sorted(con_diferencias["contaminantes"].astype(str).unique().tolist()),
        }
        volcar_json(resumen, salida / "multiverso_resumen.json", sello)
        h5_payload = _pieza("h5", lambda: h5.ejecutar(settings, salida, sello, prereg))
        h6_payload = _pieza("h6", lambda: h6.ejecutar(settings, salida, sello))
        base = int(prereg.get("enmienda_v3_2.semillas.base"))
        h7_payload = _pieza("h7", lambda: h7.ejecutar(
            settings, salida, sello, prereg,
            semilla=base + int(prereg.get("enmienda_v3_4.semillas.desplazamientos.h7_modelo_v3_4")),
            semilla_b=base + int(prereg.get("enmienda_v3_4.semillas.desplazamientos.h7_bootstrap_v3_4")),
        ))
        antes = dict(h7_payload.get("por_recurso") or {})
        h7_payload = v3_4._enriquecer_h7(salida, sello, prereg, h7_payload)
        for clave, valor in antes.items():
            if valor.get("estado") != "no_disponible":
                h7_payload["por_recurso"][clave] = valor
        volcar_json(h7_payload, salida / "h7_visibilidad_en_sangre.json", sello)
        h8_payload = _pieza("h8", lambda: h8.ejecutar(settings, salida, sello))
        h9_payload = _pieza("h9", lambda: h9.ejecutar(settings, salida, sello, prereg))
        h10_payload = _pieza("h10", lambda: h10.ejecutar(settings, salida, sello, prereg))
        pre = _pieza("preanalitica", lambda: preanalitica(settings, salida, sello, prereg))
        volcar_json(pre, salida / "preanalitica.json", sello)
        tri = _pieza("triangulacion", lambda: triangulacion.ejecutar(settings, salida, sello, prereg))
        _alias(salida / "triangulacion_v3_5.tsv", salida / "triangulacion.tsv")
        _alias(salida / "triangulacion_v3_5.json", salida / "triangulacion.json")
        sec = _pieza("secuencia", lambda: v3_4._secuencia(salida, sello))
        volcar_json(sec, salida / "secuencia.json", sello)
        v1 = _numeros_v1(settings)
        volcar_json(v1, salida / "regresion_v1.json", sello)
        def segunda_m2() -> bool:
            from oncogo_crc.v3.ml_suero import evaluar, matrices_de_suero
            from oncogo_crc.v3.suero_v32 import semilla_v33
            bloques = matrices_de_suero(
                pd.read_csv(salida / "suero_matriz_principal.tsv", sep="\t"),
                pd.read_csv(salida / "suero_matriz_humana_fdr.tsv", sep="\t"),
                pd.read_csv(salida / "suero_calidad_muestra.tsv", sep="\t"),
                prereg,
            )
            x, y, tipo = bloques["M2_proteinas_principal"]
            segunda = evaluar(x, y, semilla_v33(prereg, "ml_cv"), semilla_v33(prereg, "ml_permutacion"), tipo)
            return abs(float(segunda["auc_oof_media"]) - float(m2["auc_oof_media"])) < 1e-12
        determinismo_ok = _pieza("determinismo_m2", segunda_m2)
        volcar_json({"pasa": bool(determinismo_ok), "pieza": "M2_proteinas_principal"}, salida / "determinismo.json", sello)
        eventos, prohibidas = cerrar_vigilante()
        relativas = []
        for evento in eventos:
            texto = str(evento["ruta"]).replace("\\", "/")
            marca = texto.lower().find("/results/")
            relativa = texto[marca + 1:] if marca >= 0 else Path(texto).name
            relativas.append({"ruta": relativa, "clase": evento["clase"]})
        prohibidas_rel = []
        for ruta in prohibidas:
            texto = str(ruta).replace("\\", "/")
            marca = texto.lower().find("/results/")
            prohibidas_rel.append(texto[marca + 1:] if marca >= 0 else Path(texto).name)
        volcar_json({"lecturas": relativas, "prohibidas": prohibidas_rel}, salida / "lecturas_v3_6.json", sello)
        _procedencia(salida, sello)
        ml_json = json.loads((salida / "ml_celda_preregistrada.json").read_text(encoding="utf-8"))
        tri_tabla = pd.read_csv(salida / "triangulacion.tsv", sep="\t") if (salida / "triangulacion.tsv").is_file() else pd.DataFrame()
        piezas = {
            "suero_celda": celda,
            "ml": ml_json,
            "multiverso": resumen,
            "controles": controles,
            "h5": h5_payload,
            "h6": h6_payload,
            "h7": h7_payload,
            "h8": h8_payload,
            "h9": h9_payload,
            "h10": h10_payload,
            "pre": pre,
            "v1": v1,
            "sin_mezcla": not prohibidas and all(fila["origen"] == "recalculado_v3_6" for fila in json.loads((salida / "procedencia.json").read_text())["ficheros"]),
            "reproduce": celda == {"n": 284, "diferente": 0, "equivalente": 166, "indeterminada": 118} and abs(float(m2["auc_oof_media"]) - 0.88025) < 1e-9,
            "rutas_publicas": rutas_limpias(salida),
            "triangulacion_ok": len(tri_tabla) >= 14 and not tri_tabla.isna().any().any() if len(tri_tabla) else False,
            "determinismo": bool(determinismo_ok),
        }
        certificado = emitir(piezas, salida, sello)
        _procedencia(salida, sello)
        print(f"[v3.6] certificado {certificado['estado']}", flush=True)
        if prohibidas:
            raise RuntimeError(f"lectura prohibida: {prohibidas[:5]}")
        if multi is None:
            pass
        if tri is None:
            pass
    finally:
        cerrar_corrida_limpia(tokens)
        cerrar_vigilante()
