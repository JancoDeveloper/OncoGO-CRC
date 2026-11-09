"""K2 y K3. No fusiona cohortes. Deming usa varianza de error igual en los dos ejes."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.v3.salida import volcar_json, volcar_tsv


def pendiente_deming(x: np.ndarray, y: np.ndarray, razon_de_varianzas: float = 1.0) -> dict:
    """Regresión de Deming. La razón de varianzas queda en 1 antes de ver la pendiente."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3 or not np.isfinite(razon_de_varianzas) or razon_de_varianzas <= 0:
        return {"n": int(x.size), "pendiente": float("nan"), "intercepto": float("nan"), "razon_de_varianzas": razon_de_varianzas}
    mx, my = float(x.mean()), float(y.mean())
    sxx = float(np.sum((x - mx) ** 2))
    syy = float(np.sum((y - my) ** 2))
    sxy = float(np.sum((x - mx) * (y - my)))
    if sxy == 0:
        pendiente = float("nan")
    else:
        disco = (syy - razon_de_varianzas * sxx) ** 2 + 4.0 * razon_de_varianzas * sxy ** 2
        pendiente = (syy - razon_de_varianzas * sxx + float(np.sqrt(disco))) / (2.0 * sxy)
    return {
        "n": int(x.size),
        "pendiente": float(pendiente),
        "intercepto": float(my - pendiente * mx) if np.isfinite(pendiente) else float("nan"),
        "razon_de_varianzas": float(razon_de_varianzas),
    }


def filas_de_independencia(estado: dict) -> list[dict]:
    """Un par con hash compartido, o sin correlación de nivel B, no se fusiona."""
    filas = []
    for bloque in estado.get("nivel_a_hashes_de_pride", []):
        n_hash = int(bloque.get("n_hashes_compartidos") or 0)
        par = str(bloque["par"])
        if n_hash > 0:
            decision = "una_sola_cohorte"
            motivo = "hay ficheros con el mismo hash"
        else:
            decision = "hash_distinto"
            motivo = "ningún fichero compartido; el nivel B decide si existe"
        filas.append({
            "par": par,
            "nivel": "A_hash",
            "observado": n_hash,
            "decision": decision,
            "motivo": motivo,
        })
    nivel_b = estado.get("nivel_b_correlacion_entre_muestras", {})
    vistos = {fila["par"] for fila in filas}
    for par in vistos:
        bloque = nivel_b.get(par)
        if not bloque:
            filas.append({
                "par": par,
                "nivel": "B_correlacion",
                "observado": None,
                "decision": "una_sola_cohorte",
                "motivo": "no se descarta el solapamiento: falta la correlación muestra a muestra",
            })
            continue
        sospechosos = int(bloque.get("n_pares_sospechosos") or 0)
        filas.append({
            "par": par,
            "nivel": "B_correlacion",
            "observado": bloque.get("r_maximo"),
            "decision": "una_sola_cohorte" if sospechosos > 0 else "independientes",
            "motivo": str(bloque.get("veredicto")),
        })
    return filas


def escribir_k2(settings: Settings, salida: Path, sello: dict) -> None:
    ruta = settings.resolve("external/independencia/estado.json")
    if not ruta.is_file():
        volcar_json({"estado": "no_disponible", "motivo": str(ruta)}, salida / "k2_independencia.json", sello)
        return
    estado = json.loads(ruta.read_text(encoding="utf-8"))
    filas = filas_de_independencia(estado)
    volcar_tsv(pd.DataFrame(filas), salida / "k2_independencia.tsv", sello)
    volcar_json({
        "estado": "calculado",
        "fuente": str(ruta),
        "fecha_de_la_fuente": estado.get("fecha_de_acceso"),
        "cohorte_del_contraste_de_tejido": "PXD034575",
        "no_se_fusiona": sorted({fila["par"] for fila in filas if fila["decision"] == "una_sola_cohorte"}),
    }, salida / "k2_independencia.json", sello)


def escribir_k3(cruzado: pd.DataFrame, salida: Path, sello: dict) -> None:
    ajuste = pendiente_deming(
        cruzado["log2fc_cptac"].to_numpy(),
        cruzado["log2fc_replica"].to_numpy(),
    )
    ajuste["estado"] = "calculado"
    ajuste["ejes"] = "log2fc_cptac frente a log2fc_tanaka_primario_mucosa"
    ajuste["linea_base_otros_carcinomas"] = "pendiente_de_h5"
    ajuste["supuesto"] = "La varianza de error de los dos log2FC se toma igual. No se estimó con los datos."
    volcar_json(ajuste, salida / "k3_deming.json", sello)
