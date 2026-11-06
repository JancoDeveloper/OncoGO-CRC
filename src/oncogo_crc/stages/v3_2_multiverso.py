"""Multiverso del suero. La celda preregistrada no se elige después de ver el resto."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.stages.v3_2_controles import _abrir
from oncogo_crc.v3.suero_v32 import _contraste, _contraste_covariable, _grupos
from oncogo_crc.v3.salida_v32 import volcar_json, volcar_tsv

LFC = (0.585, 1.0)
ALFA = (0.05, 0.10)
PRESENCIA = (None, 0.50, 0.70)


def _recortar(mat: pd.DataFrame, cols_a: list[str], cols_b: list[str], presencia: float | None) -> pd.DataFrame:
    if presencia is None:
        return mat
    a = np.isfinite(mat[cols_a].to_numpy(dtype=float)).mean(axis=1)
    b = np.isfinite(mat[cols_b].to_numpy(dtype=float)).mean(axis=1)
    return mat.loc[(a >= presencia) & (b >= presencia)]


def _conclusion(tabla: pd.DataFrame) -> str:
    return "sin_diferencia" if int((tabla["clase"] == "diferente").sum()) == 0 else "hay_diferencia"


def run(settings: Settings) -> None:
    _prereg, salida, sello = _abrir(settings, "v3_2_multiverso")
    principal_path = salida / "suero_matriz_principal.tsv"
    humana_path = salida / "suero_matriz_humana_fdr.tsv"
    if not principal_path.is_file() or not humana_path.is_file():
        raise FileNotFoundError("faltan las matrices del suero; corre v3-2-suero antes")
    principal = pd.read_csv(principal_path, sep="\t").set_index("protein_id")
    humana = pd.read_csv(humana_path, sep="\t").set_index("protein_id")
    grupos = _grupos([c for c in principal.columns])
    cols_a = grupos["crc"]
    cols_b = grupos["colonoscopy_negative"]
    calidad = pd.read_csv(salida / "suero_calidad_muestra.tsv", sep="\t")
    coag = calidad.set_index("sample_id")["coagulacion"]
    filas = []
    for lfc in LFC:
        for alpha in ALFA:
            for presencia in PRESENCIA:
                for con_contaminantes in (False, True):
                    for con_coagulacion in (False, True):
                        base = humana if con_contaminantes else principal
                        mat = _recortar(base, cols_a, cols_b, presencia)
                        if con_coagulacion:
                            tabla = _contraste_covariable(mat, cols_a, cols_b, coag, lfc=lfc, alpha=alpha, nombre="multiverso")
                        else:
                            tabla = _contraste(mat, cols_a, cols_b, lfc=lfc, alpha=alpha, nombre="multiverso")
                        oficial = (
                            lfc == 0.585 and alpha == 0.05 and presencia is None
                            and not con_contaminantes and not con_coagulacion
                        )
                        filas.append({
                            "analisis": "suero_k4",
                            "lfc": lfc,
                            "alfa": alpha,
                            "presencia": "min_n_3" if presencia is None else presencia,
                            "con_contaminantes": con_contaminantes,
                            "ajustado_por_coagulacion": con_coagulacion,
                            "cuantificacion": "PG.MaxLFQ",
                            "n": int(len(tabla)),
                            "diferente": int((tabla["clase"] == "diferente").sum()),
                            "equivalente": int((tabla["clase"] == "equivalente").sum()),
                            "indeterminada": int((tabla["clase"] == "indeterminada").sum()),
                            "conclusion": _conclusion(tabla),
                            "es_la_celda_preregistrada": oficial,
                        })
    tabla = pd.DataFrame(filas)
    oficial = tabla.loc[tabla["es_la_celda_preregistrada"], "conclusion"].iloc[0]
    hipotesis = {}
    ruta = salida / "hipotesis_v3_2.tsv"
    if ruta.is_file():
        hip = pd.read_csv(ruta, sep="\t")
        for codigo in ("H5", "H7", "H10"):
            bloque = hip.loc[hip["hipotesis"].eq(codigo)]
            hipotesis[codigo] = None if bloque.empty else str(bloque["estado"].iloc[0])
    volcar_tsv(tabla, salida / "multiverso.tsv", sello)
    payload = {
        "estado": "calculado",
        "conclusion_preregistrada": oficial,
        "fraccion_de_celdas_con_esa_conclusion": float((tabla["conclusion"] == oficial).mean()),
        "n_celdas": int(len(tabla)),
        "h5": hipotesis.get("H5", "no_calculada_en_esta_tabla"),
        "h7": hipotesis.get("H7", "no_calculada_en_esta_tabla"),
        "h10": hipotesis.get("H10", "no_calculada_en_esta_tabla"),
        "nota": "H5, H7 y H10 no dependen de estas celdas del suero. No se recalculan en cada una.",
    }
    volcar_json(payload, salida / "multiverso.json", sello)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 3))
        orden = tabla["conclusion"].value_counts()
        ax.bar(orden.index.astype(str), orden.to_numpy())
        ax.set_ylabel("celdas")
        ax.set_title("Multiverso del suero")
        fig.tight_layout()
        destino = salida / "report" / "figuras"
        destino.mkdir(parents=True, exist_ok=True)
        fig.savefig(destino / "multiverso.png", dpi=300)
        fig.savefig(destino / "multiverso.svg")
        plt.close(fig)
    except Exception as exc:
        payload["figura"] = f"no_disponible: {exc}"
        volcar_json(payload, salida / "multiverso.json", sello)
    print(f"[v3.2] multiverso {payload['fraccion_de_celdas_con_esa_conclusion']:.3f} mantienen {oficial}")
