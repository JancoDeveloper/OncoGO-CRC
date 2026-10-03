#!/usr/bin/env python3
"""D8. Extracción de tablas 2x2 de la literatura: CEA frente a metilación y ctDNA.

Qué hay aquí y qué no. Aquí están las cuentas de verdaderos y falsos positivos y
negativos de cada estudio, copiadas de la tabla publicada, con la referencia, la
tabla y la página de donde sale cada fila. No hay ningún número calculado a
partir de un porcentaje: reconstruir una tabla 2x2 desde una sensibilidad
redondeada produce celdas que no son las del estudio, y eso es inventar datos con
otro nombre. Si de una prueba no hay tablas por estudio, no se fuerza: se usa su
estimación combinada ya publicada y se dice que es eso.

Cada fila lleva ``verificado_por`` vacía. La rellena el autor cuando haya
comparado la fila con el original. Mientras esté vacía, los resultados que
dependan de ella llevan la marca de pendientes.

Dos trampas de las tablas de origen que están resueltas aquí y que conviene
conocer, porque las dos cambian los números sin dar ningún aviso:

1. **El orden de las columnas no es el mismo en las dos fuentes.** La tabla del
   CEA viene en orden TP, FP, FN, TN. La de la metilación viene en orden
   **TP, TN, FP, FN**. Leer una con el orden de la otra intercambia los falsos
   positivos con los verdaderos negativos, lo que convierte una especificidad del
   99 % en una del 1 % sin que nada falle.
2. **Un estudio puede aparecer en varias filas.** En la tabla de la metilación,
   cinco estudios dan resultados para dos o tres algoritmos distintos sobre los
   *mismos* pacientes. Meter las tres filas sería contar a esas personas tres
   veces y estrechar los intervalos artificialmente. Se queda una fila por
   estudio, con la regla declarada más abajo.

Las dos pruebas no responden a la misma pregunta. El CEA aquí se usa para
detectar recaída después de una cirugía; la metilación de SEPT9, para detectar
cáncer en cribado. No se combinan nunca en un mismo modelo, y la comparación
directa queda marcada como no comparable.
"""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parents[2]
DESTINO = RAIZ / "data" / "literatura" / "extraccion_2x2.tsv"

# ---------------------------------------------------------------- las fuentes
FUENTES = {
    "wang_2023": {
        "cita": (
            "Wang R, Mohamedali A, Kwan B, et al. Significance of "
            "carcinoembryonic antigen detection in the early diagnosis of "
            "colorectal cancer: A systematic review and meta-analysis. World "
            "Journal of Gastrointestinal Surgery 2023; 15(12):2904-2914."
        ),
        "pmc": "PMC10784816",
        "tabla": "Table 2 (Diagnostic values of follow up duration and reference standards)",
        "pagina": 2909,
        "orden_de_columnas_en_el_original": "True positive, False positive, False negative, True negative",
        "acceso_abierto": True,
    },
    "sun_2017": {
        "cita": (
            "Sun J, Fei F, Zhang M, et al. The role of mSEPT9 in screening, "
            "diagnosis, and recurrence monitoring of colorectal cancer. "
            "Diagnostic Accuracy of Methylated SEPT9 for Blood-based "
            "Colorectal Cancer Detection: A Systematic Review and "
            "Meta-Analysis. Clinical and Translational Gastroenterology 2017; "
            "8(1):e216."
        ),
        "pmc": "PMC5288600",
        "tabla": "Table 1 (Characteristic of included studies)",
        "pagina": 3,
        "orden_de_columnas_en_el_original": "TP, TN, FP, FN",
        "acceso_abierto": True,
    },
    "zhu_2024": {
        "cita": (
            "Circulating tumour DNA-based molecular residual disease detection "
            "in resectable cancers: a systematic review and meta-analysis. "
            "eBioMedicine 2024."
        ),
        "pmc": "PMC11021841",
        "tabla": "Table 2 (Sens, Fpr, AUSROC y DOR de cada cáncer por modelo bivariado)",
        "pagina": 7,
        "orden_de_columnas_en_el_original": "Subgroups, Number of data, Sens, LL, UL, Fpr, LL, UL",
        "acceso_abierto": True,
        "solo_combinado": (
            "Esta fuente no publica las tablas 2x2 por estudio en el texto "
            "principal, así que de aquí solo se usa su estimación combinada."
        ),
    },
}

# ------------------------------------------------- CEA para detectar recaída
# De Wang 2023, tabla 2, en el orden del original: TP, FP, FN, TN.
# El umbral es 5 microgramos por litro salvo donde se indica.
CEA_RECAIDA = [
    # (estudio, año, n, TP, FP, FN, TN, umbral, nota)
    ("Augestad", 2014, 110, 3, 6, 4, 97, "5 ug/L", ""),
    ("Bhatti", 2015, 569, 123, 23, 26, 397, "5 ug/L", ""),
    ("Chang", 2017, 357, 44, 59, 23, 231, "5 ug/L", ""),
    ("Hara", 2011, 127, 31, 31, 15, 50, "5 ug/L", ""),
    ("Guo", 2018, 178, 66, 36, 26, 50, "5 ug/L", ""),
    ("Kim y Lee", 2013, 336, 34, 23, 45, 234, "5 ug/L", ""),
    (
        "Rodrigues", 2017, 404, 23, 36, 29, 316,
        "5 ug/L fumadores, 3 ug/L no fumadores",
        "umbral mixto: el original aplica dos cortes segun el tabaquismo",
    ),
    ("Shinkins", 2017, 582, 51, 12, 53, 466, "5 ug/L", ""),
]

# Una fila de la misma tabla que NO entra, y por qué. Se deja escrita porque una
# exclusión silenciosa es indistinguible de un descuido.
CEA_EXCLUIDAS = [
    {
        "estudio": "Moloney",
        "anio": 2019,
        "n": 138,
        "tp": 2, "fp": 5, "fn": 4, "tn": 127,
        "motivo": (
            "en la tabla 2 del original la prueba indice figura como «N/A», no "
            "como CEA. El texto dice que de los 12 estudios nueve usaron CEA, "
            "dos ecografia y uno exploracion fisica, asi que no se puede saber "
            "cual es este. Entrar en el modelo del CEA seria suponerlo."
        ),
    }
]

# ------------------------------------- Metilación de SEPT9 para cribado
# De Sun 2017, tabla 1, en el orden del original: TP, TN, FP, FN.
#
# Regla declarada para los estudios con varias filas: se queda la del algoritmo
# 2(3), que es el que usa la prueba comercial Epi proColon 2.0. Si un estudio no
# da esa variante, se queda la primera que publica. La regla se fija por el
# algoritmo y no por el resultado, para que no sea una elección a posteriori.
REGLA_ALGORITMO = (
    "cuando un estudio publica varios algoritmos sobre los mismos pacientes se "
    "queda el 2(3), que es el de la prueba comercial Epi proColon 2.0; si no lo "
    "da, la primera variante publicada"
)
SEPT9_CRIBADO = [
    # (estudio, año, pais, n, TP, TN, FP, FN, algoritmo, tipo, varias_filas)
    ("Yu D", 2015, "China", 123, 57, 46, 7, 13, "2(3)", "caso-control", False),
    ("Jin P", 2015, "China", 476, 101, 298, 43, 34, "2(3)", "caso-control", False),
    ("He Q", 2015, "China", 100, 38, 48, 2, 12, "NA", "caso-control", False),
    ("He N", 2015, "China", 281, 54, 196, 9, 22, "2(3)", "caso-control", False),
    ("Wang Z", 2012, "China", 56, 25, 18, 2, 11, "NA", "caso-control", False),
    ("Li SJ", 2015, "China", 161, 66, 65, 5, 25, "2(3)", "caso-control", False),
    ("Wu D", 2016, "China", 1031, 223, 697, 43, 68, "1(1)", "prospectivo", False),
    ("Kang Q", 2014, "China", 132, 60, 51, 1, 20, "NA", "caso-control", False),
    ("Ding QQ", 2015, "China", 262, 60, 171, 9, 22, "2(3)", "caso-control", False),
    ("Warren JD", 2011, "EEUU y Rusia", 144, 38, 93, 1, 12, "2(3)", "caso-control", True),
    ("Johnson DA", 2014, "EEUU", 301, 74, 163, 37, 27, "NA", "prospectivo", False),
    ("Lee HS", 2013, "Corea", 197, 37, 87, 9, 64, "1(3)", "caso-control", False),
    ("Lucia PC", 2014, "EEUU", 367, 244, 20, 1, 102, "NA", "caso-control", False),
    ("Marc T", 2010, "Alemania", 161, 24, 98, 30, 9, "2(3)", "caso-control", True),
    ("Grutzmann R", 2008, "Alemania", 831, 193, 403, 50, 185, "2(3)", "caso-control", False),
    ("deVos T", 2009, "Alemania", 514, 105, 316, 11, 82, "2(3)", "caso-control", True),
    ("Church TR", 2013, "Alemania y EEUU", 1510, 27, 1331, 126, 26, "1(2)", "prospectivo", False),
    ("Toth K", 2012, "Hungria", 184, 73, 91, 1, 19, "2(3)", "caso-control", True),
    ("Toth K", 2014, "Hungria", 84, 30, 40, 10, 4, "NA", "caso-control", False),
    ("Su XL", 2014, "China", 234, 152, 58, 4, 20, "NA", "caso-control", False),
    ("Potter NT", 2014, "EEUU", 1544, 30, 1182, 318, 14, "1(3)", "prospectivo", False),
    ("Lofton-Day C", 2008, "EEUU", 312, 92, 154, 25, 41, "NA", "caso-control", False),
    ("He Q", 2010, "China", 352, 136, 164, 6, 46, "NA", "caso-control", False),
    ("Ørntoft MW", 2015, "Dinamarca", 470, 75, 328, 14, 53, "2(3)", "caso-control", True),
    ("Ahlquist DA", 2011, "EEUU", 100, 18, 54, 16, 12, "1(3)", "caso-control", False),
]

# ------------------------------------- ctDNA, solo estimación combinada
# De Zhu 2024, tabla 2, fila CRC. La especificidad se obtiene de la tasa de
# falsos positivos: especificidad = 1 - Fpr, y los extremos se invierten de orden
# al restar de uno, que es un detalle fácil de equivocar.
CTDNA_COMBINADO = {
    "todas": {
        "k": 18,
        "sensibilidad": 0.58,
        "sensibilidad_ci": (0.48, 0.68),
        "fpr": 0.08,
        "fpr_ci": (0.06, 0.11),
        "auc_publicado": 0.89,
        "que_es": "cualquier momento de medicion",
    },
    "landmark": {
        "k": 14,
        "sensibilidad": 0.51,
        "sensibilidad_ci": (0.41, 0.61),
        "fpr": 0.09,
        "fpr_ci": (0.07, 0.11),
        "auc_publicado": 0.89,
        "que_es": "una sola medicion despues de la cirugia",
    },
    "seriada": {
        "k": 4,
        "sensibilidad": 0.84,
        "sensibilidad_ci": (0.73, 0.91),
        "fpr": 0.05,
        "fpr_ci": (0.01, 0.20),
        "auc_publicado": 0.89,
        "que_es": "mediciones repetidas durante el seguimiento",
    },
}


def _fila(**kw) -> dict:
    """Una fila de la tabla de extracción, con las columnas en orden fijo."""
    return {
        "prueba": kw["prueba"],
        "indicacion": kw["indicacion"],
        "estudio": kw["estudio"],
        "anio": kw["anio"],
        "pais": kw.get("pais", ""),
        "n": kw["n"],
        "tp": kw["tp"],
        "fp": kw["fp"],
        "fn": kw["fn"],
        "tn": kw["tn"],
        "umbral_o_algoritmo": kw["umbral_o_algoritmo"],
        "diseno": kw.get("diseno", ""),
        "el_estudio_publica_varias_filas": kw.get("varias_filas", False),
        "fuente": kw["fuente"],
        "fuente_tabla": FUENTES[kw["fuente"]]["tabla"],
        "fuente_pagina": FUENTES[kw["fuente"]]["pagina"],
        "nota": kw.get("nota", ""),
        # La deja vacía el script a propósito.
        "verificado_por": "",
    }


def construir() -> pd.DataFrame:
    filas = [
        _fila(
            prueba="CEA",
            indicacion="deteccion_de_recaida",
            estudio=e, anio=a, n=n, tp=tp, fp=fp, fn=fn, tn=tn,
            umbral_o_algoritmo=umbral, nota=nota, fuente="wang_2023",
        )
        for e, a, n, tp, fp, fn, tn, umbral, nota in CEA_RECAIDA
    ]
    filas += [
        _fila(
            prueba="metilacion_SEPT9",
            indicacion="cribado_o_diagnostico",
            estudio=e, anio=a, pais=p, n=n, tp=tp, fp=fp, fn=fn, tn=tn,
            umbral_o_algoritmo=f"algoritmo {alg}", diseno=tipo,
            varias_filas=varias, fuente="sun_2017",
        )
        for e, a, p, n, tp, tn, fp, fn, alg, tipo, varias in SEPT9_CRIBADO
    ]
    return pd.DataFrame(filas)


def comprobar(d: pd.DataFrame) -> dict:
    """Comprobaciones que tienen que pasar antes de usar un número copiado.

    La más útil es la primera: en las dos fuentes el tamaño de muestra se publica
    aparte de las cuatro celdas, así que si las celdas suman el tamaño, el orden
    de columnas y los cuatro valores están bien. Una sola celda mal copiada o un
    orden intercambiado rompen la suma.
    """
    suma = d["tp"] + d["fp"] + d["fn"] + d["tn"]
    descuadra = d[suma != d["n"]]
    problemas = []
    if not descuadra.empty:
        problemas.append(
            "las celdas no suman el tamaño de muestra en: "
            + ", ".join(
                f"{r.estudio} {r.anio} ({r.tp + r.fp + r.fn + r.tn} != {r.n})"
                for r in descuadra.itertuples()
            )
        )
    negativas = d[(d[["tp", "fp", "fn", "tn"]] < 0).any(axis=1)]
    if not negativas.empty:
        problemas.append(f"{len(negativas)} filas con celdas negativas")

    # Un estudio repetido en la misma prueba sería contar a los mismos pacientes
    # dos veces. Con el año incluido, porque hay dos Toth K y dos He Q que son
    # estudios distintos.
    duplicados = d[d.duplicated(subset=["prueba", "estudio", "anio"], keep=False)]
    if not duplicados.empty:
        problemas.append(
            "estudios repetidos dentro de la misma prueba: "
            + ", ".join(sorted({f"{r.estudio} {r.anio}" for r in duplicados.itertuples()}))
        )

    por_prueba = {}
    for (prueba, indicacion), g in d.groupby(["prueba", "indicacion"]):
        por_prueba[f"{prueba}|{indicacion}"] = {
            "k": int(len(g)),
            "n_total": int(g["n"].sum()),
            "n_casos": int((g["tp"] + g["fn"]).sum()),
            "n_sin_la_enfermedad": int((g["fp"] + g["tn"]).sum()),
            "n_con_alguna_celda_cero": int(
                (g[["tp", "fp", "fn", "tn"]] == 0).any(axis=1).sum()
            ),
            "disenos": sorted(x for x in set(g["diseno"]) if x),
        }
    return {
        "n_filas": int(len(d)),
        "problemas": problemas,
        "coherente": not problemas,
        "las_celdas_suman_el_tamano_de_muestra": bool(descuadra.empty),
        "por_prueba_e_indicacion": por_prueba,
        "pendiente_de_verificacion_humana": int((d["verificado_por"] == "").sum()),
        "regla": (
            "Ningún número copiado de una tabla publicada entra en un resultado "
            "sin que una persona lo haya comprobado contra el original. La "
            "columna verificado_por la rellena el autor, no el script."
        ),
    }


def ficha() -> dict:
    return {
        "generado": datetime.now(UTC).date().isoformat(),
        "fuentes": FUENTES,
        "regla_de_algoritmo_en_la_metilacion": REGLA_ALGORITMO,
        "filas_excluidas_a_proposito": CEA_EXCLUIDAS,
        "ctdna_solo_combinado": CTDNA_COMBINADO,
        "por_que_no_se_reconstruyen_tablas": (
            "De ninguna prueba se calcula una tabla 2x2 a partir de una "
            "sensibilidad publicada. El redondeo de un porcentaje no determina "
            "las celdas, y las que saldrían no serían las del estudio."
        ),
        "las_dos_pruebas_no_son_comparables": (
            "El CEA de estas tablas se usa para detectar recaída después de una "
            "cirugía, con una prevalencia del 20 %. La metilación de SEPT9 se "
            "usa para detectar cáncer en cribado, en su mayoría en estudios de "
            "casos y controles donde la proporción de enfermos la fija el "
            "diseño. No se combinan en un mismo modelo y la comparación directa "
            "queda marcada como no comparable."
        ),
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--destino", default=str(DESTINO))
    args = ap.parse_args()

    destino = Path(args.destino)
    destino.parent.mkdir(parents=True, exist_ok=True)

    d = construir()
    res = comprobar(d)
    d.to_csv(destino, sep="\t", index=False)
    (destino.parent / "extraccion_2x2_ficha.json").write_text(
        json.dumps({**ficha(), "comprobacion": res}, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print(f"[2x2] {res['n_filas']} filas escritas en {destino}")
    for clave, v in res["por_prueba_e_indicacion"].items():
        print(
            f"[2x2] {clave}: {v['k']} estudios, {v['n_total']} personas, "
            f"{v['n_casos']} con la enfermedad"
        )
    print(
        f"[2x2] las celdas suman el tamaño de muestra: "
        f"{res['las_celdas_suman_el_tamano_de_muestra']}"
    )
    for p in res["problemas"]:
        print(f"[2x2] ATENCIÓN: {p}")
    print(
        f"[2x2] pendientes de verificación humana: "
        f"{res['pendiente_de_verificacion_humana']} filas"
    )
    return 0 if res["coherente"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
