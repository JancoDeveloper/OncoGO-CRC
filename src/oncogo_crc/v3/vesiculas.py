"""Presencia en vesículas de tumor, PXD033236.

La variable es 1 solo si el resultado de Proteome Discoverer trae un q-valor
de proteína y el umbral de FDR que el propio archivo dejó escrito, y el Excel
publicado del mismo estudio tiene abundancia finita en un canal de tumor.
"""

from __future__ import annotations

import re
import sqlite3
import struct
from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.v3.rutas_publicas import ruta_publica

ARCHIVOS = (
    ("2815_11plex_0miss_SPhuman_norm_final.pdResult", "Proteomic_profiling_PDresults2815_Study1.xlsx", 2275368960),
    ("2991_11plex_0miss_SPhuman_norm_final.pdResult", "Proteomic_profiling_PDresults2991_Study2.xlsx", 2409832448),
)
GN = re.compile(r"GN=([A-Za-z0-9]+)")


def etiqueta(q_valor: float, umbral: float, cot_finito: bool) -> float:
    """1 si el q-valor no supera el umbral guardado en el estudio y hay tumor."""
    if not np.isfinite(q_valor) or not np.isfinite(umbral):
        return float("nan")
    if float(q_valor) <= float(umbral) and cot_finito:
        return 1.0
    return 0.0


def por_gen(filas: pd.DataFrame) -> pd.Series:
    """Una proteína entra si alguno de sus accesos entra en alguno de los dos estudios."""
    if filas.empty:
        return pd.Series(dtype=float)
    return filas.groupby("gen")["etiqueta"].max()


def _completo(path: Path, nbytes: int) -> bool:
    return path.is_file() and path.stat().st_size == nbytes


def _tablas(con: sqlite3.Connection) -> list[str]:
    return [fila[0] for fila in con.execute("SELECT name FROM sqlite_master WHERE type='table'")]


def umbral_y_valores(con: sqlite3.Connection) -> tuple[float | None, pd.DataFrame, str]:
    """Lee el umbral escrito en el log del workflow y el q-valor de cada proteína maestra."""
    if "TargetProteins" not in _tablas(con) or "WorkflowMessages" not in _tablas(con):
        return None, pd.DataFrame(), "faltan TargetProteins o WorkflowMessages"
    umbral = _umbral_del_log(con)
    if umbral is None:
        return None, pd.DataFrame(), "el log del workflow no deja un único fdr=número"
    filas = []
    for acceso, q_blob, descripcion in con.execute(
        "SELECT Accession, Expqvalue, Description FROM TargetProteins WHERE IsMasterProtein = 0"
    ):
        filas.append({"acceso": _acceso(acceso), "q": _q_blob(q_blob), "gen": _gen(descripcion)})
    valores = pd.DataFrame(filas)
    if valores.empty or not valores["q"].notna().any():
        return None, pd.DataFrame(), "Expqvalue no se pudo leer como número"
    return umbral, valores.dropna(subset=["acceso"]), "Expqvalue"


def _umbral_del_log(con: sqlite3.Connection) -> float | None:
    """El número sale de las líneas que el propio Proteome Discoverer escribió, como fdr=0.01."""
    encontrados: set[float] = set()
    for (mensaje,) in con.execute("SELECT Message FROM WorkflowMessages WHERE Message LIKE '%fdr=%'"):
        for hallado in re.finditer(r"fdr=([0-9]*\.[0-9]+|[0-9]+)", str(mensaje), flags=re.IGNORECASE):
            numero = float(hallado.group(1))
            if 0 < numero < 1:
                encontrados.add(numero)
    if len(encontrados) == 1:
        return encontrados.pop()
    return None


def _q_blob(blob: object) -> float:
    if not isinstance(blob, (bytes, bytearray)) or len(blob) < 8:
        return float("nan")
    return float(struct.unpack_from("<d", blob)[0])


def _acceso(texto: object) -> str:
    partes = str(texto).split("|")
    nucleo = partes[1] if len(partes) > 1 else partes[0]
    return nucleo.split("-")[0].strip()


def _gen(texto: object) -> str:
    hallado = GN.search(str(texto))
    return hallado.group(1).upper() if hallado else ""


def _el_excel_es_ese_umbral(excel: Path, valores: pd.DataFrame, umbral: float) -> bool:
    """El Excel publicado tiene que ser exactamente las proteínas con q-valor dentro del umbral."""
    publicados = set(pd.read_excel(excel, usecols=["Accession"])["Accession"].map(_acceso))
    pasan = set(valores.loc[valores["q"] <= umbral, "acceso"])
    return publicados == pasan


def _cot_finito(excel: Path) -> set[str]:
    tabla = pd.read_excel(excel)
    canales = [c for c in tabla.columns if str(c).startswith("CoT")]
    if not canales or "Accession" not in tabla.columns:
        return set()
    numeros = tabla[canales].apply(pd.to_numeric, errors="coerce")
    presentes = numeros.notna().any(axis=1)
    return set(tabla.loc[presentes, "Accession"].map(_acceso))


def anotar(settings: Settings, tabla: pd.DataFrame, variables: list[dict]) -> bool:
    """Escribe la columna. Devuelve False si el FDR del archivo no se puede leer."""
    base = settings.resolve("external/pxd033236")
    bloques = []
    umbrales = []
    for pd_nombre, xlsx_nombre, nbytes in ARCHIVOS:
        pd_path = base / pd_nombre
        xlsx = base / xlsx_nombre
        if not _completo(pd_path, nbytes) or not xlsx.is_file():
            return False
        with sqlite3.connect(f"file:{pd_path}?mode=ro", uri=True) as con:
            umbral, valores, _columna = umbral_y_valores(con)
        if umbral is None or valores.empty:
            return False
        tumor = _cot_finito(xlsx)
        if not _el_excel_es_ese_umbral(xlsx, valores, umbral):
            return False
        valores = valores.loc[valores["gen"].ne("")].copy()
        valores["etiqueta"] = [
            etiqueta(q, umbral, acceso in tumor)
            for q, acceso in zip(valores["q"], valores["acceso"], strict=True)
        ]
        bloques.append(valores[["gen", "etiqueta"]])
        umbrales.append(umbral)
    if len(set(umbrales)) != 1:
        return False
    serie = por_gen(pd.concat(bloques, ignore_index=True))
    tabla["presencia_en_vesiculas_pxd033236"] = serie.reindex(tabla.index).fillna(0.0)
    variables.append({
        "variable": "presencia_en_vesiculas_pxd033236",
        "fuente": ruta_publica(base, settings)["ruta"],
        "definicion": (
            "1 si el q-valor de proteína del pdResult no supera el umbral de FDR "
            f"guardado en el archivo ({umbrales[0]}) y hay abundancia finita en un canal CoT"
        ),
    })
    return True
