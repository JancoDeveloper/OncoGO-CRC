"""ESM-2 y TMbed sobre las secuencias de UniProt. No entran en el modelo principal."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.v3.rutas_publicas import ruta_publica

GN = re.compile(r"GN=([A-Za-z0-9]+)")
MODELO_ESM = (
    "external/modelos/models--facebook--esm2_t12_35M_UR50D/snapshots/"
    "6fbf070e65b0b7291e7bbcd451118c216cff79d8"
)
TMBED_PRED = "external/modelos/h7_tmbed.pred"
FASTA_H7 = "external/modelos/h7_secuencias.fasta"
FASTA_UNIPROT = "external/uniprot/human_reviewed_isoforms.fasta.gz"
ESM_TSV = "external/modelos/h7_esm2_embeddings.tsv"
ESM_OMITIDAS = "external/modelos/h7_esm2_omitidas.tsv"
VENTANA_ESM = 1022


def helices_y_senal(etiqueta: str) -> tuple[float, float]:
    """Cuenta segmentos de hélice (H/h) y si hay algún residuo de péptido señal (S)."""
    if not etiqueta:
        return float("nan"), float("nan")
    helices = 0
    dentro = False
    for letra in etiqueta:
        if letra in {"H", "h"}:
            if not dentro:
                helices += 1
                dentro = True
        else:
            dentro = False
    return float(helices), float("S" in etiqueta)


def secuencias_por_gen(fasta: Path, genes: set[str]) -> dict[str, str]:
    """Una secuencia revisada por gen. Si hay isoformas, se queda la canónica, sin guion."""
    mejor: dict[str, tuple[bool, str]] = {}
    encabezado = ""
    trozos: list[str] = []

    def cerrar() -> None:
        if not encabezado:
            return
        hallado = GN.search(encabezado)
        if not hallado:
            return
        gen = hallado.group(1).upper()
        if gen not in genes:
            return
        acceso = encabezado.split("|")[1] if "|" in encabezado else ""
        isoforma = "-" in acceso
        secuencia = "".join(trozos)
        previa = mejor.get(gen)
        if previa is None or (previa[0] and not isoforma):
            mejor[gen] = (isoforma, secuencia)

    opener = gzip_open(fasta)
    with opener as entrada:
        for linea in entrada:
            if linea.startswith(">"):
                cerrar()
                encabezado = linea[1:].strip()
                trozos = []
            else:
                trozos.append(linea.strip())
        cerrar()
    return {gen: secuencia for gen, (_iso, secuencia) in mejor.items()}


def gzip_open(path: Path):
    import gzip
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def escribir_fasta(destino: Path, secuencias: dict[str, str]) -> None:
    destino.parent.mkdir(parents=True, exist_ok=True)
    lineas = []
    for gen in sorted(secuencias):
        lineas.append(f">{gen}")
        lineas.append(secuencias[gen])
    destino.write_text("\n".join(lineas) + "\n", encoding="utf-8")


def leer_tmbed(pred: Path) -> pd.DataFrame:
    """Formato de tres líneas de TMbed: encabezado, secuencia, etiquetas."""
    if not pred.is_file():
        return pd.DataFrame(columns=["gen", "tmbed_helices_tm", "tmbed_peptido_senal"])
    texto = pred.read_text(encoding="utf-8").splitlines()
    filas = []
    for i in range(0, len(texto) - 2, 3):
        if not texto[i].startswith(">"):
            continue
        helices, senal = helices_y_senal(texto[i + 2].strip())
        filas.append({"gen": texto[i][1:].split()[0].upper(), "tmbed_helices_tm": helices, "tmbed_peptido_senal": senal})
    return pd.DataFrame(filas)


def sha256_fichero(path: Path) -> str:
    resumen = hashlib.sha256()
    with path.open("rb") as entrada:
        while True:
            bloque = entrada.read(1 << 20)
            if not bloque:
                break
            resumen.update(bloque)
    return resumen.hexdigest()


def sha256_directorio(path: Path) -> str:
    """SHA-256 de los ficheros del directorio, ordenados por su ruta relativa."""
    lineas = []
    ficheros = sorted(
        (item for item in path.rglob("*") if item.is_file() and not item.name.startswith("._")),
        key=lambda item: item.relative_to(path).as_posix(),
    )
    for fichero in ficheros:
        lineas.append(f"{fichero.relative_to(path).as_posix()}\t{sha256_fichero(fichero)}")
    return hashlib.sha256("\n".join(lineas).encode()).hexdigest()


def ventanas(secuencia: str, maxima: int = VENTANA_ESM) -> list[str]:
    """Trozos consecutivos, sin solaparse. 1022 es el contexto de ESM-2 menos los dos tokens especiales."""
    if maxima < 1:
        raise ValueError("la ventana tiene que ser de al menos un residuo")
    if not secuencia:
        return []
    return [secuencia[i:i + maxima] for i in range(0, len(secuencia), maxima)]


def media_ponderada(vectores: list[np.ndarray], longitudes: list[int]) -> np.ndarray:
    pesos = np.asarray(longitudes, dtype=np.float64)
    matriz = np.stack([np.asarray(vector, dtype=np.float64) for vector in vectores])
    return (matriz * pesos[:, None]).sum(axis=0) / float(pesos.sum())


def leer_esm2(tsv: Path, sello: Path, omitidas: Path, fasta: Path, modelo_dir: Path,
              genes: set[str], secuencias: dict[str, str]) -> tuple[pd.DataFrame | None, dict]:
    """Lee el TSV calculado fuera de la corrida. Si no cuadra, no hay vector y no se calcula nada."""
    info: dict = {
        "esm2": "no_disponible", "n_con_esm2": None, "n_omitidas_esm2": None, "esm2_sha256": None,
    }
    if not tsv.is_file() or not sello.is_file() or not omitidas.is_file():
        info["esm2_motivo"] = "falta el TSV, sus omisiones o su sello"
        return None, info
    if not fasta.is_file() or not modelo_dir.is_dir():
        info["esm2_motivo"] = "falta el FASTA o el modelo local"
        return None, info
    sello_datos = json.loads(sello.read_text(encoding="utf-8"))
    hash_fasta = sha256_fichero(fasta)
    hash_modelo = sha256_directorio(modelo_dir)
    hash_tsv = sha256_fichero(tsv)
    info["esm2_sha256"] = hash_modelo
    if sello_datos.get("sha256_fasta") != hash_fasta or sello_datos.get("sha256_modelo") != hash_modelo:
        info["esm2_motivo"] = "el sello no coincide con el FASTA o con el modelo"
        return None, info
    if sello_datos.get("sha256_tsv") != hash_tsv:
        info["esm2_motivo"] = "el sello no coincide con el TSV"
        return None, info
    tabla = pd.read_csv(tsv, sep="\t")
    columnas = [f"esm2_{i:03d}" for i in range(480)]
    if "gen" not in tabla.columns or any(columna not in tabla.columns for columna in columnas):
        info["esm2_motivo"] = "el TSV no trae las 480 columnas"
        return None, info
    tabla["gen"] = tabla["gen"].astype(str).str.upper()
    if tabla["gen"].duplicated().any():
        info["esm2_motivo"] = "el TSV repite un gen"
        return None, info
    omisiones = pd.read_csv(omitidas, sep="\t")
    if "gen" not in omisiones.columns:
        info["esm2_motivo"] = "las omisiones no traen el gen"
        return None, info
    omisiones["gen"] = omisiones["gen"].astype(str).str.upper()
    cubiertos = set(tabla["gen"]) | set(omisiones["gen"])
    if set(genes) - cubiertos:
        info["esm2_motivo"] = "hay proteínas sin vector y sin omisión"
        return None, info
    if "longitud" in tabla.columns:
        longitudes = tabla.set_index("gen")["longitud"]
        for gen, secuencia in secuencias.items():
            if gen in longitudes.index and int(longitudes.at[gen]) != len(secuencia):
                info["esm2_motivo"] = "la longitud del TSV no coincide con el FASTA"
                return None, info
    info.update({
        "esm2": "leido",
        "n_con_esm2": int(len(tabla)),
        "n_omitidas_esm2": int(omisiones["gen"].nunique()),
        "esm2_dim": 480,
        "esm2_fuente": MODELO_ESM,
    })
    return tabla.set_index("gen")[columnas], info


def anotar(settings: Settings, tabla: pd.DataFrame) -> dict:
    """Añade TMbed si la predicción está completa y lee ESM-2. La corrida no carga el modelo."""
    genes = set(tabla.index.astype(str).str.upper())
    fasta = settings.resolve(FASTA_UNIPROT)
    info = {"esm2": "no_disponible", "tmbed": "no_disponible", "n_con_secuencia": 0, "n_sin_secuencia": len(genes)}
    if not fasta.is_file():
        info["motivo"] = "no está el FASTA revisado de UniProt"
        return info
    secuencias = secuencias_por_gen(fasta, genes)
    info["n_con_secuencia"] = len(secuencias)
    info["n_sin_secuencia"] = len(genes - set(secuencias))
    pred = settings.resolve(TMBED_PRED)
    tmbed = leer_tmbed(pred)
    if len(tmbed) >= len(secuencias) and len(secuencias) > 0:
        tmbed = tmbed.drop_duplicates("gen").set_index("gen")
        tabla["tmbed_helices_tm"] = tmbed["tmbed_helices_tm"].reindex(tabla.index)
        tabla["tmbed_peptido_senal"] = tmbed["tmbed_peptido_senal"].reindex(tabla.index)
        info["tmbed"] = "calculado"
        info["tmbed_fuente"] = ruta_publica(pred, settings)["ruta"]
        info["tmbed_sha256"] = ruta_publica(pred, settings)["sha256"]
    else:
        info["tmbed"] = "no_disponible"
        info["tmbed_motivo"] = f"la predicción de TMbed tiene {len(tmbed)} proteínas y las secuencias son {len(secuencias)}"
    embeddings, extra = leer_esm2(
        settings.resolve(ESM_TSV),
        settings.resolve(ESM_TSV + ".sello.json"),
        settings.resolve(ESM_OMITIDAS),
        fasta,
        settings.resolve(MODELO_ESM),
        genes,
        secuencias,
    )
    info.update(extra)
    if embeddings is not None:
        for columna in embeddings.columns:
            tabla[columna] = embeddings[columna].reindex(tabla.index)
    return info
