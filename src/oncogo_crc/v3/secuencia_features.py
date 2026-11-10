"""ESM-2 y TMbed sobre las secuencias de UniProt. No entran en el modelo principal."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings

GN = re.compile(r"GN=([A-Za-z0-9]+)")
MODELO_ESM = (
    "external/modelos/models--facebook--esm2_t12_35M_UR50D/snapshots/"
    "6fbf070e65b0b7291e7bbcd451118c216cff79d8"
)
TMBED_PRED = "external/modelos/h7_tmbed.pred"
FASTA_H7 = "external/modelos/h7_secuencias.fasta"


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


def embeddings_esm2(settings: Settings, secuencias: dict[str, str]) -> pd.DataFrame:
    """Media de la última capa de ESM-2 35M, sin los tokens especiales."""
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    ruta = settings.resolve(MODELO_ESM)
    tokenizer = AutoTokenizer.from_pretrained(ruta, local_files_only=True)
    modelo = AutoModelForMaskedLM.from_pretrained(ruta, local_files_only=True)
    modelo.eval()
    dispositivo = torch.device("mps") if torch.backends.mps.is_available() else torch.device("cpu")
    modelo.to(dispositivo)
    filas = {}
    with torch.no_grad():
        for gen, secuencia in secuencias.items():
            try:
                tokens = tokenizer(list(secuencia), is_split_into_words=True, return_tensors="pt")
                tokens = {k: v.to(dispositivo) for k, v in tokens.items()}
                salida = modelo(**tokens, output_hidden_states=True)
                oculto = salida.hidden_states[-1][0, 1:-1]
            except RuntimeError:
                continue
            if oculto.shape[0] != len(secuencia):
                continue
            filas[gen] = oculto.mean(dim=0).cpu().numpy()
    if not filas:
        return pd.DataFrame()
    ancho = next(iter(filas.values())).shape[0]
    return pd.DataFrame.from_dict(filas, orient="index", columns=[f"esm2_{i:03d}" for i in range(ancho)])


def anotar(settings: Settings, tabla: pd.DataFrame) -> dict:
    """Añade TMbed si la predicción está completa y calcula ESM-2. No son variables del principal."""
    genes = set(tabla.index.astype(str).str.upper())
    fasta = settings.resolve("external/uniprot/human_reviewed_isoforms.fasta.gz")
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
        info["tmbed_fuente"] = str(pred)
    else:
        info["tmbed"] = "no_disponible"
        info["tmbed_motivo"] = f"la predicción de TMbed tiene {len(tmbed)} proteínas y las secuencias son {len(secuencias)}"
    try:
        import torch  # noqa: F401
    except ImportError:
        info["esm2"] = "no_disponible"
        info["esm2_motivo"] = "torch no está instalado"
        return info
    if not settings.resolve(MODELO_ESM).is_dir():
        info["esm2_motivo"] = "el snapshot de ESM-2 35M no está en el disco"
        return info
    emb = embeddings_esm2(settings, secuencias)
    if emb.empty or emb.shape[1] != 480:
        info["esm2"] = "no_disponible"
        info["esm2_motivo"] = "ESM-2 no devolvió las 480 dimensiones"
        return info
    for columna in emb.columns:
        tabla[columna] = emb[columna].reindex(tabla.index)
    info["esm2"] = "calculado"
    info["esm2_dim"] = 480
    info["esm2_fuente"] = str(settings.resolve(MODELO_ESM))
    return info
