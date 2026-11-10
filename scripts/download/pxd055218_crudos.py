#!/usr/bin/env python3
"""Comprueba los 60 crudos de PXD055218 contra el SHA-1 de checksum.txt.

Si el SHA-1 local no coincide, descarga el crudo a un disco local y lo
sustituye solo cuando el fichero nuevo coincide. No borra el dañado: lo mueve
a cuarentena/.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import subprocess
from pathlib import Path

URL_BASE = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/06/PXD055218"
SHA1_RE = re.compile(r"([0-9a-fA-F]{40})")
BLOQUE = 1 << 20


def sha1_de(path: Path) -> str:
    resumen = hashlib.sha1()
    with path.open("rb") as fichero:
        for bloque in iter(lambda: fichero.read(BLOQUE), b""):
            resumen.update(bloque)
    return resumen.hexdigest()


def parsear_checksum(texto: str) -> dict[str, str]:
    """Devuelve nombre de fichero → SHA-1. El depositante escribe ruta, tabulador y SHA-1."""
    salida: dict[str, str] = {}
    for linea in texto.splitlines():
        if not linea.strip() or linea.startswith("#"):
            continue
        partes = re.split(r"[\t ]+", linea.strip())
        digest = next((parte for parte in partes if SHA1_RE.fullmatch(parte)), "")
        ruta = next((parte for parte in partes if parte.lower().endswith((".raw", ".tsv"))), "")
        if not digest or not ruta:
            continue
        salida[ruta.replace("\\", "/").split("/")[-1]] = digest.lower()
    return salida


def evaluar_crudo(crudo: Path, sha1_pride: str, bytes_pride: int | None) -> dict:
    """Compara un crudo local con el SHA-1 publicado. No descarga."""
    digest = sha1_de(crudo)
    igual = digest == sha1_pride.lower()
    return {
        "crudo": crudo.name,
        "bytes_local": crudo.stat().st_size,
        "bytes_pride": bytes_pride if bytes_pride is not None else "",
        "sha1_local": digest,
        "sha1_pride": sha1_pride.lower(),
        "igual": igual,
        "pide_redescarga": not igual,
        "redescargado": False,
        "intentos": 0,
    }


def descargar(url: str, destino: Path, reanudar: bool = True) -> None:
    destino.parent.mkdir(parents=True, exist_ok=True)
    comando = ["curl", "-fL", "--retry", "5", "-o", str(destino), url]
    if reanudar:
        comando[4:4] = ["-C", "-"]
    subprocess.run(comando, check=True)


def anotar_redescarga(fila: dict, cuarentena: Path, nombre: str, fecha: str) -> dict:
    """El crudo en su sitio ya es el verificado si el viejo está en cuarentena."""
    marca = cuarentena / f"{nombre}.danado-{fecha}"
    if marca.is_file() and not marca.name.startswith("._"):
        fila["redescargado"] = True
    return fila


def colocar(temporal: Path, crudo: Path, fecha: str) -> None:
    """Mueve el crudo verificado aunque el caché y el volumen sean discos distintos."""
    cuarentena = crudo.parent / "cuarentena"
    cuarentena.mkdir(exist_ok=True)
    if crudo.exists():
        crudo.rename(cuarentena / f"{crudo.name}.danado-{fecha}")
    shutil.move(temporal, crudo)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--cache", type=Path, default=Path.home() / "Library" / "Caches" / "oncogo")
    args = parser.parse_args()
    raws = args.workspace / "data" / "raw" / "plasma" / "pxd055218_raw"
    cache = args.cache
    cache.mkdir(parents=True, exist_ok=True)
    checksum = raws / "checksum.txt"
    descargar(f"{URL_BASE}/checksum.txt", cache / "checksum.txt", reanudar=False)
    texto = (cache / "checksum.txt").read_text(encoding="utf-8", errors="replace")
    if not checksum.is_file() or checksum.read_text(encoding="utf-8", errors="replace") != texto:
        checksum.write_text(texto, encoding="utf-8")
    tabla = parsear_checksum(texto)
    fecha = "2026-11-10"
    cuarentena = raws / "cuarentena"
    if cuarentena.is_dir():
        for viejo in cuarentena.glob("*.raw.danado"):
            if viejo.name.startswith("._") or not viejo.is_file():
                continue
            viejo.rename(viejo.with_name(f"{viejo.name}-{fecha}"))
    filas = []
    for nombre in sorted(nombre for nombre in tabla if nombre.endswith(".raw")):
        crudo = raws / nombre
        esperado = tabla[nombre]
        temporal = cache / nombre
        if crudo.is_file():
            fila = evaluar_crudo(crudo, esperado, None)
        else:
            fila = {"crudo": nombre, "bytes_local": 0, "bytes_pride": "", "sha1_local": "", "sha1_pride": esperado, "igual": False, "pide_redescarga": True, "redescargado": False, "intentos": 0}
        if fila["pide_redescarga"]:
            for intento in range(1, 4):
                if not (temporal.is_file() and sha1_de(temporal) == esperado):
                    if temporal.exists():
                        temporal.unlink()
                    descargar(f"{URL_BASE}/{nombre}", temporal)
                fila["intentos"] = intento
                if sha1_de(temporal) == esperado:
                    colocar(temporal, crudo, fecha)
                    fila = evaluar_crudo(crudo, esperado, None)
                    fila["redescargado"] = True
                    fila["intentos"] = intento
                    break
            else:
                raise SystemExit(f"{nombre} no coincide con PRIDE tras 3 descargas")
        anotar_redescarga(fila, cuarentena, nombre, fecha)
        filas.append(fila)
    salida = Path("results") / "v3_7"
    salida.mkdir(parents=True, exist_ok=True)
    columnas = ["crudo", "bytes_local", "bytes_pride", "sha1_local", "sha1_pride", "igual", "redescargado", "intentos"]
    lineas = ["\t".join(columnas)]
    for fila in filas:
        lineas.append("\t".join(str(fila[columna]).lower() if isinstance(fila[columna], bool) else str(fila[columna]) for columna in columnas))
    (salida / "crudos_integridad.tsv").write_text("\n".join(lineas) + "\n", encoding="utf-8")
    print(f"{len(filas)} crudos; distintos {sum(1 for fila in filas if not fila['igual'])}")


if __name__ == "__main__":
    main()
