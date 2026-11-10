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
    """Devuelve nombre de fichero → SHA-1. Quita la ruta de Windows del depositante."""
    salida: dict[str, str] = {}
    for linea in texto.splitlines():
        encontrado = SHA1_RE.search(linea)
        if not encontrado:
            continue
        nombre = linea[encontrado.end():].strip().replace("\\", "/").split("/")[-1]
        if nombre:
            salida[nombre.lower()] = encontrado.group(1).lower()
            salida[nombre] = encontrado.group(1).lower()
    return {nombre: digest for nombre, digest in salida.items() if nombre.endswith(".raw") or nombre.endswith(".tsv")}


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


def descargar(url: str, destino: Path) -> None:
    destino.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["curl", "-fL", "--retry", "5", "-C", "-", "-o", str(destino), url],
        check=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--cache", type=Path, default=Path.home() / "Library" / "Caches" / "oncogo")
    args = parser.parse_args()
    raws = args.workspace / "data" / "raw" / "plasma" / "pxd055218_raw"
    cache = args.cache
    cache.mkdir(parents=True, exist_ok=True)
    checksum = raws / "checksum.txt"
    descargar(f"{URL_BASE}/checksum.txt", cache / "checksum.txt")
    texto = (cache / "checksum.txt").read_text(encoding="utf-8", errors="replace")
    if not checksum.is_file() or checksum.read_text(encoding="utf-8", errors="replace") != texto:
        checksum.write_text(texto, encoding="utf-8")
    tabla = parsear_checksum(texto)
    filas = []
    for crudo in sorted(raws.glob("[CH][NPC]_*.raw")):
        if crudo.name.startswith("._"):
            continue
        esperado = tabla.get(crudo.name)
        if not esperado:
            raise SystemExit(f"PRIDE no publica SHA-1 para {crudo.name}")
        fila = evaluar_crudo(crudo, esperado, None)
        if fila["pide_redescarga"]:
            for intento in range(1, 4):
                temporal = cache / crudo.name
                if temporal.exists():
                    temporal.unlink()
                descargar(f"{URL_BASE}/{crudo.name}", temporal)
                fila["intentos"] = intento
                if sha1_de(temporal) == esperado:
                    cuarentena = raws / "cuarentena"
                    cuarentena.mkdir(exist_ok=True)
                    crudo.rename(cuarentena / f"{crudo.name}.danado")
                    temporal.replace(crudo)
                    fila = evaluar_crudo(crudo, esperado, None)
                    fila["redescargado"] = True
                    fila["intentos"] = intento
                    break
            else:
                raise SystemExit(f"{crudo.name} no coincide con PRIDE tras 3 descargas")
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
