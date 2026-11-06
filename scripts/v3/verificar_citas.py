#!/usr/bin/env python3
"""Comprueba los DOI de la tesis y de docs/ contra Crossref.

Lee ``tesis-ucaldas/latex/references.bib``, ``overleaf/referencias.bib``,
``docs/borradores/referencias.bib`` y todos los ``docs/**/*.md``. Para cada DOI
consulta ``https://api.crossref.org/works/<doi>`` y compara título, primer autor
y año con lo que dice el fichero. Escribe ``results/v3/auditoria_citas.tsv``.

Sale con código 1 si un DOI no resuelve o si el título, el primer autor o el
año no coinciden. Un DOI que el texto no acompaña de cita local y que sí
resuelve no es un error: no hay nada local que comparar.

La caché de respuestas queda en ``results/v3/crossref_cache/``. No se inventa
un registro si la API falla.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from oncogo_crc.v3.citas import (
    archivos_a_revisar,
    cita_local_markdown,
    comparar,
    declarado_como_ausente,
    dois_en_texto,
    entradas_bib,
)

REPO = Path(__file__).resolve().parents[2]
WORKSPACE = Path("/Volumes/JANCODATA/OncoGO-CRC")
UA = "OncoGO-CRC/3 (https://github.com/JancoDeveloper/OncoGO-CRC)"
COLUMNAS = (
    "fuente",
    "clave",
    "doi",
    "estado",
    "titulo_local",
    "titulo_crossref",
    "autor_local",
    "autor_crossref",
    "anio_local",
    "anio_crossref",
    "volumen",
    "numero",
    "pagina",
    "revista",
    "tipo",
    "http",
)


def _anio(message: dict) -> str:
    issued = (message.get("issued") or {}).get("date-parts") or []
    if issued and issued[0]:
        return str(issued[0][0])
    return ""


def _primer_autor(message: dict) -> str:
    autores = message.get("author") or []
    if not autores:
        return ""
    primero = autores[0]
    return str(primero.get("family") or primero.get("name") or "")


def _titulo(message: dict) -> str:
    titulos = message.get("title") or []
    return str(titulos[0]) if titulos else ""


def fetch(doi: str, cache: Path, *, pausa_s: float) -> tuple[dict | None, int]:
    """Devuelve el mensaje de Crossref, o None si el DOI no resuelve."""
    destino = cache / (doi.replace("/", "_") + ".json")
    if destino.is_file():
        guardado = json.loads(destino.read_text(encoding="utf-8"))
        return guardado.get("message"), int(guardado.get("http", 0))
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    pedido = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    codigo = 0
    cuerpo: dict | None = None
    for intento in range(4):
        try:
            with urllib.request.urlopen(pedido, timeout=30) as resp:
                codigo = int(resp.status)
                cuerpo = json.loads(resp.read().decode("utf-8"))
            break
        except urllib.error.HTTPError as exc:
            codigo = int(exc.code)
            if codigo == 404:
                cuerpo = None
                break
            time.sleep(1.5 * (intento + 1))
        except (urllib.error.URLError, TimeoutError):
            codigo = 0
            time.sleep(1.5 * (intento + 1))
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_text(
        json.dumps({"http": codigo, "message": (cuerpo or {}).get("message")}, ensure_ascii=False),
        encoding="utf-8",
    )
    time.sleep(pausa_s)
    if not cuerpo:
        return None, codigo
    return cuerpo.get("message"), codigo


def filas_de(path: Path, repo: Path) -> list[dict[str, str]]:
    texto = path.read_text(encoding="utf-8", errors="replace")
    rel = str(path.relative_to(repo))
    if path.suffix == ".bib":
        return [
            {
                "fuente": rel,
                "clave": e["clave"],
                "doi": e["doi"],
                "titulo": e["titulo"],
                "autor": e["autor"],
                "anio": e["anio"],
                "ausente": "",
            }
            for e in entradas_bib(texto)
        ]
    filas = []
    for doi in dois_en_texto(texto):
        local = cita_local_markdown(texto, doi)
        filas.append(
            {
                "fuente": rel,
                "clave": "",
                "doi": doi,
                "ausente": "1" if declarado_como_ausente(texto, doi) else "",
                **local,
            }
        )
    return filas


ERRORES = {
    "doi_no_resuelve",
    "titulo_no_coincide",
    "autor_no_coincide",
    "anio_no_coincide",
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", type=Path, default=WORKSPACE)
    ap.add_argument("--pausa", type=float, default=0.15)
    args = ap.parse_args(argv)
    if not args.workspace.is_dir():
        print(f"workspace no accesible: {args.workspace}", file=sys.stderr)
        return 2

    salida = args.workspace / "results/v3"
    cache = salida / "crossref_cache"
    salida.mkdir(parents=True, exist_ok=True)

    pendientes: list[dict[str, str]] = []
    for path in archivos_a_revisar(REPO):
        pendientes.extend(filas_de(path, REPO))

    unicos = sorted({f["doi"].lower(): f["doi"] for f in pendientes}.values())
    registros: dict[str, tuple[dict | None, int]] = {}
    print(f"[citas] {len(pendientes)} citas, {len(unicos)} DOI distintos.")
    for i, doi in enumerate(unicos, start=1):
        registros[doi.lower()] = fetch(doi, cache, pausa_s=args.pausa)
        if i % 25 == 0 or i == len(unicos):
            print(f"[citas] {i}/{len(unicos)}")

    lineas = ["\t".join(COLUMNAS)]
    errores = 0
    for fila in pendientes:
        mensaje, codigo = registros[fila["doi"].lower()]
        if mensaje:
            titulo_r = _titulo(mensaje)
            autor_r = _primer_autor(mensaje)
            anio_r = _anio(mensaje)
            volumen = str(mensaje.get("volume") or "")
            numero = str(mensaje.get("issue") or mensaje.get("article-number") or "")
            pagina = str(mensaje.get("page") or mensaje.get("article-number") or "")
            revista = ""
            contenedor = mensaje.get("container-title") or []
            if contenedor:
                revista = str(contenedor[0])
            tipo = str(mensaje.get("type") or "")
            resuelve = True
        else:
            titulo_r = autor_r = anio_r = volumen = numero = pagina = revista = tipo = ""
            resuelve = False
        estado = comparar(
            titulo_local=fila.get("titulo", ""),
            autor_local=fila.get("autor", ""),
            anio_local=fila.get("anio", ""),
            titulo_remoto=titulo_r,
            autor_remoto=autor_r,
            anio_remoto=anio_r,
            resuelve=resuelve,
        )
        if fila.get("ausente") and not resuelve:
            estado = "declarado_no_resuelve"
        elif fila["doi"].lower().startswith("10.7303/"):
            estado = "identificador_synapse"
        if any(parte in ERRORES for parte in estado.split("+")):
            errores += 1
        valores = [
            fila["fuente"],
            fila.get("clave", ""),
            fila["doi"],
            estado,
            fila.get("titulo", ""),
            titulo_r,
            fila.get("autor", ""),
            autor_r,
            fila.get("anio", ""),
            anio_r,
            volumen,
            numero,
            pagina,
            revista,
            tipo,
            str(codigo),
        ]
        lineas.append("\t".join(v.replace("\t", " ").replace("\n", " ") for v in valores))

    destino = salida / "auditoria_citas.tsv"
    destino.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    print(f"[citas] {destino} errores={errores}")
    return 1 if errores else 0


if __name__ == "__main__":
    raise SystemExit(main())
