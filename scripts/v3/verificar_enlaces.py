"""Comprueba los enlaces externos. Escribe en docs/, no en results/v3_7/."""

from __future__ import annotations

import re
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ENTRADA = REPO / "config" / "enlaces_externos.tsv"
SALIDA = REPO / "docs" / "web_enlaces.tsv"
TITULO = re.compile(r"<title[^>]*>(.*?)</title>", re.I | re.S)


def leer() -> list[dict]:
    lineas = ENTRADA.read_text(encoding="utf-8").splitlines()
    cabecera = lineas[0].split("\t")
    filas = []
    for linea in lineas[1:]:
        if not linea.strip():
            continue
        partes = linea.split("\t")
        filas.append(dict(zip(cabecera, partes)))
    return filas


def pedir(url: str) -> tuple[int, str, str]:
    peticion = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.0 Safari/605.1.15"},
    )
    try:
        with urllib.request.urlopen(peticion, timeout=30) as respuesta:
            codigo = respuesta.getcode()
            final = respuesta.geturl()
            cuerpo = respuesta.read(200_000).decode("utf-8", errors="replace")
    except urllib.error.HTTPError as error:
        return error.code, str(error.url or url), ""
    except Exception as error:
        return 0, url, str(error)
    hallado = TITULO.search(cuerpo)
    titulo = re.sub(r"\s+", " ", hallado.group(1)).strip() if hallado else ""
    return codigo, final, titulo


def pedir_con_reintento(url: str) -> tuple[int, str, str]:
    codigo, final, titulo = pedir(url)
    if codigo == 0:
        codigo, final, titulo = pedir(url)
    return codigo, final, titulo


def main() -> int:
    filas = leer()
    fecha = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    salida = ["id\turl\tcodigo_final\ttitulo\turl_final\tfecha_verificacion\testado_http"]
    fallos = []
    for fila in filas:
        codigo, final, titulo = pedir_con_reintento(fila["url"])
        esperado = fila.get("titulo_esperado", "")
        texto = f"{titulo} {final}"
        bien = codigo == 200 and esperado.lower() in texto.lower()
        estado = "200" if bien else str(codigo or "fallo")
        fila["fecha_verificacion"] = fecha
        fila["estado_http"] = "200" if bien else "no_200"
        if not bien:
            fallos.append(f"{fila['id']}: HTTP {codigo}; titulo={titulo[:80]!r}")
        salida.append("\t".join([
            fila["id"], fila["url"], str(codigo), titulo.replace("\t", " "),
            final, fecha, "200" if bien else "no_200",
        ]))
        print(fila["id"], codigo, titulo[:70])
    SALIDA.write_text("\n".join(salida) + "\n", encoding="utf-8")
    cabecera = ENTRADA.read_text(encoding="utf-8").splitlines()[0].split("\t")
    nuevas = ["\t".join(cabecera)]
    for fila in filas:
        nuevas.append("\t".join(fila.get(col, "") for col in cabecera))
    ENTRADA.write_text("\n".join(nuevas) + "\n", encoding="utf-8")
    if fallos:
        print("FALLOS")
        print("\n".join(fallos))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
