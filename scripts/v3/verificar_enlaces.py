"""Comprueba los enlaces externos. Escribe en docs/, no en results/v3_7/.

Un DOI se comprueba en la agencia de registro y en Crossref.
Europe PMC se comprueba en su API REST.
Una página se comprueba con un GET. No se pide la página del editor
para un DOI ni para Europe PMC.
"""

from __future__ import annotations

import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ENTRADA = REPO / "config" / "enlaces_externos.tsv"
SALIDA = REPO / "docs" / "web_enlaces.tsv"
TITULO_HTML = re.compile(r"<title[^>]*>(.*?)</title>", re.I | re.S)
AGENTE = "OncoGO-CRC/3.7 (mailto:jancod@jancod.co)"
EDITORES = ("pubs.acs.org", "linkinghub.elsevier.com", "europepmc.org", "www.europepmc.org")


def leer() -> list[dict]:
    lineas = ENTRADA.read_text(encoding="utf-8").splitlines()
    cabecera = lineas[0].split("\t")
    filas = []
    for linea in lineas[1:]:
        if not linea.strip():
            continue
        partes = linea.split("\t")
        filas.append(dict(zip(cabecera, partes + [""] * (len(cabecera) - len(partes)))))
    return filas


def titulo_comparable(texto: str) -> str:
    limpio = " ".join(texto.split())
    if limpio.endswith("."):
        limpio = limpio[:-1].rstrip()
    return limpio.casefold()


def obtener_red(url: str) -> tuple[int, str, str]:
    peticion = urllib.request.Request(url, headers={"User-Agent": AGENTE})
    try:
        with urllib.request.urlopen(peticion, timeout=45) as respuesta:
            codigo = respuesta.getcode()
            final = respuesta.geturl()
            cuerpo = respuesta.read(400_000).decode("utf-8", errors="replace")
    except urllib.error.HTTPError as error:
        cuerpo = error.read(400_000).decode("utf-8", errors="replace")
        return error.code, str(error.url or url), cuerpo
    except Exception as error:
        return 0, url, str(error)
    return codigo, final, cuerpo


def pedir(url: str, obtener) -> tuple[int, str, str]:
    codigo, final, cuerpo = obtener(url)
    if codigo == 0:
        codigo, final, cuerpo = obtener(url)
    return codigo, final, cuerpo


def _json(cuerpo: str) -> dict:
    try:
        dato = json.loads(cuerpo)
    except json.JSONDecodeError:
        return {}
    return dato if isinstance(dato, dict) else {}


def _handle_tiene_url(dato: dict) -> str:
    if dato.get("responseCode") != 1:
        return ""
    for valor in dato.get("values") or []:
        if valor.get("type") == "URL":
            destino = str((valor.get("data") or {}).get("value") or "")
            if destino:
                return destino
    return ""


def _es_editor(url: str) -> bool:
    return any(host in url for host in EDITORES)


def comprobar_fila(fila: dict, obtener=obtener_red) -> dict:
    metodo = fila["metodo"]
    pedidas: list[str] = []
    codigo = 0
    titulo = ""
    url_final = fila["url"]
    bien = False

    def una(url: str) -> tuple[int, str, str]:
        pedidas.append(url)
        return pedir(url, obtener)

    if metodo == "doi":
        doi = fila["url"].removeprefix("https://doi.org/").strip()
        codigo_h, final_h, cuerpo_h = una(f"https://doi.org/api/handles/{doi}")
        codigo_c, final_c, cuerpo_c = una(f"https://api.crossref.org/works/{doi}")
        destino = _handle_tiene_url(_json(cuerpo_h))
        titulos = (_json(cuerpo_c).get("message") or {}).get("title") or []
        titulo = str(titulos[0]) if titulos else ""
        codigo = codigo_c if codigo_h == 200 else codigo_h
        url_final = destino or final_c
        bien = (
            codigo_h == 200
            and codigo_c == 200
            and bool(destino)
            and titulo_comparable(titulo) == titulo_comparable(fila.get("titulo_esperado", ""))
            and not _es_editor(final_h)
            and not _es_editor(final_c)
        )
    elif metodo == "europepmc":
        med = fila["url"].split("/MED/", 1)[1].split("?", 1)[0]
        consulta = (
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
            f"?query=EXT_ID:{med}%20AND%20SRC:MED&format=json"
        )
        codigo, final, cuerpo = una(consulta)
        dato = _json(cuerpo)
        resultados = ((dato.get("resultList") or {}).get("result")) or []
        primero = resultados[0] if resultados else {}
        titulo = str(primero.get("title") or "")
        doi = str(primero.get("doi") or "")
        url_final = final
        bien = (
            codigo == 200
            and int(dato.get("hitCount") or 0) >= 1
            and titulo_comparable(titulo) == titulo_comparable(fila.get("titulo_esperado", ""))
            and doi.casefold() == fila.get("doi_esperado", "").casefold()
            and not _es_editor(final)
        )
    elif metodo == "pagina":
        codigo, final, cuerpo = una(fila["url"])
        hallado = TITULO_HTML.search(cuerpo)
        titulo = re.sub(r"\s+", " ", hallado.group(1)).strip() if hallado else ""
        url_final = final
        esperado = fila.get("titulo_esperado", "")
        bien = codigo == 200 and esperado.lower() in f"{titulo} {final}".lower()
    elif metodo == "doi_datacite":
        doi = fila["url"].removeprefix("https://doi.org/").strip()
        codigo, final, cuerpo = una(
            "https://api.datacite.org/dois/" + urllib.parse.quote(doi, safe="/")
        )
        titulos = (((_json(cuerpo).get("data") or {}).get("attributes") or {}).get("titles")) or []
        primero = titulos[0] if titulos else {}
        titulo = str(primero.get("title") or "") if isinstance(primero, dict) else ""
        url_final = final
        bien = (
            codigo == 200
            and titulo_comparable(titulo) == titulo_comparable(fila.get("titulo_esperado", ""))
            and not _es_editor(final)
        )
    else:
        titulo = f"metodo desconocido: {metodo}"

    return {
        "id": fila["id"],
        "url": fila["url"],
        "metodo": metodo,
        "codigo_final": str(codigo),
        "titulo": titulo.replace("\t", " "),
        "url_final": url_final,
        "estado_http": "200" if bien else "no_200",
        "pedidas": pedidas,
    }


def main() -> int:
    filas = leer()
    fecha = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    salida = ["id\turl\tmetodo\tcodigo_final\ttitulo\turl_final\tfecha_verificacion\testado_http"]
    fallos = []
    for fila in filas:
        resultado = comprobar_fila(fila)
        fila["fecha_verificacion"] = fecha
        fila["estado_http"] = resultado["estado_http"]
        if resultado["estado_http"] != "200":
            fallos.append(
                f"{fila['id']}: HTTP {resultado['codigo_final']}; titulo={resultado['titulo'][:120]!r}"
            )
        salida.append("\t".join([
            resultado["id"], resultado["url"], resultado["metodo"], resultado["codigo_final"],
            resultado["titulo"], resultado["url_final"], fecha, resultado["estado_http"],
        ]))
        print(fila["id"], resultado["metodo"], resultado["estado_http"], " ".join(resultado["pedidas"]))
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
