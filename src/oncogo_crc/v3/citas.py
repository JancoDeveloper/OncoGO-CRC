"""Comparación de citas locales con un registro bibliográfico.

No consulta la red. La descarga de Crossref vive en
``scripts/v3/verificar_citas.py``. Aquí solo se extraen DOI y se decide si
título, primer autor y año coinciden.
"""

from __future__ import annotations

import re
import unicodedata
from pathlib import Path

DOI_RE = re.compile(
    r"(?i)(?:doi[:\s]*|https?://(?:dx\.)?doi\.org/)(10\.\d{4,9}/[^\s\]\)>\"'<,;]+)"
)
BIB_HEAD = re.compile(r"@\w+\s*\{\s*([^,\s]+)\s*,", re.IGNORECASE)


def normalizar_texto(texto: str) -> str:
    """Minúsculas, sin acentos ni puntuación, espacios simples."""
    base = unicodedata.normalize("NFKD", texto or "")
    base = "".join(c for c in base if not unicodedata.combining(c))
    base = base.lower().replace("&", " and ")
    base = re.sub(r"[^a-z0-9]+", " ", base)
    return " ".join(base.split())


def tokens(texto: str) -> list[str]:
    return normalizar_texto(texto).split()


def titulos_coinciden(local: str, remoto: str) -> bool:
    """Coinciden si, ya normalizados, son iguales o uno es prefijo del otro.

    Un prefijo cubre el caso de un subtítulo omitido. Una palabra sustituida
    (por ejemplo «proteograph» donde el registro dice otra cosa) no es prefijo
    y no coincide.
    """
    a = tokens(local)
    b = tokens(remoto)
    if not a or not b:
        return False
    if a == b:
        return True
    if len(a) >= 6 and b[: len(a)] == a:
        return True
    if len(b) >= 6 and a[: len(b)] == b:
        return True
    return False


def apellido(autor: str) -> str:
    """Apellido del primer autor de una cita o de un campo BibTeX.

    «Wannberg F» y «Wannberg, Frida and ...» devuelven «wannberg».
    Una etiqueta «Cita:» delante no cuenta como apellido.
    """
    texto = (autor or "").strip()
    if "|" in texto:
        texto = texto.split("|")[-1]
    texto = re.sub(r"(?i)^(?:[-*]\s*)?(?:\*\*)?cita(?:\*\*)?\s*:\s*", "", texto)
    texto = re.split(r"\band\b|\bet al\b|&", texto, maxsplit=1, flags=re.IGNORECASE)[0]
    texto = texto.strip(" ,.;:")
    if "," in texto:
        familia = texto.split(",", 1)[0]
    else:
        familia = texto
    partes = [p for p in re.split(r"\s+", familia) if p]
    largos = [p for p in partes if len(p.strip(".")) > 2 and not any(ch.isdigit() for ch in p)]
    elegido = largos[0] if largos else (partes[-1] if partes else "")
    return normalizar_texto(elegido)


def autores_coinciden(local: str, remoto: str) -> bool:
    """El apellido basta. Un autor colectivo coincide si el registro lo nombra."""
    ln = normalizar_texto(local)
    rn = normalizar_texto(remoto)
    if not ln or not rn:
        return False
    if ln.startswith(rn) or rn.startswith(ln):
        return True
    # Un autor colectivo («the NCI CPTAC») coincide si la cita lo nombra.
    # No se hace con apellidos cortos: «li» aparece dentro de otras palabras.
    if " " in rn and rn in ln:
        return True
    a = apellido(local)
    b = apellido(remoto)
    return bool(a) and a == b


def limpiar_doi(doi: str) -> str:
    return doi.strip().rstrip(".,;:)")


def dois_en_texto(texto: str) -> list[str]:
    vistos: list[str] = []
    ya: set[str] = set()
    for m in DOI_RE.finditer(texto):
        doi = limpiar_doi(m.group(1))
        clave = doi.lower()
        if clave not in ya:
            ya.add(clave)
            vistos.append(doi)
    return vistos


def _campo_bib(entrada: str, nombre: str) -> str:
    """Valor de un campo BibTeX con llaves balanceadas."""
    m = re.search(rf"(?im)^\s*{re.escape(nombre)}\s*=\s*\{{", entrada)
    if not m:
        m = re.search(rf"(?im)\b{re.escape(nombre)}\s*=\s*\{{", entrada)
    if not m:
        return ""
    i = m.end()
    profundidad = 1
    buf: list[str] = []
    while i < len(entrada) and profundidad:
        c = entrada[i]
        if c == "{":
            profundidad += 1
            buf.append(c)
        elif c == "}":
            profundidad -= 1
            if profundidad:
                buf.append(c)
        else:
            buf.append(c)
        i += 1
    return "".join(buf).strip()


def entradas_bib(texto: str) -> list[dict[str, str]]:
    """Una ficha por entrada con DOI. Sin DOI no hay nada que verificar."""
    partes = re.split(r"\n(?=@)", texto)
    salida: list[dict[str, str]] = []
    for parte in partes:
        if "@" not in parte:
            continue
        doi = limpiar_doi(_campo_bib(parte, "doi"))
        if not doi.startswith("10."):
            continue
        clave = BIB_HEAD.search(parte)
        salida.append(
            {
                "clave": clave.group(1) if clave else "",
                "doi": doi,
                "titulo": _campo_bib(parte, "title"),
                "autor": _campo_bib(parte, "author"),
                "anio": _campo_bib(parte, "year"),
            }
        )
    return salida


_REVISTAS = (
    "Nature|Cell|Science|Journal|Genome|PLOS|Nucleic|EMBO|Bioinformatics|"
    "Radiology|Biometrics|Molecular|Cancer|Proteomics"
)
_TITULO_SUCIO = re.compile(r"doi|http|estado|script|portal|www\.|\d{4}", re.IGNORECASE)


def _titulo_limpio(candidato: str) -> str:
    """Un título de una cita, no el párrafo entero que lo rodea."""
    titulo = " ".join(candidato.split()).strip(" .")
    palabras = tokens(titulo)
    if not (6 <= len(palabras) <= 40):
        return ""
    if _TITULO_SUCIO.search(titulo):
        return ""
    return titulo


def _bloque_del_doi(texto: str, pos: int) -> str:
    """La viñeta o el párrafo que contiene el DOI, sin la viñeta anterior."""
    inicio = 0
    cursor = pos
    while cursor > 0:
        salto = texto.rfind("\n", 0, cursor)
        if salto < 0:
            inicio = 0
            break
        linea = texto[salto + 1 : cursor]
        if linea.startswith(("- ", "* ", "#", "|")) or not linea.strip():
            inicio = salto + 1
            break
        cursor = salto
    fin = len(texto)
    cursor = texto.find("\n", pos)
    while 0 <= cursor < len(texto):
        siguiente = texto.find("\n", cursor + 1)
        linea = texto[cursor + 1 : siguiente if siguiente >= 0 else len(texto)]
        if linea.startswith(("- ", "* ", "#", "|")) or not linea.strip():
            fin = cursor
            break
        if siguiente < 0:
            break
        cursor = siguiente
    return texto[inicio:fin]


def cita_local_markdown(texto: str, doi: str) -> dict[str, str]:
    """Saca título, autor y año de la viñeta que contiene el DOI.

    No mira la viñeta anterior: si no, el autor de esa cita se le atribuye a
    esta.
    """
    patron = re.compile(re.escape(doi), re.IGNORECASE)
    m = patron.search(texto)
    if not m:
        return {"titulo": "", "autor": "", "anio": ""}
    bloque = _bloque_del_doi(texto, m.start())
    if len(dois_en_texto(bloque)) > 1:
        return {"titulo": "", "autor": "", "anio": ""}
    plana = re.sub(r"[*_`>#]", " ", bloque)
    plana = re.sub(r"\s+", " ", plana).strip()
    anio = ""
    anios = list(re.finditer(r"\b(?:19|20)\d{2}\b", plana))
    if anios:
        anio = anios[-1].group(0)
    autor = ""
    titulo = ""
    cortes = list(re.finditer(r"et al\.?,?\s+", plana, flags=re.IGNORECASE))
    if cortes:
        corte = cortes[-1]
        autor = plana[: corte.start()].split(".")[-1].strip(" -—:;")
        resto = plana[corte.end() :]
        tm = re.match(rf"(.+?)\.\s+(?:{_REVISTAS})", resto)
        if tm:
            titulo = _titulo_limpio(tm.group(1))
    return {"titulo": titulo, "autor": autor, "anio": anio}


def declarado_como_ausente(texto: str, doi: str) -> bool:
    """El propio texto dice que ese DOI no resolvió y no se usa como cita."""
    m = re.search(re.escape(doi), texto, flags=re.IGNORECASE)
    if not m:
        return False
    ventana = texto[max(0, m.start() - 180) : m.end() + 80].lower()
    return "404" in ventana and ("no se enlaza" in ventana or "no entra" in ventana)


def comparar(
    *,
    titulo_local: str,
    autor_local: str,
    anio_local: str,
    titulo_remoto: str,
    autor_remoto: str,
    anio_remoto: str,
    resuelve: bool,
) -> str:
    """Estado de una fila. El primero que falle decide el código de salida."""
    if not resuelve:
        return "doi_no_resuelve"
    fallos: list[str] = []
    if titulo_local and titulo_remoto and not titulos_coinciden(titulo_local, titulo_remoto):
        fallos.append("titulo_no_coincide")
    if autor_local and autor_remoto and not autores_coinciden(autor_local, autor_remoto):
        fallos.append("autor_no_coincide")
    if anio_local and anio_remoto and str(anio_local).strip() != str(anio_remoto).strip():
        fallos.append("anio_no_coincide")
    if fallos:
        return "+".join(fallos)
    if not titulo_local and not autor_local and not anio_local:
        return "ok_sin_cita_local"
    return "ok"


def archivos_a_revisar(repo: Path) -> list[Path]:
    fijos = [
        repo / "tesis-ucaldas/latex/references.bib",
        repo / "overleaf/referencias.bib",
        repo / "docs/borradores/referencias.bib",
    ]
    md = sorted((repo / "docs").rglob("*.md")) if (repo / "docs").is_dir() else []
    return [p for p in fijos + md if p.is_file()]
