#!/usr/bin/env python3
"""Verifica que todo número de los informes de v2 salga de una tabla (C14).

El problema que resuelve. En v1 los números del PDF y del Markdown se escribían
a mano. Cuando el pipeline cambiaba, el texto se quedaba con el valor viejo y no
había forma de detectarlo salvo leerlo todo. Así es como R4 llegó a figurar
``calculado`` y a la vez con una definición que decía «la pureza tumoral aún no
es covariable».

Cómo funciona. Se recogen todos los valores numéricos de los JSON y los TSV de
``results/v2/``. Después se extraen los números de los informes y se comprueba
que cada uno aparezca en ese conjunto, con la tolerancia del redondeo con el que
se imprimió. Si alguno no aparece, el script falla y dice cuál.

Qué no cuenta como número a verificar. Lo que no es un resultado: el 95 de
«intervalo del 95 %», los años de las citas, los identificadores de proyecto y
los umbrales que vienen del preregistro, que también es una fuente legítima.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import yaml

# Números que son estructura del texto, no resultados.
PERMITIDOS_SIEMPRE = {
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 50.0, 95.0, 100.0,
    # Años de citas y referencias.
    1995.0, 2002.0, 2004.0, 2009.0, 2012.0, 2015.0, 2019.0, 2020.0,
    2022.0, 2023.0, 2024.0, 2025.0, 2026.0,
}

# Un número dentro de un identificador no es una medida.
PATRONES_IGNORADOS = (
    re.compile(r"\bPXD\d+\b"),
    re.compile(r"\bPDC\d+\b"),
    re.compile(r"\bGO:\d+\b"),
    re.compile(r"\bR-HSA-\d+\b"),
    re.compile(r"\bR[1-5]_\w+\b"),
    re.compile(r"\bH[1-4]\b"),
    re.compile(r"\bC\d{1,2}\b"),
    re.compile(r"\bE\d\b"),
    re.compile(r"\blog2\w*\b", re.IGNORECASE),
    re.compile(r"\bI²\b"),
    re.compile(r"`[^`]*`"),          # código en línea
    re.compile(r"\[[^\]]*\]\([^)]*\)"),  # enlaces
    re.compile(r"^\s{0,3}#{1,6}\s.*$", re.MULTILINE),  # títulos
    # Primera celda de una fila de tabla cuando es un entero suelto: es el orden
    # de los pasos de un procedimiento, no una medida.
    re.compile(r"^\|\s*\d+\s*\|", re.MULTILINE),
)

# El signo menos de los informes es U+2212, no un guion. Sin capturarlo, un
# −0,228 se leía como +0,228 y no coincidía con su fuente: todo número negativo
# del texto parecía escrito a mano. El guion normal se deja fuera a propósito,
# porque en «0,492–0,708» separa los extremos de un intervalo y no niega nada.
NUMERO = re.compile(r"(?<![\w.,\u2212])(\u2212?\d+(?:[.,]\d+)?)(?![\w])")


def valores_de_json(obj, salida: set[float]) -> None:
    if isinstance(obj, bool):
        return
    if isinstance(obj, (int, float)):
        salida.add(float(obj))
    elif isinstance(obj, dict):
        for v in obj.values():
            valores_de_json(v, salida)
    elif isinstance(obj, list):
        for v in obj:
            valores_de_json(v, salida)
    elif isinstance(obj, str):
        # Un número dentro de una cadena de un JSON también es una fuente: los
        # motivos y las notas los escribe el código a partir de los datos.
        for m in NUMERO.finditer(obj):
            try:
                salida.add(float(m.group(1).replace(",", ".")))
            except ValueError:
                pass


def valores_de_tsv(path: Path, salida: set[float]) -> None:
    n_filas = 0
    with path.open(encoding="utf-8", errors="replace") as fh:
        for i, linea in enumerate(fh):
            if i == 0:
                continue
            n_filas += 1
            for campo in linea.rstrip("\n").split("\t"):
                try:
                    salida.add(float(campo))
                except ValueError:
                    # Un número dentro de un campo de texto también es fuente: los
                    # motivos del changelog los escribe el código desde los datos.
                    for m in NUMERO.finditer(campo):
                        try:
                            salida.add(float(m.group(1).replace(",", ".")))
                        except ValueError:
                            continue
    # Cuántas filas tiene la tabla es un dato de la tabla. El informe dice «de los
    # 11 números del changelog», y ese 11 es exactamente este recuento.
    salida.add(float(n_filas))


# Solo las tablas de resumen cuentan como fuente. Las tablas por proteína traen
# cientos de miles de valores, y con tantos cualquier número del texto coincidiría
# con alguno por azar: el verificador dejaría de verificar nada.
TSV_DE_RESUMEN = (
    "indice_R_v2.tsv",
    "changelog_numeros.tsv",
    "tamano_necesario_cohorte_cea.tsv",
)


def recoger_fuentes(
    raiz: Path, extra: list[Path], fichas: Path | None = None
) -> tuple[set[float], list[str]]:
    valores: set[float] = set()
    usadas: list[str] = []
    rutas = sorted(raiz.rglob("*.json"))
    if fichas is not None and fichas.exists():
        # Las fichas de las fuentes externas también son artefactos calculados: las
        # escriben los scripts de descarga con su hash y su fecha de acceso. Un
        # número citado de una ficha no está escrito a mano, pero vive fuera de
        # `results/` y sin esto el verificador no podía verlo.
        rutas += sorted(fichas.glob("*/estado.json"))
    for path in rutas:
        if path.name.startswith("._"):
            continue
        try:
            valores_de_json(json.loads(path.read_text(encoding="utf-8")), valores)
            usadas.append(path.name if raiz not in path.parents else str(path.relative_to(raiz)))
        except (json.JSONDecodeError, UnicodeDecodeError):
            print(f"[verificar] aviso: {path} no es JSON válido", file=sys.stderr)
    for nombre in TSV_DE_RESUMEN:
        for path in raiz.rglob(nombre):
            if path.name.startswith("._"):
                continue
            valores_de_tsv(path, valores)
            usadas.append(str(path.relative_to(raiz)))
    for path in extra:
        if path.exists():
            valores_de_json(yaml.safe_load(path.read_text(encoding="utf-8")), valores)
            usadas.append(path.name)
    return valores, usadas


def limpiar(texto: str) -> str:
    for patron in PATRONES_IGNORADOS:
        texto = patron.sub(" ", texto)
    return texto


def coincide(valor: float, fuentes: set[float], decimales: int) -> bool:
    if valor in PERMITIDOS_SIEMPRE:
        return True
    # La tolerancia es media unidad del último decimal impreso.
    tol = 0.5 * 10.0 ** (-decimales) if decimales else 0.5
    for f in fuentes:
        if abs(f - valor) <= tol:
            return True
        # Las fracciones se imprimen a veces como porcentaje.
        if abs(f * 100.0 - valor) <= tol:
            return True
    return False


def verificar(path: Path, fuentes: set[float]) -> list[dict]:
    texto = limpiar(path.read_text(encoding="utf-8"))
    fallos = []
    for m in NUMERO.finditer(texto):
        bruto = m.group(1)
        decimales = len(bruto.split(",")[1]) if "," in bruto else (
            len(bruto.split(".")[1]) if "." in bruto else 0
        )
        try:
            valor = float(bruto.replace(",", ".").replace("\u2212", "-"))
        except ValueError:
            continue
        if not coincide(valor, fuentes, decimales):
            inicio = max(0, m.start() - 60)
            fallos.append(
                {
                    "fichero": path.name,
                    "numero": bruto,
                    "contexto": " ".join(texto[inicio : m.end() + 40].split()),
                }
            )
    return fallos


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    ap.add_argument("--repo", default=str(Path(__file__).resolve().parents[2]))
    ap.add_argument("informes", nargs="*", help="Markdown a verificar. Por defecto, los de v2.")
    args = ap.parse_args()

    raiz = Path(args.workspace) / "results" / "v2"
    if not raiz.exists():
        print(f"[verificar] no existe {raiz}", file=sys.stderr)
        return 1
    repo = Path(args.repo)
    fuentes, usadas = recoger_fuentes(
        raiz,
        [
            repo / "config" / "preregistro_v2.yaml",
            repo / "config" / "registro_subindices.yaml",
        ],
        fichas=Path(args.workspace) / "external",
    )
    print(
        f"[verificar] {len(fuentes)} valores numéricos distintos "
        f"en {len(usadas)} tablas de resumen."
    )

    informes = [Path(p) for p in args.informes]
    if not informes:
        # El informe de brecha se copia a report/ para el despliegue, así que
        # aparece dos veces. Verificarlo una basta.
        vistos: set[str] = set()
        for p in sorted(raiz.rglob("*.md")):
            if p.name.startswith("._") or p.name in vistos:
                continue
            vistos.add(p.name)
            informes.append(p)
    if not informes:
        print("[verificar] no hay informes que verificar", file=sys.stderr)
        return 1

    todos: list[dict] = []
    for informe in informes:
        fallos = verificar(informe, fuentes)
        estado = "sin números huérfanos" if not fallos else f"{len(fallos)} sin fuente"
        print(f"[verificar] {informe.name}: {estado}")
        todos.extend(fallos)

    if todos:
        print("\nNúmeros que no salen de ninguna tabla:", file=sys.stderr)
        for f in todos[:40]:
            print(f"  {f['fichero']}: {f['numero']}  …{f['contexto']}…", file=sys.stderr)
        return 1
    print("[verificar] todos los números de los informes salen de una tabla.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
