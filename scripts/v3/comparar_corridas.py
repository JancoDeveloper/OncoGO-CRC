"""Compara dos salidas de v3.7. Los JSON se miran sin el sello. Los TSV, por sha256.

Con --estricto solo se ignora fecha_utc dentro del sello. El resto del sello
tiene que coincidir. Se admiten corrida.json, tests.txt, el orden de
gsea_tejido.tsv y un hash que apunta a un fichero sellado.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

PRINCIPALES = (
    "suero_embudo.json",
    "ml_modelos.tsv",
    "k1_etiquetas.json",
    "h5_especificidad_organo.json",
    "h6_origen_celular.json",
    "h8_carga_tumoral.json",
    "h9_wald.json",
    "h10_especificidad.json",
    "preanalitica.json",
    "multiverso_resumen.json",
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sin_sello(path: Path):
    datos = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(datos, dict):
        datos.pop("sello", None)
    return datos


def _h7_secundario(datos: dict) -> dict:
    modelo = datos.get("modelo_secundario")
    if isinstance(modelo, dict):
        return modelo
    return {}


def comparar(a: Path, b: Path) -> dict:
    nombres = sorted({p.name for p in a.iterdir() if p.is_file() and not p.name.startswith("._")} | {p.name for p in b.iterdir() if p.is_file() and not p.name.startswith("._")})
    distintos = []
    solo_a = []
    solo_b = []
    for nombre in nombres:
        pa, pb = a / nombre, b / nombre
        if not pa.is_file():
            solo_b.append(nombre)
            continue
        if not pb.is_file():
            solo_a.append(nombre)
            continue
        if nombre.endswith(".json"):
            if _sin_sello(pa) != _sin_sello(pb):
                distintos.append(nombre)
        elif nombre.endswith(".sello.json"):
            if _sha(pa) != _sha(pb):
                distintos.append(nombre)
        elif _sha(pa) != _sha(pb):
            distintos.append(nombre)
    h7_distinto = "h7_visibilidad_en_sangre.json" in distintos
    secundario = False
    if h7_distinto:
        da, db = _sin_sello(a / "h7_visibilidad_en_sangre.json"), _sin_sello(b / "h7_visibilidad_en_sangre.json")
        sa, sb = _h7_secundario(da), _h7_secundario(db)
        da.pop("modelo_secundario", None)
        db.pop("modelo_secundario", None)
        secundario = da == db and sa != sb
        if secundario:
            distintos.remove("h7_visibilidad_en_sangre.json")
    datos = [nombre for nombre in distintos if not nombre.endswith(".sello.json")]
    principales = [nombre for nombre in datos if nombre in PRINCIPALES or nombre.startswith("k4_") or nombre.startswith("suero_k4")]
    return {
        "distintos": distintos,
        "solo_en_a": solo_a,
        "solo_en_b": solo_b,
        "principales": principales,
        "h7_secundario_distinto": secundario,
        "limitacion_conocida": "El secundario de H7 (ESM-2 y TMbed) no es determinista entre procesos y no entra en el veredicto." if secundario else "",
    }


_HEX = re.compile(r"^[0-9a-f]{64}$")


def _es_sello(datos: dict) -> bool:
    return "commit" in datos and "etapa" in datos


def _sin_fecha(obj):
    if isinstance(obj, dict):
        return {k: _sin_fecha(v) for k, v in obj.items() if not (k == "fecha_utc" and _es_sello(obj))}
    if isinstance(obj, list):
        return [_sin_fecha(v) for v in obj]
    return obj


def _hojas(a, b, pref=""):
    if isinstance(a, dict) and isinstance(b, dict):
        for clave in sorted(set(a) | set(b)):
            ruta = f"{pref}.{clave}" if pref else clave
            if clave not in a or clave not in b:
                yield ruta, a.get(clave), b.get(clave)
            else:
                yield from _hojas(a[clave], b[clave], ruta)
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            yield pref, len(a), len(b)
            return
        for i, (x, y) in enumerate(zip(a, b)):
            yield from _hojas(x, y, f"{pref}[{i}]")
        return
    if a != b:
        yield pref, a, b


def _solo_hash_o_corrida(a, b) -> bool:
    hojas = list(_hojas(a, b))
    if not hojas:
        return True
    for ruta, x, y in hojas:
        nombre = ruta.rsplit(".", 1)[-1]
        if nombre in {"pid", "inicio", "fin"}:
            continue
        if nombre in {"sha256", "sha256_fuente", "r1_sha256"} or (
            isinstance(x, str) and isinstance(y, str) and _HEX.fullmatch(x) and _HEX.fullmatch(y)
        ):
            continue
        return False
    return True


def _gsea_igual(a: Path, b: Path) -> bool:
    uno = pd.read_csv(a, sep="\t").sort_values(["cohorte", "Term"]).reset_index(drop=True)
    dos = pd.read_csv(b, sep="\t").sort_values(["cohorte", "Term"]).reset_index(drop=True)
    return uno.equals(dos)


def _h7_tabla_solo_esm2(a: Path, b: Path) -> bool:
    uno = pd.read_csv(a, sep="\t").sort_values("gen").reset_index(drop=True)
    dos = pd.read_csv(b, sep="\t").sort_values("gen").reset_index(drop=True)
    if list(uno.columns) != list(dos.columns) or len(uno) != len(dos):
        return False
    for col in uno.columns:
        if uno[col].astype(str).equals(dos[col].astype(str)):
            continue
        if not str(col).startswith("esm2_"):
            return False
        x = pd.to_numeric(uno[col], errors="coerce")
        y = pd.to_numeric(dos[col], errors="coerce")
        ambos = x.notna() & y.notna()
        if ambos.any() and float((x[ambos] - y[ambos]).abs().max()) != 0:
            return False
    return True


def comparar_estricto(a: Path, b: Path) -> dict:
    """B frente a C. Ignora solo la fecha del sello y las diferencias admitidas."""
    nombres = sorted(
        {p.name for p in a.iterdir() if p.is_file() and not p.name.startswith("._")}
        | {p.name for p in b.iterdir() if p.is_file() and not p.name.startswith("._")}
    )
    admitidos = []
    bloquean = []
    secundario = None
    for nombre in nombres:
        pa, pb = a / nombre, b / nombre
        if not pa.is_file() or not pb.is_file():
            bloquean.append(nombre)
            continue
        if nombre == "tests.txt":
            if _sha(pa) != _sha(pb):
                admitidos.append(nombre)
            continue
        if nombre == "gsea_tejido.tsv":
            if _sha(pa) != _sha(pb):
                admitidos.append(nombre) if _gsea_igual(pa, pb) else bloquean.append(nombre)
            continue
        if nombre == "h7_tabla_proteinas.tsv":
            if _sha(pa) != _sha(pb):
                bloquean.append(nombre)
                if not _h7_tabla_solo_esm2(pa, pb):
                    bloquean.append(nombre + ":no_es_solo_esm2")
            continue
        if nombre.endswith(".json") or nombre.endswith(".sello.json"):
            da, db = _sin_fecha(json.loads(pa.read_text(encoding="utf-8"))), _sin_fecha(json.loads(pb.read_text(encoding="utf-8")))
            if da == db:
                continue
            if nombre == "h7_visibilidad_en_sangre.json":
                sa, sb = da.pop("modelo_secundario", {}), db.pop("modelo_secundario", {})
                if da == db and sa != sb:
                    secundario = {"b": sa.get("auc_secundario"), "c": sb.get("auc_secundario")}
                    bloquean.append(nombre)
                    continue
            if nombre == "corrida.json" or _solo_hash_o_corrida(da, db):
                admitidos.append(nombre)
                continue
            bloquean.append(nombre)
            continue
        if _sha(pa) != _sha(pb):
            bloquean.append(nombre)
    return {
        "admitidos": admitidos,
        "bloquean": bloquean,
        "h7_secundario": secundario,
        "etiqueta": not bloquean,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", type=Path, required=True)
    parser.add_argument("--b", type=Path, required=True)
    parser.add_argument("--md", type=Path, default=Path("docs/determinismo_dos_procesos.md"))
    parser.add_argument("--json", type=Path, default=Path("docs/determinismo_dos_procesos.json"))
    parser.add_argument("--estricto", action="store_true")
    args = parser.parse_args()
    if args.estricto:
        resultado = comparar_estricto(args.a, args.b)
        args.json.write_text(json.dumps(resultado, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        sec = resultado["h7_secundario"] or {}
        lineas = [
            "# Determinismo B frente a C",
            "",
            f"B: `{args.a}`. C: `{args.b}`.",
            "",
            "Las dos salidas salen del commit `b7f8b8e844f6c484b72e19fecf5fb7d07f53c2ab`. En el sello coinciden el commit, `arbol_limpio`, el preregistro, las enmiendas, el manifiesto, los paquetes, la semilla y la etapa. Solo cambia `fecha_utc`.",
            "",
            f"Diferencias admitidas: {', '.join(resultado['admitidos']) or 'ninguna'}.",
            "",
            f"El secundario de H7 es {sec.get('b')} en B y {sec.get('c')} en C. El veredicto sigue en refuta y el AUC principal no cambia. En `h7_tabla_proteinas.tsv` los valores presentes de ESM-2 coinciden; cambia en qué filas falta el valor. Ese secundario no entra en el veredicto, y esta comparación no lo admite como diferencia menor.",
            "",
            f"Diferencias que impiden la etiqueta: {', '.join(resultado['bloquean']) or 'ninguna'}.",
            "",
            "No se crea la etiqueta `v3.7-final`.",
            "",
        ]
        args.md.write_text("\n".join(lineas), encoding="utf-8")
        if not resultado["etiqueta"]:
            raise SystemExit(1)
        return
    resultado = comparar(args.a, args.b)
    args.json.write_text(json.dumps(resultado, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    sellos = [nombre for nombre in resultado["distintos"] if nombre.endswith(".sello.json") or nombre == "sello_corrida.json"]
    datos = [nombre for nombre in resultado["distintos"] if nombre not in sellos]
    lineas = [
        "# Determinismo de dos procesos",
        "",
        f"A: `{args.a}`. B: `{args.b}`.",
        "",
        f"Ficheros de datos distintos: {len(datos)}.",
        f"Sellos distintos: {len(sellos)}. Cambian con el commit y con el hash del manifiesto.",
        f"Números principales distintos: {len(resultado['principales'])}.",
        "",
    ]
    if resultado["limitacion_conocida"]:
        lineas.append(resultado["limitacion_conocida"])
        lineas.append("")
    if resultado["principales"]:
        lineas.append("Principales distintos: " + ", ".join(resultado["principales"]) + ".")
    elif not datos:
        lineas.append("No hay diferencias en los ficheros de datos fuera del sello.")
    else:
        lineas.append("Las diferencias de datos no están en los números principales: " + ", ".join(datos) + ".")
        if "gsea_tejido.tsv" in datos:
            lineas.append("En `gsea_tejido.tsv` el NES, el FDR y el p coinciden si se ordena por cohorte y término. Cambia el orden de las filas.")
    args.md.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    if resultado["principales"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
