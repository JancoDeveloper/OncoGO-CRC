"""Compara D y E sin dar por iguales dos sha256 distintos.

No modifica comparar_corridas.py. Quita el bloque sello y las claves
fecha_utc, pid, inicio y fin. El resto, sha256 incluidos, tiene que coincidir.
Si un sha256 difiere, el fichero al que apunta solo puede diferir en esas claves
o en otros sha256 que se expliquen igual.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CLAVES = {"fecha_utc", "pid", "inicio", "fin"}
HEX = re.compile(r"^[0-9a-f]{64}$")
HASH_A_RUTA = {
    "sha256_fuente": ("fuente",),
    "r1_sha256": ("r1_fuente",),
    "sha256": ("fichero", "fuente", "ruta", "archivo"),
}


def sha256(path: Path) -> str:
    resumen = hashlib.sha256()
    with path.open("rb") as entrada:
        while True:
            bloque = entrada.read(1 << 20)
            if not bloque:
                break
            resumen.update(bloque)
    return resumen.hexdigest()


def quitar(obj):
    if isinstance(obj, dict):
        return {k: quitar(v) for k, v in obj.items() if k != "sello" and k not in CLAVES}
    if isinstance(obj, list):
        return [quitar(v) for v in obj]
    return obj


def mismo(a, b) -> bool:
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    return a == b


def hojas_con_padre(a, b, pref=(), padre_a=None, padre_b=None):
    if isinstance(a, dict) and isinstance(b, dict):
        for clave in sorted(set(a) | set(b)):
            ruta = pref + (clave,)
            if clave not in a or clave not in b:
                yield ruta, a.get(clave), b.get(clave), a, b
            else:
                yield from hojas_con_padre(a[clave], b[clave], ruta, a, b)
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            yield pref, len(a), len(b), padre_a, padre_b
            return
        for i, (x, y) in enumerate(zip(a, b)):
            yield from hojas_con_padre(x, y, pref + (i,), a, b)
        return
    if not mismo(a, b):
        yield pref, a, b, padre_a, padre_b


def resolver(corrida: Path, referencia: str) -> Path | None:
    texto = str(referencia).replace("\\", "/")
    if texto.startswith("results/v3_7/"):
        candidato = corrida / texto.split("results/v3_7/", 1)[1]
        if candidato.is_file():
            return candidato
    if (corrida / texto).is_file():
        return corrida / texto
    if (REPO / texto).is_file():
        return REPO / texto
    return None


def es_hash(clave, valor_a, valor_b) -> bool:
    return (
        isinstance(clave, str)
        and "sha256" in clave
        and isinstance(valor_a, str)
        and isinstance(valor_b, str)
        and HEX.fullmatch(valor_a) is not None
        and HEX.fullmatch(valor_b) is not None
    )


def normalizar_tests(texto: str) -> str:
    texto = re.sub(r"\d+\.\d+s", "TIEMPO", texto)
    texto = re.sub(r"\(0:\d+:\d+\)", "(TIEMPO)", texto)
    return re.sub(r"pid=\d+", "pid=PID", texto)


def gsea_igual(a: Path, b: Path) -> bool:
    import pandas as pd
    uno = pd.read_csv(a, sep="\t").sort_values(["cohorte", "Term"]).reset_index(drop=True)
    dos = pd.read_csv(b, sep="\t").sort_values(["cohorte", "Term"]).reset_index(drop=True)
    return uno.equals(dos)


def fichero_explicado(a: Path, b: Path, vistos: set[tuple[str, str]]) -> tuple[bool, str]:
    marca = (str(a), str(b))
    if marca in vistos:
        return True, "ciclo"
    vistos.add(marca)
    if not a.is_file() or not b.is_file():
        return False, "falta el fichero apuntado"
    if a.name == "tests.txt":
        if normalizar_tests(a.read_text(encoding="utf-8")) == normalizar_tests(b.read_text(encoding="utf-8")):
            return True, "tests.txt, tiempos y pid"
        return False, "tests.txt cambia fuera de los tiempos y los pid"
    if a.name == "gsea_tejido.tsv":
        if gsea_igual(a, b):
            return True, "gsea_tejido.tsv, orden de filas"
        return False, "gsea_tejido.tsv cambia el contenido"
    if a.suffix == ".json" or a.name.endswith(".sello.json") or b.suffix == ".json":
        try:
            da = quitar(json.loads(a.read_text(encoding="utf-8")))
            db = quitar(json.loads(b.read_text(encoding="utf-8")))
        except json.JSONDecodeError:
            return a.read_bytes() == b.read_bytes(), "no es JSON"
        pendientes = []
        for ruta, x, y, padre_a, padre_b in hojas_con_padre(da, db):
            clave = ruta[-1] if ruta else ""
            if es_hash(clave, x, y) and isinstance(padre_a, dict) and isinstance(padre_b, dict):
                pendientes.append((ruta, x, y, padre_a, padre_b))
            elif isinstance(x, float) and isinstance(y, float) and math.isnan(x) and math.isnan(y):
                continue
            else:
                return False, f"cambia {'.'.join(str(p) for p in ruta)}"
        for ruta, x, y, padre_a, padre_b in pendientes:
            clave = str(ruta[-1])
            nombres = HASH_A_RUTA.get(clave, ("fichero", "fuente", "ruta"))
            ref_a = next((padre_a[n] for n in nombres if n in padre_a), None)
            ref_b = next((padre_b[n] for n in nombres if n in padre_b), None)
            if ref_a is None or ref_a != ref_b:
                return False, f"sha256 sin fichero en {'.'.join(str(p) for p in ruta)}"
            corrida_a = _corrida_de(a)
            corrida_b = _corrida_de(b)
            destino_a = resolver(corrida_a, str(ref_a))
            destino_b = resolver(corrida_b, str(ref_b))
            if destino_a is None or destino_b is None:
                return False, f"no se encuentra {ref_a}"
            ok, motivo = fichero_explicado(destino_a, destino_b, vistos)
            if not ok:
                return False, f"{a.name} → {ref_a}: {motivo}"
        return True, "solo fecha, proceso o un hash que los resume"
    iguales = a.read_bytes() == b.read_bytes()
    return iguales, "bytes iguales" if iguales else f"{a.name} cambia los bytes"


def _corrida_de(path: Path) -> Path:
    for padre in [path, *path.parents]:
        if (padre / "sello_corrida.json").is_file() and (padre / "certificado.json").is_file():
            return padre
    return path.parent


def comparar(a: Path, b: Path) -> dict:
    nombres = sorted(
        {p.name for p in a.iterdir() if p.is_file() and not p.name.startswith("._")}
        | {p.name for p in b.iterdir() if p.is_file() and not p.name.startswith("._")}
    )
    admitidos = []
    bloquean = []
    for nombre in nombres:
        pa, pb = a / nombre, b / nombre
        if not pa.is_file() or not pb.is_file():
            bloquean.append(f"{nombre}: falta en una corrida")
            continue
        if nombre == "tests.txt":
            if normalizar_tests(pa.read_text(encoding="utf-8")) == normalizar_tests(pb.read_text(encoding="utf-8")):
                if sha256(pa) != sha256(pb):
                    admitidos.append(nombre)
            else:
                bloquean.append(nombre)
            continue
        if nombre == "gsea_tejido.tsv":
            if sha256(pa) == sha256(pb):
                continue
            if gsea_igual(pa, pb):
                admitidos.append(nombre)
            else:
                bloquean.append(nombre)
            continue
        if nombre.endswith(".json") or nombre.endswith(".sello.json"):
            ok, motivo = fichero_explicado(pa, pb, set())
            if ok:
                if sha256(pa) != sha256(pb):
                    admitidos.append(f"{nombre} ({motivo})")
            else:
                bloquean.append(f"{nombre}: {motivo}")
            continue
        if sha256(pa) != sha256(pb):
            bloquean.append(nombre)
    return {"admitidos": admitidos, "bloquean": bloquean, "identicas": not bloquean}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", type=Path, required=True)
    parser.add_argument("--b", type=Path, required=True)
    parser.add_argument("--md", type=Path, default=REPO / "docs" / "determinismo_d_e_estricto.md")
    args = parser.parse_args()
    resultado = comparar(args.a, args.b)
    lineas = [
        "# Determinismo D frente a E, regla estricta",
        "",
        f"D: `{args.a}`. E: `{args.b}`.",
        "",
        "Se quita el bloque `sello` y las claves `fecha_utc`, `pid`, `inicio` y `fin`. El resto se compara tal cual, sha256 incluidos. Si un sha256 difiere, el fichero al que apunta solo puede diferir en esas claves o en otro sha256 explicado de la misma forma. `gsea_tejido.tsv` se compara ordenado por cohorte y término. En `tests.txt` solo se normalizan los tiempos y los pid.",
        "",
        f"Diferencias explicadas por fecha o proceso: {len(resultado['admitidos'])}.",
        "",
        f"Diferencias sin explicar: {', '.join(resultado['bloquean']) or 'ninguna'}.",
        "",
        "D y E coinciden con esta regla." if resultado["identicas"] else "Hay una diferencia que no es de fecha ni de proceso. No se publica la web hasta explicarla.",
        "",
    ]
    args.md.write_text("\n".join(lineas), encoding="utf-8")
    print(json.dumps({"identicas": resultado["identicas"], "bloquean": resultado["bloquean"][:12], "n_admitidos": len(resultado["admitidos"])}, ensure_ascii=False))
    return 0 if resultado["identicas"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
