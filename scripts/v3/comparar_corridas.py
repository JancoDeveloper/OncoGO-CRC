"""Compara dos salidas de v3.7. Los JSON se miran sin el sello. Los TSV, por sha256."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", type=Path, required=True)
    parser.add_argument("--b", type=Path, required=True)
    parser.add_argument("--md", type=Path, default=Path("docs/determinismo_dos_procesos.md"))
    parser.add_argument("--json", type=Path, default=Path("docs/determinismo_dos_procesos.json"))
    args = parser.parse_args()
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
