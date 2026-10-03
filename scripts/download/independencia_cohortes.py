#!/usr/bin/env python3
"""Auditoría de independencia entre cohortes de PRIDE (H0 nivel A y nivel B).

Por qué hace falta. Un meta-análisis que incluya dos veces la misma cohorte da
una homogeneidad falsamente alta y un intervalo falsamente estrecho. El
repositorio ya documenta un caso: PXD019504 y PXD031705 comparten 44 ficheros
con el mismo SHA-1. Y PXD019504 y PXD034575 son del mismo grupo, así que podrían
compartir pacientes aunque los ficheros sean distintos.

Dos niveles de comprobación.

- **Nivel A, los ficheros**. Se piden a PRIDE las sumas de verificación de todos
  los ficheros de cada proyecto y se buscan coincidencias. Un hash repetido
  significa el mismo fichero depositado dos veces.
- **Nivel B, los números**. Dos depósitos pueden ser ficheros distintos y aun así
  venir del mismo paciente. Se correlaciona cada muestra de una cohorte contra
  cada muestra de la otra sobre las proteínas comunes. Una correlación muy alta
  entre dos muestras concretas es la firma de un espécimen repetido.

El nivel B es el que decide, porque es el que mira el dato y no la etiqueta.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from datetime import UTC, datetime
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

API = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{}/files?pageSize=500"
UMBRAL_CORRELACION = 0.98
MIN_PROTEINAS_COMUNES = 200


def _checksums(acc: str) -> tuple[dict[str, str], str]:
    """Hash → nombre de fichero, tal como lo publica PRIDE."""
    try:
        pedido = urllib.request.Request(
            API.format(acc), headers={"User-Agent": "OncoGO-CRC/2.0"}
        )
        with urllib.request.urlopen(pedido, timeout=180) as r:
            datos = json.loads(r.read())
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ValueError) as e:
        return {}, f"PRIDE no responde para {acc}: {e}"
    ficheros = datos if isinstance(datos, list) else datos.get("_embedded", {}).get("files", [])
    salida = {}
    for f in ficheros:
        suma = str(f.get("checksum") or "").strip().lower()
        if suma:
            salida[suma] = str(f.get("fileName", ""))
    return salida, ""


def _matriz_lfq(path: Path) -> pd.DataFrame | None:
    """Log2 de las columnas LFQ, indexado por el primer UniProt del grupo."""
    if not path.exists():
        return None
    columnas = pd.read_csv(path, sep="\t", nrows=0, low_memory=False).columns
    lfq = [c for c in columnas if str(c).startswith("LFQ intensity ")]
    clave = "Majority protein IDs" if "Majority protein IDs" in columnas else "Protein IDs"
    if not lfq or clave not in columnas:
        return None
    tabla = pd.read_csv(path, sep="\t", usecols=[clave, *lfq], low_memory=False)
    # Las banderas de MaxQuant se combinan antes de filtrar: aplicarlas una a una
    # sobre la tabla ya recortada desalinea las máscaras.
    presentes = [
        b
        for b in ("Reverse", "Potential contaminant", "Only identified by site")
        if b in columnas
    ]
    if presentes:
        marcas = pd.read_csv(path, sep="\t", usecols=presentes, low_memory=False)
        limpio = marcas.isna().all(axis=1).to_numpy()
        tabla = tabla.loc[limpio]
    ids = tabla[clave].astype(str).str.split(";").str[0].str.split("|").str[-1].str.upper()
    valores = tabla[lfq].apply(pd.to_numeric, errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        log2 = np.where(valores > 0, np.log2(valores), np.nan)
    salida = pd.DataFrame(log2, columns=[c.replace("LFQ intensity ", "") for c in lfq])
    salida.index = ids
    return salida.loc[~salida.index.duplicated(keep="first")]


def _correlaciones_cruzadas(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    """Cada muestra de A contra cada muestra de B, sobre proteínas comunes."""
    comunes = a.index.intersection(b.index)
    if len(comunes) < MIN_PROTEINAS_COMUNES:
        return {
            "estado": "no_evaluable",
            "motivo": f"solo {len(comunes)} proteínas comunes, hacen falta {MIN_PROTEINAS_COMUNES}",
            "n_proteinas_comunes": int(len(comunes)),
        }
    x, y = a.loc[comunes], b.loc[comunes]
    filas = []
    for ca in x.columns:
        va = x[ca].to_numpy(dtype=float)
        for cb in y.columns:
            vb = y[cb].to_numpy(dtype=float)
            ok = np.isfinite(va) & np.isfinite(vb)
            if ok.sum() < MIN_PROTEINAS_COMUNES:
                continue
            r = float(np.corrcoef(va[ok], vb[ok])[0, 1])
            filas.append({"muestra_a": ca, "muestra_b": cb, "r": r, "n": int(ok.sum())})
    if not filas:
        return {"estado": "no_evaluable", "motivo": "ningún par de muestras con solape suficiente"}
    detalle = pd.DataFrame(filas).sort_values("r", ascending=False, kind="stable")
    sospechosos = detalle.loc[detalle["r"] >= UMBRAL_CORRELACION]
    return {
        "estado": "evaluado",
        "n_proteinas_comunes": int(len(comunes)),
        "n_pares_de_muestras": int(len(detalle)),
        "r_maximo": float(detalle["r"].iloc[0]),
        "r_mediano": float(detalle["r"].median()),
        "umbral": UMBRAL_CORRELACION,
        "n_pares_sospechosos": int(len(sospechosos)),
        "pares_sospechosos": sospechosos.head(20).to_dict("records"),
        "veredicto": (
            "cohortes independientes: ninguna pareja de muestras se parece lo "
            "suficiente para ser el mismo espécimen"
            if sospechosos.empty
            else "ATENCIÓN: hay parejas de muestras casi idénticas; no son cohortes "
            "independientes y no pueden entrar las dos en el mismo meta-análisis"
        ),
        "detalle": detalle,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accessions", nargs="+", default=["PXD019504", "PXD031705", "PXD034575"])
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    args = ap.parse_args()

    ws = Path(args.workspace)
    salida = ws / "external" / "independencia"
    salida.mkdir(parents=True, exist_ok=True)

    # Nivel A: los hashes que publica PRIDE.
    hashes: dict[str, dict[str, str]] = {}
    fallos: dict[str, str] = {}
    for acc in args.accessions:
        h, motivo = _checksums(acc)
        if motivo:
            fallos[acc] = motivo
            print(f"[independencia] {acc}: {motivo}", file=sys.stderr)
        hashes[acc] = h
        print(f"[independencia] {acc}: {len(h)} ficheros con suma de verificación.")

    nivel_a = []
    for a, b in combinations(args.accessions, 2):
        compartidos = sorted(set(hashes.get(a, {})) & set(hashes.get(b, {})))
        nivel_a.append(
            {
                "par": f"{a}_vs_{b}",
                "n_ficheros_a": len(hashes.get(a, {})),
                "n_ficheros_b": len(hashes.get(b, {})),
                "n_hashes_compartidos": len(compartidos),
                "ejemplos": [
                    {"hash": h, "nombre_en_a": hashes[a][h], "nombre_en_b": hashes[b][h]}
                    for h in compartidos[:5]
                ],
                "veredicto": (
                    "sin ficheros compartidos"
                    if not compartidos
                    else f"ATENCIÓN: {len(compartidos)} ficheros idénticos; es el mismo dato"
                ),
            }
        )
        print(f"[independencia] nivel A {a} vs {b}: {nivel_a[-1]['veredicto']}")

    # Nivel B: los números. Solo entre las cohortes cuya matriz esté en disco.
    rutas = {
        "PXD019504": ws / "external" / "pride" / "PXD019504" / "proteinGroups.txt",
        "PXD034575": ws / "data" / "raw" / "pride" / "Cohort1_proteinGroups.txt",
    }
    for acc in args.accessions:
        if acc not in rutas:
            candidatos = list((ws / "external" / "pride" / acc).glob("proteinGroups.txt"))
            if candidatos:
                rutas[acc] = candidatos[0]
    matrices = {}
    for acc, path in rutas.items():
        m = _matriz_lfq(path)
        if m is not None:
            matrices[acc] = m
            print(f"[independencia] {acc}: matriz de {m.shape[0]} proteínas y {m.shape[1]} muestras.")

    nivel_b = {}
    for a, b in combinations(sorted(matrices), 2):
        res = _correlaciones_cruzadas(matrices[a], matrices[b])
        detalle = res.pop("detalle", None)
        if detalle is not None:
            detalle.to_csv(salida / f"correlaciones_{a}_vs_{b}.tsv", sep="\t", index=False)
        nivel_b[f"{a}_vs_{b}"] = res
        print(f"[independencia] nivel B {a} vs {b}: {res.get('veredicto', res.get('motivo'))}")

    ficha = {
        "fecha": datetime.now(UTC).isoformat(timespec="seconds"),
        "accessions": args.accessions,
        "fallos_de_api": fallos,
        "nivel_a_hashes_de_pride": nivel_a,
        "nivel_b_correlacion_entre_muestras": nivel_b,
        "umbral_correlacion": UMBRAL_CORRELACION,
        "nota": (
            "El nivel B es el que decide. Dos depósitos pueden tener ficheros "
            "distintos y venir del mismo paciente; lo que delata eso es que dos "
            "muestras concretas se parezcan casi perfectamente."
        ),
    }
    (salida / "estado.json").write_text(
        json.dumps(ficha, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[independencia] ficha en {salida / 'estado.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
