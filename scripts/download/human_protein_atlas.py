#!/usr/bin/env python3
"""Descarga idempotente del Human Protein Atlas (fuente D4, necesaria para C8).

Para qué. C8 pregunta por la ventana de detección: una proteína sobreexpresada
en tejido solo puede verse en sangre si su concentración en circulación está por
encima de lo que el método alcanza. Sin concentraciones, «no se detectó en
plasma» y «no está en plasma» se confunden.

Qué se usa. ``proteinatlas.tsv.zip`` trae, por gen, la concentración en sangre
medida por inmunoensayo y por espectrometría de masas, y la clasificación de
proteína secretada o de membrana. Se extraen solo esas columnas a una tabla
limpia; el fichero completo no entra en el repositorio.

Si la descarga falla, la fuente queda ``no_disponible`` con el motivo y C8 se
marca sin evaluar. No se inventan concentraciones.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import urllib.error
import urllib.request
import zipfile
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

URL = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
TIEMPO_LIMITE = 600
CITA = "Uhlén M et al. Tissue-based map of the human proteome. Science 2015. proteinatlas.org"
LICENCIA = "CC BY-SA 4.0"

# Nombres tal como vienen en el TSV. Se buscan por subcadena porque el HPA
# cambia la redacción entre versiones.
COLUMNAS = {
    "gene": ["Gene"],
    "uniprot": ["Uniprot"],
    "clase_proteica": ["Protein class"],
    "localizacion_secretoma": ["Secretome location"],
    "conc_sangre_inmunoensayo_pg_l": ["Blood concentration - Conc. blood IM [pg/L]"],
    "conc_sangre_ms_pg_l": ["Blood concentration - Conc. blood MS [pg/L]"],
}


def sha256_bytes(datos: bytes) -> str:
    return hashlib.sha256(datos).hexdigest()


def _columna(df: pd.DataFrame, candidatos: list[str]) -> str | None:
    for c in candidatos:
        if c in df.columns:
            return c
    bajas = {str(c).lower(): c for c in df.columns}
    for c in candidatos:
        clave = c.lower()
        for baja, original in bajas.items():
            if clave in baja:
                return original
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    ap.add_argument("--forzar", action="store_true")
    args = ap.parse_args()

    base = Path(args.workspace) / "external" / "hpa"
    base.mkdir(parents=True, exist_ok=True)
    limpio = base / "hpa_concentracion_sangre.tsv"
    ficha_path = base / "estado.json"
    ficha = json.loads(ficha_path.read_text(encoding="utf-8")) if ficha_path.exists() else {}

    if limpio.exists() and not args.forzar and ficha.get("hpa", {}).get("estado") == "usada":
        print("[hpa] ya está la tabla limpia, no se descarga.")
        return 0

    try:
        pedido = urllib.request.Request(URL, headers={"User-Agent": "OncoGO-CRC/2.0"})
        with urllib.request.urlopen(pedido, timeout=TIEMPO_LIMITE) as r:
            crudo = r.read()
        with zipfile.ZipFile(io.BytesIO(crudo)) as z:
            nombre = next(n for n in z.namelist() if n.endswith(".tsv"))
            tabla = pd.read_csv(io.BytesIO(z.read(nombre)), sep="\t", low_memory=False)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, StopIteration, OSError) as e:
        print(f"[hpa] NO DISPONIBLE: {e}", file=sys.stderr)
        ficha["hpa"] = {
            "estado": "no_disponible",
            "motivo": str(e)[:500],
            "url": URL,
            "fecha": datetime.now(UTC).isoformat(timespec="seconds"),
            "cita": CITA,
            "licencia": LICENCIA,
            "consecuencia": "C8 queda sin evaluar: no hay concentraciones en sangre.",
        }
        ficha_path.write_text(json.dumps(ficha, indent=2, ensure_ascii=False), encoding="utf-8")
        return 0

    renombra: dict[str, str] = {}
    faltan: list[str] = []
    for destino, candidatos in COLUMNAS.items():
        origen = _columna(tabla, candidatos)
        if origen is None:
            faltan.append(destino)
        else:
            renombra[origen] = destino
    salida = tabla[list(renombra)].rename(columns=renombra)
    salida.to_csv(limpio, sep="\t", index=False)

    con_im = int(salida.get("conc_sangre_inmunoensayo_pg_l", pd.Series(dtype=float)).notna().sum())
    con_ms = int(salida.get("conc_sangre_ms_pg_l", pd.Series(dtype=float)).notna().sum())
    ficha["hpa"] = {
        "estado": "usada",
        "url": URL,
        "miembro_del_zip": nombre,
        "sha256_del_zip": sha256_bytes(crudo),
        "sha256_de_la_tabla_limpia": hashlib.sha256(limpio.read_bytes()).hexdigest(),
        "n_genes": int(len(salida)),
        "n_con_concentracion_inmunoensayo": con_im,
        "n_con_concentracion_ms": con_ms,
        "columnas_no_encontradas": faltan,
        "fecha": datetime.now(UTC).isoformat(timespec="seconds"),
        "cita": CITA,
        "licencia": LICENCIA,
        "nota": (
            "Solo se guardan las columnas necesarias. El TSV completo del HPA no "
            "entra en el repositorio ni en el workspace."
        ),
    }
    ficha_path.write_text(json.dumps(ficha, indent=2, ensure_ascii=False), encoding="utf-8")
    print(
        f"[hpa] {len(salida)} genes; {con_im} con concentración por inmunoensayo, "
        f"{con_ms} por espectrometría. Columnas no encontradas: {faltan or 'ninguna'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
