#!/usr/bin/env python3
"""Descarga idempotente de las colecciones de vías de v2 (C4).

Tres colecciones. GO Biological Process ya está en el repositorio desde v1, así
que aquí solo se traen las dos que faltan.

- **Reactome**. `ReactomePathways.gmt` es descarga libre.
- **MSigDB Hallmark**. Históricamente pide registro. Si la descarga falla, la
  colección queda marcada `no_disponible` con el motivo y el análisis sigue con
  las otras dos. No se inventa el fichero.

Idempotente: si el fichero ya está y el hash coincide con el registrado, no se
vuelve a descargar. Se escribe una ficha por colección en
`external/vias/estado.json`.
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

TIEMPO_LIMITE = 120

FUENTES: dict[str, dict] = {
    "reactome": {
        "url": "https://reactome.org/download/current/ReactomePathways.gmt.zip",
        "fichero": "reactome_pathways.gmt",
        "dentro_del_zip": "ReactomePathways.gmt",
        "licencia": "CC BY 4.0",
        "cita": "Milacic M et al. The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res.",
    },
    "hallmark": {
        # Varios espejos. Se prueban en orden; el primero que responda gana.
        "urls": [
            "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt",
            "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt",
        ],
        "fichero": "msigdb_hallmark.gmt",
        "licencia": "CC BY 4.0 para Hallmark; el resto de MSigDB tiene condiciones propias",
        "cita": "Liberzon A et al. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Syst 2015.",
    },
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for trozo in iter(lambda: fh.read(1 << 20), b""):
            h.update(trozo)
    return h.hexdigest()


def _bajar(url: str) -> bytes:
    pedido = urllib.request.Request(url, headers={"User-Agent": "OncoGO-CRC/2.0"})
    with urllib.request.urlopen(pedido, timeout=TIEMPO_LIMITE) as r:
        return r.read()


def _descargar_reactome(destino: Path, cfg: dict) -> dict:
    crudo = _bajar(cfg["url"])
    with zipfile.ZipFile(io.BytesIO(crudo)) as z:
        nombre = cfg["dentro_del_zip"]
        if nombre not in z.namelist():
            candidatos = [n for n in z.namelist() if n.endswith(".gmt")]
            if not candidatos:
                raise RuntimeError(f"el zip no trae ningún .gmt: {z.namelist()}")
            nombre = candidatos[0]
        destino.write_bytes(z.read(nombre))
    return {"url": cfg["url"], "miembro_del_zip": nombre}


def _descargar_hallmark(destino: Path, cfg: dict) -> dict:
    errores = []
    for url in cfg["urls"]:
        try:
            datos = _bajar(url)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            errores.append(f"{url}: {e}")
            continue
        if not datos.lstrip().startswith(b"HALLMARK"):
            errores.append(f"{url}: la respuesta no parece un GMT de Hallmark")
            continue
        destino.write_bytes(datos)
        return {"url": url}
    raise RuntimeError("; ".join(errores))


def _lineas_y_terminos(path: Path) -> tuple[int, int]:
    n_term = 0
    n_gen = 0
    for linea in path.read_text(encoding="utf-8", errors="replace").splitlines():
        campos = linea.rstrip("\n").split("\t")
        if len(campos) < 3:
            continue
        n_term += 1
        n_gen += len([g for g in campos[2:] if g.strip()])
    return n_term, n_gen


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    ap.add_argument("--forzar", action="store_true", help="vuelve a descargar aunque esté")
    args = ap.parse_args()

    base = Path(args.workspace) / "external" / "vias"
    base.mkdir(parents=True, exist_ok=True)
    ficha_path = base / "estado.json"
    ficha = json.loads(ficha_path.read_text(encoding="utf-8")) if ficha_path.exists() else {}

    descargadores = {"reactome": _descargar_reactome, "hallmark": _descargar_hallmark}
    for nombre, cfg in FUENTES.items():
        destino = base / cfg["fichero"]
        previo = ficha.get(nombre, {})
        if destino.exists() and not args.forzar and previo.get("sha256") == sha256(destino):
            print(f"[vias] {nombre}: ya está y el hash coincide, no se descarga.")
            continue
        try:
            detalle = descargadores[nombre](destino, cfg)
        except Exception as e:  # noqa: BLE001 - cualquier fallo deja la fuente no disponible
            print(f"[vias] {nombre}: NO DISPONIBLE. {e}", file=sys.stderr)
            ficha[nombre] = {
                "estado": "no_disponible",
                "motivo": str(e)[:500],
                "fecha": datetime.now(UTC).isoformat(timespec="seconds"),
                "licencia": cfg["licencia"],
                "cita": cfg["cita"],
            }
            continue
        n_term, n_gen = _lineas_y_terminos(destino)
        ficha[nombre] = {
            "estado": "usada",
            "fichero": str(destino),
            "sha256": sha256(destino),
            "bytes": destino.stat().st_size,
            "n_terminos": n_term,
            "n_pares_termino_gen": n_gen,
            "fecha": datetime.now(UTC).isoformat(timespec="seconds"),
            "licencia": cfg["licencia"],
            "cita": cfg["cita"],
            **detalle,
        }
        print(f"[vias] {nombre}: {n_term} términos, {n_gen} pares término-gen.")

    ficha_path.write_text(json.dumps(ficha, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[vias] ficha en {ficha_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
