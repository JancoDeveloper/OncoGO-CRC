#!/usr/bin/env python3
"""Descarga idempotente de una cohorte extra de PRIDE (C11 y D3).

Para qué. El meta-análisis de efectos aleatorios de C10 necesita al menos tres
cohortes de tejido para que la heterogeneidad (R3) sea calculable; con dos no lo
es. v1 tiene dos. Esta tercera cohorte es la que lo desbloquea.

Qué se descarga. Solo la tabla de proteínas de la búsqueda ya hecha por los
autores, no los ficheros crudos. En MaxQuant eso es ``proteinGroups.txt``. No se
reprocesa nada: reanalizar espectros con otros parámetros introduciría una
diferencia de método que luego se confundiría con heterogeneidad biológica.

Comprobación de nivel A de H0. Si dos proyectos distintos resultan ser el mismo
dato con otro identificador, el meta-análisis contaría dos veces la misma
cohorte. Se compara el hash del fichero descargado con el de las cohortes que ya
están y se avisa si coinciden.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

# El endpoint correcto es projects/{acc}/files. `files/byProject` responde 200
# con cuerpo vacío, que es peor que un error porque parece funcionar.
API = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{}/files?pageSize=500"
TIEMPO_LIMITE = 900

# Nombres que interesan, en orden de preferencia.
PREFERIDOS = ("proteingroups.txt", "proteins.txt", "proteins.tsv")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for trozo in iter(lambda: fh.read(1 << 20), b""):
            h.update(trozo)
    return h.hexdigest()


def _json(url: str):
    pedido = urllib.request.Request(url, headers={"User-Agent": "OncoGO-CRC/2.0"})
    with urllib.request.urlopen(pedido, timeout=120) as r:
        return json.loads(r.read())


def _ftp_a_https(enlace: str) -> str:
    if enlace.startswith("ftp://"):
        return "https://" + enlace[len("ftp://") :]
    return enlace


def _elegir(ficheros: list[dict]) -> dict | None:
    for preferido in PREFERIDOS:
        for f in ficheros:
            if str(f.get("fileName", "")).lower() == preferido:
                return f
    # Si no hay un nombre exacto, cualquier tabla de resultados con proteínas.
    for f in ficheros:
        nombre = str(f.get("fileName", "")).lower()
        if "protein" in nombre and nombre.endswith((".txt", ".tsv", ".csv")):
            return f
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accession", default="PXD019504")
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    ap.add_argument("--forzar", action="store_true")
    args = ap.parse_args()

    acc = args.accession.upper()
    base = Path(args.workspace) / "external" / "pride" / acc
    base.mkdir(parents=True, exist_ok=True)
    ficha_path = base / "estado.json"
    ficha = json.loads(ficha_path.read_text(encoding="utf-8")) if ficha_path.exists() else {}

    if ficha.get("estado") == "usada" and not args.forzar:
        destino = Path(ficha["fichero"])
        if destino.exists() and sha256(destino) == ficha.get("sha256"):
            print(f"[{acc}] ya está y el hash coincide, no se descarga.")
            return 0

    def _fallo(motivo: str) -> int:
        print(f"[{acc}] NO DISPONIBLE: {motivo}", file=sys.stderr)
        ficha_path.write_text(
            json.dumps(
                {
                    "accession": acc,
                    "estado": "no_disponible",
                    "motivo": motivo[:800],
                    "fecha": datetime.now(UTC).isoformat(timespec="seconds"),
                    "consecuencia": (
                        "C10 se queda con dos cohortes de tejido, así que R3 "
                        "(heterogeneidad) sigue sin ser calculable."
                    ),
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        return 0

    try:
        listado = _json(API.format(acc))
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ValueError) as e:
        return _fallo(f"el listado de PRIDE no responde: {e}")

    ficheros = listado if isinstance(listado, list) else listado.get("_embedded", {}).get("files", [])
    if not ficheros:
        return _fallo("PRIDE no devuelve ningún fichero para este accession")

    plano = []
    for f in ficheros:
        enlaces = f.get("publicFileLocations") or []
        url = ""
        for e in enlaces:
            valor = e.get("value") or ""
            if valor.startswith(("http", "ftp")):
                url = _ftp_a_https(valor)
                break
        plano.append({"fileName": f.get("fileName", ""), "url": url, "bytes": f.get("fileSizeBytes")})

    elegido = _elegir(plano)
    if elegido is None or not elegido["url"]:
        nombres = sorted({p["fileName"] for p in plano})[:40]
        return _fallo(
            "no hay tabla de proteínas de la búsqueda de los autores. "
            f"Ficheros disponibles (primeros 40): {nombres}"
        )

    destino = base / elegido["fileName"]
    try:
        pedido = urllib.request.Request(
            elegido["url"], headers={"User-Agent": "OncoGO-CRC/2.0"}
        )
        with urllib.request.urlopen(pedido, timeout=TIEMPO_LIMITE) as r, destino.open("wb") as fh:
            while trozo := r.read(1 << 20):
                fh.write(trozo)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
        return _fallo(f"la descarga de {elegido['fileName']} falló: {e}")

    hash_nuevo = sha256(destino)

    # H0 nivel A: ¿es el mismo fichero que alguna cohorte que ya tengo?
    duplicados = []
    for otro in (Path(args.workspace) / "data" / "raw").rglob("*"):
        if otro.is_file() and otro.stat().st_size == destino.stat().st_size:
            if sha256(otro) == hash_nuevo:
                duplicados.append(str(otro))

    ficha_path.write_text(
        json.dumps(
            {
                "accession": acc,
                "estado": "usada",
                "fichero": str(destino),
                "fichero_elegido": elegido["fileName"],
                "url": elegido["url"],
                "bytes": destino.stat().st_size,
                "sha256": hash_nuevo,
                "fecha": datetime.now(UTC).isoformat(timespec="seconds"),
                "n_ficheros_en_el_proyecto": len(plano),
                "h0_nivel_a_duplicados_por_hash": duplicados,
                "h0_nivel_a_veredicto": (
                    "cohorte independiente: ningún fichero existente tiene el mismo hash"
                    if not duplicados
                    else "ATENCIÓN: hash idéntico a un fichero ya presente; no es una cohorte nueva"
                ),
                "nota": (
                    "Se usa la tabla de proteínas de la búsqueda de los autores. No se "
                    "reprocesan espectros: hacerlo introduciría una diferencia de método "
                    "que después se confundiría con heterogeneidad biológica."
                ),
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    print(
        f"[{acc}] {elegido['fileName']}: {destino.stat().st_size} bytes. "
        f"Duplicados por hash: {duplicados or 'ninguno'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
