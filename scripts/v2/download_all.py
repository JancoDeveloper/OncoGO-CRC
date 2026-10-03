#!/usr/bin/env python3
"""Orquestador de las descargas de v2.

Una regla de diseño: el fallo de una fuente no detiene a las demás. Si Reactome
no responde, el resto sigue y Reactome queda marcado ``no_disponible`` con el
motivo. Lo contrario convertiría cada caída de red en una tarde perdida y, peor,
invitaría a rellenar el hueco con algo.

Este script no descarga nada por sí mismo. Llama a los de ``scripts/download/``,
que son idempotentes: si el fichero ya está y el hash cuadra, no vuelven a
pedirlo. Así que ejecutarlo dos veces seguidas es barato y seguro.

Al final escribe ``external/estado_de_las_fuentes.json``, que es el resumen que
lee la sección 6 de la validación: cada fuente con su estado, y ninguna sin él.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[2]
DESCARGAS = RAIZ / "scripts" / "download"

# Orden de ejecución. Las colecciones de vías van primero porque C4 depende de
# ellas y son las más rápidas; la auditoría de independencia va después de las
# cohortes porque las necesita descargadas para comparar.
#
# La tercera columna es la carpeta donde la fuente deja su ``estado.json``, que no
# siempre coincide con su nombre: PXD019504 vive bajo ``pride/`` junto a las demás
# cohortes. Se declara en lugar de deducirse para que un cambio de ruta salte
# aquí y no se traduzca en un «sin declarar» silencioso.
FUENTES: list[tuple[str, str, str, str]] = [
    (
        "vias",
        "colecciones_vias.py",
        "vias",
        "colecciones de vías: Reactome y MSigDB Hallmark (C4)",
    ),
    (
        "hpa",
        "human_protein_atlas.py",
        "hpa",
        "concentraciones en sangre del Human Protein Atlas (C8)",
    ),
    (
        "pride_cohorte_extra",
        "pride_cohorte_extra.py",
        "pride/PXD019504",
        "tercera cohorte de tejido pareado, PXD019504 (C10)",
    ),
    (
        "independencia_cohortes",
        "independencia_cohortes.py",
        "independencia",
        "auditoría de independencia entre cohortes (C11)",
    ),
    (
        "tcga_coadread",
        "tcga_coadread.py",
        "tcga_coadread",
        "CEA sérico, pureza ABSOLUTE y CEACAM5 en TCGA (D1)",
    ),
    (
        "alvez_2023",
        "alvez_2023.py",
        "alvez_2023",
        "recurso de sangre pan-cáncer, Olink Explore (D4b)",
    ),
]

# Tiempo máximo por fuente. El HPA descarga un zip grande; el resto es rápido.
TIEMPO_MAXIMO_S = 1800


def _ejecutar(script: Path, workspace: Path) -> dict:
    if not script.exists():
        return {
            "resultado": "script_ausente",
            "motivo": f"no existe {script.relative_to(RAIZ)}",
        }
    try:
        proc = subprocess.run(
            [sys.executable, str(script), "--workspace", str(workspace)],
            capture_output=True,
            text=True,
            timeout=TIEMPO_MAXIMO_S,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {
            "resultado": "tiempo_agotado",
            "motivo": f"el script pasó de {TIEMPO_MAXIMO_S} s",
        }
    salida = (proc.stdout or "").strip().splitlines()
    error = (proc.stderr or "").strip().splitlines()
    return {
        "resultado": "ok" if proc.returncode == 0 else "fallo",
        "codigo_de_salida": proc.returncode,
        # Solo las últimas líneas: el resto es ruido de progreso.
        "ultimas_lineas": salida[-6:],
        "ultimas_lineas_de_error": error[-6:],
    }


def _estado_declarado(workspace: Path, carpeta: str) -> dict:
    """Lee el ``estado.json`` que deja cada descarga.

    El estado que vale es el que escribe la propia fuente, no el código de salida
    del proceso: un script puede terminar en 0 y haber marcado la fuente
    ``requiere_acceso``, que es información y no un error.
    """
    path = workspace / "external" / carpeta / "estado.json"
    if not path.exists():
        return {"estado": "sin_declarar", "motivo": f"no existe {path.name}"}
    try:
        d = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        return {"estado": "sin_declarar", "motivo": f"{path.name} ilegible: {e}"}
    if "estado" in d:
        return {
            "estado": d["estado"],
            "motivo": d.get("motivo"),
            "fecha_de_acceso": d.get("fecha_de_acceso"),
        }

    # Algunas descargas traen varias subfuentes en el mismo fichero: `vias` lleva
    # Reactome y Hallmark por separado. El estado del conjunto es el peor de los
    # suyos, porque decir «usada» cuando una mitad no está sería falso.
    orden = ["usada", "requiere_acceso", "no_disponible"]
    subfuentes = {
        k: v["estado"] for k, v in d.items() if isinstance(v, dict) and "estado" in v
    }
    if not subfuentes:
        return {
            "estado": "sin_declarar",
            "motivo": f"{path.name} no trae ninguna clave 'estado'",
        }
    peor = max(subfuentes.values(), key=lambda e: orden.index(e) if e in orden else 99)
    return {
        "estado": peor,
        "subfuentes": subfuentes,
        "motivo": None if peor == "usada" else f"el estado peor es de {peor}",
        "fecha_de_acceso": next(
            (
                v.get("fecha_de_acceso")
                for v in d.values()
                if isinstance(v, dict) and v.get("fecha_de_acceso")
            ),
            None,
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    ap.add_argument(
        "--solo",
        default=None,
        help="ejecuta una sola fuente, por su nombre",
    )
    args = ap.parse_args()
    workspace = Path(args.workspace)

    fuentes = [f for f in FUENTES if args.solo is None or f[0] == args.solo]
    if not fuentes:
        nombres = ", ".join(f[0] for f in FUENTES)
        print(f"[descargas] {args.solo!r} no es una fuente. Hay: {nombres}")
        return 2

    informe = {}
    for nombre, script, carpeta, descripcion in fuentes:
        print(f"[descargas] {nombre}: {descripcion}")
        ejecucion = _ejecutar(DESCARGAS / script, workspace)
        declarado = _estado_declarado(workspace, carpeta)
        informe[nombre] = {
            "descripcion": descripcion,
            "script": f"scripts/download/{script}",
            "carpeta": f"external/{carpeta}",
            "ejecucion": ejecucion,
            **declarado,
        }
        marca = declarado["estado"]
        if ejecucion["resultado"] != "ok" and marca == "sin_declarar":
            marca = f"{marca} ({ejecucion['resultado']})"
        print(f"[descargas] {nombre}: {marca}")

    resumen = {
        "fecha": datetime.now(UTC).isoformat(timespec="seconds"),
        "fuentes": informe,
        "recuento": {
            estado: sum(1 for v in informe.values() if v["estado"] == estado)
            for estado in sorted({v["estado"] for v in informe.values()})
        },
        "regla": (
            "El fallo de una fuente no detiene a las demás. Ninguna fuente "
            "marcada no_disponible o requiere_acceso se sustituye por datos "
            "simulados."
        ),
    }
    destino = workspace / "external" / "estado_de_las_fuentes.json"
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_text(
        json.dumps(resumen, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[descargas] recuento: {resumen['recuento']}")
    print(f"[descargas] resumen en {destino}")

    # Solo es un fallo que una fuente no sepa decir en qué estado está. Que una
    # fuente requiera cuenta o no responda es información, no un error del
    # orquestador.
    sin_declarar = [k for k, v in informe.items() if v["estado"] == "sin_declarar"]
    if sin_declarar:
        print(
            f"[descargas] sin estado declarado: {', '.join(sin_declarar)}",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
