"""Sello de una corrida de v3.2. Si el árbol de git está sucio, la etapa no sigue."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path


class ArbolSucio(RuntimeError):
    """La corrida no se firma si hay cambios sin commitear."""


def _git(root: Path, *args: str) -> str:
    resultado = subprocess.run(
        ["git", *args],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    return resultado.stdout.strip()


def _sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    resumen = hashlib.sha256()
    with path.open("rb") as entrada:
        while True:
            bloque = entrada.read(1 << 20)
            if not bloque:
                break
            resumen.update(bloque)
    return resumen.hexdigest()


def _version(modulo: str) -> str:
    try:
        importado = __import__(modulo)
    except ImportError:
        return "no_instalado"
    return str(getattr(importado, "__version__", "desconocida"))


def hacer_sello(root: Path, *, etapa: str, manifest: Path | None = None) -> dict:
    sucio = _git(root, "status", "--porcelain")
    paquetes = {
        "numpy": _version("numpy"),
        "pandas": _version("pandas"),
        "scipy": _version("scipy"),
        "scikit-learn": _version("sklearn"),
    }
    return {
        "commit": _git(root, "rev-parse", "HEAD"),
        "arbol_limpio": sucio == "",
        "fecha_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sha256_preregistro": _sha256(root / "config" / "preregistro_v3.yaml"),
        "sha256_enmiendas": _sha256(root / "docs" / "enmiendas_v3.md"),
        "sha256_manifest_entradas": _sha256(manifest) if manifest else None,
        "python": sys.version.split()[0],
        "paquetes": paquetes,
        "semilla": 20260914,
        "etapa": etapa,
    }


def exigir_arbol_limpio(sello: dict) -> None:
    if not sello.get("arbol_limpio"):
        raise ArbolSucio("el árbol de git no está limpio; no se firma la corrida")


def escribir_sello(path: Path, sello: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(sello, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
