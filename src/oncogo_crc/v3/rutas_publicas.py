"""Rutas que pueden publicarse. Nunca una ruta del disco local."""

from __future__ import annotations

import hashlib
from pathlib import Path


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


def ruta_publica(path: Path, settings) -> dict:
    rel = settings.relative_to_project(path)
    if rel.startswith("/") or rel.startswith("\\") or (len(rel) > 2 and rel[1] == ":"):
        rel = Path(path).name
    return {"ruta": rel.replace("\\", "/"), "sha256": _sha256(Path(path))}


def texto_con_ruta_local(texto: str) -> bool:
    return "/Volumes/" in texto or "/Users/" in texto or "C:\\" in texto or "C:/" in texto
