"""Registra cada open() durante la corrida limpia y rechaza lecturas de corridas anteriores."""

from __future__ import annotations

import sys
from pathlib import Path

_ACTIVO = False
_ENGANCHADO = False
_EXCEPCIONES: tuple[str, ...] = ()
_EVENTOS: list[dict] = []
_PROHIBIDAS: list[str] = []
_VERSION = "v3_6"


def clasificar(ruta: str, excepciones: tuple[str, ...] = (), version: str = "v3_6") -> str:
    texto = str(ruta).replace("\\", "/")
    marca = "/results/"
    indice = texto.lower().find(marca)
    if indice < 0:
        return "fuera"
    resto = texto[indice + 1:]
    if resto.startswith(f"results/{version}/") or resto == f"results/{version}":
        return "permitida"
    if resto.startswith("results/certified/"):
        return "permitida"
    for excepcion in excepciones:
        if resto == excepcion or resto.startswith(excepcion.rstrip("/") + "/"):
            return "excepcion"
    return "prohibida"


def _gancho(evento: str, args) -> None:
    if not _ACTIVO or evento != "open" or not args:
        return
    ruta = args[0]
    if isinstance(ruta, int):
        return
    clase = clasificar(str(ruta), _EXCEPCIONES, _VERSION)
    if clase == "fuera":
        return
    modo = args[1] if len(args) > 1 and isinstance(args[1], str) else ""
    _EVENTOS.append({"ruta": str(ruta), "clase": clase, "modo": modo})
    if clase == "prohibida":
        _PROHIBIDAS.append(str(ruta))


def activar(excepciones: tuple[str, ...], version: str = "v3_6") -> None:
    global _ACTIVO, _EXCEPCIONES, _ENGANCHADO, _VERSION
    _VERSION = version
    _EXCEPCIONES = excepciones
    _EVENTOS.clear()
    _PROHIBIDAS.clear()
    if not _ENGANCHADO:
        sys.addaudithook(_gancho)
        _ENGANCHADO = True
    _ACTIVO = True


def cerrar() -> tuple[list[dict], list[str]]:
    global _ACTIVO
    _ACTIVO = False
    return list(_EVENTOS), list(_PROHIBIDAS)


def excepciones_de(ruta: Path) -> tuple[str, ...]:
    if not ruta.is_file():
        return ()
    import yaml
    datos = yaml.safe_load(ruta.read_text(encoding="utf-8")) or {}
    filas = datos.get("entradas") or datos.get("excepciones") or []
    return tuple(str(item["ruta"]) for item in filas if isinstance(item, dict) and "ruta" in item)
