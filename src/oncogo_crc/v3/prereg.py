"""Lee config/preregistro_v3.yaml. No inventa umbrales."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


class FaltaEnPreregistro(KeyError):
    """Se pidió un valor que el preregistro de v3 no declara."""


class PreregistroV3:
    def __init__(self, root: Path):
        path = Path(root) / "config" / "preregistro_v3.yaml"
        if not path.is_file():
            raise FileNotFoundError(path)
        self._data = yaml.safe_load(path.read_text(encoding="utf-8"))

    def get(self, ruta: str) -> Any:
        actual: Any = self._data
        for parte in ruta.split("."):
            if not isinstance(actual, dict) or parte not in actual:
                raise FaltaEnPreregistro(ruta)
            actual = actual[parte]
        return actual

    def semilla(self, nombre: str) -> int:
        base = int(self.get("semillas.base"))
        return base + int(self.get(f"semillas.desplazamientos.{nombre}"))
