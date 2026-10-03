"""Lee el preregistro de v2. Ningún umbral se escribe en el código.

Si un análisis necesita un número que no está en ``config/preregistro_v2.yaml``,
el acceso falla en vez de inventar un valor por defecto. Así un umbral no puede
aparecer a mitad del análisis sin quedar registrado.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

PREREG_PATH = "config/preregistro_v2.yaml"
REGISTRO_PATH = "config/registro_subindices.yaml"


class FaltaEnPreregistro(KeyError):
    """Se pidió un umbral que el preregistro no declara."""


@lru_cache(maxsize=4)
def _load(root: str, rel: str) -> dict[str, Any]:
    path = Path(root) / rel
    if not path.is_file():
        raise FileNotFoundError(f"no existe {path}: el preregistro es obligatorio en v2")
    return yaml.safe_load(path.read_text(encoding="utf-8"))


class Preregistro:
    """Acceso por ruta con puntos: ``prereg.get("c1_permutaciones.cptac_n_permutaciones")``."""

    def __init__(self, root: Path):
        self._root = Path(root)
        self._data = _load(str(self._root), PREREG_PATH)
        self._subindices = _load(str(self._root), REGISTRO_PATH)

    @property
    def data(self) -> dict[str, Any]:
        return self._data

    def get(self, dotted: str) -> Any:
        node: Any = self._data
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                raise FaltaEnPreregistro(
                    f"{dotted!r} no está en {PREREG_PATH}. "
                    "Añádelo al preregistro y regístralo en docs/enmiendas_v2.md."
                )
            node = node[part]
        return node

    def semilla(self, nombre: str) -> int:
        """Semilla derivada de la base con el desplazamiento declarado para ese análisis."""
        base = int(self.get("semillas.base"))
        offsets = self.get("semillas.desplazamientos")
        if nombre not in offsets:
            raise FaltaEnPreregistro(
                f"el análisis {nombre!r} no tiene desplazamiento de semilla en {PREREG_PATH}"
            )
        return base + int(offsets[nombre])

    def subindices(self) -> list[dict[str, Any]]:
        return list(self._subindices["subindices"])

    def subindice(self, codigo: str) -> dict[str, Any]:
        for sub in self.subindices():
            if sub["codigo"] == codigo:
                return dict(sub)
        raise FaltaEnPreregistro(f"{codigo!r} no está en {REGISTRO_PATH}")


def cargar(root: Path) -> Preregistro:
    return Preregistro(root)
