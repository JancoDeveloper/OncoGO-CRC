"""Etapas que usa la corrida v3.7 y la regeneración de sus entradas v1.

Las etapas de v2 y de v3.2 a v3.6 siguen en la etiqueta v3.7-final.
No se registran aquí.
"""

from __future__ import annotations

from collections.abc import Callable

from oncogo_crc.settings import Settings
from oncogo_crc.stages import certification, certified_inference, v3_7

StageFn = Callable[[Settings], None]

STAGES: dict[str, StageFn] = {
    "certified_inference": certified_inference.run,
    "certification": certification.run,
    "v3_7": v3_7.run,
}

STAGE_ORDER = ["certified_inference", "certification", "v3_7"]

# `oncogo-crc run --stage all` ejecuta solo v3_7. Las dos etapas anteriores
# regeneran results/certified y results/certification y se lanzan con su nombre
# o con `make certified`, para no reescribir esas entradas en cada corrida.
EXCLUDED_FROM_ALL: frozenset[str] = frozenset({"certified_inference", "certification"})
