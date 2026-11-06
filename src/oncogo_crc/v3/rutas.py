"""v3 no lee el pipeline antiguo ni escribe fuera de results/v3."""

from __future__ import annotations

from pathlib import Path

LECTURAS_PROHIBIDAS = (
    "results/differential/",
    "results/go/",
    "results/gsea/",
    "results/ml/",
)


def _texto(path: Path) -> str:
    return str(path).replace("\\", "/")


def rechazar_si_es_pipeline_antiguo(path: Path) -> None:
    texto = _texto(path)
    for relativo in LECTURAS_PROHIBIDAS:
        if relativo in texto:
            raise RuntimeError(
                f"v3 no lee {relativo}: ese directorio es del pipeline antiguo"
            )


def exigir_salida_v3(path: Path) -> None:
    texto = _texto(path)
    rechazar_si_es_pipeline_antiguo(path)
    if "/results/v3/" not in texto and not texto.endswith("/results/v3"):
        raise RuntimeError(f"v3 solo escribe bajo results/v3: {path}")
