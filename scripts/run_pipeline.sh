#!/usr/bin/env bash
# Equivalente operativo a la secuencia documentada en Informacion/Manual/manual/03_comandos_por_fase.txt
# usando la CLI moderna (los módulos `python -m src.01_*` no son válidos en Python por nombres numéricos).
#
# Volumen externo: si está definido ONCOGO_CRC_DATA_WORKSPACE (ruta absoluta recomendada), se pasa a
# `oncogo-crc run --data-workspace …` para que descargas, data/, results/, reports/, metadata/ y docs/web
# vivan fuera del clon. Ver GUIA_CONFIGURACION.md § almacenamiento externo.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
if [[ ! -d .venv ]]; then
  python3 -m venv .venv
fi
# shellcheck source=/dev/null
source .venv/bin/activate
pip install -q -e .
DW="${ONCOGO_CRC_DATA_WORKSPACE:-}"
if [[ -n "$DW" ]]; then
  exec oncogo-crc run --data-workspace "$DW" --stage all "$@"
else
  exec oncogo-crc run --stage all "$@"
fi
