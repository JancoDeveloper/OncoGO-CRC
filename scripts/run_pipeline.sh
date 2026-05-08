#!/usr/bin/env bash
# Equivalente operativo a la secuencia documentada en Informacion/Manual/manual/03_comandos_por_fase.txt
# usando la CLI moderna (los módulos `python -m src.01_*` no son válidos en Python por nombres numéricos).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
if [[ ! -d .venv ]]; then
  python3 -m venv .venv
fi
# shellcheck source=/dev/null
source .venv/bin/activate
pip install -q -e .
exec oncogo-crc run --stage all "$@"
