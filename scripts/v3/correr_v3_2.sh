#!/bin/bash
# Corre v3.2 de punta a punta. No arranca el suero si el report no coincide.
set -euo pipefail
cd "$(dirname "$0")/../.."
REPORT="/Volumes/JANCODATA/OncoGO-CRC/data/raw/plasma/DIA-NN_report.tsv"
ESPERADO=213313241
if [[ ! -f "$REPORT" ]]; then
  echo "No está $REPORT"
  exit 1
fi
BYTES=$(wc -c < "$REPORT" | tr -d ' ')
if [[ "$BYTES" != "$ESPERADO" ]]; then
  echo "El report va en $BYTES de $ESPERADO bytes. Cuando llegue, vuelve a lanzar este script."
  exit 1
fi
make all-v3-2
echo "Listo. Lee results/v3_2/suero_resumen.json y results/v3_2/report/informe.html"
