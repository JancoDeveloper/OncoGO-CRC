#!/bin/bash
# Corre el cálculo de v3.3. No arranca si el report no coincide.
# Los controles del tejido pueden fallar si una proteína no está en MaxQuant:
# eso se anota y el resto de las etapas sigue, porque el suero no depende de ellas.
set -u
cd "$(dirname "$0")/../.."
REPORT="/Volumes/JANCODATA/OncoGO-CRC/data/raw/plasma/DIA-NN_report.tsv"
ESPERADO=213313241
if [[ ! -f "$REPORT" ]]; then
  echo "No está $REPORT"
  exit 1
fi
BYTES=$(wc -c < "$REPORT" | tr -d ' ')
if [[ "$BYTES" != "$ESPERADO" ]]; then
  echo "El report va en $BYTES de $ESPERADO bytes."
  exit 1
fi
HASH=$(shasum -a 256 "$REPORT" | awk '{print toupper($1)}')
if [[ "$HASH" != "A8A5AADD74391AFEF68CF4D6C3A7D3860023101E5C92590E2E87C152A39E59E8" ]]; then
  echo "SHA-256 distinto: $HASH"
  exit 1
fi
fallo=0
for objetivo in v3-3-suero v3-3-tejido; do
  make "$objetivo" || fallo=1
done
make v3-3-controles || echo "CONTROLES_CON_FALLO"
for objetivo in v3-3-hipotesis v3-3-multiverso v3-3-ml v3-3-report; do
  make "$objetivo" || fallo=1
done
exit "$fallo"
