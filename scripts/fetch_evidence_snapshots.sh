#!/usr/bin/env bash
# Orquestador de snapshots de evidencia externa para la Capa C (priorización convergente).
# Genera los TSV en <workspace>/metadata/evidence/ que consume el stage `prioritization`.
#
# Requiere que `differential` ya se haya ejecutado (usa la firma como universo de genes),
# o pásale --genes-tsv / --differential-tsv a cada script.
#
# Uso:
#   bash scripts/fetch_evidence_snapshots.sh                 # APIs ligeras (DGIdb, OpenTargets, cBioPortal, GTEx)
#   WITH_GTEX=0 bash scripts/fetch_evidence_snapshots.sh     # omitir GTEx (descarga ~50 MB)
#   OT_LIMIT=300 bash scripts/fetch_evidence_snapshots.sh    # limitar OpenTargets (lento, per-gen)
#
# DepMap y HPA requieren archivos/descargas grandes → ejecútalos aparte (ver más abajo).
#
# Variables:
#   ONCOGO_CRC_DATA_WORKSPACE   → workspace externo (metadata/evidence, results/…)
#   DIFF_TSV                    → firma diferencial usada como universo de genes

set -euo pipefail
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

PY="${PYTHON:-python}"
WS="${ONCOGO_CRC_DATA_WORKSPACE:-$REPO_ROOT}"
EV="$WS/metadata/evidence"
DIFF_TSV="${DIFF_TSV:-$WS/results/differential/cptac_tumor_vs_normal_diff.tsv}"
mkdir -p "$EV"

echo "==> Workspace: $WS"
echo "==> Evidence dir: $EV"
echo "==> Universo de genes (firma): $DIFF_TSV"
if [[ ! -s "$DIFF_TSV" ]]; then
  echo "!! No existe $DIFF_TSV. Ejecuta antes: oncogo-crc run --stage differential"
  echo "   (o define DIFF_TSV=/ruta/otra_firma.tsv)"
fi

echo "==> [1/4] DGIdb (druggability) …"
"$PY" scripts/fetch_evidence_dgidb.py --differential-tsv "$DIFF_TSV" \
  --out "$EV/dgidb_gene_drugs.tsv" || echo "   (DGIdb falló; continúo)"

echo "==> [2/4] cBioPortal (frecuencia de mutación COADREAD) …"
"$PY" scripts/fetch_evidence_cbioportal.py --differential-tsv "$DIFF_TSV" \
  --out "$EV/cbioportal_alteration_freq.tsv" || echo "   (cBioPortal falló; continúo)"

echo "==> [3/4] Open Targets (tractabilidad) — puede tardar (per-gen) …"
"$PY" scripts/fetch_evidence_opentargets.py --differential-tsv "$DIFF_TSV" \
  --limit "${OT_LIMIT:-0}" --out "$EV/opentargets_tractability.tsv" || echo "   (OpenTargets falló; continúo)"

if [[ "${WITH_GTEX:-1}" == "1" ]]; then
  echo "==> [4/4] GTEx (seguridad: expresión en tejidos normales) …"
  "$PY" scripts/fetch_evidence_gtex.py --differential-tsv "$DIFF_TSV" \
    --cache "$EV/.gtex_cache/gene_median_tpm.gct.gz" \
    --out "$EV/gtex_normal_max.tsv" || echo "   (GTEx falló; continúo)"
else
  echo "==> [4/4] GTEx omitido (WITH_GTEX=0)"
fi

echo
echo "==> DepMap y HPA (descargas grandes / manuales):"
echo "    DepMap: descarga CRISPRGeneEffect.csv y Model.csv de https://depmap.org/portal/download/ y corre:"
echo "      $PY scripts/fetch_evidence_depmap.py --gene-effect-csv <..>/CRISPRGeneEffect.csv --model-csv <..>/Model.csv --out $EV/depmap_crc_gene_effect.tsv"
echo "    HPA:"
echo "      $PY scripts/fetch_evidence_hpa.py --differential-tsv \"$DIFF_TSV\" --out $EV/hpa_tumor_specificity.tsv"
echo
echo "==> Snapshots en: $EV"
echo "    Ahora ejecuta la Capa C:  oncogo-crc run --stage prioritization"
