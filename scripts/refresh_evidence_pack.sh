#!/usr/bin/env bash
# Orquestador end-to-end de los TSV auxiliares que potencian la corrida real.
# Ejecuta UNA vez (o cuando quieras actualizar) los descargadores caros antes del pipeline.
#
# Uso:
#   bash scripts/refresh_evidence_pack.sh           # todo (uniprot + TCGA log2FC + Cox + CMS)
#   SKIP_SURVIVAL=1 bash scripts/refresh_evidence_pack.sh   # omitir Cox (lento)
#   SKIP_CMS=1 bash scripts/refresh_evidence_pack.sh        # omitir CMS predictor (requiere matriz armonizada previa)
#
# Variables de entorno respetadas:
#   ONCOGO_CRC_DATA_WORKSPACE  → directorio externo para datasets/resultados
#   COX_ONLY_PROTEIN_GENES=1   → restringe Cox a los ~6 000 genes presentes en la matriz CPTAC log2norm
#                                (muy recomendado: reduce de ~10 min a ~3 min)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

PY="${PYTHON:-python}"
WS="${ONCOGO_CRC_DATA_WORKSPACE:-$REPO_ROOT}"
META_TCGA="$WS/metadata/tcga_auxiliary"
META_CLIN="$WS/metadata/clinical_layers"
META_MAPS="$WS/metadata/protein_id_maps"

echo "==> Workspace efectivo: $WS"
mkdir -p "$META_TCGA" "$META_CLIN" "$META_MAPS"

# 1) UniProt bulk (gen + alias + UniProt canónico para ~20 k entradas humanas reviewed)
if [[ ! -s "$META_MAPS/gene_uniprot_human_reviewed_bulk.tsv" ]] || [[ "${REFRESH_UNIPROT:-0}" == "1" ]]; then
  echo "==> [1/4] Descargando UniProt bulk (human reviewed) …"
  "$PY" "$REPO_ROOT/scripts/fetch_uniprot_human_reviewed_gene_map.py" \
    --out "$META_MAPS/gene_uniprot_human_reviewed_bulk.tsv"
else
  echo "==> [1/4] UniProt bulk ya presente: $META_MAPS/gene_uniprot_human_reviewed_bulk.tsv (export REFRESH_UNIPROT=1 para regenerar)"
fi

# 2) TCGA RNA-seq log2FC (UCSC Xena)
if [[ ! -s "$META_TCGA/tcga_coadread_log2fc.tsv" ]] || [[ "${REFRESH_TCGA:-0}" == "1" ]]; then
  echo "==> [2/4] Calculando TCGA COAD+READ log2FC desde Xena …"
  "$PY" "$REPO_ROOT/scripts/fetch_tcga_coadread_log2fc.py" \
    --out "$META_TCGA/tcga_coadread_log2fc.tsv" \
    --cache "$META_TCGA/.xena_cache"
else
  echo "==> [2/4] TCGA log2FC ya presente (export REFRESH_TCGA=1 para regenerar)"
fi

# 3) TCGA Survival Cox (opcional, lento)
if [[ "${SKIP_SURVIVAL:-0}" != "1" ]]; then
  if [[ ! -s "$META_TCGA/tcga_coadread_cox_os.tsv" ]] || [[ "${REFRESH_SURVIVAL:-0}" == "1" ]]; then
    EXTRA=""
    if [[ "${COX_ONLY_PROTEIN_GENES:-0}" == "1" ]] && [[ -s "$WS/data/processed/cptac_colorectal_log2norm.tsv" ]]; then
      EXTRA="--only-genes-tsv $WS/data/processed/cptac_colorectal_log2norm.tsv"
      echo "    (restringiendo Cox a genes de matriz proteómica CPTAC)"
    fi
    echo "==> [3/4] Cox HR por gen (TCGA OS) — puede tardar ~5-15 min …"
    "$PY" -m pip install --quiet "lifelines>=0.27" || true
    "$PY" "$REPO_ROOT/scripts/fetch_tcga_coadread_survival.py" \
      --out "$META_TCGA/tcga_coadread_cox_os.tsv" \
      --cache "$META_TCGA/.xena_cache" \
      $EXTRA
  else
    echo "==> [3/4] TCGA Cox ya presente (export REFRESH_SURVIVAL=1 para regenerar)"
  fi
else
  echo "==> [3/4] Cox omitido (SKIP_SURVIVAL=1)"
fi

# 4) CMS predictor proteómico (requiere matriz CPTAC log2norm ya generada por harmonize/normalize)
if [[ "${SKIP_CMS:-0}" != "1" ]]; then
  MAT="$WS/data/processed/cptac_colorectal_log2norm.tsv"
  MD="$WS/data/processed/cptac_colorectal_metadata.tsv"
  if [[ -s "$MAT" ]] && [[ -s "$MD" ]]; then
    if [[ ! -s "$META_CLIN/cptac_cms_predicted.tsv" ]] || [[ "${REFRESH_CMS:-0}" == "1" ]]; then
      echo "==> [4/4] Asignando CMS por firma proteómica simplificada …"
      "$PY" "$REPO_ROOT/scripts/predict_cms_from_proteome.py" \
        --matrix "$MAT" \
        --metadata "$MD" \
        --out "$META_CLIN/cptac_cms_predicted.tsv"
    else
      echo "==> [4/4] CMS predicted ya presente (export REFRESH_CMS=1 para regenerar)"
    fi
  else
    echo "==> [4/4] CMS omitido: falta matriz/metadata CPTAC log2norm. Correr pipeline harmonize+normalize primero."
  fi
else
  echo "==> [4/4] CMS predictor omitido (SKIP_CMS=1)"
fi

echo "==> Evidence pack listo en $WS/metadata/"
echo "    Verificar config.yaml apunta a:"
echo "      external_evidence.tcga_rnaseq_log2fc_tsv  = metadata/tcga_auxiliary/tcga_coadread_log2fc.tsv"
echo "      external_evidence.survival_gene_cox_tsv   = metadata/tcga_auxiliary/tcga_coadread_cox_os.tsv"
echo "      harmonize_public_metadata.cptac.cms_subtype_supplement_tsv = metadata/clinical_layers/cptac_cms_predicted.tsv"
