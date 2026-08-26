# Snapshots de evidencia externa (Capa C — priorización convergente)

Cada archivo es un **snapshot local versionado** que un proveedor pluggable
(`src/oncogo_crc/utils/evidence_providers.py`) lee para anotar los candidatos. Si un archivo **no existe**,
esa fuente se marca `available=false` y **no penaliza** el score de convergencia (diseño no-promisorio).

Formato: TSV con encabezado. La columna de gen puede llamarse `gene_symbol`, `gene`, `symbol`, etc.
(se detecta por heurística). Opcionalmente añade un sidecar `<archivo>.version` con la versión/fecha de la fuente.

| Archivo | Columnas mínimas | Interpretación | Fuente / cómo obtener |
|---|---|---|---|
| `depmap_crc_gene_effect.tsv` | `gene_symbol`, `chronos_mean` (opc. `essential_fraction`) | Chronos más **negativo** = gen más esencial en líneas CCR | [DepMap](https://depmap.org/portal/download/) → `CRISPRGeneEffect.csv`, promediar líneas CCR |
| `opentargets_tractability.tsv` | `gene_symbol`, `tractability` (0–1) | Mayor = más tratable/accionable | [Open Targets](https://platform.opentargets.org/downloads) (tractability) |
| `dgidb_gene_drugs.tsv` | `gene_symbol`, `n_drugs` | Nº de fármacos que interactúan (druggability) | [DGIdb](https://www.dgidb.org/downloads) `interactions.tsv`, agrupar por gen |
| `hpa_tumor_specificity.tsv` | `gene_symbol`, `tumor_specificity` | Mayor = más específico de tumor | [Human Protein Atlas](https://www.proteinatlas.org/about/download) |
| `gtex_normal_max.tsv` | `gene_symbol`, `max_normal_tpm` | **Menor** expresión normal = mejor ventana de seguridad | [GTEx](https://gtexportal.org/home/downloads/adult-gtex) (medianas por tejido → máx) |
| `cbioportal_alteration_freq.tsv` | `gene_symbol`, `alteration_frequency` (0–1) | Mayor frecuencia de alteración en CCR | [cBioPortal](https://www.cbioportal.org/) (COADREAD, export por gen) |

## Generación automática (conectores)

En lugar de armar los TSV a mano, usa los conectores incluidos (ver `docs/TUTORIAL_EJECUCION.md` §5):

```bash
# APIs ligeras (DGIdb, cBioPortal, Open Targets, GTEx) usando la firma como universo de genes
bash scripts/fetch_evidence_snapshots.sh

# Individuales
python scripts/fetch_evidence_dgidb.py       --differential-tsv results/differential/cptac_tumor_vs_normal_diff.tsv
python scripts/fetch_evidence_cbioportal.py  --differential-tsv results/differential/cptac_tumor_vs_normal_diff.tsv
python scripts/fetch_evidence_opentargets.py --differential-tsv results/differential/cptac_tumor_vs_normal_diff.tsv --limit 300
python scripts/fetch_evidence_gtex.py        --differential-tsv results/differential/cptac_tumor_vs_normal_diff.tsv
# Descargas grandes / manuales:
python scripts/fetch_evidence_depmap.py --gene-effect-csv CRISPRGeneEffect.csv --model-csv Model.csv
python scripts/fetch_evidence_hpa.py    --differential-tsv results/differential/cptac_tumor_vs_normal_diff.tsv
```

Cada conector escribe también un sidecar `<archivo>.version` con la procedencia.

## Uso

1. Genera los TSV con los conectores (arriba) o colócalos a mano aquí (o en el `data_workspace`).
2. Ejecuta: `oncogo-crc run --stage prioritization` (o el pipeline completo).
3. Salidas en `results/prioritization/`:
   - `hypothesis_priority.tsv` — tabla priorizada.
   - `convergence_evidence_matrix.tsv` — detalle por señal (transparencia).
   - `prioritization_provenance.json` — versión/ruta de cada fuente + pesos.

## Notas de reproducibilidad y ética

- Cada snapshot debe registrar su **versión** (sidecar `.version`) para trazabilidad FAIR.
- Los resultados son **hipótesis exploratorias**, no dianas terapéuticas validadas.
- Respeta las licencias/términos de uso de cada base de datos al versionar snapshots.
