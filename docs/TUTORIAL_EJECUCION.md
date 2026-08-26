# Tutorial end-to-end — ejecutar OncoGO-CRC y obtener resultados

Guía práctica para instalar, ejecutar el pipeline, incorporar las **bases de datos externas** (Capa C) y
explorar los resultados en la **web interactiva**. Al final tendrás una **tabla de hipótesis priorizadas**
y un dashboard consultable.

> Recuerda: los resultados son **hipótesis exploratorias con evidencia convergente**, no dianas validadas.

---

## 0. Requisitos

- Python 3.10–3.14, `git`, ~2–5 GB libres (más si descargas GTEx/DepMap/HPA).
- Conexión a internet (descargas CPTAC/PRIDE y APIs de evidencia).

---

## 1. Instalación

```bash
git clone <repo> OncoGO-CRC && cd OncoGO-CRC
python -m venv .venv && source .venv/bin/activate
pip install -e ".[go-remote,explain,survival,web]"
```

Extras: `go-remote` (GSEA/KEGG), `explain` (SHAP), `survival` (Cox), `web` (Streamlit+Plotly), `docs` (MkDocs).

---

## 2. Configurar el workspace de datos (recomendado)

Los datos pesados (descargas, resultados, web) van a un directorio aparte del código:

```bash
export ONCOGO_CRC_DATA_WORKSPACE=/ruta/con/espacio/OncoGO-CRC
```

También puedes fijar `paths.data_workspace` en `config/config.yaml`. Prioridad:
`--data-workspace` (CLI) > variable de entorno > `config.yaml`.

---

## 3. Prueba rápida en modo demo (sin descargar nada)

Valida que todo el flujo corre con datos sintéticos:

```bash
oncogo-crc run --stage all --demo
```

> En `--demo` el workspace se ignora (no contamina resultados reales). Sirve para verificar la instalación.

---

## 4. Corrida real del pipeline

### 4.1 (Opcional) Paquete auxiliar TCGA/UniProt/CMS
```bash
bash scripts/refresh_evidence_pack.sh          # UniProt bulk + TCGA log2FC + Cox + CMS
```

### 4.2 Pipeline completo
```bash
oncogo-crc run --stage all
```

O por etapas (útil para depurar). Orden principal:

| Etapa | Qué hace |
|---|---|
| `download_real` | Descarga CPTAC/PRIDE/plasma + snapshots de BD |
| `harmonize` | Proteína→gen canónico (UniProt + alias) |
| `qc` | log2, KNN, z-score, PCA |
| `differential` | t-test Welch/pareado + BH-FDR (**genera la firma**) |
| `cohort_overlap` | Firma reproducible CPTAC∩PRIDE |
| `go` / `gsea` / `kegg` | Enriquecimiento funcional |
| `network` | Red de coexpresión + hubs |
| `ml` | L1/EN/RF + nested CV + SHAP |
| `external_evidence` | Concordancia TCGA / Cox |
| **`prioritization`** | **Capa C: score de convergencia** |
| `publish_results` | Documentación web (MkDocs) |

Al terminar tendrás, entre otros:
`results/differential/*_diff.tsv`, `results/reproducibility/…`, `results/run_manifest.json`.

---

## 5. Incorporar las bases de datos externas (Capa C)

La Capa C anota la firma con evidencia pública. Cada fuente es un **snapshot local** en
`metadata/evidence/` (ver `metadata/evidence/README.md`). Genera los snapshots así:

### 5.1 Conectores automáticos (APIs ligeras)
Usan la firma diferencial como universo de genes:

```bash
bash scripts/fetch_evidence_snapshots.sh
```

Esto ejecuta:

| Fuente | Script | Aporta | Notas |
|---|---|---|---|
| **DGIdb** | `fetch_evidence_dgidb.py` | druggability (nº fármacos) | API GraphQL, rápido |
| **cBioPortal** | `fetch_evidence_cbioportal.py` | frecuencia de mutación COADREAD | API REST, rápido |
| **Open Targets** | `fetch_evidence_opentargets.py` | tractabilidad (0–1) | per-gen; usa `OT_LIMIT=300` para acotar |
| **GTEx** | `fetch_evidence_gtex.py` | seguridad (expresión normal) | descarga ~50 MB (`WITH_GTEX=0` para omitir) |

Ejemplos de uso individual:
```bash
python scripts/fetch_evidence_dgidb.py --differential-tsv results/differential/cptac_tumor_vs_normal_diff.tsv
python scripts/fetch_evidence_cbioportal.py --genes KRAS TP53 APC EGFR   # prueba con pocos genes
OT_LIMIT=200 bash scripts/fetch_evidence_snapshots.sh
```

### 5.2 Conectores con descarga grande (manual)
- **DepMap** (esencialidad CRISPR): descarga `CRISPRGeneEffect.csv` y `Model.csv` de
  <https://depmap.org/portal/download/> y ejecuta:
  ```bash
  python scripts/fetch_evidence_depmap.py \
      --gene-effect-csv ~/Downloads/CRISPRGeneEffect.csv \
      --model-csv ~/Downloads/Model.csv \
      --out metadata/evidence/depmap_crc_gene_effect.tsv
  ```
- **HPA** (especificidad tumoral):
  ```bash
  python scripts/fetch_evidence_hpa.py \
      --differential-tsv results/differential/cptac_tumor_vs_normal_diff.tsv \
      --out metadata/evidence/hpa_tumor_specificity.tsv
  ```

> Si falta una fuente, **no penaliza** el score (solo baja `evidence_count`). Puedes correr la Capa C
> con las fuentes que tengas.

### 5.3 Ejecutar la Capa C
```bash
oncogo-crc run --stage prioritization
```

Salidas en `results/prioritization/`:

| Archivo | Contenido |
|---|---|
| `hypothesis_priority.tsv` | Tabla priorizada (rank, gen, `convergence_score`, tier, evidencias) |
| `convergence_evidence_matrix.tsv` | Detalle por señal (valor, favorabilidad, peso) |
| `prioritization_provenance.json` | Versión/ruta de cada fuente + pesos + parámetros |

Ajusta pesos y candidatos en `config.yaml → prioritization`.

---

## 6. Explorar los resultados en la web

### 6.1 Dashboard interactivo (Streamlit)
```bash
pip install -r webapp/requirements.txt
streamlit run webapp/streamlit_app.py
```
Si usas workspace externo:
```bash
export ONCOGO_CRC_RESULTS=$ONCOGO_CRC_DATA_WORKSPACE/results
streamlit run webapp/streamlit_app.py
```
Pestañas: **Resumen · Hipótesis (filtros + descarga) · Volcano · Red · Procedencia**.
Sin artefactos reales arranca en **modo demostración**.

### 6.2 Documentación estática (MkDocs)
```bash
pip install -e ".[docs]"
mkdocs serve      # http://127.0.0.1:8000
mkdocs build      # genera el sitio en site/
```

---

## 7. Cómo interpretar la tabla de hipótesis

- **`convergence_score` (0–1):** media ponderada de las señales **disponibles** (no penaliza ausencias).
- **`evidence_tier` (A–D):** cuartil del score (A = mayor convergencia).
- **`evidence_count` / `external_evidence_count`:** cuántas señales / cuántas fuentes externas respaldan.
- **Señales:** `diff_effect`, `diff_significance`, `reproducible`, `network_hub`,
  `depmap_essentiality`, `opentargets_tractability`, `dgidb_druggability`, `hpa_tumor_specificity`,
  `gtex_safety`, `cbioportal_genomic`.

**Criterio sugerido de lectura (no clínico):** prioriza candidatos con **alto score**, **`evidence_count`
alto** y **reproducibles** entre cohortes. Un candidato con score alto pero pocas evidencias es más
especulativo. Contrasta siempre con la literatura antes de proponer seguimiento experimental.

---

## 8. Reproducibilidad

- `results/run_manifest.json`: qué etapas corrieron, hashes y estado (evita mezclar corridas).
- `metadata/evidence/*.version`: versión/fecha de cada snapshot externo.
- `prioritization_provenance.json`: pesos y procedencia del score.

Para una corrida limpia nueva: reejecuta `oncogo-crc run --stage all` (el manifest se reinicia) y
regenera los snapshots si quieres actualizarlos.
