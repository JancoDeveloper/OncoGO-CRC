# Guía de configuración OncoGO-CRC

Este documento explica **qué edita el usuario** y **qué efecto tiene** cada bloque de los YAML en `config/`.

## `config/config.yaml`

| Clave | Propósito |
|-------|-----------|
| `project.name` | Nombre legible del proyecto (aparece en informes). |
| `project.seed` | Semilla global para operaciones estocásticas (ML, demo). |
| `project.primary_dataset` | Clave en `datasets.yaml` usada para red, ML y traslación núcleo. |
| `paths.*` | Directorios de datos, resultados e informes (rutas relativas al árbol de trabajo; ver abajo). |
| `paths.data_workspace` | Directorio base opcional (ruta absoluta recomendada) donde residen `data/`, `results/`, `reports/`, **`metadata/`** y la salida **`docs/web`** cuando quieres todo el volumen de trabajo fuera del clon. Vacío = esas rutas bajo la raíz del repo. Cada arranque de `oncogo-crc` con workspace activo **sincroniza** desde el repo hacia `<workspace>/metadata/` las carpetas `protein_id_maps/`, `gene_ontology/` y `tcga_auxiliary/` (copia si falta o si el archivo del repo es **más reciente** por fecha de modificación). `sample_metadata/` y `clinical_layers/` solo se rellenan si **faltan** en el workspace (no se sobrescriben harmonizaciones ya guardadas allí). |
| `preprocessing.missing_threshold_*` | Umbral máximo de NA por proteína / muestra antes de excluir. |
| `preprocessing.log2_transform` | Aplicar transformación log2 antes de imputación. |
| `preprocessing.imputation_method` | En esta versión: `knn` (vecinos en `knn_neighbors`). |
| `preprocessing.normalize_method` | `zscore_protein`: estandariza cada fila (proteína) a media 0 y varianza 1 entre muestras. |
| `preprocessing.batch_correction` | Reservado (`none` por defecto); COMBAT/pyCombat puede añadirse en extensiones. |
| `differential.alpha` | Nivel α para significancia tras corrección FDR. |
| `differential.fdr_method` | Método `statsmodels` (por defecto Benjamini–Hochberg). |
| `go.fdr_threshold` | Umbral de enriquecimiento GO a nivel de término. |
| `go.ontologies` | Subconjuntos BP/MF/CC en la tabla larga `go_gene_long.tsv`. |
| `go.semantic_reduction` | Activa agrupación de términos redundantes por solapamiento de genes. |
| `go.ora_foreground_max_genes` | Tope de genes significativos enviados a ORA (ordenados por `padj`). Si casi todo el proteoma es significativo, ORA se diluye; valores típicos 500–1200. `0` = sin tope. |
| `go.min_genes_per_term` | Tamaño mínimo del término en el universo (genes anotados al término) para entrar al test de Fisher. |
| `network.correlation_threshold` | Umbral absoluto de Pearson para crear aristas. |
| `network.min_degree_for_hub` | Grado mínimo para marcar un nodo como hub. |
| `ml.models` | Lista: `lasso` (SGD L1), `elastic_net` (SGD elasticnet), `random_forest`. |
| `ml.linear_alpha` | Fuerza de regularización por defecto para L1/EN en hold-out y transferencia. |
| `ml.feature_stability_bootstrap_n` | Si ≥10, escribe frecuencia bootstrap de genes top‑k en L1. |
| `ml.nested_cv.only_rf` | Si `true`, nested CV solo sobre **RandomForest** (más rápido). Ignorado si `nested_cv.models` no está vacío. |
| `ml.nested_cv.models` | Lista opcional de modelos nested: alias `l1`, `elastic_net`, `rf` (o nombres largos). |
| `ml.nested_cv.rf_n_estimators` | Número de árboles **solo** en el grid nested (por defecto igual a `random_forest_n_estimators` si se omite). |
| `ml.external_transfer.enabled` | Si `true`, entrena en `train_dataset` y evalúa métricas en `test_dataset` (genes alineados). |
| `ml.external_transfer.calibration_method` | `sigmoid` (Platt), `isotonic`, o vacío/`none` para desactivar calibración en transferencia. |
| `ml.external_transfer.calibration_cv` | Pliegues internos de `CalibratedClassifierCV` sobre CPTAC (se recorta si hay pocas muestras/clase). |
| `external_evidence.tcga_rnaseq_log2fc_tsv` | Ruta a tu TSV (vacío = no ejecutar concordancia RNA-seq). Plantilla de columnas: `metadata/tcga_auxiliary/gene_rnaseq_log2fc_tcga_template.tsv`. |
| `external_evidence.survival_gene_cox_tsv` | Ruta a tu TSV Cox por gen (vacío = omitir). Plantilla: `metadata/tcga_auxiliary/gene_survival_cox_template.tsv`. |
| `external_evidence.enabled` / `merge_translation` | Matriz `hit_*` desde `crc_context_layers.yaml`; fusión con traslación clínica. |
| `translation.evidence_tiers` | Nombres de los tres niveles de evidencia traslacional. |
| `demo.*` | Parámetros del generador sintético (`--demo`). |

### Almacenamiento en otro disco

<span id="data-workspace"></span>

Configuración YAML: `paths.data_workspace`. Si el volumen del repo no tiene espacio, define un directorio externo con la **misma estructura lógica** que bajo el repo (`data/raw`, …, `results`, `reports`, **`metadata/protein_id_maps`**, **`metadata/gene_ontology`**, **`metadata/tcga_auxiliary`**, `metadata/sample_metadata`, `metadata/clinical_layers`, `docs/web`, etc.). El pipeline creará subcarpetas al escribir.

**Sincronización al iniciar el CLI:** con workspace activo, el código copia o **actualiza** (si el repo lleva una versión más nueva) los archivos bajo `metadata/protein_id_maps/`, `metadata/gene_ontology/` y `metadata/tcga_auxiliary/`. Así, tras `git pull` o `python scripts/fetch_uniprot_human_reviewed_gene_map.py` en el clon, la siguiente corrida arrastra el mapa bulk al volumen externo. Las rutas `sample_metadata/` y `clinical_layers/` del repo solo se copian cuando **no existen** en el workspace (evita machacar `cptac_metadata.tsv` u otros TSV ya harmonizados en el disco grande).

**Qué no se mueve:** `config/`, el código en `src/`, y el resto de `docs/` (p. ej. `docs/pipeline.md`) siguen en la raíz del repositorio. Se redirigen a `data_workspace` las rutas bajo `data/`, `results/`, `reports/`, **`metadata/`** y **`docs/web`**. Tras `publish_results`, el proyecto **copia** `docs/web` de vuelta al clon para que `mkdocs build` siga encontrando las páginas del bloque *Proyecto Web* sin cambiar `mkdocs.yml`.

**Checklist — otra corrida en el mismo workspace:** no hace falta borrar el volumen; vuelve a ejecutar `oncogo-crc run --stage all` (se regenera `results/run_manifest.json` al inicio). Si cambiaste solo YAML en el repo, no toques el workspace salvo que quieras forzar re-descarga (`download_real.overwrite`, etc.).

**Prioridad** (la primera no vacía gana): argumento CLI `--data-workspace` → variable de entorno `ONCOGO_CRC_DATA_WORKSPACE` → `paths.data_workspace` en `config.yaml`.

Ejemplos:

```bash
export ONCOGO_CRC_DATA_WORKSPACE=/mnt/datadisk/oncogo_crc_workspace
oncogo-crc run --stage all
```

```bash
oncogo-crc run --data-workspace /mnt/datadisk/oncogo_crc_workspace --stage all
```

## `config/datasets.yaml`

Cada entrada bajo `datasets` declara:

- `matrix_path`: TSV con identificadores en la primera columna (`protein_id` por defecto).
- `sample_metadata`: tabla con `sample_id` y columnas fenotípicas usadas en `comparisons.yaml`.
- `id_column` / `protein_column`: nombres de columnas si se desvían del estándar.

`protein_id_map` apunta al mapa maestro de identificadores.

### Modo `--demo`

No es necesario editar los YAML: `download` crea matrices y escribe `data/interim/_demo_<dataset>_samples.tsv`. `settings.sample_metadata_path` prioriza esos archivos cuando `demo=True`.

## `config/comparisons.yaml`

Cada comparación define:

- `dataset`: clave de cohorte.
- `condition_column`: columna en la metadata (p. ej. `condition`).
- `group_a` / `group_b`: valores categóricos contrastados.
- `optional: true`: se omite por defecto (evita contrastes sin estratos reales).

## Archivos de metadatos

- `metadata/sample_metadata/*.tsv`: plantillas con columnas mínimas; amplíelas según CPTAC/PRIDE reales.
- `metadata/protein_id_maps/protein_mapping_master.tsv`: columnas `protein_id`, `gene_symbol` (requerido para GO), `uniprot_accession` (opcional).
- `metadata/gene_ontology/go_gene_long.tsv`: anotación larga **offline** para ORA sin depender de APIs externas.

## Variables de entorno

No son obligatorias salvo tokens para proveedores remotos (p. ej. `PDC_TOKEN`). El proyecto resuelve rutas relativas al directorio que contiene `config/config.yaml` (normalmente la raíz del repo). Con `data_workspace` activo, las rutas pesadas se resuelven bajo ese directorio y, al arrancar, se sincroniza metadata desde el repo según la sección *Almacenamiento en otro disco*.

| Variable | Uso |
|----------|-----|
| `ONCOGO_CRC_DATA_WORKSPACE` | Ruta absoluta (o relativa al repo) del directorio de trabajo pesado; equivale a `paths.data_workspace` y puede anularse con `--data-workspace` en la CLI. |
