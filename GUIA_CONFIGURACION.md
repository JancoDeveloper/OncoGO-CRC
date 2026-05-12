# Guía de configuración OncoGO-CRC

Este documento explica **qué edita el usuario** y **qué efecto tiene** cada bloque de los YAML en `config/`.

## `config/config.yaml`

| Clave | Propósito |
|-------|-----------|
| `project.name` | Nombre legible del proyecto (aparece en informes). |
| `project.seed` | Semilla global para operaciones estocásticas (ML, demo). |
| `project.primary_dataset` | Clave en `datasets.yaml` usada para red, ML y traslación núcleo. |
| `paths.*` | Directorios de datos, resultados e informes (rutas relativas al árbol de trabajo; ver abajo). |
| `paths.data_workspace` | Directorio base opcional (ruta absoluta recomendada) donde residen `data/`, `results/`, `reports/`, **`metadata/`** y la salida **`docs/web`** cuando quieres todo el volumen de trabajo fuera del clon. Vacío = esas rutas bajo la raíz del repo. |
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
| `ml.*` | Fracción de test, hiperparámetros de Random Forest y banderas de explicabilidad extendida. |
| `translation.evidence_tiers` | Nombres de los tres niveles de evidencia traslacional. |
| `demo.*` | Parámetros del generador sintético (`--demo`). |

### Almacenamiento en otro disco

<span id="data-workspace"></span>

Configuración YAML: `paths.data_workspace`. Si el volumen del repo no tiene espacio, define un directorio externo con la **misma estructura lógica** que bajo el repo (`data/raw`, …, `results`, `reports`, **`metadata/sample_metadata`**, **`metadata/protein_id_maps`**, **`metadata/gene_ontology`**, `docs/web`, etc.). El pipeline creará subcarpetas al escribir. Al iniciar el pipeline con `data_workspace` activo, si falta algún archivo bajo `metadata/protein_id_maps`, `metadata/gene_ontology` o `metadata/sample_metadata` en el volumen externo, se **copian desde el repo solo los que faltan** (no se pisan tablas ya generadas, p. ej. tras `harmonize_public_metadata`). También puedes sincronizar a mano con `rsync` si lo prefieres.

**Prioridad** (la primera no vacía gana): argumento CLI `--data-workspace` → variable de entorno `ONCOGO_CRC_DATA_WORKSPACE` → `paths.data_workspace` en `config.yaml`.

**Qué no se mueve:** `config/`, el código en `src/`, y el resto de `docs/` (p. ej. `docs/pipeline.md`) siguen en la raíz del repositorio. Se redirigen a `data_workspace` las rutas bajo `data/`, `results/`, `reports/`, **`metadata/`** y **`docs/web`**. Tras `publish_results`, el proyecto **copia** `docs/web` de vuelta al clon para que `mkdocs build` siga encontrando las páginas del bloque *Proyecto Web* sin cambiar `mkdocs.yml`.

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

No son obligatorias salvo tokens para proveedores remotos (p. ej. `PDC_TOKEN`). El proyecto resuelve rutas relativas al directorio que contiene `config/config.yaml` (normalmente la raíz del repo).

| Variable | Uso |
|----------|-----|
| `ONCOGO_CRC_DATA_WORKSPACE` | Ruta absoluta (o relativa al repo) del directorio de trabajo pesado; equivale a `paths.data_workspace` y puede anularse con `--data-workspace` en la CLI. |
