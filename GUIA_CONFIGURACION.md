# Guía de configuración OncoGO-CRC

Este documento explica **qué edita el usuario** y **qué efecto tiene** cada bloque de los YAML en `config/`.

## `config/config.yaml`

| Clave | Propósito |
|-------|-----------|
| `project.name` | Nombre legible del proyecto (aparece en informes). |
| `project.seed` | Semilla global para operaciones estocásticas (ML, demo). |
| `project.primary_dataset` | Clave en `datasets.yaml` usada para red, ML y traslación núcleo. |
| `paths.*` | Directorios de datos, resultados e informes (rutas relativas al repo). |
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
| `network.correlation_threshold` | Umbral absoluto de Pearson para crear aristas. |
| `network.min_degree_for_hub` | Grado mínimo para marcar un nodo como hub. |
| `ml.*` | Fracción de test, hiperparámetros de Random Forest y banderas de explicabilidad extendida. |
| `translation.evidence_tiers` | Nombres de los tres niveles de evidencia traslacional. |
| `demo.*` | Parámetros del generador sintético (`--demo`). |

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

No son obligatorias. El proyecto resuelve rutas relativas al directorio que contiene `config/config.yaml` (normalmente la raíz del repo).
