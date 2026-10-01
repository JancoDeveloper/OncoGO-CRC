# OncoGO-CRC

**OncoGO-CRC** es un entorno reproducible de bioinformática para integrar proteómica pública de cáncer colorrectal (CCR), usar **Gene Ontology como eje central** de interpretación (no solo como capa descriptiva final) y priorizar hallazgos con **modelos de aprendizaje automático interpretables**, alineado con el anteproyecto y el reajuste técnico descritos en la carpeta `Informacion/`.

## Documentación (raíz del repositorio)

| Archivo | Contenido |
|-----------|-----------|
| [docs/GUIA_ESTUDIANTE_CERTIFICACION.md](docs/GUIA_ESTUDIANTE_CERTIFICACION.md) | **Empieza aquí si vas a escribir la tesis.** Qué se busca, conceptos desde cero, diagramas, qué se actualizó y el orden de escritura de los capítulos. |
| [GUIA_EJECUCION_Y_REDACTION_TESIS.md](GUIA_EJECUCION_Y_REDACTION_TESIS.md) | Checklist operativo + mapa **pipeline → capítulos de la tesis**; redacción en paralelo al análisis. |
| [PLAN_EJECUCION_TESIS.md](PLAN_EJECUCION_TESIS.md) | Plan por fases (técnico + científico + material de tesis); ejecución incremental del pipeline. |
| [CATALOGO_DATASETS_Y_ESTADO.md](CATALOGO_DATASETS_Y_ESTADO.md) | Fuente, contenido, utilidad y estado funcional (remoto) de cada dataset. |
| [DOCUMENTACION_TESIS.md](DOCUMENTACION_TESIS.md) | Marco científico, fases, archivos generados, limitaciones y cierre para tesis. |
| [GUIA_CONFIGURACION.md](GUIA_CONFIGURACION.md) | Parámetros YAML y convención de rutas. |
| [ESTRUCTURA_DEL_PROYECTO.md](ESTRUCTURA_DEL_PROYECTO.md) | Árbol de carpetas y responsabilidad de cada módulo Python. |

## Requisitos

- Python 3.10+ (recomendado 3.11; Conda opcional vía `environment.yml`).

## Instalación y demo reproducible

```bash
python3 -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -e .
oncogo-crc run --demo --stage all
```

Con `--demo` se generan matrices sintéticas tipo CPTAC/plasma, metadatos coherentes y anotaciones GO mínimas para validar el flujo **sin descargar datos reales**.

## Ejecución con datos propios

1. Coloque las matrices TSV en las rutas declaradas en `config/datasets.yaml` (primera columna `protein_id`, resto `sample_id` numéricos).
2. Complete `metadata/sample_metadata/*.tsv` y el mapa `metadata/protein_id_maps/protein_mapping_master.tsv`.
3. Opcional: amplíe `metadata/gene_ontology/go_gene_long.tsv` (formato largo: `gene_symbol`, `go_id`, `namespace`).
4. Ejecute:

```bash
oncogo-crc run --stage all
```

Si el disco del repo es pequeño, redirija `data/`, `results/`, `reports/`, `metadata/` y la salida `docs/web` a otro volumen (misma jerarquía de carpetas). Prioridad: `--data-workspace` → `ONCOGO_CRC_DATA_WORKSPACE` → `paths.data_workspace` en `config/config.yaml`. Detalle en [GUIA_CONFIGURACION.md — almacenamiento externo](GUIA_CONFIGURACION.md#data-workspace).

```bash
oncogo-crc run --data-workspace /ruta/al/volumen/grande --stage all
```

Equivalente shell: `bash scripts/run_pipeline.sh` (respeta `ONCOGO_CRC_DATA_WORKSPACE` si está exportada antes de llamar al script).

## Equivalencia con la reconstrucción original

Los documentos en `Informacion/` proponían `python -m src.01_download_data`, etc. En Python **no es válido** importar módulos cuyo nombre comienza con dígitos; por ello este repositorio expone la CLI `oncogo-crc` y el módulo `python -m oncogo_crc`. La correspondencia de etapas está en [DOCUMENTACION_TESIS.md](DOCUMENTACION_TESIS.md#equivalencia-con-scripts-01_-).

## Descarga real vía API/URL (opcional)

La etapa `download_real` permite descargar matrices/metadatos desde endpoints HTTP definidos por dataset en `config/datasets.yaml` (bloque `download_real`), con `timeout`, `retries` y `overwrite` configurables en `config/config.yaml`.
Si `matrix_url` o `metadata_url` están vacíos, intenta **autodescubrimiento** desde:
- PRIDE (`provider: pride` + `source_accession`)
- `source_catalog_url` (scraping best-effort de enlaces)

También soporta placeholders de entorno en headers/query, por ejemplo:
- `Authorization: "Bearer ${PDC_TOKEN}"`

```bash
# 1) habilitar en config/config.yaml: download_real.enabled: true
# 2) CPTAC primario: PDC000116 vía GraphQL en config/datasets.yaml (`pdc_study_submitter_id`), o bien:
#    pip install 'oncogo-crc[ingest]' && oncogo-crc run --stage ingest_cptac
# 3) PRIDE / plasma: accesiones PXD en datasets.yaml (sin URLs landing).
oncogo-crc run --stage download_real
oncogo-crc run --stage harmonize_public_metadata
oncogo-crc run --stage version_inputs
oncogo-crc run --stage download
```

Tras `download_real`, **`harmonize_public_metadata`** cruza el mapping TMT10 + biospecimen PDC para CPTAC, recorta la matriz PRIDE a columnas LFQ con `protein_id`, y solo aplica el harmonizador péptido/ratios DIA si la capa líquida incluye **PXD036594**. La capa `plasma_crc_2024` declara **PXD015163** (cohorte plasma Kopylov 2020) y descarga la matriz DIA-NN de **PXD055218** con pivote automático a ancho + `plasma_metadata.tsv` (CP\_* = crc, Normal\_* = healthy); ver comentarios en `config/datasets.yaml`.

La etapa **ML** puede usar **SHAP** (`TreeExplainer` sobre el Random Forest) y **validación cruzada** ROC-AUC 5-fold (`config.yaml` → `ml.use_shap`, `ml.cv_folds`; SHAP requiere `pip install -e ".[explain]"`).

Para reproducibilidad CPTAC ∩ PRIDE: etapa `cohort_overlap` (tras `differential`) o el script `python3 scripts/differential_overlap_cptac_pride.py`.

`ingest_cptac` no se ejecuta con `--stage all`; convócalo explícitamente si quieres la vía rápida `cptac.Colon()` (documenta versión en `data/raw/cptac/cptac_ingest_provenance.json`).

La traza de esta etapa queda en:
- `data/raw/download_real_manifest.tsv`
- `data/raw/download_real_manifest.json`
- `data/versions/matrix_versions.tsv`
- `reports/executive_summary/matrix_inventory.md`

## Validación y publicación web

```bash
oncogo-crc run --stage validate_inputs
oncogo-crc run --stage publish_results
mkdocs build
```

Salidas principales:
- `results/validation/input_validation_report.json`
- `results/validation/input_validation_report.tsv`
- `reports/web/*.md` y copia en `docs/web/*.md` (MkDocs), incluyendo **`thesis_synthesis.md`** (guía de síntesis), **`results_base.md`** (tablas de proteínas, GO, ML) e **`results_inventory.md`**
- `reports/executive_summary/results_summary.md` (Markdown con extractos + **resumen cuantitativo** FDR)
- Prompt para redactar conclusiones con un LLM: **`docs/PROMPT_SINTESIS_TESIS.md`**
- Integración con enriquecimiento visual tipo **ProteinsLanguage**: guía **`INTEGRACION_PROTEINSLANGUAGE.md`** y script **`scripts/export_protein_context_bundle.py`**

## Sitio MkDocs

```bash
pip install mkdocs mkdocs-material
mkdocs serve
```

## Licencia

MIT (ver `LICENSE`).
