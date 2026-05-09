# OncoGO-CRC

**OncoGO-CRC** es un entorno reproducible de bioinformática para integrar proteómica pública de cáncer colorrectal (CCR), usar **Gene Ontology como eje central** de interpretación (no solo como capa descriptiva final) y priorizar hallazgos con **modelos de aprendizaje automático interpretables**, alineado con el anteproyecto y el reajuste técnico descritos en la carpeta `Informacion/`.

## Documentación (raíz del repositorio)

| Archivo | Contenido |
|-----------|-----------|
| [GUIA_EJECUCION_Y_REDACTION_TESIS.md](GUIA_EJECUCION_Y_REDACTION_TESIS.md) | Checklist operativo + mapa **pipeline → capítulos de la tesis**; redacción en paralelo al análisis. |
| [PLAN_EJECUCION_TESIS.md](PLAN_EJECUCION_TESIS.md) | Plan por fases (técnico + científico + material de tesis); ejecución incremental del pipeline. |
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

Equivalente shell: `bash scripts/run_pipeline.sh`.

## Equivalencia con la reconstrucción original

Los documentos en `Informacion/` proponían `python -m src.01_download_data`, etc. En Python **no es válido** importar módulos cuyo nombre comienza con dígitos; por ello este repositorio expone la CLI `oncogo-crc` y el módulo `python -m oncogo_crc`. La correspondencia de etapas está en [DOCUMENTACION_TESIS.md](DOCUMENTACION_TESIS.md#equivalencia-con-scripts-01_-).

## Sitio MkDocs

```bash
pip install mkdocs mkdocs-material
mkdocs serve
```

## Licencia

MIT (ver `LICENSE`).
