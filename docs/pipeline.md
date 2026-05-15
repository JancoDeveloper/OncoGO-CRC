# Pipeline por fases

| Fase | CLI (`oncogo-crc run --stage`) | Rol |
|------|----------------------------------|-----|
| A | `download` | Registro de rutas / manifest; en `--demo`, datos sintéticos. |
| B | `harmonize` | Armonización de IDs y colapso por gen. |
| C | `qc` | Missingness, log2, KNN, z-score por proteína. |
| D | `differential` | t-test + FDR por comparación declarada. |
| E | `go` | ORA (Fisher) + reducción semántica por solapamiento. |
| F | `network` | Red por correlación y hubs. |
| G | `ml` | L1 / Elastic Net (SGD) + RF; nested CV opcional; transfer CPTAC→PRIDE; bootstrap L1 opcional. |
| H | `clinical` | Evidencia traslacional por niveles. |
| H2 | `external_evidence` | Matriz `hit_*` + TCGA/supervivencia (TSV opcionales). |
| I | `report` | Resumen ejecutivo en Markdown. |

Ejecutar todo en orden: `oncogo-crc run --stage all` (o `make demo` con datos sintéticos).

Guía ampliada con **datasets**, **lineage de datos** y **diagramas Mermaid** (secuencia + flujo): [Proceso y datos (detallado)](proceso_datos_y_pipeline.md).
