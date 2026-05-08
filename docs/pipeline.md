# Pipeline por fases

| Fase | CLI (`oncogo-crc run --stage`) | Rol |
|------|----------------------------------|-----|
| A | `download` | Registro de rutas / manifest; en `--demo`, datos sintéticos. |
| B | `harmonize` | Armonización de IDs y colapso por gen. |
| C | `qc` | Missingness, log2, KNN, z-score por proteína. |
| D | `differential` | t-test + FDR por comparación declarada. |
| E | `go` | ORA (Fisher) + reducción semántica por solapamiento. |
| F | `network` | Red por correlación y hubs. |
| G | `ml` | L1-logit + Random Forest y ranking integrado. |
| H | `clinical` | Evidencia traslacional por niveles. |
| I | `report` | Resumen ejecutivo en Markdown. |

Ejecutar todo en orden: `oncogo-crc run --stage all` (o `make demo` con datos sintéticos).
