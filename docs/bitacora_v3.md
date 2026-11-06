# Bitácora v3.2

## 2026-11-06 16:33 (UTC−5). Fase 0, comprobaciones

Comando: `git checkout main`, `git pull --ff-only origin main`, `make manifest-verify`, `make regresion-v1`.

Resultado. `main` estaba un commit por delante de `origin/main` (`1c8a128`, el banco de preguntas). El pull no trajo nada. Manifiestos de entradas y de v1 coinciden. La regresión reproduce 688, 140, 28, 104/119, 226 y 307.

Decisiones, con el valor por defecto del encargo porque no hubo otra respuesta:

- No se descargan los 60 crudos de PXD055218 (65,6 GB). AlphaDIA queda como extensión documentada.
- No se despliega Firebase. La web se construye solo en local.

El entorno es Python 3.14.3, no 3.11. Queda anotado en `docs/auditoria_v3.md`.

## 2026-11-06 16:34 (UTC−5). Fase 0, enmienda

Comando: commit solo de `config/preregistro_v3.yaml` y `docs/enmiendas_v3.md`.

Resultado. Commit `cdce000`. Mensaje: `prereg: enmienda v3.2 (FDR de proteína, especie, contaminantes, paneles de calidad)`. No había ningún fichero en `results/v3_2/`.

Decisión. `min_pares_observados` de v3.2 es 8, el del cierre. La librería de contaminantes queda fijada por commit y SHA-256 antes de mirar el suero filtrado. La tabla de Geyer no se descargó (404) y los paneles son la lista del texto.
