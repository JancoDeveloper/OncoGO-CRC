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

## 2026-11-06 16:36 (UTC−5). Fase 0, sello

Comando: `pytest tests/test_v3_2_sello.py` y `scripts/v3/congelar_manifiestos.py`.

Resultado. Tres tests pasan. `results/v3_manifest.sha256` resume 66 ficheros de `results/v3/` en el disco de datos. Commit `862502f`.

## 2026-11-06 16:37 (UTC−5). Fase 1, descarga del report

Comando: `curl` de `https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/06/PXD055218/DIA-NN_report.tsv` hacia `data/raw/plasma/DIA-NN_report.tsv`.

Resultado. La primera descarga se cortó a los 14 minutos con `curl` código 56 (tiempo de espera), en unos 50 MB de 203. Se reanudó con `curl -C -`. El parser nuevo está en `src/oncogo_crc/utils/dia_nn_plasma.py`. Los tests sintéticos pasan. El test con el report real se salta mientras el fichero no tenga 213 313 241 bytes. `download_real.py` ya no escribe encima del report. `plasma_crc_2024` tiene `overwrite: false` para no pisar la matriz ancha histórica. La etapa `v3_2_suero` está en el commit `590412e` y no se ha corrido: sin el report completo no hay cifra.

## 2026-11-06 17:00 (UTC−5). Fase 1 y 2, código antes de la cifra

Comando: `pytest tests/test_v3_2_filtros.py tests/test_v3_2_sello.py tests/test_v3_2_suero.py`.

Resultado. 8 pasan y 1 se salta (el report aún no tiene 213 313 241 bytes). El análisis S2–S9 está en `src/oncogo_crc/v3/suero_v32.py` y no se ejecuta hasta que el SHA-256 coincida. El tejido, que ya está en disco, se corre en la etapa `v3_2_tejido` con el filtro CON__/REV__ y 8 pares.

Decisión. No se interpretan queratinas ni el ranking exploratorio. La frase del suero sale de las clases de treat y TOST cuando el report esté completo.
