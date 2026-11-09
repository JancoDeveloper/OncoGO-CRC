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

## 2026-11-06 17:05 (UTC−5). Fase 2, tejido

Comando: `oncogo-crc run --stage v3_2_tejido --data-workspace /Volumes/JANCODATA/OncoGO-CRC`.

Resultado. `results/v3_2/k1_etiquetas.json`. Se excluyen 237 filas (CON__, REV__ y las marcas de MaxQuant). Con mínimo de 8 pares, primario frente a mucosa es 89/90, IC 0,940–1,000, Spearman 0,713, Lin 0,704. H1 sigue en `rechaza_h1`. Metástasis frente a mucosa: 69/82. No pareado: 95/98. Hígado: 13 pasan, 2 no evaluables, 1 falla. Mucosa: 16 no evaluables (pocos genes observados en el par). Cromosoma Y: no detectado en la matriz. Nadie se excluye en la sensibilidad, porque nadie falla dos de tres pruebas.

Decisión. No se reetiqueta ninguna muestra. El mínimo de 8 pares es el de la enmienda, no se movió al ver 89/90.

## 2026-11-06 17:15 (UTC−5). Fases 3 a 6, código listo para correr

Comando: `pytest tests/test_v3_2_resto.py tests/test_v3_2_filtros.py tests/test_v3_2_sello.py tests/test_v3_2_suero.py`.

Resultado. Quedan escritas las etapas `v3_2_controles`, `v3_2_hipotesis`, `v3_2_multiverso`, `v3_2_ml` y `v3_2_report`. H5 descarga LinkedOmicsKB y HGNC al correr. H7, H8 y H9 quedan declaradas si el fichero no está. No se ha lanzado la corrida: el report de PXD055218 todavía no tiene 213 313 241 bytes.

Decisión. La corrida la lanza el autor con `bash scripts/v3/correr_v3_2.sh` cuando el report esté completo. El sello exige el árbol limpio, así que este código se commitea antes.

## 2026-11-06 18:45 (UTC−5). El suero se detuvo en el control 307

Comando: `bash scripts/v3/correr_v3_2.sh`.

Resultado. Manifiestos y `regresion-v1` pasaron. El report tiene 213 313 241 bytes y el SHA-256 de la enmienda. La etapa `v3_2_suero` se detuvo: humanos con los dos q-valores son 312. De esos, 307 tienen gen. Los cinco sin gen son P0DOX2, P0DOX3, P0DOX5, P0DOX7 y P0DOX8.

Decisión. Gana el fichero. El control acepta 312 y guarda 307 como el subconjunto con gen. No se sacan esos grupos: el índice del análisis es `Protein.Group`.

## 2026-11-06 18:50 (UTC−5). El suero se cayó al contar precursores

Comando: `bash scripts/v3/correr_v3_2.sh`.

Resultado. K4 y el nulo de K5 terminaron. Después, `analizar_suero` buscó la columna `Run` en la tabla del último contraste, porque el bucle había reutilizado el nombre del report. No se escribió K7 ni el resumen. La sensibilidad por gen tampoco vio el report.

Decisión. El report y las matrices no se editan. Se separan los dos nombres y la etapa se vuelve a correr entera, para que el sello corresponda al código que produce las tablas.

## 2026-11-06. La corrida v3.2 terminó

Comando: `bash scripts/v3/correr_v3_2.sh`.

Resultado. El informe quedó en el disco, con `faltan []`. La frase de `results/v3_2/suero_resumen.json` es «No se detecta diferencia; 105 proteínas son equivalentes». Matriz principal, 205 proteínas: 0 diferentes, 105 equivalentes, 100 indeterminadas entre cáncer y colonoscopia negativa. Humanos con los dos q-valores: 312. Tejido: 89/90, H1 en `rechaza_h1`. Sello del suero: commit `293e6ee`.

Decisión. Esas tablas se copiaron al repositorio en el commit `35fad60`. Los números de `results/v3/` no se reescriben.

## 2026-11-09 14:25 (UTC−5). v3.3, Fase 0, congelar v3.2

Comando: `git checkout main`, `git pull --ff-only origin main`, `make regresion-v1`, `python scripts/v3/congelar_manifiestos.py --carpeta results/v3_2 --salida results/v3_2_manifest.sha256`, `git tag -a v3.2-final`.

Resultado. `main` coincide con `origin/main` en `c2be8e4`. El árbol estaba limpio. La regresión pasó: 688, 140, 28, 104/119, 226 y 307. El SHA-256 del report es `A8A5AADD74391AFEF68CF4D6C3A7D3860023101E5C92590E2E87C152A39E59E8`. El disco del sistema tiene 14 GB libres y el externo, 464 GB. El manifiesto de `results/v3_2/` queda en `results/v3_2_manifest.sha256`. La etiqueta `v3.2-final` apunta a `c2be8e4`, antes de este commit. La cronología está en `docs/cronologia_v3.md`.

Decisión del autor, el 9 de noviembre de 2026:

- Sí se despliega Firebase, solo si el checklist de la sección 16 queda en verde.
- Sí se descargan los pesos de ProtT5-XL, en el disco externo. El disco del sistema no tiene 20 GB libres.
- Sí se descargan los 60 crudos de PXD055218 al disco externo. AlphaDIA sigue sin ejecutarse en esta pasada: DIA-NN no se corre y los crudos no entran en el cálculo de v3.3.

El mensaje `35fad60` empieza por `feat:`. El encargo lo contó sin prefijo. Gana el mensaje real y queda anotado en `docs/auditoria_v3.md` y en la cronología. `50caef8` y `c2be8e4` sí se publicaron sin prefijo y no se reescriben.

## 2026-11-09 14:40 (UTC−5). v3.3, Fase 1, enmienda

Comando: edición de `config/preregistro_v3.yaml` y `docs/enmiendas_v3.md`, y `yaml.safe_load`.

Resultado. El YAML parsea. Commit `0edf54a`. No hay ficheros en `results/v3_3/`.

Decisión. La enmienda queda antes del código de cálculo. Los umbrales no se mueven.

## 2026-11-09 14:50 (UTC−5). v3.3, renombrar etapas

Comando: `git mv` de `stages/v3_2_*.py` y `tests/test_v3_2_*.py` a `v3_3`, y `pytest tests/test_v3_3_suero.py tests/test_v3_3_filtros.py tests/test_v3_3_sello.py tests/test_v3_3_resto.py tests/test_v3_3_congelado.py tests/test_cronologia.py`.

Resultado. 19 tests pasan. El registro ya no tiene etapas `v3_2_*`. La salida vigente es `results/v3_3/`. `salida_v32.py` se queda, sin uso. Los 112 ficheros de `results/v3_2/` coinciden con el manifiesto.

Decisión. El Makefile de v3.3 se añade después de los parches, como pide el encargo. La descarga de los 60 crudos sigue en el disco externo y no entra en este commit.

## 2026-11-09 15:05 (UTC−5). v3.3, contaminantes por acceso

Comando: `pytest tests/test_v3_3_contaminantes.py` antes y después del parche.

Resultado. Antes, el test no se recoge: `ImportError`, no existe `acceso_base`. Después, 9 tests pasan, incluido el de la librería real: 381 entradas, 151 humanas (`OX=9606`) y ningún gen de plasma clásico en esa lista.

Decisión. El cruce por símbolo queda retirado. KRT7, KRT8, KRT18, KRT19 y KRT20 siguen fuera del principal y vuelven en `con_queratinas_epiteliales`. Los grupos mixtos se quedan y salen en `sin_grupos_mixtos`. Si ALB, APOA1, C3, TF, SERPINA1, A2M, APOE, FGA, FGB, FGG o LRG1 salen como contaminante, la etapa se detiene.

## 2026-11-09 15:15 (UTC−5). v3.3, ML sin muestras duplicadas

Comando: `pytest tests/test_v3_3_ml.py::test_cada_muestra_entra_una_vez` antes del parche, y `pytest tests/test_v3_3_ml.py tests/test_v3_3_resto.py` después.

Resultado. Antes, el test falla: `TypeError`, `matrices_de_suero` exigía `prereg` y recorría las columnas de las dos matrices, así que cada `CP_` y cada `CN_` entraba dos veces. Después, 6 tests pasan. Cada modelo tiene 40 filas y 20 de la clase cáncer. La dispersión entre pliegues ya no se llama intervalo. El intervalo que se informa es el de LeDell sobre las predicciones fuera del pliegue.

Decisión. El titular no se escribe aquí. Sale de la regla, cuando corra la etapa.

## 2026-11-09 15:25 (UTC−5). v3.3, mapa de genes del tejido

Comando: `pytest tests/test_v3_3_genes_tanaka.py` antes de cablear `_genes_bajos`, y otra vez después.

Resultado. Antes, `test_control_positivo_no_evaluable_es_falla` falla: el código viejo escribía `pasa=no_evaluable` y la etapa seguía. El mapa sintético ya encontraba CA1, AQP8 por `GN=` y SLC26A3 por el acceso. Después, 7 tests pasan, con `tests/test_v3_3_resto.py`.

Decisión. AQP8 y GUCA2A no están en el proteinGroups. Queda escrito en `docs/auditoria_v3.md`. No se inventa una fila para que el control pase.

## 2026-11-09 15:40 (UTC−5). v3.3, normalización del tejido

Comando: `pytest tests/test_v3_3_normalizacion.py`.

Resultado. Con un desplazamiento constante de +1 en las columnas P, la mediana del log2FC pareado queda en 0 después de centrar. El análisis principal no se centra. La sensibilidad escribe `k1_sensibilidad_normalizacion.json`. Si R1 cruza 0,70, H1 queda frágil a la normalización.

Decisión. La mediana de las diferencias no es la diferencia de las medianas. El test usa un desplazamiento constante, que es el caso que describe la enmienda.

## 2026-11-09 15:50 (UTC−5). v3.3, H5 y H10

Comando: `pytest tests/test_v3_3_h5_h10.py tests/test_v3_3_resto.py`.

Resultado. 10 tests pasan. Con 7 carcinomas la clase coincide con v3. Con 6 y 4 concordantes sale pan_carcinoma. Con 5 sale no_clasificable. H10 lee `Healthy`, `Class` y `All other diseases`. El test viejo, que daba por específica la primera fila, ahora da `no_medida` porque faltan los tres controles.

Decisión. UCEC no se reetiqueta: sus normales se llaman NX1, NX2 y no comparten identificador con el tumor. H10 queda descriptiva.

## 2026-11-09 16:10 (UTC−5). v3.3, H7

Comando: `pytest tests/test_v3_3_h7.py`.

Resultado. 4 tests pasan. Con etiquetas aleatorias el veredicto no es sostiene. Con una variable que determina la etiqueta, sí. La PCA se ajusta dentro del tubo y no ve las proteínas de prueba. Con menos de 10 por clase, el intervalo queda no_evaluable.

Decisión. Papier 2024 no está en el disco: H7 usa Álvez 2023 y el atlas del HPA, y la regla de dos recursos exige los dos. ESM-2 y TMbed no entran en el modelo principal. La rejilla de C empieza en 0,1: el entero 10 de scikit-learn empieza en 10⁻⁴ y, si todas empatan en AUC, deja las probabilidades en 0,5.

## 2026-11-09 16:20 (UTC−5). v3.3, Makefile y CI

Comando: `ruff check src tests`.

Resultado. Pasa, después de quitar imports sin uso. La CI corre en macOS 14 y Ubuntu 24.04, con Python 3.12 y 3.14, y con la historia completa de git para el test de cronología.

Decisión. `correr_v3_3.sh` sigue si los controles del tejido fallan, y lo anota. AQP8 y GUCA2A no están en el proteinGroups: eso no bloquea el suero.

## 2026-11-09 15:40 (UTC−5). v3.3, fuentes y certificado

Comando: `pytest tests/test_v3_3_cierre.py`.

Resultado. 3 tests pasan. H1, H5, H7 y H10 quedan con fuente `results/v3_3/`. El certificado copia el titular, H5, H7 y H10. El changelog anota solo la cifra que cambia y le pone la causa 2.1 a 2.7.

Decisión. La corrida de las 15:25 escribió la fuente como `results/v3_2/` aunque el número ya era el de v3.3. Se vuelve a correr hipótesis e informe con el árbol limpio. No se reescriben a mano los JSON.

## 2026-11-09 15:45 (UTC−5). v3.3, resultados

Comando: `make v3-3-hipotesis`, `make v3-3-report` y `python scripts/v3/demostrar_fuga_ml.py`.

Resultado. `results/v3_3/certificado_v3_3.json` y `changelog_numeros.tsv`. Suero principal: 284 proteínas, 0 diferentes, 166 equivalentes, 118 indeterminadas. Titular: `separacion_explicada_por_calidad`. H1 `rechaza_h1`, 147/147, no frágil a la normalización. H5 `refuta`, fracción pan-carcinoma 0,1186. H7 `refuta`. H10 descriptiva, 5 específicas. CA1 pasa. AQP8, SLC26A3 y GUCA2A fallan. El ruido sin duplicar da AUC 0,475 y, duplicado, 0,93125.

Decisión. Esos ficheros entran en el repositorio. No se publican: los tres controles del tejido no pasan, el multiverso completo no está calculado y ESM-2 y TMbed no corrieron.

## 2026-11-09 15:55 (UTC−5). v3.3, guion de figuras

Comando: `ruff check scripts/v3/figuras_v3_3.py`.

Resultado. El guion lee `results/v3_3/` y escribe F1 a F14 en PNG y SVG. El título de F1 es el texto del titular. F7 dice que AQP8, SLC26A3 y GUCA2A no tienen log2FC.

Decisión. El guion se commitea antes de firmar las figuras, porque el sello exige árbol limpio.

## 2026-11-09 16:05 (UTC−5). v3.3, figuras

Comando: `python scripts/v3/figuras_v3_3.py`.

Resultado. Catorce figuras en `results/v3_3/report/figuras/`, `sitio/figuras/` y `tesis-ucaldas/latex/figures/`. F7 no dibuja un log2FC en cero cuando el control no lo tiene. F1 usa el texto del titular.

Decisión. Entran las figuras. El título no dice que las etiquetas del tejido quedaron confirmadas.
