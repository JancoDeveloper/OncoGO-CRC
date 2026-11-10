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

## 2026-11-09 16:15 (UTC−5). v3.3, página

Comando: `python scripts/v3/exportar_sitio.py`.

Resultado. `sitio/index.html` y `sitio/datos/v3_3.json` salen del certificado. La frase vigente es la de 166 equivalentes. La tabla de cambios conserva el 105 y el 0,996 en la columna de v3.2. El simulador aplica Hanley y McNeil a números que escribe el lector.

Decisión. No se despliega. La página local ya no presenta la corrida v3.2 como vigente.

## 2026-11-09 16:20 (UTC−5). v3.3, cifras de la tesis

Comando: `python scripts/v3/exportar_cifras_tesis.py`.

Resultado. `tesis-ucaldas/latex/cifras.tex`, `chapters/09-v33.tex`, `word/cifras.json` y `docs/resumen_director_v3_3.md`. El capítulo de resultados deja la corrida de 205 proteínas como v3.2.

Decisión. No se regenera el Word ni el PDF. El resumen y el abstract siguen por encima de 250 palabras. No hay 18 diapositivas nuevas ni el banco de 50 preguntas.

## 2026-11-09 16:45 (UTC−5). v3.3, resumen, diapositivas y preguntas

Comando: `python scripts/v3/escribir_documentacion_v3_3.py`.

Resultado. El resumen queda en 243 palabras y el abstract en 221. Hay 18 diapositivas y 52 preguntas. Las cifras salen de `results/v3_3/`. El Word no se reescribe: el constructor parte de la plantilla y reemplaza el documento entero.

Decisión. Esta documentación se puede cerrar. El checklist de publicación no: siguen abiertos los tres controles del tejido, el multiverso fuera de la rejilla K4 y los modelos de secuencia.

## 2026-11-09 16:50 (UTC−5). v3.3, Word

Comando: `python tesis-ucaldas/build_word.py`.

Resultado. El documento se reconstruye desde la plantilla y los capítulos. La sección «Corrida v3.3» entra con 284, 166 y 147/147. Las macros de `cifras.tex` se sustituyen antes de quitar los comandos. No hay LuaLaTeX en el equipo, así que el PDF no se compiló.

Decisión. El Word se regenera. El PDF/A y Firebase siguen esperando el checklist.

## 2026-11-09 16:22 (UTC−5). v3.4, Fase 0

Comando: `git checkout main`, `git fetch origin`, `git pull --ff-only origin main`, `python scripts/v3/congelar_manifiestos.py --carpeta results/v3_3 --salida results/v3_3_manifest.sha256`, `pytest tests/test_cronologia.py tests/test_v3_3_congelado.py`.

Resultado. `main` en `4f85911`, limpio y al día con `origin/main`. El SHA-256 del report DIA-NN es `A8A5AADD74391AFEF68CF4D6C3A7D3860023101E5C92590E2E87C152A39E59E8`. El manifiesto resume 189 ficheros en `results/v3_3_manifest.sha256`. La etiqueta `v3.3-final` apunta a `4f85911`. En el disco externo hay 460 GB libres y 4 crudos de PXD055218, no los 60.

Decisiones, con el valor por defecto del encargo:

- Firebase sí, solo si el checklist de la sección 18 queda en verde.
- El PDF se entrega como zip de Overleaf. No se instala MacTeX.
- ESM-2 y TMbed sí: el disco externo tiene más de 20 GB libres.
- El orden de adquisición se intenta solo si están los 60 crudos. Con 4, la etapa queda `no_disponible`.
- H9 no se intenta. Queda `requiere_acceso`.

Tests. `test_enmienda_v3_4_es_anterior_a_sus_resultados` y `test_el_manifiesto_de_v3_3_sigue_igual` no existían y fallaban por ausencia; pasan al añadir el manifiesto y la guarda de cronología. Commit `b904b5e`.

## 2026-11-09 16:35 (UTC−5). v3.4, Fase 1

Comando: commit solo de `config/preregistro_v3.yaml` y `docs/enmiendas_v3.md`.

Resultado. Commit `a70d46b`. La fecha de la enmienda es 2026-11-09. No hay ficheros en `results/v3_4/`.

Decisión. Firebase solo con el checklist en verde. PDF por Overleaf. ESM-2 y TMbed sí. H9 no se intenta. La matriz de PXD064655 no se abre en este commit.

## 2026-11-09 16:45 (UTC−5). v3.4, H5

Comando: `pytest tests/test_h5_v34.py`.

Resultado antes del parche: `ImportError: cannot import name 'concordancias_mismo_universo'`. 1 error en la recolección.

Resultado después: 4 passed en 2,04 s. `test_el_bucle_de_v3_3_no_completa_ninguna_repeticion` documenta 0 de 200. `test_submuestreo_completa_200_y_n_sub_es_10` da `n_sub` 10 y 200 repeticiones. Commit `9c7725d`.

## 2026-11-09 16:55 (UTC−5). v3.4, controles del tejido

Comando: `pytest tests/test_controles_v34.py`.

Resultado antes del parche: `ModuleNotFoundError: No module named 'oncogo_crc.stages.v3_4_controles'`.

Resultado después: 5 passed. Un gen ausente es `no_detectado`. `n_obs` 4 es `no_evaluable_por_pares`. CA1 con −2,31 es `pasa`. El bloque con QC 16/16 pasa y con QC incompleto falla. Commit `a9e9c8b`.

## 2026-11-09 17:05 (UTC−5). v3.4, certificado

Comando: el certificado de v3.3 sobre una H7 con variables faltantes, y `pytest tests/test_certificado_v34.py`.

Resultado antes: `armar_certificado` de v3.3 devuelve estado `calculado` y `faltan` `[]`.

Resultado después: 1 passed. El certificado queda `incompleto` y nombra `h7_variables` y `h7_constructos`.

## 2026-11-09 21:58 (UTC−5). v3.4, publicación

Comando: `make all-v3-4`, luego `python scripts/v3/exportar_sitio_v3_4.py` y `python scripts/v3/exportar_tesis_v3_4.py`.

Resultado. La corrida firmada es el commit `3c3f666`, sello `2026-11-10T02:42:55Z`, árbol limpio. El certificado queda `publicable`. Las casillas verdes son H5, controles, H7, H6, multiverso (288, 72 y 2) y H9 (`exploratorio`). H8 queda `limitacion_declarada`: `report.pr_matrix.tsv` es una matriz de precursores y directLFQ no está instalado. H5 `sostiene`. H6 `no_sostiene`. H7 `refuta`, con el modelo principal completo. H10 `descriptiva`, 5 proteínas específicas. La preanalítica queda `no_disponible`. ESM-2 y TMbed entran solo en el modelo secundario: 668 genes con secuencia y 20 sin ella.

Decisión. La página y `results/v3_4/` publican esos JSON. El Word añade la sección generada desde los mismos ficheros. Ninguna cifra de esta entrada se escribió a mano.

## 2026-11-09 22:15 (UTC−5). v3.4, página de lectura

Comando: `python scripts/v3/exportar_sitio_v3_4.py results/v3_4` y `pytest tests/test_v3_4_pendientes.py`.

Resultado. 9 passed. La página deja de ser solo el checklist. Explica la pregunta, el orden de los ocho pasos, la razón de cada fuente, el resultado de H5, controles, H6, H7, la celda de suero, H8, H9 y H10, las catorce figuras de v3.3 y la conclusión. Las cifras salen de `results/v3_4/`. H7 sigue en refuta. H8 y la preanalítica siguen con el motivo de la corrida firmada `3c3f666`.

## 2026-11-10 06:41 (UTC−5). v3.5, Fase 0

Comando: `git checkout main`, `git pull --ff-only origin main`, `make regresion-v1`, comparación de manifiestos, `stat` de los seis crudos y lectura de `external/pxd055218_metadata/CN_02_DIA.json`.

Resultado. Cabeza `bff0c5e`, árbol limpio, disco montado. Regresión v1: 4 passed en 9,90 s (688, 140, 28, 104/119, 226, 307). Manifiestos v3.2 (112), v3.3 (189) y v3.4 (50) coinciden. SHA-256 del report DIA-NN coincide. Los seis crudos rechazados tienen el tamaño de PRIDE. El campo de fecha es `Content Creation Date`. Hay 60 AppleDouble. Etiquetas `v3.2-final`, `v3.3-final` y `v3.4-final` (esta sobre `bff0c5e`) subidas a origin.

Decisión. No se vuelven a descargar los seis crudos. La enmienda usa el campo de fecha leído en el disco.

## 2026-11-10 06:48 (UTC−5). v3.5, H8

Comando antes del parche: `pytest tests/test_v3_5_h8.py`.

Resultado antes: 1 error. `ImportError: cannot import name 'agregar_directlfq'`.

Comando después: el mismo pytest.

Resultado después: 4 passed. El mapa une 2 y excluye 1 estándar. directLFQ 0.3.3 escribe 2 proteínas y 4 columnas. Pandas entrega arrays de solo lectura; el cálculo copia el array antes de que directLFQ escriba encima. No se suman precursores.

## 2026-11-10 06:52 (UTC−5). v3.5, orden de adquisición

Comando antes: `pytest tests/test_v3_5_preanalitica.py`.

Resultado antes: 1 error. `ModuleNotFoundError: No module named 'oncogo_crc.stages.v3_5_preanalitica'`.

Comando después: el mismo pytest.

Resultado después: 3 passed. Los `._*` no cuentan. El campo `Content Creation Date` se lee como fecha. Con 54 de 60 el estado es `no_evaluable` y el AUC queda en `calculado_con_legibles`.

## 2026-11-10 07:05 (UTC−5). v3.5, parches

Certificado. Antes: `ModuleNotFoundError: oncogo_crc.stages.v3_5_report`. Después: 1 passed. H8 `no_evaluable` con solo la regla de H6 queda en rojo.

H7. `partir_recursos` no existía en el commit anterior. Después: 1 passed. Un recurso 1/125 no entra en `entrenables`.

H6. Después: 1 passed. La tabla `[[20, 0], [98, 0]]` marca `tabla_degenerada` y el estado sigue en `no_sostiene`.

Controles. Después: 1 passed. `n_fallos` 3 da `falla` y 0 da `pasa`.

Triangulación. Después: 1 passed. 14 genes y ninguna celda vacía: lo que falta dice `no_medida`.

El conjunto de estos tests, con la cronología, queda en 19 passed.

## 2026-11-10 07:06 (UTC−5). v3.5, corrida

Comando: `make v3-5`. Duración 135613 ms. Salida `0`.

Resultado. `results/v3_5/certificado_v3_5.json` queda `publicable` y `puede_publicar` es verdadero. H5 sigue en sostiene, copiada de v3.4. H6 queda `no_sostiene` con `tabla_degenerada` verdadero (`h6_origen_celular.json`): no es evidencia ni a favor ni en contra. H7 sigue en `refuta`. Papier queda `no_evaluable` (1 y 125) y no entra en el entrenamiento. H8 une 100 de 100, excluye 4 estándares y directLFQ 0.3.3 escribe `pr_matrix_SH.aq_reformat.tsv.protein_intensities.tsv`. El estado es `no_evaluable`: menos de 10 tumorales epiteliales en el universo (`h8_carga_tumoral.json`). No demuestra que el estadio IV no lleve proteína del tumor. La preanalítica cuenta 60 crudos y 54 con fecha; el estado es `no_evaluable` (`preanalitica_v3_5.json`). Los crudos legibles salen por bloques: CP, luego CN, luego HC. La celda de ML reproduce 0,88025 (`ml_celda_preregistrada.json`). El valor 0,888875 de v3.4 queda como variación de semilla. La triangulación tiene 14 filas y 0 celdas vacías. La queratina epitelial de la traza, `A5A6N0;P08729;Q9DCV7`, no está en la matriz de proteína, así que el eje no cambia el contraste (`multiverso_resumen_v3_5.json`).

Decisión. Se copian los ficheros al repositorio y se congela `results/v3_5_manifest.sha256`. Las limitaciones declaradas son H8, la preanalítica y el recurso de H7. Están en la lista de la enmienda.

## 2026-11-10 07:15 (UTC−5). v3.5, página

Comando: `python scripts/v3/exportar_sitio_v3_5.py`. Salida: las páginas de `sitio/`.

Prueba. Sin `sitio/index.html`, `pytest tests/test_sitio_v3_5.py` da 2 failed y 1 passed. Con la página regenerada, da 3 passed.

Lighthouse sobre `http://127.0.0.1:8765/index.html`: rendimiento 100, accesibilidad 100, buenas prácticas 96. El informe quedó en `/tmp/lh-v35.json`.

Decisión. La página conserva `v3_3.html`, `lectura-31-oct.html` y copia la página anterior como `v3_4.html`. Cada cifra lleva `data-fuente`.

## 2026-11-10 07:46 (UTC−5). v3.6, fase 0

Comandos: `git checkout main`, `git fetch origin --tags`, `git pull --ff-only origin main`, `git status`, `pytest tests/test_regresion_v1.py` con `ONCOGO_REGRESION_V1=1`, y la comprobación de los cuatro manifiestos.

Resultado. La cabeza es `ef7262f` y el árbol está limpio. El disco `/Volumes/JANCODATA/OncoGO-CRC` está montado. La regresión v1 da 4 passed en 9,69 s (`tests/test_regresion_v1.py`). Los manifiestos de `results/v3_2_manifest.sha256` (112), `results/v3_3_manifest.sha256` (189), `results/v3_4_manifest.sha256` (50) y `results/v3_5_manifest.sha256` (70) coinciden, con 0 hashes distintos. El SHA-256 de `data/raw/plasma/DIA-NN_report.tsv` en el disco es `a8a5aadd74391afef68cf4d6c3a7d3860023101e5c92590e2e87c152a39e59e8`. Las etiquetas locales y de GitHub son v3.2-final@c2be8e4, v3.3-final@4f85911, v3.4-final@bff0c5e y v3.5-final@ef7262f. ThermoRawFileParser en `external/bin/trfp/ThermoRawFileParser` informa `2.0.0.0`. No hay otra copia más reciente en el Mac.

Decisión. Se sigue con la enmienda v3.6. El inventario de lecturas de corridas anteriores queda en `docs/auditoria_v3.md`.

## 2026-11-10 08:20 (UTC−5). v3.6, código de la corrida limpia

Comandos: `pytest tests/test_v3_6_preanalitica.py tests/test_v3_5_preanalitica.py` y `pytest tests/test_v3_6_sin_mezcla.py tests/test_v3_6_preanalitica.py tests/test_formato_v3_6.py tests/test_v3_6_entradas.py`.

Resultado. La primera prueba da 5 passed. La segunda da 7 passed. El formato había fallado en la colección (salida 2) antes de existir el módulo y después dio 3 passed, anotado en la fase del módulo.

Decisión. Los commits `e33c1f7`, `76d14ef`, `1b5e0a7`, `80fb952` y `8eb0f4c` quedan en main antes de `make all-v3-6`. La corrida todavía no empieza en esta entrada. La hora de esta entrada dice 08:20; el commit de la bitácora fue a las 07:57 (UTC−5).

## 2026-11-10 08:46 (UTC−5). v3.6, tesis y página

Comandos: `python scripts/v3/exportar_cifras_tesis.py`, `python scripts/v3/figuras_tesis_v3_6.py`, `python tesis-ucaldas/build_word.py`, `python scripts/v3/exportar_sitio_v3_6.py` y `pytest tests/test_tesis_v3_6.py tests/test_sitio_v3_6.py`.

Resultado. El exportador escribe 49 macros. La prueba de la tesis da 3 passed y la del sitio da 2 passed, 5 passed en total. Chrome sin interfaz escribe `docs/capturas_v3_6/inicio-1440.png`, `inicio-390.png` y `resultados-1440.png`. El volcado del DOM de la portada sale con código 0. No se ejecutó Lighthouse: no está en las herramientas de esta corrida.

Decisión. La tesis, el Word, el zip de Overleaf y la página quedan en el repositorio. La discusión y la introducción siguen con `% REDACTAR-AUTOR`. Los objetivos llevan `% REVISAR-AUTOR`.

## 2026-11-10 09:20 (UTC−5). v3.7, fase 0

Comandos: `git checkout main`, `git fetch origin --tags`, `git pull --ff-only origin main`, `git status`, `git tag -l`, `git ls-remote --tags origin`, la comprobación de los cinco manifiestos, `ONCOGO_REGRESION_V1=1 pytest tests/test_regresion_v1.py -q` y el SHA-256 y el SHA-1 de `data/raw/plasma/DIA-NN_report.tsv`.

Resultado. La cabeza es `859b86b` y el árbol estaba limpio. El disco está montado. Las etiquetas locales y de GitHub son v3.2-final@c2be8e4, v3.3-final@4f85911, v3.4-final@bff0c5e, v3.5-final@ef7262f y v3.6-final@859b86b. Los manifiestos coinciden: v3.2 112, v3.3 189, v3.4 50, v3.5 70 y v3.6 210, con 0 hashes distintos. La regresión v1 da 4 passed en 11,54 s. El report pesa 213313241 bytes. Su SHA-256 es `a8a5aadd74391afef68cf4d6c3a7d3860023101e5c92590e2e87c152a39e59e8` y su SHA-1 es `ddf7f2a5d615a319577d0299754f186836881ad1`. Hay 60 crudos en `data/raw/plasma/pxd055218_raw/`. Python 3.14.3, numpy 2.4.4, pandas 3.0.2, scipy 1.17.1, scikit-learn 1.8.0, statsmodels 0.14.6, matplotlib 3.10.9, gseapy 1.2.1, directlfq 0.3.3, Pillow 12.2.0, Jinja2 3.1.6 y openpyxl 3.1.5. Chrome está instalado. No están Playwright, Lighthouse, latexmk, lualatex, conda ni pyteomics.

Decisión. Se sigue. El inventario de lecturas queda en `docs/auditoria_v3.md`. Playwright y Chromium se instalan cuando llegue la comprobación visual. No se instala MacTeX. pyteomics queda POR VERIFICAR hasta el control del PRTC.

## 2026-11-10 09:22 (UTC−5). v3.7, enmienda

Comando, antes de escribir la clave: `pytest tests/test_cronologia.py::test_enmienda_v3_7_es_anterior_a_sus_resultados`.

Resultado. 1 failed. `AssertionError: no aparece enmienda_v3_7 en el preregistro`. El texto de trabajo fechaba la revisión el 11-nov-2026. El commit es del 10-nov-2026 y así queda escrito en `config/preregistro_v3.yaml`. La diferencia está en `docs/auditoria_v3.md` no: queda en esta entrada.

Decisión. La enmienda se escribe antes de cualquier fichero de `results/v3_7/`. No cambia umbrales. Después del commit `2b3c93d`, la misma prueba da 1 passed.

## 2026-11-10 09:24 (UTC−5). v3.7, integridad de crudos

Comando, antes del módulo: `pytest tests/test_v3_7_crudos.py`.

Resultado. Salida 2 en la colección: `ModuleNotFoundError: No module named 'scripts.download.pxd055218_crudos'`. Con el módulo y la carga por ruta, 1 passed.

Decisión. El comprobador compara el SHA-1 y pide redescarga si no coincide. La pasada sobre los 60 crudos reales empieza después de este commit.

## 2026-11-10 09:25 (UTC−5). v3.7, el checksum trae la ruta delante

Comando: `python scripts/download/pxd055218_crudos.py --workspace /Volumes/JANCODATA/OncoGO-CRC`.

Resultado. Sale con código 1: `PRIDE no publica SHA-1 para CN_01_DIA.raw`. `checksum.txt` tiene 4838 bytes y cada línea es `D:\Akhila\PhD\DIA\PRIDE\<crudo>`, tabulador y el SHA-1. El lector tomaba el texto que va detrás del hash, y ahí no hay nombre. La prueba nueva falla con `KeyError: CN_01_DIA.raw` y, con el lector corregido, `pytest tests/test_v3_7_crudos.py` da 2 passed.

Decisión. Se relanza la pasada de los 60.

## 2026-11-10 09:28 (UTC−5). v3.7, el crudo nuevo no cruza de disco

Comando: la misma pasada. `CN_01_DIA.raw` se descargó entero (1001 MB) y su SHA-1 coincidió. Al moverlo del caché de APFS al volumen exFAT, `Path.replace` salió con `OSError: [Errno 18] Cross-device link`. El crudo viejo ya estaba en `cuarentena/` y el nuevo seguía en `~/Library/Caches/oncogo/`.

Prueba. Antes del arreglo, `test_colocar_el_crudo_nuevo_deja_el_viejo_en_cuarentena` da `AttributeError`. Después, `pytest tests/test_v3_7_crudos.py` da 3 passed.

Decisión. `shutil.move` hace la copia entre discos. Si el fichero del caché ya tiene el SHA-1 de PRIDE, no se vuelve a descargar.

## 2026-11-10 09:29 (UTC−5). v3.7, AppleDouble en la cuarentena

Comando: la pasada relanzada. Sale con `FileNotFoundError` al renombrar `cuarentena/._CN_01_DIA.raw.danado`. El crudo real ya se llama `CN_01_DIA.raw.danado-2026-11-10`. El fichero `._` es un AppleDouble y desaparece al renombrar el crudo.

Decisión. Esos nombres se saltan. La pasada se relanza.

## 2026-11-10 11:28 (UTC−5). v3.7, los 60 crudos coinciden con PRIDE

Comandos: `ascp` 4.4.3 de IBM Aspera Connect, usuario `fasp-public`, puerto 33001, y `python scripts/download/pxd055218_crudos.py --workspace /Volumes/JANCODATA/OncoGO-CRC`.

Resultado. HTTPS no entrega `CN_16_DIA.raw`: el objeto publicado pesa 4096 bytes. Aspera sí entrega el crudo de 1205043346 bytes y el SHA-1 `5ee92fddbb3603f60c4e07819031d06566c560d4`. `CP_17_DIA.raw` se cortó por HTTPS (curl, código 18) y se completó por Aspera. El primer envío de `HC_17_DIA.raw` pesaba 1065175806 bytes, el tamaño de PRIDE, pero el SHA-1 era `4c27517c8675f035a480d52bb42a0a905a4af493`. El segundo y el tercero dieron `0f40e8209011141ef37d6033f681b8a5620c4c6c`, el de `checksum.txt`. Se instaló el segundo. La pasada final imprime `60 crudos; distintos 0`. `results/v3_7/crudos_integridad.tsv` tiene 60 filas con `igual` verdadero. Doce quedaron `redescargado` verdadero porque el crudo viejo está en `cuarentena/`: CN_01, CN_08, CN_13, CN_15, CN_16, CP_10, CP_17, CP_18, HC_10, HC_11, HC_17 y HC_20.

Prueba. Antes de `anotar_redescarga`, `pytest tests/test_v3_7_crudos.py::test_la_cuarentena_marca_el_crudo_como_redescargado` da 1 failed (`AttributeError`). Después, `pytest tests/test_v3_7_crudos.py` da 4 passed.

Decisión. No es un fichero dañado en origen: el SHA-1 del que quedó en el volumen coincide con PRIDE. No se redacta el correo. `make all-v3-7` puede seguir. La clave de Aspera no entra al repositorio.

## 2026-11-10 11:50 (UTC−5). v3.7, etapa que recalcula en su carpeta

Comando, antes del parche: `pytest tests/test_v3_7_sin_mezcla.py`.

Resultado. 2 failed. `clasificar()` no acepta `version` y `h8.py` todavía asigna el gen de Tanaka con `uni.get(acceso_base`. Después del parche, `pytest tests/test_v3_7_sin_mezcla.py tests/test_v3_6_sin_mezcla.py tests/test_v3_7_crudos.py tests/test_v3_5_triangulacion.py` da 10 passed.

Decisión. La etapa `v3_7` escribe en `results/v3_7/`. H8 y la triangulación usan `genes_tanaka.mapa_genes`. El report trae `Run` igual a `CP_` y `CN_`, y a los sanos los llama `Normal_` en un `.mzML`, no `HC_` en un `.raw`. Esa unión no se inventa: queda `verificada` falso. ThermoRawFileParser 2.0.0 ya escribió la fecha de los seis crudos que no tenían JSON.

## 2026-11-10 14:20 (UTC−5). v3.7, certificado de la recálculo

Comandos: `make v3-7` y, al cortarse, el tramo desde H9 con el sello ya escrito.

Resultado. `make v3-7` salió con código 2 al imprimir `[v3.7] h9`, sin traza de Python. El mismo corte, también con código 2 y sin traza, había ocurrido antes al imprimir `[v3.7] h10`. H9, H10, la unión, la preanalítica, la triangulación, la repetición de M2 y el certificado se calcularon después con el sello `2026-11-10T18:31:55Z` del commit `c1c2289`. `results/v3_7/certificado.json` queda `publicable` y `rojas` está vacío. La limitación aceptada es H8. La celda sigue en 284, 0, 166 y 118. M2 sigue en 0,88025. H8 sigue `no_evaluable` y el universo pasa de 33 a 53. La preanalítica pasa de `no_evaluable` con 54 fechas a `calculado` con 60. `union_crudos_calidad.tsv` deja verificadas las 40 corridas CP y CN, y falsas las 20 `Normal_`, con la evidencia «el report no nombra ese crudo». `lecturas_v3_7.json` solo cubre el segundo proceso. `sha256_manifest_entradas` queda nulo. `determinismo.json` solo repite M2.

Decisión. Se publican los ficheros de esta recálculo y la documentación dice qué cambió y qué no se calculó. No se etiqueta `v3.7-final` y no se despliega Firebase: faltan `vias_suero.json`, `h1_vias.json`, `go_ora_tejido.tsv`, `catalogo_auditoria.tsv`, `neoantigenos.tsv`, `prtc_por_corrida.tsv` y `manifest_entradas.sha256`. El detalle está en `docs/estado_v3_7.md` y `docs/changelog_v3_6_frente_a_v3_7.md`.
