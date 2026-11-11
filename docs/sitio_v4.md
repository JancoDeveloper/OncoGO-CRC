# Web documental

La página pública sale de `sitio_fuente/` con `scripts/sitio/construir.py`. El commit desplegado es `4b525f4`. La corrida que cita es la v3.7-final, etiqueta en `595118ae09ea6006c744c86b28cb0a10d42e0442`. No se volvió a correr la ciencia ni se recompiló la tesis.

El generador escribió 39 páginas de markdown. Hay 284 fichas, una por `protein_id` de `results/v3_7/suero_k4_crc_cn.tsv`. Los 65 grupos que llevan `;` enlazan esa carpeta sin codificar el punto y coma: el emulador devolvía 404 o 500 si la ruta iba con `%3B`.

La priorización deja 20 genes en el paso 5 y ninguno en el paso 6. La etiqueta de candidata sale en 21 fichas porque KNG1 está en dos grupos. El panorama de `/proteinas/candidatas` entra cada gen una vez.

## Anotación guardada

UniProt `2026_04`. Gene Ontology `releases/2026-07-26`.

AlphaFold: 352 respuestas guardadas. 348 traen `latestVersion` 6. Cuatro respuestas están vacías: `P04114`, `P22352`, `P49908` y `Q9Y6R7`. Los BinaryCIF alojados son 348 ficheros, 62080630 bytes, 59.205 MB, por debajo del corte de 83886080. La ficha de P02749 usa el modelo `AF-P02749-F1-model_v6` (`toolUsed`: AlphaFold Monomer v2.0 pipeline, `modelCreatedDate` 2025-08-01T00:00:00Z). Mol* 4.18 abre ese fichero con `loadStructureFromUrl(url, "mmcif", true)`: en este paquete el nombre registrado es `mmcif` y el tercer argumento marca el binario.

STRING quedó en estado `ok`, 59 aristas. La cabecera de versión está vacía.

El glosario tiene 96 términos.

## Revisión, antes de publicar

`pytest`: 598 passed, 1 skipped, 88 warnings, 96.02 s. El salto es `tests/test_regresion_v1.py::test_reejecutar_v1_reproduce_los_numeros`, que solo corre con `ONCOGO_REGRESION_V1=1`.

axe-core, en el emulador `http://127.0.0.1:5002`, sobre 328 HTML fuera de `archivo`: 0 hallazgos serious o critical. El informe es `docs/sitio_v4_axe.json`.

Lighthouse 13.6.0, las mismas tres páginas, rendimiento, accesibilidad y buenas prácticas:

| Página | Escritorio 1440 | Móvil 390 |
|---|---|---|
| `/` | 97, 100, 96 | 98, 100, 96 |
| `/metodos/canalizacion` | 99, 100, 96 | 100, 100, 96 |
| `/proteinas/P02749` | 96, 100, 96 | 100, 100, 96 |

Los JSON están en `docs/sitio_v4_lighthouse/`. Las capturas de 15 vistas, a 1440 y a 390, en claro y en oscuro, están en `docs/sitio_v4_capturas/`. Las cuatro de referencia de UniProt no se tocaron.

No se añadieron citas nuevas en este cierre, así que no se volvió a ejecutar `scripts/v3/verificar_enlaces.py`. La pasada que terminó en cero es la del commit `1aa32de`.

La portada de la tesis no se llenó: Director, Codirector, Grupo, Línea y Leyenda siguen vacíos. Lo demás que esta web no afirma está en `docs/sitio_v4_omitido.md`.
