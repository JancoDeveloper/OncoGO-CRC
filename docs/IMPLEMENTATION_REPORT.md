# Informe de la implementación

Fecha: 2026-11-11. Publicado. Hosting `161f72104d6d9530`, liberada el 2026-11-11T16:19:47.041Z, sitio `e77d0e1`. La comprobación en vivo está en `docs/despliegue_sitio_v4.md`.

## Qué se modificó

La portada (`sitio_fuente/contenido/inicio.md`) cuenta la pregunta de la tesis, la celda, las dos matrices, el diagrama D1, el embudo de priorización y el alcance. El menú principal pasa a cinco entradas. El resto de secciones queda en el pie. Los tokens de color claro cambian a texto `#142d4e` sobre fondo `#f5f7fa`, comprobados a 4,5:1. El título `h1` usa `clamp(1.8rem, 4vw, 2.7rem)`.

`scripts/sitio/construir.py` añade estas piezas, todas leídas de ficheros o de una negativa explícita:

- `{{celda}}`: los cuatro enteros de `suero_embudo.json#celda_preregistrada`.
- `{{embudo:priorizacion}}`: nombre y n de cada paso de `priorizacion_embudo.json`. La anchura es `n / n_mayor`. El paso 6 y el 7 miden 0. Cada paso se abre y repite el nombre guardado en el JSON.
- `{{matrices}}`: la concordancia de tejido, 147 de 147, y los diferentes del suero, 0 de 284. Las dos pistas no se unen.
- `{{anatomia:colon}}`: esquema con seis nombres. Cada botón dice que la tesis no asigna muestras a un segmento y repite las cifras del estudio, no del segmento.
- `{{llamadas}}` y `{{explorar}}`: enlaces a páginas que ya existían.

En cada ficha, el suero y el tejido quedan marcados como resultado calculado. La función es el texto de UniProt. La etiqueta de candidata, cuando sale, sigue siendo la hipótesis del paso 5.

`/ayuda` añade tres comprobaciones: CN no se describe como grupo sin enfermedad, el paso 6 y el estado de H7.

La hoja usa `font-display: optional` y la plantilla precarga las tres fuentes ya alojadas. La barra del embudo y la de las matrices cambian de escala sin mover el texto.

## Qué se reutilizó

El generador Jinja, las fichas, el comparador, las candidatas, Mol*, los diagramas D1–D11, los gráficos, el glosario, la bibliografía, el archivo v3.7, las redirecciones y el aviso en dos sitios. No hay `package.json`. No se añadieron fuentes descargadas: siguen Source Sans 3 y Source Serif 4. Firebase Hosting sigue sirviendo `sitio/`.

## Fuentes

`results/v3_7/suero_embudo.json`, `priorizacion_embudo.json`, `k1_etiquetas.json`, `h7_visibilidad_en_sangre.json` y `config/preregistro_v3.yaml`, más las citas ya comprobadas de Rai 2025, Vasaikar 2019, Tanaka 2024, Zhang 2014 y García-Foncillas 2026.

## Qué no estaba disponible

Instrumento, tamaño y grupos de varios depósitos, listados en `docs/sitio_v4_omitido.md`. Segmento del colon de las muestras. Versión de STRING. Modelos de P04114, P22352, P49908 y Q9Y6R7. Por eso no hay una red nueva, ni un quiz clínico, ni una ilustración de proteínas pasando del tumor a la sangre.

## Pruebas

La suite completa, con la portada ya generada y el lector de `test_el_estado_v3_7_coincide_con_sus_ficheros` todavía sin índices de lista, dio 1 failed, 598 passed, 1 skipped, en 83,09 s. El fallo era `nadie_pasa_el_paso_6`: el HTML guarda `true` y `str(True)` es `True`. El lector tampoco bajaba a `pasos.0.n`.

El lector quedó alineado con `caminar`: si el nodo es una lista, el tramo es un entero; si es booleano, se compara en minúsculas. La suite completa, ya con ese lector y con la hoja `v4.d6da1f8c47.css`, dio 599 passed, 1 skipped, 88 warnings, en 79,18 s. El salto sigue siendo `test_reejecutar_v1_reproduce_los_numeros`. Son 39 páginas de markdown.

Después de añadir las matrices, el embudo en `details` y la hoja `v4.a1fc0bd2e5.css`, `tests/test_sitio_v4.py` y `tests/test_sitio_v3_6.py` dieron 18 passed. La suite completa no se volvió a correr sobre esa hoja.

## Publicación

Hosting `161f72104d6d9530`. La comprobación con `?v=4` pasó: las páginas nuevas en 200 con el aviso dos veces, las 16 rutas viejas en 301, el archivo v3.7 y v3.6 en 200, y las cinco fichas con su `log2fc`. El detalle está en `docs/despliegue_sitio_v4.md`.
