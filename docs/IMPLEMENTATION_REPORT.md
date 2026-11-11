# Informe de la implementación

Fecha: 2026-11-11. No se publicó una versión nueva de Hosting. La página en vivo sigue siendo `0fe6e84d75d1c6df`.

## Qué se modificó

La portada (`sitio_fuente/contenido/inicio.md`) cuenta la pregunta de la tesis, la celda, las dos matrices, el diagrama D1, el embudo de priorización y el alcance. El menú principal pasa a cinco entradas. El resto de secciones queda en el pie. Los tokens de color claro cambian a texto `#142d4e` sobre fondo `#f5f7fa`, comprobados a 4,5:1. El título `h1` usa `clamp(1.8rem, 4vw, 2.7rem)`.

`scripts/sitio/construir.py` añade cuatro piezas, todas leídas de ficheros o de una negativa explícita:

- `{{celda}}`: los cuatro enteros de `suero_embudo.json#celda_preregistrada`.
- `{{embudo:priorizacion}}`: nombre y n de cada paso de `priorizacion_embudo.json`. La anchura es `100 * n / n_mayor`. El paso 6 y el 7 miden 0.
- `{{anatomia:colon}}`: esquema con seis nombres. Cada botón dice que la tesis no asigna muestras a un segmento.
- `{{llamadas}}` y `{{explorar}}`: enlaces a páginas que ya existían.

`/ayuda` añade tres comprobaciones: CN no se describe como grupo sin enfermedad, el paso 6 y el estado de H7.

## Qué se reutilizó

El generador Jinja, las fichas, el comparador, las candidatas, Mol*, los diagramas D1–D11, los gráficos, el glosario, la bibliografía, el archivo v3.7, las redirecciones y el aviso en dos sitios. No hay `package.json`. No se añadieron fuentes descargadas: siguen Source Sans 3 y Source Serif 4.

## Fuentes

`results/v3_7/suero_embudo.json`, `priorizacion_embudo.json`, `k1_etiquetas.json`, `h7_visibilidad_en_sangre.json` y `config/preregistro_v3.yaml`, más las citas ya comprobadas de Rai 2025, Vasaikar 2019, Tanaka 2024, Zhang 2014 y García-Foncillas 2026.

## Qué no estaba disponible

Instrumento, tamaño y grupos de varios depósitos, listados en `docs/sitio_v4_omitido.md`. Segmento del colon de las muestras. Versión de STRING. Modelos de P04114, P22352, P49908 y Q9Y6R7. Por eso no hay una red nueva, ni un quiz clínico, ni una ilustración de proteínas pasando del tumor a la sangre.

## Pruebas

La suite completa, con la portada ya generada y el lector de `test_el_estado_v3_7_coincide_con_sus_ficheros` todavía sin índices de lista, dio 1 failed, 598 passed, 1 skipped, en 83,09 s. El fallo era `nadie_pasa_el_paso_6`: el HTML guarda `true` y `str(True)` es `True`. El lector tampoco bajaba a `pasos.0.n`.

El lector quedó alineado con `caminar`: si el nodo es una lista, el tramo es un entero; si es booleano, se compara en minúsculas. La suite completa, ya con ese lector y con la hoja `v4.d6da1f8c47.css`, dio 599 passed, 1 skipped, 88 warnings, en 79,18 s. El salto sigue siendo `test_reejecutar_v1_reproduce_los_numeros`. Son 39 páginas de markdown.

## Limitación

Este cambio está en el árbol de trabajo. No hay commit ni `firebase deploy` de esta pasada.
