# Auditoría UX/UI

Fecha: 2026-11-11. Objeto: el sitio publicado en https://oncogo-crc-tesis.web.app antes de este rediseño (Hosting `0fe6e84d75d1c6df`, commit `e4b8c81`).

## Qué es el proyecto

No hay `package.json` ni aplicación React. La web se genera con `scripts/sitio/construir.py` (Jinja2) desde `sitio_fuente/` y `results/v3_7/`. Firebase Hosting sirve la carpeta `sitio/`. Los gráficos son SVG propios. El visor de estructura es Mol* 4.18, local, y solo se pide al pulsar el botón.

No se reemplaza ese generador: ya cita la corrida, archiva la v3.7 y pasa las pruebas de contraste, aviso y citas.

## Problemas de UX

El menú listaba una entrada por sección: Investigación, Datos, Métodos, Resultados, Proteínas, Discusión, Conclusiones, Entregables, Archivo, Ayuda, Glosario y Bibliografía. En un móvil esa barra se parte en muchas filas antes de llegar al contenido.

La portada abría con la celda (0 diferentes, 166 equivalentes, 118 indeterminadas, de 284) y enseguida con tarjetas de sección. No había un recorrido que separara la pregunta, las dos matrices, el embudo y el alcance.

La búsqueda existe en `/buscar` y lista páginas, fichas y términos. No explica, en la portada, que el índice no inventa proteínas.

## Problemas de UI

La página se leía como un documento continuo: una columna, títulos en serif y cifras en el mismo cuerpo que el párrafo. El cero de la celda no tenía una jerarquía distinta de una frase secundaria.

Las tarjetas de la portada repetían el nombre de la sección y no el resumen de la página de destino.

## Navegación

Cada sección seguía en su ruta. Las rutas viejas de la v3.7 redirigen con 301. El problema no era la falta de páginas, sino que el menú principal y el secundario eran el mismo listado.

## Legibilidad

Los tokens de entonces pasaban WCAG AA (4,5:1) en claro y en oscuro. El cuerpo usa Source Sans 3 y los títulos Source Serif 4, alojados en el sitio. No hay fuentes de un CDN.

El aviso de que ningún resultado sirve para decidir sobre un paciente está dos veces en cada página. Esa repetición es una regla del sitio, no un defecto a quitar.

## Responsive

Por debajo de 800 px la disposición pasa a una columna y el índice deja de ser fijo. Los gráficos SVG usan `width: 100%`. El visor de estructura mide 24 rem de alto.

## Visualización de datos

Los gráficos r01–r08 y las figuras de la corrida salen de `results/v3_7/`. Cada diagrama D1–D11 tiene tabla y sello. El embudo de priorización estaba en `/resultados/priorizacion` como figura, no como una barra cuya anchura se lee del mismo JSON.

No hay un mapa del colon en las fuentes. Dibujar lesiones o asignar muestras a un segmento sería inventar.

## Accesibilidad

axe-core, el 2026-11-11, sobre 328 HTML fuera de `archivo`: 0 hallazgos serious o critical (`docs/sitio_v4_axe.json`). Los nodos de los diagramas no son enlaces dentro de un SVG con `role="img"`: el enlace está en la tabla.

Lighthouse 13.6.0, en `/`, `/metodos/canalizacion` y `/proteinas/P02749`, a 1440 y a 390: rendimiento 96–100, accesibilidad 100, buenas prácticas 96.

## Inventario de funcionalidades

Portada, investigación, datos, métodos, resultados, 284 fichas, comparador, candidatas, discusión, conclusiones, entregables, archivo de la v3.7, ayuda, glosario (96 términos), bibliografía, búsqueda, tema claro/oscuro, visor Mol* bajo demanda.

## Inventario de datos

El detalle está en `docs/SCIENTIFIC_DATA_INVENTORY.md`. La corrida sellada es `results/v3_7/`, sello `71b310a3fcae`, fecha 2026-11-11T00:48:45Z.

## Riesgos de integridad

Llamar a una proteína un hallazgo de sangre cuando el paso 6 tiene n 0. Dibujar el paso del tumor a la sangre. Completar instrumento o n de un depósito que no tiene respuesta de API guardada. Usar STRING como resultado de la tesis. Poner «sanos» junto a CN. Añadir el aviso una tercera vez con la misma frase, o quitar una de las dos.

## Arquitectura visual propuesta

Cinco entradas en el menú principal y el resto en el pie. Portada en este orden: pregunta, celda, problema, hipótesis de visibilidad, recorrido D1, embudo de priorización, tejido frente a suero, esquema del colon que declara que no hay localización, accesos a páginas ya existentes, alcance. Paleta clara con texto `#142d4e` sobre `#f5f7fa`, ambos por encima de 4,5:1. Sin CDN y sin cifras nuevas.
