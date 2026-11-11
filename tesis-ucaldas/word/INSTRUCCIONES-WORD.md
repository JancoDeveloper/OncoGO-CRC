# Word

El archivo `tesis-maestria-ucaldas.docx` es el andamiaje del 11 de noviembre de 2026, no la prosa final. Cada párrafo gris que empieza por «TODO Juan Diego» es lo que hay que completar. El orden y las cifras están en `../GUIA-DE-REDACCION.md`.

El archivo es `tesis-maestria-ucaldas.docx`. Se partió de la guía oficial:

https://biblio.ucaldas.edu.co/wp-content/uploads/2022/07/guia-presentacion-tesis-trabajos-grado.docx

Se conservó el papel carta y se aplicaron los estilos que la propia guía escribe y que sus ejemplos muestran en Arial. Los estilos de fábrica del archivo venían en Cambria 13; los párrafos de ejemplo ya estaban en Arial. El documento de trabajo usa Arial.

## Al abrirlo

1. Si Word pregunta si actualiza los campos, acepte.
2. Lea la hoja «Qué completar antes de exportar el PDF», justo después de la portada.
3. Sustituya los cuatro corchetes grises: director, codirector, grupo y línea. El repositorio no tiene esos nombres.
4. Clic derecho en la tabla de contenido, en la lista de tablas y en la lista de figuras, y elija **Actualizar campos → Actualizar toda la tabla**.
5. Los textos en gris y cursiva que queden son huecos. La dedicatoria y los agradecimientos no están: la guía de 2022 los deja como opcionales. Si se escriben, usan Título 1 y empiezan en página impar.

No exporte el PDF mientras un corchete de la portada siga visible. Las cifras del apartado 6.7 son las de `results/v3_7/`, sello `71b310a3fcaebe9267fac971156df3aab491e009`. No las cambie por las de `results/v3_3/` ni `results/v3_4/`.

## Estilos

| Estilo | Uso | Formato aplicado |
|---|---|---|
| Normal | Cuerpo | Arial 11, interlineado 1,5, espacio posterior de 16 pt (equivale al Enter adicional de la guía, sin párrafos vacíos) |
| Título 1 | Capítulos y preliminares | Arial 20, negrita, 100 pt antes, 24 pt después, interlineado sencillo |
| Título 2 | Subtítulos | Arial 16, 20 pt antes, 24 pt después |
| Título 3 | Subtítulos de segundo nivel | Arial 14, negrita |

Los capítulos nuevos deben usar **Título 1** y empezar en página impar: Diseño de página → Saltos → Página impar. No numere 13A ni «bis».

## Página

Carta, 21,59 cm × 27,94 cm. Márgenes simétricos: interior 3,60 cm; superior, inferior y exterior 2,54 cm. El número va en el encabezado, Arial 11: a la derecha en impares y a la izquierda en pares. Cubierta y portada no muestran número y sí se cuentan. Los preliminares van en romanos mayúsculos. La Introducción empieza en 1.

La guía lo formula como sugerencia. El archivo de ejemplo no traía activados los márgenes simétricos ni los encabezados pares e impares; este documento sí, para cumplir esa prosa.

## Figuras, tablas y ecuaciones

La figura lleva la leyenda debajo, con el formato `Figura 1: …`, y la fuente en la línea siguiente, no como nota al pie. La tabla lleva el título arriba: `Tabla 1: …`. La primera fila está marcada para repetirse si la tabla pasa de página. La ecuación (1) está centrada y el número queda a la derecha.

## Bibliografía

La lista está en APA 7, orden alfabético, sangría francesa. La guía dice que puede usarse cualquier norma, aplicada con rigor, y recomienda un gestor. Para rehacerla:

1. Instale Zotero o Mendeley y el complemento de Word.
2. Importe `latex/references.bib`.
3. Elija el estilo APA 7.ª edición.
4. Sustituya la lista pegada por las citas del gestor.

Si el programa pide Vancouver o NLM, cambie el estilo del gestor. No hay un estilo obligatorio publicado.

## Exportar PDF/A

Siga el instructivo de la Biblioteca:

https://biblio.ucaldas.edu.co/wp-content/uploads/2022/07/guia-creacion-archivos-pdf.pdf

1. Archivo → Guardar como → PDF.
2. Nombre: `JuanDiego_PolancoFranco_2026.pdf` (sin tildes).
3. En Opciones, active el formato que el instructivo marca para PDF/A y deje libres copiar, imprimir y seleccionar.
4. Metadatos, en el cuadro de Word antes de guardar: Autor `Polanco Franco, Juan Diego`; palabras clave separadas por punto y coma; Asunto `Tesis para optar al título de Magíster en Bioinformática y Biología Computacional`. Esos tres campos ya están en las propiedades del archivo.

El manual de autoarchivo del 12 de febrero de 2026 muestra otro nombre (`NombresApellidos_Programa_Año`). Confirme cuál pide el formulario el día del depósito.

## Opción de capítulo-artículo

La guía admite que un capítulo sea un artículo. Si se usa, hay que mirar la carta de autorización del repositorio: un artículo ya enviado a revista puede requerir embargo. No es un requisito de este archivo.
