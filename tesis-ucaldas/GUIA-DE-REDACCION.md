# Guía de redacción

Este archivo dice qué tienes que escribir tú. El PDF `tesis-ucaldas/latex/main.pdf` y el Word `tesis-ucaldas/word/tesis-maestria-ucaldas.docx` son el andamiaje, no la tesis. Cada `% TODO Juan Diego:` del LaTeX está repetido en el Word como un párrafo en gris y cursiva.

La rama es `tesis-andamiaje-v4`. No se modificaron `results/` ni el código que calcula la ciencia.

## Por qué la estructura no es la del encargo, palabra por palabra

El encargo pedía seis capítulos y meter planteamiento, justificación y objetivos en la introducción. La guía de 2022 y las tesis en español del programa hacen lo mismo. Eso se respetó.

El repo tenía un capítulo aparte, «Desviaciones del preregistro y trabajo futuro». Esas desviaciones siguen escritas en `latex/chapters/11-desviaciones.tex`, pero **ya no entran en `main.tex`**. En el andamiaje son una sección de la Discusión, para no pasar de seis capítulos. No borres ese archivo: es la prosa anterior y no se debe pegar sin reescribirla.

No hay codirector. Esa línea salió de la portada.

## Orden sugerido para escribir

1. Introducción: problema, justificación, pregunta y objetivos. De 2 a 4 páginas.
2. Resumen, máximo 250 palabras, y después el Abstract como traducción.
3. Marco teórico y estado del arte.
4. Métodos: reescribe con tus palabras los parámetros que ya están. No cambies un umbral.
5. Resultados: la interpretación de cada figura. El caption técnico se queda.
6. Discusión, con las limitaciones y las desviaciones.
7. Conclusiones, una por objetivo, y las recomendaciones.
8. Anexos: elige qué tablas suplementarias pegas.
9. En Word: clic derecho en la tabla de contenido, la lista de tablas y la lista de figuras, y actualiza toda la tabla.
10. Exporta el PDF/A como `JuanDiego_PolancoFranco_2026.pdf` cuando ya no quede un pendiente que tú debas resolver.

## Qué hay en cada sección

### Portada

Director: Doctor en Ingeniería, Jorge Alberto Jaramillo Garzón. Grupo: GITIR. Línea: Modelos Computacionales y Bioinformática. No hay codirector. No hay dedicatoria ni agradecimientos.

### Declaración de IA

Cursor está en la versión 3.28.7. SciSpace figura como consultado entre 2026-03 y 2026-11. El párrafo borrador sigue en Métodos y en esa declaración para que lo ajustes. Cursor y SciSpace no se citan como fuente de un hecho.

### Resumen y Abstract

Objetivo: 200 a 250 palabras, con objetivo, diseño, métodos, resultados y conclusión. Usa la celda 0 / 166 / 118 de 284, la concordancia 147 de 147, los AUC 0,880, 0,915 y 0,647 (p 0,12), y los estados H5 a H10. No conviertas el cero en «no hay diferencia biológica». Palabras clave: deja entre 3 y 7.

### 1. Introducción

Objetivo: plantear el problema sin prometer un examen. Preguntas: qué se conoce en tejido, qué parte aborda la tesis, cuál es la pregunta y cuáles son los objetivos en infinitivo. Extensión: 2 a 4 páginas. Citas para tener a mano: Vasaikar 2019, Zhang 2014, Tanaka 2024, Rai 2025, García-Foncillas 2026, Pepe 2001.

### 2. Marco teórico y estado del arte

Objetivo: definir solo lo que luego lee un resultado. Preguntas: por qué las matrices no se fusionan; qué significan equivalente, diferente e indeterminada; qué trabajo sí cambia la lectura. Extensión: 8 a 15 páginas.

Estas claves están en `references.bib` y **todavía no están citadas en una sección**. Decides si entran aquí o se quedan solo en la bibliografía, que hoy las incluye todas con `\nocite{*}`: `seo2026`, `dhami2026`, `malmstrom2025`, `lesur2023`, `li2024`, `jumper2021`, `uhlen2015`, `verathamjamras2023`, `yu2026`, `ogren2026`, `robles2024`, `hulstaert2020`, `hermoso2026`, `locker2006`, `ley1384`.

### 3. Materiales y métodos

Objetivo: que otra persona repita la corrida. Lo técnico ya está, marcado para que lo reescribas: alfa 0,05, Benjamini-Hochberg, margen 0,585, sin imputar, sin fusionar, celda preregistrada, modelos M1 a M6, embudo 688, 135, 135, 118, 20, 0 y 0. Entorno: Python 3.14.3, numpy 2.4.4, pandas 3.0.2, scipy 1.17.1, scikit-learn 1.8.0, semilla registrada 20260914. Extensión: 8 a 12 páginas.

### 4. Resultados

Objetivo: cifra, figura y qué no demuestra, en la misma página. No redactes el caption otra vez. Escribe la interpretación debajo de cada `% TODO`. Extensión: la que pidan las nueve figuras, sin repetir Métodos.

### 5. Discusión

Objetivo: interpretar y limitar. Incluye que es trabajo *in silico* y que no propone un examen. La sección de desviaciones reemplaza al antiguo capítulo 11. Extensión: 6 a 10 páginas. No afirmes que no hay población latinoamericana en los depósitos: eso quedó `[PENDIENTE]` porque el repo no lo dice en los ficheros leídos para este andamiaje.

### 6. Conclusiones y recomendaciones

Una frase por objetivo, con la cifra y con lo que no demuestra. En el trabajo futuro no anticipes un AUC ni un tamaño de cohorte.

### Anexos

La tabla de correspondencia sale de `config/correspondencia_propuesta.tsv`. La palabra «cumplido» es la etiqueta de ese fichero. No la conviertas en un biomarcador.

## Bibliografía

Estilo APA 7 (`biblatex`, `style=apa`). El 11 de noviembre de 2026 se consultó `api.crossref.org` para cada DOI. Los 33 DOI devolvieron un título que coincide con la entrada. No hubo fallos.

Sin DOI, a propósito: `ley1384` (la ley no tiene DOI), `cursor2026` (versión 3.28.7) y `scispace2026` (consultado entre 2026-03 y 2026-11).

`\nocite{*}` mete en el PDF toda la biblioteca verificada. Cuando cites en la prosa cada trabajo que sí uses, puedes quitar el `\nocite{*}` y volver a compilar. Hasta entonces no hay citas huérfanas en el PDF: las que no están en una sección están en la bibliografía por esa orden.
