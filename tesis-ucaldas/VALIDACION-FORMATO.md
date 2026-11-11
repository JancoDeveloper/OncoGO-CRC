# Validación de formato

## Chequeo del andamiaje, 11 de noviembre de 2026

La guía oficial sigue siendo `tesis-ucaldas/assets/guia-presentacion-tesis.docx`, descargada de la Biblioteca. La tabla de abajo, del 6 de noviembre, no se borró. Este bloque dice qué se volvió a comprobar hoy sobre el andamiaje.

| Requisito | Estado en este andamiaje |
|---|---|
| Papel carta 21,59 × 27,94 cm | LaTeX: `letterpaper` y geometry. Word: las secciones del .docx miden 21,59 × 27,94 cm |
| Márgenes 2,54 cm y 3,60 cm interior | LaTeX: `inner=3.60cm`, `outer=2.54cm`, `top` y `bottom` 2,54 cm. Word: 3,60 cm a un lado y 2,54 cm al otro |
| Arial 11, interlineado 1,5 | Word: estilo Normal del archivo que salió de la guía. LaTeX: TeX Gyre Heros 11 pt, porque Arial no está en TeX Live. Desviación ya anotada |
| Título de capítulo en 20 pt, página impar | LaTeX: `titlesec` 20 pt y `openright`. Word: Título 1 y salto `ODD_PAGE` antes de cada capítulo numerado |
| Preliminares en romanos; arábigos desde la Introducción | LaTeX: romanos desde `\AtBeginDocument` y `\mainmatter` justo antes de la Introducción |
| Cubierta y portada sin número, contadas | LaTeX: `\thispagestyle{empty}` |
| Resumen ≤ 250 palabras y 3 a 7 palabras clave | No hay prosa de resumen. El hueco está marcado para que lo escriba el autor |
| Seis capítulos como máximo | Introducción, Marco, Métodos, Resultados, Discusión, Conclusiones. Las desviaciones pasaron a una sección de la Discusión |
| APA 7 | `biblatex` `style=apa`. Los 33 DOI se consultaron en Crossref el 11 de noviembre de 2026 y el título coincidió |
| PDF/A `JuanDiego_PolancoFranco_2026.pdf` | No se exportó. El PDF compilado es `latex/main.pdf`, 53 páginas, LuaLaTeX y Biber, sin citas indefinidas. No es el PDF/A de entrega |
| Codirector | Eliminado, porque no hay |
| Director, grupo y línea | Director: Jorge Alberto Jaramillo Garzón. Grupo: GITIR. Línea: Modelos Computacionales y Bioinformática. Sin nota de confirmación |
| Dedicatoria y agradecimientos | Eliminados del `main.tex` y del Word |
| PNG a 300 dpi de f01–f09 | No hay cairosvg ni Inkscape. LaTeX incluye los PDF vectoriales. Los SVG quedan en `tesis-ucaldas/figures/`. El Word no incrusta esas figuras |

El detalle del 6 de noviembre sigue debajo.

Se leyeron, el 6 de noviembre de 2026, los archivos guardados en `tesis-ucaldas/assets/`. Si la prosa de la guía y el estilo de fábrica del mismo .docx no coincidían, se siguió lo que el ejemplo muestra en pantalla y lo que la guía escribe como sugerencia. La guía dice «se sugiere»: no es un reglamento rígido.

| Requisito | Fuente | Word | LaTeX | Notas |
|---|---|---|---|---|
| Papel carta 21,59 × 27,94 cm | Guía 2022, secciones del .docx | Sí | Sí, `letterpaper` y geometry | No se usa A4 |
| Margen superior, inferior y exterior 2,54 cm; interior 3,60 cm | Guía 2022, prosa y `pgMar` (1440 y 2041 twips) | Sí, márgenes simétricos | Sí, `inner=3.60cm`, `outer=2.54cm` | El .docx de ejemplo deja el margen ancho siempre a la izquierda y no activa `mirrorMargins`. Se activó porque la misma guía habla de margen interna y de número distinto en pares e impares |
| Fuente Arial | Guía 2022. Los ejemplos visibles están en Arial 11 y los capítulos de ejemplo en Arial 20 | Sí | TeX Gyre Heros | Los estilos de fábrica del .docx dicen Cambria 13. Los párrafos de ejemplo ya están forzados a Arial. Arial propietario no está en TeX Live |
| Interlineado 1,5 y otro Enter después del punto aparte | Guía 2022. El ejemplo separa párrafos con un párrafo vacío y espaciado 0 | Espacio posterior de 16 pt, sin párrafos vacíos | `\parskip` de 0,9 líneas e interlineado 1,5 | Misma separación, más fácil de editar |
| Título de capítulo Arial 20, negrita, 100 pt antes, 24 pt después, interlineado sencillo | Guía 2022. El ejemplo «Introducción» está en Arial 20 negrita; el espaciado 100/24 está en la prosa, no en ese párrafo (tiene antes 0) | Estilo Título 1 | `titlesec`, 20 pt, 100 pt / 24 pt | Se aplicó la instrucción escrita |
| Subtítulo 16 pt; el siguiente 14 pt negrita | Ejemplo del .docx: Subtítulo Arial 16, antes 20 pt, después 24 pt; el de 14 pt va en negrita | Título 2 y Título 3 | `\section` 16 pt; `\subsection` 14 pt negrita | |
| Número de página en el encabezado, Arial 11; impar a la derecha, par a la izquierda | Guía 2022 | Sí, encabezados pares e impares | `fancyhdr` | El .docx de ejemplo no tiene `evenAndOddHeaders`. Se activó para cumplir la prosa |
| Cubierta y portada sin número visible, pero contadas | Guía 2022 | Sección inicial sin encabezado | `\thispagestyle{empty}` sobre numeración romana | |
| Preliminares en romanos mayúsculos; arábigos desde la Introducción, en 1 | Guía 2022 | Sí | `\pagenumbering{Roman}` y `\mainmatter` | |
| Capítulo en página impar | Guía 2022 y el encargo | Salto de página impar | `book`, `openright` | |
| Sin numeración 13A ni «bis» | Guía 2022 | Sí | Sí | Capítulos 1 a 5 |
| Resumen analítico, máximo 250 palabras, 3 a 7 palabras clave | Guía 2022, NTC 1486 / ISO 214, citadas por la guía | 241 palabras, 6 palabras clave | El mismo texto | Objetivos, diseño, métodos, resultados y conclusiones |
| Abstract y Keywords, traducción del resumen | Guía 2022 | 233 palabras, 6 keywords | `otherlanguage{english}` | |
| Orden: cubierta, portada, resumen, abstract, contenido, listas, introducción, capítulos, conclusiones, anexos, referencias | Guía 2022. Dedicatoria y agradecimientos son opcionales y se omitieron. Las tesis en español del programa meten planteamiento, justificación, objetivos y política dentro de la introducción | El mismo orden | El mismo orden | No se copió el texto de esas tesis |
| Fórmula de la portada | Guía 2022: «Tesis o trabajo de investigación presentada(o) como requisito parcial para optar al título de:» | Texto de la guía | Texto de la guía | El encargo traía una frase más corta. Ganó la guía. El título del grado lleva tilde, como la Resolución 003253: «MAGÍSTER…» |
| Figura: leyenda y fuente al pie, sin «No.» ni «#» | Guía 2022. El ejemplo usa `Figura 1:` | `Figura 1:` y línea «Fuente:» | Igual, `labelsep=colon` | El encargo escribía «Figura 1.» con punto. Ganó el ejemplo de la guía |
| Tabla: título arriba; si continúa, encabezados repetidos y «Tabla X: (Continuación)» | Guía 2022 | Fila de encabezado marcada para repetirse. El anexo tiene 25 fuentes | `longtable` con `\caption*{Tabla \thetable: (Continuación)}` | |
| Ecuación centrada, número a la derecha | Guía 2022, ejemplo (1) | Tabla sin borde: fórmula y (1) | `equation` | Es el umbral real, log2(1,5) = 0,585 |
| Redacción impersonal; otro idioma en cursiva | Guía 2022 | Sí | Sí | |
| Extensión sugerida: 3 a 6 capítulos, introducción de 2 a 4 páginas, no más de 120 | Guía 2022, «se sugiere» | 5 capítulos numerados, conclusiones y anexos incluidos | Igual | No se bloquea la compilación si crece |
| Bibliografía obligatoria; cualquier norma, con rigor | Guía 2022. Cita NTC 1486 y el Manual APA, sin fijar edición | APA 7, orden alfabético, sangría francesa | `biblatex` `style=apa` | Vancouver queda comentado en la clase. No hay estilo obligatorio del programa |
| Nombre del PDF `JuanDiego_PolancoFranco_2026` | Instructivo PDF/A, 2022: `Nombres_Apellidos_Año`, sin tildes ni ñ | Documentado | Documentado | El manual del 12 de febrero de 2026 muestra `NombresApellidos_Programa_Año`. Confirmar cuál pide el formulario |
| Metadatos: autor «Apellidos, Nombres»; palabras clave con punto y coma; asunto | Instructivo PDF/A, 2022 | Propiedades del .docx | `hyperref` | El encargo escribía «;;». El instructivo dice punto y coma |
| PDF/A, sin bloquear copiar ni imprimir | Instructivo PDF/A, 2022 | Pasos en `word/INSTRUCCIONES-WORD.md` | Overleaf no garantiza PDF/A. El PDF/A se exporta desde Word | |
| Entrega solo digital, sin ejemplar físico | Circular 008. El PDF descargado de la URL indicada está fechado el 17 de marzo de 2021, no en 2025 | Documentado en el README | Igual | El archivo de esa URL es la circular de 2021. El contenido coincide: el repositorio es el canal y no hay recepción física |
| Autoarchivo, carta de SABio, correo personal, 5 a 8 días, tipología | Manual de autoarchivo, 12 de febrero de 2026, y preguntas frecuentes de 2025 | README, sin automatizar el trámite | Igual | La URL del manual devolvió un PDF, no un .docx |
| 60 créditos, investigación, convenio, Manizales, título de magíster | Resolución MEN 003253 del 2 de marzo de 2021, PDF de 3 páginas | Una línea en la portada | La misma línea | |
| SNIES 103917 y facultad actual | Página del programa, Facultad de Inteligencia Artificial e Ingenierías, consultada el 6 de noviembre de 2026 | Portada | Portada | La resolución descargada no trae el código SNIES. Las tesis de 2023 dicen «Facultad de Ingenierías». Se usó el nombre de la página vigente |
| Tesis de Maestría I = 9 créditos y Tesis II = 10 | Encargo. No está en la Resolución 003253 descargada | No se escribió | No se escribió | [NO ENCONTRADO EN NORMA OFICIAL — confirmar con coordinación] |
| Logos del convenio | No están en la guía ni en la resolución | No se inventaron | No se inventaron | [NO ENCONTRADO EN NORMA OFICIAL — confirmar con coordinación] |
| Plantilla LaTeX oficial | No existe en las fuentes | — | Clase propia que imita el Word | |
| Acta, hoja de vida, ISBN, Turnitin, número de jurados, nota mínima, plazo interno, inglés de Caldas | No están en las fuentes leídas | No se exigen | No se exigen | El Reglamento Estudiantil (Acuerdo 16 de 2007) no se volvió a descargar en esta sesión. El encargo cita los artículos 24 y 27; no se convirtieron en requisitos del archivo |

## Fuentes guardadas

- `assets/guia-presentacion-tesis.docx` — https://biblio.ucaldas.edu.co/wp-content/uploads/2022/07/guia-presentacion-tesis-trabajos-grado.docx
- `assets/guia-creacion-pdf.pdf` — https://biblio.ucaldas.edu.co/wp-content/uploads/2022/07/guia-creacion-archivos-pdf.pdf
- `assets/circular-008.pdf` — https://biblio.ucaldas.edu.co/wp-content/uploads/2025/03/Circular-008-V.-Acad-Entrega-y-publicacion-de-tesis-y-trabajos-de-grado-en-Repositorio-Inst.pdf
- `assets/manual-autoarchivo.pdf` — https://biblio.ucaldas.edu.co/wp-content/uploads/2026/02/Manual-autoarchivo-Repositorio-institucional-febrero-12_2026.docx-1.pdf
- `assets/preguntas-frecuentes.pdf` — https://biblio.ucaldas.edu.co/wp-content/uploads/2025/04/PREGUNTAS-FRECUENTES-REPOSITORIO-INSTITUCIONAL-UNIVERSIDAD-DE-CALDAS.pdf
- `assets/resolucion-003253.pdf` — https://viceacademica.ucaldas.edu.co/wp-content/uploads/2025/07/Maestria-en-Bioinformatica-y-Biologia-Computacional-Res.-003253-02-03-2021.pdf
- Página del programa: https://ingenierias.ucaldas.edu.co/oferta-academica/maestria-en-bioinformatica-y-biologia-computacional/
- Autoarchivo: https://repositorio.ucaldas.edu.co/autoarchivo
- Guía de publicación: https://biblio.ucaldas.edu.co/guia-publicacion-repositorio-institucional/

## Tesis del programa usadas solo como diseño

Se miró la osamenta, no el contenido, de las tesis en español de Willington Londoño Ramos, Marcela Orjuela Rodríguez y, para el cierre en inglés, Simón Gaviria Orrego. Orjuela numera la introducción y cierra cada capítulo con referencias. La guía oficial deja la introducción fuera de la numeración 1, 2, 3 y junta las referencias al final. Se siguió la guía. De Orjuela y de Londoño se conservó, dentro de la introducción, el planteamiento, la justificación, la política pública y los objetivos.
