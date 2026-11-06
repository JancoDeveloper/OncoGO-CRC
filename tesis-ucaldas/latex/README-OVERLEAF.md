# Overleaf

El archivo principal es `main.tex`. El compilador es **LuaLaTeX** y la bibliografía es **Biber**. Overleaf lee esta línea al inicio de `main.tex`:

```
% !TeX program = lualatex
% !BIB program = biber
```

## Subir el proyecto

1. Comprima esta carpeta (`main.tex`, `ucaldas-tesis.cls`, `chapters/`, `references.bib`, `figures/`, `latexmkrc`). El zip ya armado está en `tesis-ucaldas/oncogo-crc-ucaldas-overleaf.zip`.
2. En [overleaf.com](https://www.overleaf.com): **New Project → Upload Project** y elija ese zip.
3. Si el menú no queda en LuaLaTeX, ábralo y seleccione **LuaLaTeX**.
4. Pulse **Recompile**. La tabla de contenido queda bien a la segunda compilación.

No suba el repositorio de análisis ni el disco de datos. Overleaf solo necesita este documento.

## Qué imita y qué no

La universidad no publica plantilla LaTeX. Esta clase imita la guía Word de la Biblioteca (2022): carta, margen interno 3,60 cm, margen externo 2,54 cm, interlineado 1,5, capítulos en página impar, romanos en preliminares y arábigos desde la Introducción.

La fuente es **TeX Gyre Heros**, la familia libre de métrica afín a Arial que trae Overleaf. Arial propietario no está en TeX Live. Esa desviación está en `VALIDACION-FORMATO.md`.

## Bibliografía

El valor por defecto es APA 7 (`style=apa` en `ucaldas-tesis.cls`). La guía de 2022 no impone estilo. Para Vancouver, cambie `style=apa` por `style=vancouver`, borre la línea `\DeclareLanguageMapping{spanish}{spanish-apa}` y vuelva a compilar.

`references.bib` trae solo las fuentes cuyo DOI se confirmó en Crossref el 31 de octubre de 2026, más la Ley 1384 de 2010.

## Nombre del PDF

Al descargar el PDF desde Overleaf, renómbrelo según el instructivo de 2022:

`JuanDiego_PolancoFranco_2026.pdf`

Sin tildes. La ñ no aparece en este nombre. El manual de autoarchivo del 12 de febrero de 2026 muestra otro patrón (`NombresApellidos_Programa_Año`). Confirme con el repositorio cuál exige el formulario vigente. El PDF/A se produce desde Word con el instructivo de la Biblioteca; Overleaf no garantiza PDF/A.
