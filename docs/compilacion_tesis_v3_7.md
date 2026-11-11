# Compilación de la tesis v3.7

El ZIP para Overleaf es `tesis-ucaldas/overleaf_v3_7.zip`, armado desde `tesis-ucaldas/latex/`. Overleaf pidió inicio de sesión y en este entorno no hay una sesión abierta, así que no se subió el proyecto a overleaf.com. La comprobación se hizo con ese mismo ZIP, con LuaLaTeX y Biber, que es el par que declaran `main.tex` y `latexmkrc` (`pdf_mode = 4`).

Motor: LuaHBTeX 1.24.0 (TeX Live 2026) y Biber 2.22, en TinyTeX del Mac del autor. Tras la enmienda de ESM-2 se recompiló igual. `latexmk -lualatex` terminó con código 0. El log tiene 0 líneas de error (`^!`), 0 referencias indefinidas y 0 citas indefinidas. El PDF tiene 58 páginas.

En el PDF están la tabla de correspondencia (con «Frente a HC salen 3 (CP) y 1 (CN)»), el capítulo «Desviaciones del preregistro y trabajo futuro» y la figura 3.3. Esa figura titula «60 muestras en tres bloques»: CP en 1–20, CN en 21–40 y HC en 41–60. El pie cita `results/v3_7/preanalitica_orden.tsv`.

Dos arreglos hicieron falta para llegar a ese log. `CON__` y `REV__` en métodos se leían como subíndice. En `references.bib`, la entrada `zhang2014` cerraba el autor después de Wang, Jing, y Biber rechazaba el fichero entero. La entrada de García-Foncillas 2026 no se tocó.

Últimas líneas de `main.log`:

```
Package rerunfilecheck Info: File `main.out' has not changed.
(rerunfilecheck)             Checksum: 1A083925553C37824DF38B69115BFCD6;9846.
Package logreq Info: Writing requests to 'main.run.xml'.

Here is how much of LuaTeX's memory you used:
 26043 strings out of 475522
 125081,1977958 words of node,token memory allocated
 791 words of node memory still in use:
   6 hlist, 2 vlist, 2 rule, 2 glue, 4 kern, 1 glyph, 45 attribute, 73 glue_spec, 23 attribute_list, 1 write, 15 pdf_action nodes
   avail lists: 1:2,2:1639,3:615,4:199,5:780,6:151,7:7670,8:50,9:1127,10:27,11:507
 48334 multiletter control sequences out of 65536+600000
 70 fonts using 7542287 bytes
 75i,13n,105p,1070b,4980s stack positions out of 10000i,1000n,20000p,200000b,200000s
Output written on main.pdf (58 pages, 305172 bytes).

PDF statistics: 790 PDF objects out of 1000 (max. 8388607)
 697 compressed objects within 7 object streams
 168 named destinations out of 1000 (max. 1048576)
 489 words of extra memory for PDF output out of 10000 (max. 100000000)
```
