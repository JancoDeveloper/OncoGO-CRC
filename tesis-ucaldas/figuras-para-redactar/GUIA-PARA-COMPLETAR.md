# Guía para completar la tesis OncoGO-CRC

Carpeta de figuras: `/Users/jancod/GitHub/OncoGO-CRC/tesis-ucaldas/figuras-para-redactar/`. Rutas de datos relativas a `/Volumes/JANCODATA/OncoGO-CRC/` salvo `.tex`/tests/scripts (repositorio GitHub). Citas con `\textcite{clave}` como en el documento actual.

## Cómo escribir para que suene a ti
- Frases cortas, en pasado: se observó, en este análisis, la tabla X muestra.
- Primero el resultado con número y figura; luego la lectura; el límite en el mismo párrafo.
- Nada de adjetivos apilados, preguntas retóricas, 'en el vasto campo de', sinónimos de relleno.
- Reescribe cada frase guía con tus palabras y más corta. No pegues las guías tal cual.

## Introducción (cap. 1, apertura)
Extensión: 2–4 págs. Hoy ~2,6. Máximo +1 pág.

Idea a desarrollar:
- El problema: CRC se detecta tarde; la colonoscopia es el estándar y una prueba en sangre sería útil
- La pregunta de la tesis tal como está en el título: qué parte del tumor se deja ver en la sangre
- Por qué datos públicos: no hay reclutamiento, se reanalizan tres cohortes sin fusionarlas

Qué decir (guías, reescríbelas):
- El cáncer colorrectal … (dato de incidencia en Colombia, con cita).
- La pregunta no es si existe un marcador, sino cuánto del cambio del tejido llega al suero.
- Se usaron tres conjuntos públicos: PDC000116, PXD034575 y PXD055218.
- El control que importa es la colonoscopia negativa, no solo la persona sana.

Dónde buscar:
- `tesis-ucaldas/latex/chapters/06-introduccion.tex`
- `docs/RESULTADOS_Y_LIMITACIONES_TESIS.md`

Bibliografía: garciafoncillas2026, ley1384, rai2025, malmstrom2025

## 1.1 Proteómica por espectrometría de masas: TMT, LFQ y DIA
Extensión: 3–4 págs.

Idea a desarrollar:
- Qué mide cada modo de adquisición y en qué escala queda el dato
- CPTAC es TMT con log-ratio; Tanaka es LFQ de MaxQuant; Rai es DIA de suero
- Por qué esas escalas no se pueden poner en una sola matriz (justifica no fusionar)

Qué decir (guías, reescríbelas):
- En TMT cada canal es una muestra dentro de un plex; el valor es una razón frente a la referencia.
- En LFQ la intensidad es por corrida; el cero de MaxQuant significa no detectado.
- En DIA se fragmenta todo en ventanas; el suero tiene rango dinámico amplio y pocas proteínas de tejido.
- Por eso en este trabajo las cohortes se comparan por signo y no se fusionan.

Dónde buscar:
- `data/manifests/*.json (JANCODATA)`
- `CONCLUSIONES_CERTIFICADAS.md (párrafo de escalas)`

Bibliografía: zhang2014, vasaikar2019, cox2008, tanaka2024, rai2025, lesur2023

Figuras/tablas: Tabla 2.1 ya existente (datasets)

## 1.2 Contraste diferencial: log2FC, ausencias, treat y Benjamini–Hochberg
Extensión: 3–4 págs.

Idea a desarrollar:
- Qué es log2FC y por qué el umbral 0,585 (1,5 veces)
- Por qué no imputar: una ausencia no es un cero biológico
- t moderada (limma), treat contra un umbral, y FDR de BH

Qué decir (guías, reescríbelas):
- log2FC = 0,585 equivale a un cambio de 1,5 veces.
- treat prueba que el efecto supera el umbral, no solo que sea distinto de cero.
- BH controla la proporción esperada de falsos entre los significativos.
- Sin imputación, n_obs cambia por proteína; CHMP4A es el ejemplo (16 de 40).

Dónde buscar:
- `results/tesis_cierre/cierre_audit.json (alpha, treat_lfc, min_por_grupo)`
- `results/tesis_cierre/crc_vs_sano.tsv (columnas n_obs, padj)`

Bibliografía: smyth2004, mccarthy2009, benjamini1995, cox2008

Figuras/tablas: Equation del umbral ya en métodos

## 1.3 Ontología génica, GSEA y lectura de vías
Extensión: 2–3 págs.

Idea a desarrollar:
- GO como vocabulario; GSEA como prueba sobre una lista ordenada
- NES y FDR ≤ 0,25 del GSEA
- Por qué un suero con ~700 proteínas cubre pocos genes por término

Qué decir (guías, reescríbelas):
- GSEA no usa un corte; ordena todas las proteínas por su estadístico.
- Un término con menos de diez genes medidos no se puede leer en el suero.
- La cobertura del suero limita lo que se puede afirmar sobre vías del tumor.

Dónde buscar:
- `results/gsea/ (JANCODATA, archivos cptac_*)`
- `cierre_audit.json → vias_suero`

Bibliografía: ashburner2000, subramanian2005, uhlen2015

Figuras/tablas: Fig3-8 (referencia cruzada)

## 1.4 Del tejido a la sangre: secreción, fase aguda y neoantígenos
Extensión: 3–4 págs.

Idea a desarrollar:
- Qué proteínas pasan a sangre (secretadas, fuga, vesículas)
- Respuesta de fase aguda: SAA1, LBP, CRP, HP como señal del huésped
- Neoantígeno vs. proteína diferencial: péptido mutado presentado por MHC

Qué decir (guías, reescríbelas):
- Una proteína alta en tumor no implica que suba en suero.
- La respuesta de fase aguda aparece en inflamación y en cáncer; no distingue el origen.
- Rai reporta 19 péptidos; la validación PRM cubre uno.
- Xiang mide inmunopeptidoma de tejido; Hembrough, péptidos de un paciente en vesículas.

Dónde buscar:
- `results/tesis_cierre/peptidos_neoantigenos.tsv`
- `results/tesis_cierre/fichas_proteinas.tsv (secreción, GO)`

Bibliografía: malmstrom2025, uhlen2015, seo2026, rai2025, xiang2026, hembrough2026, dhami2026

Figuras/tablas: Tabla3-6 (referencia cruzada)

## 1.5 Estado del arte: biomarcadores proteicos de CRC en sangre
Extensión: 4–5 págs., una tabla larga

Idea a desarrollar:
- Tabla: autor, matriz, comparación, n, hallazgo, qué no resuelve
- Separar estudios contra sanos y contra otra enfermedad intestinal
- Cerrar con el hueco: casi nadie compara con colonoscopia negativa

Qué decir (guías, reescríbelas):
- Verathamjamras (plasma LFQ) reporta LRG, C9, FN y A1AT.
- Ögren combina proteínas y metabolitos en plasma para varios cánceres.
- Yu propone candidatas tempranas; Li sigue el plasma en el tiempo.
- Robles hace metaanálisis de datos públicos de CRC.
- Rai incluye colonoscopia negativa; es la base del diseño de tres grupos.

Dónde buscar:
- `references.bib`
- `results/tesis_cierre/candidatas_literatura.tsv (columna fuente)`

Bibliografía: verathamjamras2023, ogren2026, yu2026, li2024, robles2024, rai2025, seo2026, malmstrom2025

Figuras/tablas: Tabla 1.x nueva: estado del arte (hacerla a mano desde los artículos)

## 2.1 Conjuntos de datos y diseño
Extensión: 2–3 págs.

Idea a desarrollar:
- Una ficha por conjunto: acceso, tecnología, muestras, grupos
- Por qué no se fusionan
- Leyenda de columnas de Tanaka: P primario, M metástasis hepática, N mucosa

Qué decir (guías, reescríbelas):
- PXD055218: 20 cáncer, 20 colonoscopia negativa, 20 sanos.
- PXD034575: 16 pacientes con P, M y N; no hay pólipo.
- CPTAC: 94 pares tumor–normal del mismo participante (cuantificación Unshared Log Ratio).
- Los pólipos no entran en ningún contraste.

Dónde buscar:
- `cierre_audit.json (n_por_grupo, pxd034575.leyenda)`
- `CONCLUSIONES_CERTIFICADAS.md`

Bibliografía: vasaikar2019, tanaka2024, rai2025

Figuras/tablas: Tabla 2.1 existente; Fig2-1_flujo_de_analisis.png

## 2.2 Preprocesamiento y reglas fijadas antes del contraste
Extensión: 2–3 págs.

Idea a desarrollar:
- Lista de reglas escritas antes de ver resultados
- Cero → ausencia; log2 de intensidad positiva; sin imputación
- Mínimo 8 observaciones por grupo; semilla 20261031

Qué decir (guías, reescríbelas):
- Las reglas se fijaron en el archivo de auditoría antes de correr el contraste.
- No se imputaron ausencias.
- Las marcas de fase aguda, hemólisis e EII se escribieron en el catálogo antes del contraste.

Dónde buscar:
- `cierre_audit.json`
- `tests/test_tesis_cierre.py (GitHub)`
- `scripts/ (GitHub)`

Bibliografía: cox2008

Figuras/tablas: Anexo: tabla de reglas

## 2.3 Contraste, especificidad y fichas de proteína
Extensión: 2–3 págs.

Idea a desarrollar:
- t moderada + treat + BH
- Regla de especificidad: AUC > 0,70, IC bajo > 0,5, mayor que enfermedad vs. sano
- Cómo se arma la ficha: UniProt, GO, secreción, AlphaFold

Qué decir (guías, reescríbelas):
- El AUC se calculó por proteína, sin combinar en panel.
- El intervalo es bootstrap con 2000 remuestreos.
- El AUC reportado es el mayor entre el valor y su complemento.

Dónde buscar:
- `especificidad_candidatas.tsv`
- `fichas_proteinas.tsv`

Bibliografía: smyth2004, mccarthy2009, benjamini1995, ashburner2000, jumper2021, uhlen2015

Figuras/tablas: Equation de la regla

## 2.4 Límites operativos
Extensión: 1 pág. (hoy 37 palabras)

Idea a desarrollar:
- Qué no se hizo y por qué
- Sin pureza tumoral, sin canales TMT, sin metaanálisis de efectos aleatorios
- PXD034575 no es la réplica PXD019504; sin crudos de PRM

Qué decir (guías, reescríbelas):
- No hay pureza tumoral en la clínica pública.
- No se reanalizaron crudos de péptidos.
- El suero tiene 20 personas por grupo; un AUC alto con ese n es inestable.

Dónde buscar:
- `CONCLUSIONES_CERTIFICADAS.md → Qué queda fuera`

## 3.1 El suero no separa el cáncer de la otra enfermedad
Extensión: 4–5 págs.

Idea a desarrollar:
- Primero el conteo (0, 3, 1), luego los volcanes
- Cajas de SAA1, LBP, CHMP4A, APOA4 por los tres grupos
- CHMP4A: la mayoría de sueros sin medición en col. negativa y sanos

Qué decir (guías, reescríbelas):
- La figura 3.1 muestra que ninguna proteína pasa los dos cortes en cáncer frente a colonoscopia negativa.
- En la figura 3.2, SAA1 y LBP suben en cáncer y también en colonoscopia negativa frente a sanos.
- CHMP4A se midió en 13 de 20 cánceres, 4 de 20 colonoscopias negativas y 3 de 20 sanos.
- APOA4 baja en la enfermedad intestinal y no distingue el cáncer de esa enfermedad.

Dónde buscar:
- `results/tesis_cierre/crc_vs_*.tsv, enfermedad_vs_sano.tsv`
- `data/raw/plasma/plasma_crc_matrix.tsv`

Bibliografía: rai2025, malmstrom2025

Figuras/tablas: Fig3-1_volcanes_suero_PXD055218.png; Fig3-2_cajas_SAA1_LBP_CHMP4A_APOA4.png; Tabla3-1_conteos_suero.csv; Tabla3-2_n_y_medianas_cajas.csv

## 3.2 El área bajo la curva es más laxa que el FDR
Extensión: 4–5 págs.

Idea a desarrollar:
- Ocho candidatas pasan la regla de AUC; ninguna pasa FDR en ese contraste
- Cinco con marca de fase aguda
- Mostrar todas las curvas y el gráfico de AUC con IC

Qué decir (guías, reescríbelas):
- La figura 3.3 muestra las curvas ROC de las ocho candidatas, una por proteína.
- TPM4 tiene AUC 0,80 con solo 5 cánceres y 3 colonoscopias negativas medidas.
- En la figura 3.4 el intervalo de varias candidatas toca 0,5–0,6.
- Un AUC alto con 20 personas no reemplaza el control de falsos.

Dónde buscar:
- `especificidad_candidatas.tsv`
- `control_auc_recalculado_vs_tesis_cierre.csv (los 8 AUC se reproducen)`

Bibliografía: rai2025, verathamjamras2023, ogren2026

Figuras/tablas: Fig3-3_ROC_ocho_candidatas.png; Fig3-4_AUC_candidatas_literatura.png; Tabla3-3_ocho_candidatas_AUC.csv

## 3.3 El tejido, leído con la leyenda del artículo
Extensión: 4–5 págs.

Idea a desarrollar:
- Volcanes P vs N, M vs N, M vs P con la leyenda correcta
- Dispersión primario vs. metástasis por paciente
- Concordancia con CPTAC (ya en tabla)

Qué decir (guías, reescríbelas):
- La figura 3.5 muestra 58, 160 y 16 proteínas significativas en los tres contrastes.
- En la figura 3.6 el signo coincide en 708 de 876 proteínas.
- Solo el 46,2 % de las que coinciden tiene menor magnitud en el primario.
- No es un gradiente adenoma–carcinoma; no hay pólipos.

Dónde buscar:
- `tanaka_*_vs_*.tsv`
- `cierre_audit.json → pxd034575`

Bibliografía: tanaka2024, vasaikar2019, zhang2014

Figuras/tablas: Fig3-5_volcanes_tejido_Tanaka_PXD034575.png; Fig3-6_tanaka_primario_vs_metastasis.png; Anexo: TablaA_tanaka_*_significativas.csv

## 3.4 Las vías del tumor casi no están en este suero
Extensión: 2–3 págs.

Idea a desarrollar:
- 12–13 términos GO evaluables por contraste
- Correlación de NES negativa o cerca de cero; p nulo 0,30–0,93
- Decir cuáles términos (sacar la lista del GSEA de CPTAC)

Qué decir (guías, reescríbelas):
- La figura 3.8 muestra la correlación de NES entre tejido y suero por contraste.
- Ninguna se distingue del nulo de 200 permutaciones.
- Con este número de términos la conclusión es de ausencia de evidencia, no de evidencia de ausencia.

Dónde buscar:
- `cierre_audit.json → vias_suero`
- `results/gsea/cptac_*_gsea_go.tsv`

Bibliografía: subramanian2005, ashburner2000

Figuras/tablas: Fig3-8_vias_tejido_en_suero.png; Tabla3-4_vias_tejido_en_suero.csv

## 3.5 Catálogo de literatura, péptidos y neoantígenos
Extensión: 3–4 págs.

Idea a desarrollar:
- Mapa de veredictos del catálogo en las cuatro cohortes
- 7 replicadas en CPTAC; marcas de hemólisis y EII
- Tabla de péptidos como límite, no como hallazgo

Qué decir (guías, reescríbelas):
- La figura 3.7 resume el veredicto de las 37 candidatas en cada contraste.
- En CPTAC se detectan 33 y siete repiten el signo de la literatura.
- La matriz de grupos proteicos no conserva la secuencia mutada.

Dónde buscar:
- `auditoria_*.tsv`
- `peptidos_neoantigenos.tsv`

Bibliografía: rai2025, xiang2026, hembrough2026

Figuras/tablas: Fig3-7_auditoria_catalogo_heatmap.png; Tabla3-5_auditoria_catalogo_por_cohorte.csv; Tabla3-6_peptidos_neoantigenos.csv

## 3.6 Discusión
Extensión: 5–7 págs.

Idea a desarrollar:
- Un párrafo por estudio comparado: qué encontraron, qué se vio aquí, por qué difiere
- Señal del huésped vs. señal del tumor
- Limitaciones con calma: n=20, sin CEA, sin pureza, sin péptidos
- Qué haría un estudio siguiente

Qué decir (guías, reescríbelas):
- Verathamjamras reporta C9, LRG y A1AT; aquí pasan el AUC y llevan marca de fase aguda.
- Ögren usa CEA en el panel; esta matriz no lo mide.
- Malmström y Seo describen señal del huésped; encaja con SAA1 y LBP.
- Rai reporta AUC 0,758 con CEA; no se puede auditar sin el inmunoensayo.
- Un estudio con más colonoscopias negativas resolvería BAIAP3.

Dónde buscar:
- `09-resultados.tex (Discusión actual)`
- `CONCLUSIONES_CORRIDA.md`

Bibliografía: verathamjamras2023, ogren2026, yu2026, malmstrom2025, seo2026, rai2025, xiang2026, hembrough2026, robles2024

## Anexos 5.1–5.3
Extensión: 5–10 págs.

Idea a desarrollar:
- Catálogo completo de 37 candidatas con fuente
- Fichas de 21 proteínas
- Listas de significativas de Tanaka
- Versiones de software y comandos para reproducir

Qué decir (guías, reescríbelas):
- Tabla A.1 reúne el catálogo de candidatas con su fuente y marca.
- Tabla A.2 muestra la ficha de cada proteína evaluada.
- Los scripts gen.py y gen2.py regeneran las figuras de este capítulo.

Dónde buscar:
- `figuras-para-redactar/TablaA_*.csv`
- `pip freeze / renv del repositorio`

Figuras/tablas: TablaA_catalogo_37_candidatas_literatura.csv; TablaA_fichas_21_proteinas.csv; TablaA_tanaka_*_significativas.csv
