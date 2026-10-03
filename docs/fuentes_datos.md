# Fuentes de datos de v2

Una ficha por fuente. El estado solo puede ser uno de tres:

- **usada**: descargada, parseada y en uso en algún resultado de `results/v2/`.
- **no_disponible**: la descarga falló o el recurso no responde. Se anota el
  motivo exacto y qué análisis queda sin evaluar por su causa.
- **requiere_acceso**: existe, pero hace falta una solicitud, un acuerdo de uso
  o credenciales que este trabajo no tiene.

Ninguna fuente marcada `no_disponible` o `requiere_acceso` se sustituye por
datos simulados. Los datos sintéticos solo existen dentro de `tests/`, en
ficheros con el prefijo `SINTETICO_`.

El estado legible por máquina vive junto a cada descarga, en
`external/<fuente>/estado.json`. Esta página es la versión para leer.

---

## Cohortes del certificado (heredadas de v1, sin volver a descargar)

### PDC000116 — CPTAC cáncer de colon, cohorte de descubrimiento

- **Estado**: usada
- **Qué es**: proteoma de tumor y mucosa normal del mismo paciente, marcaje
  isobárico TMT, cuantificación del Proteomic Data Commons.
- **Qué aporta a v2**: los 94 pares tumor–normal del contraste principal.
- **Columnas**: se usa exclusivamente la familia `Unshared Log Ratio`, que son
  los péptidos no compartidos. La matriz **ya viene en log2**, así que no se
  vuelve a transformar. Cada espécimen aparece dos veces en el fichero, una por
  familia de columnas, y las dos familias nunca se mezclan.
- **Trampa conocida**: las filas `Mean`, `Median` y `StdDev` no son proteínas y
  se descartan.

### PXD034575 — Tanaka 2024, réplica independiente

- **Estado**: usada
- **Qué es**: LFQ con MaxQuant de tumor, mucosa normal y pólipo del mismo
  paciente, 16 pacientes.
- **Qué aporta a v2**: la réplica de dirección (R1), el brazo pareado del nulo
  exacto de C1, el brazo de pólipo de C5 y las 16 unidades de validación de C6.
- **Nombres de columna**: `LFQ intensity CRC<n><M|N|P>`, donde **M es el tumor**
  (maligno), N el tejido normal y P el pólipo. No es `T`.
- **Trampa conocida**: el cero de MaxQuant significa ausencia de identificación,
  no abundancia cero. No se imputa.

### PXD055218 — suero DIA, cáncer frente a sano

- **Estado**: usada
- **Qué es**: 20 sueros de cáncer colorrectal y 20 de control sano, DIA sin
  enriquecimiento.
- **Qué aporta a v2**: el lado de circulación de la brecha tejido→plasma y de la
  ventana de detección de C8.
- **Limitación estructural**: no hay ningún paciente con tejido y sangre a la
  vez. Toda comparación tejido→plasma es entre cohortes, nunca intraindividual.

---

## Fuentes nuevas de v2

### Gene Ontology — Biological Process

- **Estado**: usada
- **Fichero**: `metadata/gene_ontology/go_gene_long.tsv`, el mismo de v1
- **Qué aporta**: una de las tres colecciones de la cláusula de vías (C4).
  11 248 términos BP.
- **Resultado en v2**: de 132 términos significativos en CPTAC, solo 6 llegan a
  ser evaluables en la réplica, por debajo del mínimo preregistrado de 10. La
  colección queda **no evaluable**.

### Reactome

- **Estado**: usada
- **Script**: `scripts/download/colecciones_vias.py`
- **Fichero**: `external/vias/reactome_pathways.gmt`
- **Origen**: `https://reactome.org/download/current/ReactomePathways.gmt.zip`
- **Licencia**: CC BY 4.0
- **Cita**: Milacic M et al. The Reactome Pathway Knowledgebase 2024. *Nucleic
  Acids Research*.
- **Contenido**: 2 868 términos, 142 657 pares término-gen.
- **Resultado en v2**: la única de las tres colecciones evaluable. 449 términos
  significativos en CPTAC, 51 evaluables en la réplica, concordancia 47/51.

### MSigDB Hallmark

- **Estado**: usada
- **Script**: `scripts/download/colecciones_vias.py`
- **Fichero**: `external/vias/msigdb_hallmark.gmt`
- **Licencia**: CC BY 4.0 para la colección Hallmark. El resto de MSigDB tiene
  condiciones propias y no se usa.
- **Cita**: Liberzon A et al. The Molecular Signatures Database Hallmark Gene Set
  Collection. *Cell Systems* 2015.
- **Contenido**: 50 términos, 7 322 pares término-gen.
- **Resultado en v2**: 7 términos significativos en CPTAC, 2 evaluables en la
  réplica. **No evaluable**.

### PXD019504 — tercera cohorte de tejido pareado

- **Estado**: usada
- **Script**: `scripts/download/pride_cohorte_extra.py`
- **Fichero**: `external/pride/PXD019504/proteinGroups.txt`
- **Qué es**: LFQ con MaxQuant de tumor y mucosa normal del mismo paciente.
- **Qué aporta a v2**: la tercera cohorte del meta-análisis REML de C10. Con dos
  cohortes la heterogeneidad no se puede estimar; con tres, sí.
- **Nombres de columna**: `LFQ intensity CRC<n><N|T>`, donde **T es el tumor**.
  Es la convención contraria a PXD034575, que usa `M`. Confundirlas invertiría
  el signo de una cohorte entera, y el error no daría ningún aviso.
- **Identificadores**: las accesiones vienen sueltas, separadas por `;`. Otros
  depósitos usan el formato `sp|P12345|NOMBRE_HUMAN`, donde la accesión es el
  campo del medio. Quedarse con el último campo daría el nombre de la proteína.
- **Limitación**: es LFQ, igual que PXD034575 y a diferencia del TMT de
  PDC000116. La heterogeneidad que mide C10 podría venir de esa diferencia de
  tecnología, y por eso se calcula también el I² solo entre las dos cohortes LFQ.

### Auditoría de independencia entre cohortes

- **Estado**: usada
- **Script**: `scripts/download/independencia_cohortes.py`
- **Fichero**: `external/independencia/estado.json`, con las correlaciones cruzadas en
  `external/independencia/correlaciones_<a>_vs_<b>.tsv`
- **Qué es**: no es una fuente de datos nueva, es la comprobación de que las
  cohortes que entran en el meta-análisis son independientes de verdad.
- **Cómo decide**: dos niveles. El nivel A compara las sumas de verificación que
  publica PRIDE: ficheros compartidos entre dos depósitos significan muestras
  compartidas. El nivel B correlaciona cada muestra con cada muestra de las otras
  cohortes sobre las proteínas comunes; por encima de 0,98 con al menos 200
  proteínas comunes, se considera la misma muestra reanalizada. **El nivel B es el
  que decide**, porque un fichero puede haberse renombrado.
- **Por qué importa**: un meta-análisis que suma dos veces la misma muestra
  estrecha artificialmente los intervalos y presenta como replicación lo que es
  una repetición.
- **Regla dura**: sin auditoría no entra ninguna cohorte. La ausencia del fichero
  no se interpreta como «independientes».

### TCGA COAD/READ — CEA sérico, pureza y expresión de CEACAM5

- **Estado**: usada
- **Script**: `scripts/download/tcga_coadread.py`
- **Fichero**: `external/tcga_coadread/tcga_coadread_cea_pureza.tsv`
- **Origen**: BCR Biotab del GDC (`nationwidechildrens.org_clinical_patient_coad.txt`
  y `..._read.txt`), pureza ABSOLUTE del GDC, y la API REST de cBioPortal para
  el estudio `coadread_tcga_pan_can_atlas_2018`.
- **Licencia**: datos abiertos del GDC. No hay datos controlados en esta descarga.
- **Cita**: The Cancer Genome Atlas Network. Comprehensive molecular
  characterization of human colon and rectal cancer. *Nature* 2012. Y
  Hoadley KA et al. *Cell* 2018 para el PanCanAtlas.
- **Contenido**: 594 muestras de tumor primario. CEA pretratamiento en 380
  (64,0 %), pureza ABSOLUTE en 586, las dos en 376. Expresión RNA-seq de
  CEACAM5, CEACAM6 y CEACAM1.
- **Qué aporta**: el eje D1. Mide cuánto de la variación del CEA en sangre
  explica la expresión tisular de su propio gen, que es el caso mejor estudiado
  que existe de un marcador y su origen.
- **Columnas usadas**: `cea_level_pretreatment`, `ajcc_pathologic_tumor_stage`,
  `ajcc_metastasis_pathologic_pm`, `icd_o_3_site`, `gender`,
  `age_at_initial_pathologic_diagnosis`, la pureza y la fracción de ADN tumoral
  de ABSOLUTE, `MSI_SCORE_MANTIS` y `SUBTYPE` de cBioPortal.
- **Tres discrepancias con lo que se esperaba**, comprobadas contra el servicio
  en lugar de supuestas:
  1. El campo del CEA no se llama `preoperative_pretreatment_cea_level`. En el
     Biotab actual se llama `cea_level_pretreatment`.
  2. `HEAD` contra `api.gdc.cancer.gov/data/<id>` devuelve 400 aunque el fichero
     exista. Hay que comprobarlo con `GET`.
  3. El tarball de cBioPortal en S3 responde 403. La API REST sí funciona y es
     la que se usa.
- **CMS**: `requiere_acceso`. La clasificación consensuada vive en Synapse
  (syn2623706) y hace falta cuenta. El script queda preparado y no intenta
  autenticarse. Como sustituto declarado se usa el `SUBTYPE` del PanCanAtlas
  (CIN, MSI, GS, POLE), que **no es equivalente** a CMS1–CMS4.
- **Limitaciones**: el CEA es una determinación clínica de laboratorios
  distintos, sin un ensayo único; la expresión es RNA-seq y no proteína; y hay
  594 tumores con proteoma ausente, así que la pureza medida no sirve todavía
  para un R4 alternativo (eso depende de D2).

### Human Protein Atlas — concentraciones en sangre

- **Estado**: usada
- **Script**: `scripts/download/human_protein_atlas.py`
- **Fichero**: `external/hpa/hpa_concentracion_sangre.tsv`
- **Origen**: `https://www.proteinatlas.org/download/proteinatlas.tsv.zip`
- **Licencia**: CC BY-SA 4.0
- **Cita**: Uhlén M et al. Tissue-based map of the human proteome. *Science*
  2015.
- **Contenido**: 20 162 genes; 4 285 con concentración por espectrometría de
  masas y 457 por inmunoensayo. Del TSV completo solo se guardan las columnas
  necesarias; el resto no entra en el workspace.
- **Qué aporta**: la corrección C8. Permite estimar la probabilidad de detección
  en función de la concentración y medir la concordancia solo dentro de la
  ventana de detección.
- **Limitación**: las concentraciones son de sujetos de referencia sanos, no de
  los pacientes de PXD055218. Sirven para el orden de magnitud, no como medida
  del individuo. Además solo una fracción del proteoma tiene concentración
  publicada, y está sesgada hacia proteínas abundantes y ya estudiadas.

### Álvez et al. 2023 — proteoma de sangre pan-cáncer (Olink Explore)

- **Estado**: usada
- **Script**: `scripts/download/alvez_2023.py`
- **Ficheros**: `external/alvez_2023/alvez_2023_colorrectal_por_proteina.tsv` y
  `..._auc.tsv`
- **Origen**: información suplementaria del artículo, en
  `static-content.springer.com`. Se descargan tres suplementarios: el 3 (resumen
  de cohortes), el 6 (AUC por número de proteínas) y el 11 (datos de las
  figuras, con los resultados por proteína).
- **Licencia**: CC BY 4.0, acceso abierto.
- **Cita**: Álvez MB, Edfors F, von Feilitzen K, et al. Next generation
  pan-cancer blood proteome profiling using proteograph. *Nature Communications*
  2023; 14:4308. doi:10.1038/s41467-023-39765-y
- **Contenido**: 1 463 ensayos de Olink Explore en 221 pacientes de cáncer
  colorrectal (101 mujeres, 120 hombres; estadios I 26, II 68, III 99, IV 28)
  frente a 74 controles sanos del programa Wellness. 780 proteínas
  significativas: 328 suben y 452 bajan.
- **Qué aporta**: la comprobación externa que a v1 le faltaba. v1 encontró tres
  proteínas de plasma en 20 contra 20 sueros; aquí hay una cohorte diez veces
  mayor en la que mirar.
- **Resultado en v2**: de las tres proteínas de plasma de v1, solo LBP (P18428)
  está en el panel, y **no replica**: +0,14 en NPX con p ajustado 0,082. SAA1
  (P0DJI8) y CHMP4A (Q9BY43) **no son ensayos del panel**, así que este recurso
  no dice nada sobre ellas, ni a favor ni en contra. CEACAM5 sí sube (+0,52, p
  ajustado 4,7·10⁻⁵) y CEACAM1 baja, lo que concuerda con el eje D1.
- **Qué rinden tres proteínas**: el clasificador del artículo da AUC 0,854 con 3
  proteínas y 0,994 con las 1 463. Es la escala de lo que puede esperarse de un
  panel del tamaño del de v1.
- **Solo estadísticas publicadas**. No hay acceso a datos individuales, así que
  no se recalcula nada ni se aplica el `treat` de este trabajo.
- **Limitaciones**:
  - NPX es una escala log2 relativa de inmunoensayo, no comparable en magnitud
    con el log2 de una intensidad de espectrometría de masas. Solo se cruzan el
    signo y el hecho de ser significativa.
  - Los controles son sanos de 50 a 65 años, no pacientes con síntomas ni con
    enfermedad benigna. Un AUC contra sanos es optimista para cualquier uso
    clínico: en la consulta el problema es separar cáncer de enfermedad benigna.
    Eso es justo lo que mira D6.
  - Los grupos están desbalanceados (221 contra 74) y vienen de reclutamientos
    distintos, así que parte de la separación puede ser de cohorte.
  - Que una proteína no esté en el panel no es evidencia de que no cambie. Dos
    de las tres proteínas de v1 están en esa situación.

### Wannberg et al. 2026 — cáncer frente a enfermedad benigna

- **Estado**: usada
- **Script**: `scripts/download/wannberg_2026.py`
- **Fichero**: `external/wannberg_2026/wannberg_2026_tabla_s5.tsv`
- **Origen**: información suplementaria 1 del artículo (PDF de 9 páginas).
- **Licencia**: CC BY 4.0, acceso abierto.
- **Cita**: Wannberg F, Bueno Álvez M, Qvick A, et al. Plasma protein profiling
  predicts cancer in patients with non-specific symptoms. *Nature
  Communications* 2026. doi:10.1038/s41467-025-67688-3
- **Contenido**: Olink Explore con las mismas 1 463 proteínas que D4b. Cohorte
  de descubrimiento MEDECA: 456 pacientes derivados por síntomas no
  específicos, 160 con cáncer y 296 sin él. Réplica ALLVOS: 238 pacientes, 35
  con cáncer. AUC 0,80 en descubrimiento y 0,82 en réplica.
- **Qué aporta**: el grupo de control que a v1 y a D4b les falta. Aquí los
  controles no son personas sanas: son pacientes con síntomas que, tras el
  estudio, no tenían cáncer. De los 155 descritos en la cohorte de
  descubrimiento, 55 tenían una enfermedad autoinmune y 32 una infección.
- **Resultado en v2, la comparación que importa**: con el mismo panel de 1 463
  proteínas, el AUC es **0,994 contra personas sanas** (D4b) y **0,80 contra
  pacientes con síntomas**. La diferencia no viene del método ni de la
  plataforma: viene de contra quién se compara. Es la advertencia que hay que
  hacerle a los tres resultados de plasma de v1, obtenidos también contra
  personas sanas.
- **Segundo resultado**: las 22 proteínas de la tabla S5 **no son de respuesta de
  fase aguda**. Ocho son nucleares, cinco intracelulares, tres mitocondriales y
  dos citoqueratinas; solo una, S100A12, es claramente inflamatoria, y es la de
  menor efecto de las 22 (OR 1,27 en descubrimiento). En plasma, las proteínas
  nucleares y las citoqueratinas indican recambio y muerte celular, no
  inflamación. Ninguna de las tres proteínas de plasma de v1 aparece.
- **Columnas usadas**: de la tabla S5, `Cohort`, `Assay`, `OR`, `SE`,
  `Statistic`, `P value`, `Adjusted P value`, `CI low`, `CI high`. De las tablas
  S1 a S4, los tamaños de los grupos y los diagnósticos no malignos. Del
  resumen del artículo, los dos AUC.
- **Números extraídos de un PDF**: las 44 filas de la tabla S5 se extraen del
  texto del PDF con un patrón de fila completa. Una fila que no encaje entera no
  entra. Se comprueba además que el OR caiga dentro de su intervalo, que el p
  ajustado no sea menor que el crudo y que las dos cohortes traigan las mismas
  proteínas. Cada fila lleva `fuente_tabla`, `fuente_pagina` y
  **`verificado_por` vacía**, que rellena el autor tras comparar con el
  original.
- **Limitaciones**:
  - Es un trabajo de cáncer en general. Solo 5 de los 160 cánceres de
    descubrimiento y 1 de los 35 de la réplica eran colorrectales, así que de
    aquí no sale ninguna afirmación sobre cáncer colorrectal en particular.
  - La localización de cada proteína es la anotación de UniProt, no un resultado
    del artículo. Sirve solo para decir si la lista se parece o no a una
    respuesta de fase aguda.
  - Las dos cohortes son suecas y de dos hospitales. Un grupo de control de
    pacientes con síntomas no es el mismo en todos los sistemas sanitarios.
  - Poner un AUC de este trabajo al lado de uno de D4b no es un contraste
    estadístico: son dos números publicados puestos uno al lado del otro. Vale
    porque la plataforma y el número de ensayos coinciden, y nada más.

---

## Recordatorio

Ningún resultado de esta tesis sirve para decidir sobre un paciente.
