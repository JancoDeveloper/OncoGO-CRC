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

### PXD034575 — Tanaka et al., Cell Reports 2024

- **Estado**: usada
- **Artículo**: Tanaka A. et al. *Cell Reports* 2024; 43:113810.
  doi:10.1016/j.celrep.2024.113810
- **Qué es**: LFQ con MaxQuant de 16 tripletes. El artículo describe tumor
  primario (pCRC), metástasis hepática (mCRC) y mucosa colónica benigna.
  En las columnas, **P es el primario, M la metástasis y N la mucosa**.
  No es un pólipo. La versión 2 guardó P como `polyp` y M como `tumor`;
  esa etiqueta no se usa en el cierre de tesis.
- **Qué aporta a v2**: la réplica de dirección se calculó con M frente a N.
  Con la leyenda del artículo, eso es metástasis hepática frente a mucosa, no
  tumor primario frente a mucosa. El brazo que v2 llamó pólipo es el primario.
- **Nombres de columna**: `LFQ intensity CRC<n><M|N|P>`. **M es metástasis, N
  mucosa normal y P tumor primario.** No es `T` y no es un pólipo.
- **Trampa conocida**: el cero de MaxQuant significa ausencia de identificación,
  no abundancia cero. No se imputa.

### PXD055218 — suero DIA, tres grupos

- **Estado**: usada
- **Qué es**: suero DIA de Rai et al. 2025, 20 cáncer, 20 colonoscopia negativa
  y 20 sanos. El contraste binario histórico solo lee cáncer y sanos. El
  cierre de tesis usa los tres.
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

### PDC000111 — proteoma retrospectivo de colon de TCGA (Zhang et al. 2014)

- **Estado**: usada
- **Script**: `scripts/download/tcga_cptac_retro.py`
- **Etapa**: `make v2-pureza`
- **Portal**: <https://pdc.cancer.gov/pdc/study/PDC000111>
- **Identificador**: `PDC000111`, `TCGA COAD Proteome S016-1`, CPTAC2
  Retrospective, label free.
- **Cita**: Zhang B, Wang J, Wang X, et al. Proteogenomic characterization of
  human colon and rectal cancer. *Nature* 2014; 513(7518):382-387.
  doi:10.1038/nature13438
- **Licencia**: datos abiertos del Proteomic Data Commons del NCI. Las matrices
  derivadas del CDAP son públicas y no requieren acuerdo de uso.
- **Qué contiene**: 95 alícuotas de 90 pacientes, 5 561 genes. 64 de colon y 31
  de recto.
- **Columnas usadas**: `Gene` y las que acaban en `Unshared Area`, que son los
  péptidos que solo pueden venir de ese gen. Es el mismo criterio de «no
  compartido» que v1 usó en CPTAC.
- **Qué se descarga**: solo la matriz derivada del CDAP, 8,7 MB. Los espectros
  crudos del estudio pasan del terabyte y no hacen falta para nada de esto.
- **Para qué sirve**: es la única cohorte con proteoma y **pureza del tumor
  medida por ADN en la misma persona**, porque son pacientes de TCGA. Eso permite
  preguntar si el índice estromal que R4 usa como covariable se parece a una
  pureza de verdad. El preregistro ya contemplaba esta vía
  (`pureza_absolute_tcga` entre los métodos de sensibilidad de C7).
- **Tres discrepancias con el enunciado, encontradas al consultar la API en vez
  de suponer**:
  1. **El identificador no era adivinable.** El candidato que más se parece por
     el nombre, `PDC000289` («TCGA Colorectal Cancer CompRef Samples»), son
     muestras de referencia para control de calidad y no pacientes. Usarlo habría
     dado una cohorte de duplicados técnicos presentada como cohorte de personas.
  2. **Esta cohorte no tiene normales**: las 95 alícuotas son tumor primario. Sin
     pares no hay contraste tumor-normal, así que D2 **no recalcula R4**. Valida
     su covariable, que es una pregunta distinta y se declara como tal.
  3. **`PDC000109` no es una cohorte independiente.** Es «Prospective COAD
     Proteome S037-1», label free, y 97 de sus 100 pacientes son los mismos de
     `PDC000116`, el descubrimiento del certificado. Es la misma gente
     recuantificada por otro laboratorio. Si entrara como réplica, contaría a
     esas personas dos veces y el acuerdo entre cohortes saldría altísimo por
     construcción. El script lo comprueba cruzando identificadores de caso y deja
     el veredicto escrito.
- **El cero es ausencia**: el 36 % de la matriz son áreas de cero, que en label
  free significa que el péptido no se detectó. Se tratan como ausentes y no se
  imputa nada.
- **Las alícuotas repetidas se usan**: cinco tumores se midieron dos veces. Esas
  réplicas dan la reproducibilidad del ensayo, que es lo que decide si una
  correlación débil es biología o es el método.
- **Limitaciones**:
  - Es label free y `PDC000116` es TMT. Lo que se valida es el constructo del
    índice, no los valores concretos de R4.
  - Son tumores sin normales, así que el índice se calcula entre tumores y no
    como diferencia dentro de un par.
  - La pureza ABSOLUTE es a su vez una estimación de ADN, no un recuento de
    células al microscopio: la comparación es entre dos medidas indirectas.
  - La pieza del proteoma y la del ADN no son el mismo trozo de tumor, así que
    parte de la discrepancia es muestreo y no error de ninguna de las dos.
  - El CEA sérico solo existe en una parte de estos pacientes, así que el
    análisis del CEACAM5 corre sobre menos gente que el de la pureza. Los dos
    tamaños se publican por separado.
- **No entra en el certificado**: ningún número de esta fuente es un subíndice.
  Lo que decide es cómo hay que describir R4, y queda recogido en la enmienda E8.

### PXD036434 — ocho regiones del mismo tumor

- **Estado**: usada, **fuera del meta-análisis**
- **Script**: `scripts/download/pride_heterogeneidad.py`
- **Etapa**: `make v2-heterogeneidad`
- **Portal**: <https://www.ebi.ac.uk/pride/archive/projects/PXD036434>
- **Cita**: Herzog R, Amaral FMR, Lerner R, et al. Panomics reveals patient
  individuality as the major driver of colorectal cancer progression. *Genome
  Medicine* 2023. doi:10.1186/s13073-023-01176-5
- **Licencia**: PRIDE Archive, datos públicos.
- **Qué contiene**: 2 880 genes y 19 muestras de 4 pacientes. De cada paciente,
  tumor primario, mucosa normal y metástasis hepática; y del cuarto paciente,
  **ocho regiones distintas del mismo tumor**.
- **Columnas usadas**: las que siguen el patrón `T`, `NM` y `LM` con el número de
  paciente. El identificador se toma del símbolo de gen que viene en `GN=` dentro
  de `DESCRIPTION`, **no** de la columna `PROTEIN ID`, que trae nombres de entrada
  de UniProt (`PLEC_HUMAN`). Cruzar un nombre de entrada con un símbolo de gen no
  falla: devuelve cero coincidencias y parece que las cohortes no comparten
  proteínas.
- **No entra en el meta-análisis**: el preregistro exige al menos 10 pacientes por
  cohorte (`d3_inclusion_cohortes.min_pacientes`). Son 4. El umbral no se dobla
  para que entre.
- **Para qué sirve entonces**: ninguna cohorte del certificado tiene más de un
  trozo del mismo tumor, así que ninguna puede comprobar lo que todas suponen:
  que el trozo medido representa al tumor. Esta sí.
- **Dos trampas del análisis, las dos resueltas**:
  1. **El número de observaciones.** Comparar la desviación de ocho valores con
     la de cuatro no es justo, porque el estimador de cuatro es más ruidoso. Se
     eligen cuatro de las ocho regiones al azar, 500 veces.
  2. **A qué paciente se compara.** El paciente de las ocho regiones también es
     uno de los cuatro, y cuál de sus regiones lo represente cambia la variación
     entre pacientes. Fijar una sería escoger el resultado, así que se promedia
     sobre todas las elecciones.
- **La concordancia se publica junto a lo que daría el azar**: con cuatro pares,
  el azar ya pone de acuerdo a una proteína de cada ocho. Un porcentaje de
  concordancia sin esa referencia no se puede leer.
- **Limitaciones**:
  - Son cuatro pacientes. De aquí no sale ninguna estimación de efecto comparable
    con las cohortes del certificado.
  - Las ocho regiones son de **un** tumor de **un** paciente. Lo que se puede
    decir es que existe al menos un tumor así, no cuántos tumores son así.
  - La cuantificación es un cociente contra una referencia común y la plataforma
    no es la de ninguna cohorte del certificado. Lo que se compara aquí es
    variación contra variación dentro de la misma matriz.
  - Las regiones del mismo tumor se midieron probablemente más juntas en el
    tiempo que las muestras de pacientes distintos. Si eso sesga algo, las haría
    parecer **más** parecidas entre sí: el sesgo va en contra del resultado.

### La búsqueda de cohortes de D3, con los descartes

El preregistro fija cuatro criterios en `d3_inclusion_cohortes`: tejido humano,
tumor y normal, matriz descargable y al menos 10 pacientes. Se revisaron los 571
depósitos de PRIDE que responden a «colorectal» y los estudios colorrectales del
PDC. **Ningún depósito nuevo cumple los cuatro.** El registro queda en
`results/v2/axes/heterogeneidad/busqueda_de_cohortes_d3.tsv` porque «no
encontramos más» y «buscamos y esto es lo que hay» no son lo mismo.

| Fuente | Veredicto | Motivo |
| --- | --- | --- |
| PXD019504 | incluida | cumple los cuatro criterios; es la que C10 ya usa |
| PDC000278 | no disponible | 30 personas con epitelio de colon normal, pero el depósito no publica matriz de proteínas: solo espectros (393 GB) y PSM |
| PDC000109 | excluida | 97 de sus 100 pacientes son los de PDC000116 |
| PXD002137 | no disponible | el diseño es el ideal (mucosa, adenoma, cáncer) pero no hay matriz publicada |
| PXD001676 | no disponible | sin matriz de proteínas publicada |
| PXD060933 | excluida | tiene matriz, pero sus 12 muestras son todas de cáncer y no hay normal |
| PXD059755 | excluida | 594 proteínas y el campo de enfermedad viene como «not available» |
| PXD051882 | excluida | el diseño da el canal de TMT pero no dice cuál es tumor y cuál mucosa |
| PXD036434 | fuera del meta-análisis | 4 pacientes frente a los 10 que exige el preregistro |

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

### UK Biobank Pharma Proteomics Project — pQTL del plasma

- **Estado**: **requiere_acceso**
- **Script**: `scripts/download/ukb_ppp_pqtl.py`, preparado y con tests
- **Origen**: las estadísticas de resumen están en Synapse, proyecto
  `syn51364943` (doi:10.7303/syn51364943), bajo los términos del UKB-PPP y de
  UK Biobank. El portal público es `ukb-ppp.gwas.eu`.
- **Por qué no está**: hace falta cuenta de Synapse y aceptar las condiciones.
  El script **no intenta autenticarse**: comprueba qué páginas son públicas,
  deja el análisis escrito y marca la fuente.
- **Cita**: Sun BB, Chiou J, Traylor M, et al. Plasma proteomic associations
  with genetics and health in the UK Biobank. *Nature* 2023; 622:329-338.
- **Qué aportaría**: una tercera explicación de la brecha tejido→plasma de C8,
  junto a la ventana de detección (C8) y la elección del grupo de control (D6).
  Si el nivel de una proteína en plasma depende en buena parte de la variante
  que cada persona hereda, las diferencias entre dos individuos no hablan de la
  enfermedad: hablan del genotipo. SAA1 y LBP, las dos proteínas de fase aguda
  de v1, son de las que se sabe que tienen variantes de efecto grande.
- **Qué calcularía**: por cada proteína, su cis-pQTL de mayor varianza explicada
  —no el de menor p, que no es lo mismo: en 54 000 personas un efecto pequeño da
  un p minúsculo y explica casi nada— y si pasa del 5 % de la varianza del
  nivel. Las proteínas que se consultarían están fijadas en el script antes de
  tener los datos: SAA1, LBP, CHMP4A, CEACAM5, S100A12, KRT19 y PARP1.
- **Cómo completarlo**: quien tenga cuenta descarga el fichero y ejecuta
  `scripts/download/ukb_ppp_pqtl.py --estadisticas RUTA`. No hay nada más que
  programar: el análisis tiene 16 tests sobre estadísticas de resumen
  sintéticas.
- **Lo que no se hace**: no se sustituye por otra fuente ni por datos
  simulados. Los datos sintéticos de este análisis existen solo en `tests/`.
- **Limitación del propio recurso**: un cis-pQTL no invalida un marcador. Dice
  que un umbral fijo no vale igual para todo el mundo, que es otra cosa.

### Literatura: tablas 2×2 de CEA, metilación de SEPT9 y ctDNA

- **Estado**: usada, **pendiente de verificación humana**
- **Script**: `scripts/literatura/extraer_2x2.py`
- **Fichero**: `data/literatura/extraccion_2x2.tsv`, en el repositorio y no en el
  workspace de datos, porque es un artefacto curado a mano que necesita revisión
  e historia de cambios.
- **Qué es**: 33 tablas 2×2 copiadas de tablas publicadas, cada fila con su
  referencia, su tabla y su página.
- **Fuentes**:
  - **CEA para detectar recaída**, umbral 5 µg/L: 8 estudios de Wang R et al.,
    *World Journal of Gastrointestinal Surgery* 2023; 15(12):2904-2914, tabla 2
    (PMC10784816). 2 663 personas, 596 con recaída.
  - **Metilación de SEPT9 para cribado o diagnóstico**: 25 estudios de Sun J et
    al., *Clinical and Translational Gastroenterology* 2017; 8(1):e216, tabla 1
    (PMC5288600). 9 927 personas, 2 975 con cáncer.
  - **ctDNA para detectar recaída**: solo la estimación combinada publicada, de
    *eBioMedicine* 2024, tabla 2, fila CRC (PMC11021841). Esa fuente no publica
    las tablas por estudio en el texto principal.
- **Dos trampas de las tablas de origen, las dos resueltas en el script**:
  1. **El orden de las columnas no es el mismo.** La del CEA viene en orden TP,
     FP, FN, TN. La de la metilación viene en orden **TP, TN, FP, FN**. Leer una
     con el orden de la otra intercambia los falsos positivos con los
     verdaderos negativos, y convierte una especificidad del 99 % en una del
     1 % sin que nada falle.
  2. **Un estudio puede aparecer en varias filas.** Cinco estudios de
     metilación publican dos o tres algoritmos sobre los *mismos* pacientes.
     Meter todas las filas sería contar a esas personas varias veces. Se queda
     una fila por estudio, con la regla declarada: el algoritmo 2(3), que es el
     de la prueba comercial Epi proColon 2.0, y si no lo hay, la primera
     variante publicada. La regla va por algoritmo y no por resultado.
- **Comprobación automática**: en las dos fuentes el tamaño de muestra se
  publica aparte de las cuatro celdas, así que si las celdas suman el tamaño, el
  orden y los cuatro valores están bien. Las 33 filas suman. Además, los totales
  de la metilación reproducen exactamente los 9 927 y los 2 975 que declara el
  artículo, lo que confirma que la selección de una fila por estudio es la
  correcta.
- **Una fila excluida a propósito**: Moloney 2019. En la tabla de origen su
  prueba índice figura como «N/A» y no como CEA, y el texto dice que de los 12
  estudios nueve usaron CEA, dos ecografía y uno exploración física. Entrar en
  el modelo del CEA sería suponerlo. La exclusión queda escrita porque una
  exclusión silenciosa es indistinguible de un descuido.
- **Ninguna tabla se reconstruye desde un porcentaje.** El redondeo de una
  sensibilidad publicada no determina las celdas, y las que saldrían no serían
  las del estudio.
- **`verificado_por` vacía** en las 33 filas. La rellena el autor tras comparar
  cada fila con el original.
- **Limitación de fondo**: el CEA de estas tablas detecta recaída después de una
  cirugía; la metilación detecta cáncer en cribado. Son dos preguntas distintas,
  no se combinan nunca en un mismo modelo y la comparación directa entre las dos
  queda marcada como no comparable.

### Cohorte de CEA seriado — protocolo, sin datos

- **Estado**: `pendiente_de_aprobacion`. El estudio **no se ha hecho** y no hay
  ningún dato de ninguna persona.
- **Protocolo**: `docs/protocolo_cohorte_cea.md`
- **Análisis**: `src/oncogo_crc/v2/cohorte_cea.py`
- **Tests**: `tests/test_cohorte_cea_d9.py` (31)
- **Diana**: `make v2-cohorte-cea`
- **De dónde sale la pregunta**: del meta-análisis de la fuente anterior. El
  ctDNA medido una sola vez tras la cirugía tiene sensibilidad 0,51 y medido de
  forma seriada 0,84, con la misma prueba y la misma población. Al CEA no se le
  ha hecho esa pregunta: en los ocho estudios se evalúa como en la consulta,
  cada valor contra un umbral fijo de 5 ng/ml.
- **Qué compararía**: tres modelos sobre las mismas visitas. El umbral de
  5 ng/ml, el valor en log2 sin umbral, y la trayectoria del paciente frente a
  sus propias mediciones anteriores.
- **Datos que haría falta recoger**: ocho columnas, de las que cinco son
  obligatorias. Ni nombre, ni fecha de nacimiento, ni número de historia, ni
  ninguna fecha de calendario: todos los tiempos son meses desde la cirugía. El
  diccionario completo está en el protocolo.
- **Tamaño necesario**: 330 pacientes para detectar una diferencia de AUC de
  0,10 y 1 340 para una de 0,05, con potencia del 80 % y alfa del 5 %. Es la
  varianza de Hanley y McNeil corregida por el efecto de diseño del
  conglomerado, declarada como aproximación.
- **Cohorte sintética**: `tests/SINTETICO_cohorte_cea.csv`, generada por
  `scripts/literatura/generar_sintetico_cohorte_cea.py` con la semilla
  20260914 + 909. No corresponde a nadie. Vive solo en `tests/` y el prefijo
  `SINTETICO_` está ahí para que no se confunda nunca con datos reales. **Sus
  números no entran en ningún resultado ni en ningún informe.**
- **Qué demuestran los tests**: que el análisis detecta una diferencia cuando el
  generador la programó (el valor mejora sobre el umbral y el intervalo de la
  diferencia no incluye el cero) y que no la fabrica cuando no la programó (la
  trayectoria no mejora sobre el valor y el intervalo sí incluye el cero).
- **Limitaciones declaradas**: el estudio no diría si conviene hacer más CEA
  —detectar antes no es vivir más—, no daría un punto de corte para la consulta,
  no valdría para quien tenga una sola determinación, y no se trasladaría a otra
  población sin comprobarlo, porque el CEA basal depende del tabaquismo y de la
  función hepática.

---

### D5 — Perfiles por tipo celular del Human Protein Atlas

- **Qué es**: la expresión media de cada gen en cada tipo celular, medida con
  secuenciación de una sola célula y publicada ya agregada. No contiene ningún
  dato de ninguna persona.
- **URL**: `https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip`
- **Versión**: la vigente en la fecha de acceso; el portal documenta las
  versiones en `https://www.proteinatlas.org/about/download`.
- **Fecha de acceso**: 2026-10-03.
- **Licencia**: Creative Commons Attribution-ShareAlike 4.0 International.
  Permite reutilización citando la fuente.
- **Cita**: Karlsson M, Zhang C, Méar L, et al. A single-cell type
  transcriptomics map of human tissues. *Science Advances* 2021;
  7(31):eabh2169.
- **Columnas usadas**: `Gene name`, `Cell type`, `nCPM`.
- **Qué se construye**: los once tipos celulares del colon se fijan en
  `scripts/download/hpa_single_cell.py` antes de mirar ningún resultado, y de
  ellos se extraen 3890 genes marcadores con la regla de la enmienda E9. Cuatro
  de los once tipos son epitelio.
- **Para qué sirve**: es la referencia externa del método de sensibilidad que el
  preregistro llama `deconvolucion_referencia_externa`. Al venir de fuera del
  proteoma que se contrasta, rompe la circularidad que C7 solo pudo mitigar en el
  índice estromal.
- **Script**: `scripts/download/hpa_single_cell.py`; análisis en
  `src/oncogo_crc/v2/deconvolucion.py` y etapa `v2_deconvolucion`.

**Las tres cosas que esta fuente no puede hacer, y que conviene leer antes que
cualquier número.**

La primera: la referencia es de ARN y lo que se deconvuelve es proteína. D2 midió
el acuerdo entre el transcrito y la proteína del mismo gen en CEACAM5 y salió una
correlación de rangos de 0,24. Toda cifra que salga de aquí hereda esa
aproximación.

La segunda: los tipos celulares del atlas son de todos los tejidos juntos. Un
fibroblasto de piel y uno de colon no tienen el mismo perfil, así que
«fibroblastos» aquí es un promedio entre tejidos.

La tercera, y la que decide qué se publica: **las proporciones tipo a tipo no son
utilizables**. El residuo del ajuste es de 0,90 en una cohorte y de 0,93 en la
otra, es decir que la mezcla reproduce menos de una décima parte del perfil
medido; y cuatro o cinco tipos celulares salen ausentes en más de la mitad de las
muestras, que es lo que ocurre cuando otro tipo con perfil parecido se lleva su
señal entera. Un cero en esas tablas no significa «este tipo celular no está»:
significa que el método ha asignado su señal a otro. Lo único que se usa es el
agregado epitelio frente a no epitelio, y solo después de comprobar en los datos
que ordena los pares en el sentido que la histología espera.

### D5 — La referencia específica de colon, declarada y sin usar

- **Qué es**: `GSE178341`, el atlas unicelular de cáncer colorrectal de Pelka et
  al. Sería la referencia correcta para este tejido, porque los perfiles no
  serían un promedio entre tejidos.
- **Portal**: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178341`
- **Cita**: Pelka K, Hofree M, Chen JH, et al. Spatially organized multicellular
  immune hubs in human colorectal cancer. *Cell* 2021; 184(18):4734-4752.
- **Estado**: `no_comprobable_desde_esta_maquina`. **No** es `requiere_acceso`:
  GEO es público y el fichero no pide ninguna cuenta. Lo que falla es la red.
- **Motivo**: el servidor `ftp.ncbi.nlm.nih.gov` resuelve de forma intermitente
  en el DNS de esta red, mientras `www.ncbi.nlm.nih.gov` resuelve siempre. En un
  momento en que resolvía, el fichero respondió con código 200 y 1 203 550 558
  bytes anunciados, y la descarga se cortó a los 100 MB. El parcial se borró para
  que ninguna etapa leyera un fichero truncado que parecería válido.
- **Cómo completarlo**: descargar los tres ficheros desde una red que resuelva el
  servidor y volver a ejecutar `scripts/download/hpa_single_cell.py --pelka
  CARPETA`. El `.h5` trae la matriz de células por genes y el `cluster.csv` la
  etiqueta de tipo celular de cada célula; la firma se construye promediando por
  etiqueta.
- **Qué quedaría mejor con ella**: la segunda limitación de la fuente anterior
  —tipos celulares promediados entre tejidos— desaparecería. Las otras dos no:
  seguiría siendo ARN contra proteína, y el residuo del ajuste depende de que un
  proteoma sea una mezcla lineal de perfiles, que es lo que falla.

### D5 — Por qué la deconvolución no se corre en CPTAC

No es una decisión de conveniencia y no se afirma de memoria: la etapa lo
comprueba con los datos delante y deja el número en la auditoría.

Mezclar perfiles celulares supone que lo medido es una combinación **no
negativa** y **lineal** de ellos. CPTAC publica `Unshared Log Ratio`, que es el
logaritmo de un cociente contra un canal de referencia. **El 52,5 % de esos
valores es negativo.** Una mezcla de perfiles no negativos con pesos no negativos
solo puede producir valores no negativos, así que no existe ninguna proporción
que reproduzca ese perfil. No es que el ajuste salga malo: es que el sistema no
tiene solución con sentido.

Hay dos motivos más, independientes del primero. El logaritmo no es lineal:
aunque se recortaran los negativos, el doble de fibroblastos no da el doble de
log-intensidad. Y el cociente es contra un canal de referencia que es a su vez
una mezcla distinta en cada lote, de modo que dividir por él ya ha mezclado la
composición que se quiere estimar con la del canal.

La consecuencia es que el R4 alternativo no se calcula en la cohorte donde vive
R4. Se calcula en las cohortes LFQ, cuyas intensidades sí están en escala de
abundancia. Son cohortes distintas, con menos pares y otra tecnología: el número
no sustituye al R4 del certificado, lo acompaña.

---

## Recordatorio

Ningún resultado de esta tesis sirve para decidir sobre un paciente.
