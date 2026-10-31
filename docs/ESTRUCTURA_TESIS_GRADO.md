# Estructura de la tesis de grado — OncoGO-CRC

Programa: Bioinformática y Biología Computacional, Universidad de Caldas.
Proyecto: OncoGO-CRC. Fecha de esta estructura: 31 de octubre de 2026.

Este documento dice cómo armar la tesis con lo que el repositorio ya calcula y con el cierre que se puede ejecutar. Una referencia solo se cita con la frase que su texto sostiene. La columna de verificación distingue lo reabierto el 31 de octubre de 2026 de lo que viene de la memoria de revisión del 30 de octubre y todavía no se ha vuelto a abrir.

Regla de resultados: en toda tabla de proteínas, la columna `Origen` vale `literatura`, `recurso externo` o `resultado propio`. Una proteína de otro artículo no se presenta como hallazgo de esta tesis.

## Cómo correr el cierre

Con el disco de datos montado:

```bash
make tesis-cierre
```

Equivale a `oncogo-crc run --data-workspace /Volumes/JANCODATA/OncoGO-CRC --stage tesis_cierre`. Si el volumen no está montado, la etapa escribe en `results/tesis_cierre/` del repositorio y deja `estado_suero` en `matriz_ausente`. No reescribe el certificado de la versión 1 ni las tablas de la versión 2.

Salidas:

| Archivo | Qué es |
|---|---|
| `candidatas_literatura.tsv` | Lista fijada antes del contraste. Todas con `origen = literatura`. |
| `auditoria_cptac.tsv` | Esa lista cruzada con el contraste certificado de tejido, si la tabla existe. |
| `crc_vs_sano.tsv`, `crc_vs_enfermedad.tsv`, `enfermedad_vs_sano.tsv` | Contrastes del suero, solo si los nombres de corrida traen `CP_`, `CN_` y `Normal_`. |
| `especificidad_candidatas.tsv` | Veredicto por candidata con la regla de más abajo. |
| `peptidos_neoantigenos.tsv` | Qué de los péptidos de 2025-2026 se puede afirmar y qué no. |
| `cierre_audit.json` | Conteos, motivo si un contraste no corrió, y la regla usada. |

La regla de especificidad, fijada en el código antes de leer un resultado real, es esta: una candidata es `especifica` si el AUC de discriminación de cáncer frente a colonoscopia negativa supera 0,70, el límite inferior del intervalo bootstrap supera 0,5, y ese AUC supera al de enfermedad frente a sanos. Si la enfermedad separa de los sanos tanto o más que el cáncer separa de la enfermedad, el veredicto es `compartida_con_enfermedad`. Con menos de 30 muestras por lado el audit marca el alcance como `cribado`.

PXD036594 no entra en ese contraste. Es un panel dirigido y el nombre de la corrida no trae el diagnóstico. El archivo `metadata/sample_metadata/plasma_metadata.tsv` queda en `sin_etiqueta_clinica`. Partir esas 20 corridas por orden alfabético fabricaba un grupo de cáncer y un grupo sano.

## Qué se puede aplicar de los neoantígenos y péptidos de 2025-2026

Tres trabajos se reabrieron. Ninguno convierte esta tesis en un descubrimiento de neoantígenos.

**Rai et al., 2025** (verificado). Suero DIA, tres grupos de 20: colonoscopia positiva, colonoscopia negativa y sanos. Accesión de descubrimiento PXD055218. Validación dirigida en PXD065161, 10 por grupo. El resumen y el método PRM del suplemento dicen que el panel dirigido cubre cuatro péptidos y que uno de ellos es el neoantígeno. Reportan 19 péptidos descritos como específicos del suero de cáncer colorrectal, la proteína BAIAP3, y un panel con CEA de AUC 0,758 (sensibilidad 67 %, especificidad 64 %). Lo aplicable aquí es el diseño de tres grupos y la auditoría de BAIAP3 y CEACAM5 en la matriz de proteínas. Lo que no se puede hacer es confirmar o rechazar los 19 péptidos: una fila de grupo proteico no guarda la secuencia mutada. Rebuscar los crudos con una base que incluya las 19 secuencias queda como trabajo futuro, y solo cuando esas secuencias se copien del suplemento.

**Xiang et al., 2026** (verificado en el texto de *Cell Genomics*). Diez pares de tumor y normal. Espectrometría del inmunopeptidoma de MHC-I del tejido, no del suero. 96 neoepítopos mutados; el 80,21 % proviene de regiones no codificantes. En tumores hipermutados dominan las regiones intergénicas e intrónicas; en los no hipermutados, las variaciones de codificación y el empalme alternativo. Lo aplicable es una frontera de interpretación: un péptido presentado por MHC en el tumor no es una proteína diferencial de sangre. No se cruzan sus secuencias con CPTAC ni con PXD055218.

**Hembrough et al., 2026** (verificado como resumen de congreso AACR, abstract 7710). De vesículas de 0,5 mL de plasma recuperan péptidos de MHC-I específicos de paciente. El resumen nombra un caso de cáncer colorrectal avanzado, con 5 péptidos, y otro de osteosarcoma. Parte de esos péptidos se confirmaron en tejido y algunos tuvieron linfocitos T CD8 en sangre. Lo aplicable es el límite: la detección en plasma es posible, es de una persona, y no es un panel que esta tesis pueda replicar. No hay secuencias ni cohorte pública.

La degradómica sérica sin sesgo y la búsqueda abierta de péptidos semitripticos no se ejecutan. Piden los crudos, una base de señuelo y un control de falsos positivos por clase de péptido. En la tesis van al capítulo de trabajo futuro, con esa razón, no como un resultado pendiente de redactar.

## Estructura de capítulos

Ajustar títulos y extensión a la plantilla de la Universidad de Caldas. El orden de abajo sigue lo que el repositorio puede sostener.

### Resumen

Un párrafo con el problema (la sangre no repite el proteoma del tumor), el método (cohortes separadas, tejido CPTAC, réplica tisular PXD034575, suero de tres grupos cuando la matriz lo permite) y el aporte (concordancia, especificidad frente a otra enfermedad intestinal, auditoría de candidatas publicadas, y el límite explícito de los neoantígenos). 250 a 300 palabras. Ningún número que no esté en `cierre_audit.json` o en el informe de la versión 2.

### 1. Introducción

Cáncer colorrectal y la necesidad de una prueba de sangre. Por qué la proteómica de espectrometría de masas y por qué Gene Ontology como lectura de proceso, no como lista de proteínas. Citar Vasaikar 2019 para el proteoma de tejido que usa la tesis, y Rai 2025 para el problema clínico de separar cáncer de otra enfermedad intestinal. No abrir con la promesa de un panel diagnóstico.

### 2. Planteamiento del problema

Tres hechos, en este orden. El cambio en el tejido no se transfiere al suero: la versión 2 ya midió una concordancia de signo compatible con el azar entre CPTAC y el suero. Las candidatas publicadas se solapan con inflamación, hemólisis y proteínas del huésped. Los neoantígenos séricos publicados no están replicados a nivel de péptido en una cohorte independiente que esta tesis pueda rebuscar.

La correlación de 0,104 entre CPTAC y PXD034575 es tejido contra tejido (TMT frente a LFQ). No se usa como prueba de que la sangre refleja mal al tumor. Esa frase corresponde al contraste tejido→suero.

### 3. Justificación

El hueco que esta tesis cierra es operativo: medir, por cohorte y con una regla escrita antes, cuánto de la señal de tejido se ve en suero, si esa señal separa el cáncer de otra enfermedad intestinal, y cuántas candidatas publicadas sobreviven a ese control. Integrar tejido y sangre de cáncer colorrectal ya está publicado. No se reclama esa prioridad.

### 4. Objetivos

Objetivo general. Cuantificar la concordancia entre el proteoma del tumor colorrectal y el suero, y evaluar la especificidad y la replicación de candidatas publicadas, sin fusionar matrices de plataformas distintas.

Objetivos específicos.

1. Medir la concordancia de signo y, cuando el enriquecimiento esté calculado en las dos matrices, la concordancia de términos GO entre tejido y suero.
2. Evaluar si las candidatas publicadas separan el cáncer de la colonoscopia negativa además de separarlo de los sanos.
3. Auditar las candidatas de la literatura en CPTAC y en el suero, con origen declarado.
4. Delimitar qué parte de los neoantígenos de 2025-2026 es evaluable con las matrices de esta tesis.

### 5. Marco teórico

Espectrometría de masas: DDA, DIA, TMT y PRM, con una frase de para qué sirve cada una en los datasets de la tesis. Proteoma de plasma: la mayor parte de la señal circulante no sale del tumor. Gene Ontology: ORA sobre una lista con umbral, GSEA sobre el ranking completo; un ORA de cero términos no es ausencia de biología. Efecto de lote: no se combinan intensidades de estudios distintos. Aprendizaje automático: el bosque aleatorio de la versión 2 se valida en otra cohorte y por paciente; un AUC interno no se presenta como rendimiento diagnóstico.

### 6. Estado del arte

Una subsección por trabajo verificado, más una tabla de la memoria de revisión marcada como no reabierta. La matriz de ejes del 30 de octubre puede ir al anexo. En el cuerpo solo entra lo que cambia la interpretación de un resultado propio.

### 7. Metodología

Datasets, cada uno por separado. Tabla con matriz, adquisición, diseño y accesión (capítulo de referencias de datos, más abajo). Prerregistro de la versión 2 en `config/preregistro_v2.yaml`. Reglas nuevas de este cierre en `src/oncogo_crc/tesis/cierre.py`, con la fecha de la etapa. Preprocesamiento sin imputar ausencias. Contraste con t moderada y umbral de magnitud `treat`. Especificidad por AUC y bootstrap estratificado, semilla 20261031. Catálogo de candidatas leído de `metadata/tesis/candidatas_literatura.tsv`.

### 8. Resultados

8.1 Concordancia tejido–tejido y tejido–suero, con los números del informe de la versión 2 y, si el cierre corrió, los tres contrastes de suero.

8.2 Especificidad de las candidatas. Si `estado_suero` no es `tres_grupos`, este apartado dice que el contraste no se calculó y por qué. No se rellena con el corte alfabético de PXD036594.

8.3 Auditoría. Tabla con gen, origen, cohorte, log2FC, FDR y veredicto.

8.4 Péptidos y neoantígenos. La tabla `peptidos_neoantigenos.tsv` es el resultado: el límite del método, no una lista de péptidos descubiertos.

8.5 Modelo de la versión 2, como priorización dentro de tejido, con el veredicto de que los intervalos se solapan y no hay ganador.

### 9. Discusión

Comparar el resultado de suero con Rai 2025 (mismo diseño, otra pregunta: aquí se audita y se mide concordancia con el tejido). Comparar el inmunopeptidoma de Xiang 2026 con la imposibilidad de ver esos péptidos en una matriz de proteínas. Comparar el resumen de Hembrough 2026 con el tamaño de cohorte de esta tesis. Discutir hemólisis y fase aguda cuando una candidata caiga en `lista_ruido`.

### 10. Limitaciones

Cohortes de suero de decenas de personas. Plataformas distintas. PXD036594 sin etiqueta clínica utilizable. El grupo de colonoscopia negativa de Rai incluye otras enfermedades intestinales, no solo enfermedad inflamatoria intestinal. No hay tejido y sangre del mismo paciente en las matrices usadas. No hay población latinoamericana en estos datos públicos. Los neoantígenos no se rebuscaron en los crudos. Varias fichas de la memoria del 30 de octubre no se reabrieron para esta estructura.

### 11. Conclusiones

Una frase por objetivo. Si un contraste no corrió, la conclusión de ese objetivo es que el dato disponible no alcanza, y se cita el campo `motivo` del audit.

### 12. Trabajo futuro

Copiar las secuencias de los 19 péptidos del suplemento de Rai y rebuscarlas en un suero independiente, con q de precursor, fragmentos coeluyentes y descarte de variantes comunes. Contraste con suero de enfermedad inflamatoria intestinal de diagnóstico reciente. Concordancia en el mismo paciente cuando se incorpore un dataset pareado de tejido y vesículas. Nada de esto se anticipa con un número.

### Anexos

Catálogo de candidatas. Tabla de péptidos. Prerregistro. Fichas de la memoria del 30 de octubre, conservadas como revisión previa y no como bibliografía certificada.

## Referencias verificadas el 31 de octubre de 2026

Se abrió el DOI o el texto. La anotación de abajo es la que puede entrar en la bibliografía. El uso dice en qué capítulo va y qué frase autoriza.

### Rai AB, Codi JAK, Kotimoole CN, Ramesh P, Narayana VK, Shaji V, Prasad TSK

Data-Independent Acquisition Analysis to Identify Serum Biomarkers and Neoantigens for Differentiating Colorectal Cancer from Inflammatory Bowel Diseases. *Journal of Proteome Research*. 2025;24(8):4166-4180. doi:10.1021/acs.jproteome.5c00285. PMID: 40652393.

Datos: PXD055218 (DIA, descubrimiento, n = 20 por grupo) y PXD065161 (PRM, n = 10 por grupo).

Qué se extrajo. Tres grupos clínicos en el nombre de la corrida. Tres proteínas alteradas en colonoscopia positiva frente a negativa. BAIAP3 sobreexpresada en cáncer. 19 péptidos descritos como neoantígenos de suero. Validación dirigida de cuatro péptidos, uno de ellos neoantígeno. Panel con CEA: AUC 0,758, sensibilidad 67 %, especificidad 64 %, exactitud 65 %. Los autores piden una cohorte más grande antes de hablar de utilidad clínica.

Cómo usarla. Metodología, para el diseño de tres grupos. Resultados, para auditar BAIAP3 y CEACAM5 como literatura. Discusión, para no presentar los 19 péptidos como replicados. No usar el AUC 0,758 como listón que esta tesis deba superar: es el rendimiento de un panel elegido en el mismo estudio.

Estado: verificado el 31 de octubre de 2026 en el DOI, el resumen de ACS y la ficha PRIDE.

### Xiang H, Guan X, Wei Y, et al.

Predominant mutated non-canonical tumor-specific antigens identified by proteogenomics demonstrate immunogenicity and tumor suppression in CRC. *Cell Genomics*. 2026;6:101062. doi:10.1016/j.xgen.2025.101062. Publicado en línea el 13 de noviembre de 2025; número del 14 de enero de 2026.

Qué se extrajo. Inmunopeptidoma de MHC-I en tejido tumoral pareado, 10 pacientes. 96 neoepítopos mutados. 80,21 % de regiones no codificantes. La composición cambia con la carga mutacional. Hay validación funcional en modelo murino, que esta tesis no repite.

Cómo usarla. Marco teórico y discusión, una sola idea: la presentación por MHC en el tumor es otra medición que la abundancia de la proteína en suero. No citar los 96 péptidos como candidatos de sangre.

Estado: verificado el 31 de octubre de 2026 en el PDF del artículo (encabezado, highlights y métodos de muestra).

### Hembrough T, Ezrin AM, Opheim Z, Bandoski C, Bennett S, Fagan-Solis K, Vincent BG, Hoos W

Extracellular vesicle immunopeptidomics recovers immunogenic neoantigen peptides from cancer patient plasma. *Cancer Research*. 2026;86(7 Suppl):Abstract 7710. doi:10.1158/1538-7445.am2026-7710.

Qué se extrajo. Resumen de la reunión anual de la AACR de 2026. Vesículas de 0,5 mL de plasma, enriquecimiento de MHC-I, espectrometría dirigida con péptidos pesados. Un paciente con cáncer colorrectal avanzado, 5 péptidos; un paciente con osteosarcoma, 3. Confirmación parcial en tejido y linfocitos T CD8 para varios péptidos.

Cómo usarla. Trabajo futuro y limitaciones. Es un resumen, de dos pacientes, con un método propietario. No es una cohorte ni una replicación de Rai 2025.

Estado: verificado el 31 de octubre de 2026 en la página del DOI. No se leyó un artículo completo porque el trabajo es un abstract.

### Vasaikar S, Huang C, Wang X, et al.

Proteogenomic Analysis of Human Colon Cancer Reveals New Therapeutic Opportunities. *Cell*. 2019;177(4):1035-1049.e19. doi:10.1016/j.cell.2019.03.030.

Qué se extrajo para esta tesis. Es la fuente del proteoma de tejido PDC000116 que el pipeline usa como ancla tumoral (tumor y normal adyacente, TMT). El artículo no se releyó entero el 31 de octubre de 2026. La accesión y el papel de la matriz sí están en `config/datasets.yaml` y en la corrida certificada.

Cómo usarla. Metodología, como dataset de tejido. No atribuirle los resultados de suero.

Estado: cita estándar del recurso que el repositorio analiza. El DOI es el del artículo de CPTAC colon de 2019. Conviene confirmar páginas contra el PDF de *Cell* antes de la entrega impresa.

### Robles J, Prakash A, Vizcaíno JA, Casal JI

Integrated meta-analysis of colorectal cancer public proteomic datasets for biomarker discovery and validation. *PLOS Computational Biology*. 2024;20(1):e1011828. doi:10.1371/journal.pcbi.1011828. PMC10833860.

Qué se extrajo. Ya existe un meta-análisis de datasets públicos de cáncer colorrectal que incluye tejido y biopsia líquida. La tesis no se presenta como el primer trabajo de ese tipo. Candidatas que ese artículo destaca, entre ellas CD14 y TXNDC5, son literatura.

Cómo usarla. Estado del arte y justificación, para acotar la novedad. Si una de sus proteínas entra en una tabla, `origen = literatura`.

Estado: el identificador PMC está en el repositorio y en la memoria del 30 de octubre. El texto completo no se reabrió el 31 de octubre de 2026. Comprobar autores y paginación en el PDF antes de cerrar la bibliografía.

### Smyth GK

Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. *Statistical Applications in Genetics and Molecular Biology*. 2004;3:Article 3.

McCarthy DJ, Smyth GK. Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics*. 2009;25(6):765-771.

Qué se extrajo. Son la referencia del contraste implementado en `src/oncogo_crc/utils/moderated.py`: varianza moderada y prueba contra un cambio mínimo, no contra cero. La tesis los cita en metodología como el método estadístico, no como un resultado biológico.

Estado: implementación propia en Python, declarada en `docs/enmiendas_v2.md`. No se ejecutó limma en R.

## Datasets que la tesis puede nombrar porque el repositorio los usa

| Accesión | Matriz | Papel en la tesis | Artículo |
|---|---|---|---|
| PDC000116 | Tejido tumoral y normal, TMT | Ancla tumoral | Vasaikar et al., *Cell*, 2019 |
| PXD034575 | Tejido LFQ, 16 tripletes | Primario, metástasis hepática y mucosa normal. No es pólipo ni sangre | Tanaka et al., *Cell Reports*, 2024. doi:10.1016/j.celrep.2024.113810 |
| PXD055218 | Suero, DIA | Tres grupos, si la matriz está descargada | Rai et al., 2025 |
| PXD065161 | Suero, PRM | Validación dirigida de los autores, no reanalizada aquí | Rai et al., 2025 |
| PXD036594 | Plasma, panel dirigido PQ500 | No se contrasta. Sin etiqueta clínica en el nombre de corrida | Citar solo para explicar la exclusión |

PXD015163 figura como accesión en la configuración, pero la matriz que el archivo de configuración descarga es el informe DIA-NN de PXD055218. En la tesis se nombra el archivo que se analizó.

## Memoria de revisión del 30 de octubre de 2026

La memoria larga aportada por el autor reúne fichas de Yu, Ögren, Seo, Dhami, Franzetti, Malmström, Ramsbottom, Yarmolinsky, Lee, Zhang y Tan, todas con año 2026 en esa memoria, más trabajos complementarios. Esas fichas sirven para el anexo y para no reclamar como propio un resultado que la memoria ya atribuye a otros.

No entran en la bibliografía certificada de esta estructura. No se reabrió el PDF el 31 de octubre de 2026, y varios identificadores de esa memoria quedan por comprobar en Europe PMC antes de copiarlos al documento de grado. Hasta entonces, en el texto se dice «según la revisión compilada el 30 de octubre de 2026» y no se da el PMID como verificado.

Lo que sí se aprovecha de esa memoria, sin convertirla en cita certificada:

- No presentar como novedad la intersección tejido–sangre, el origen del huésped, el solapamiento con enfermedad inflamatoria intestinal, ni el fallo de ComBat entre estudios de suero.
- Marcar hemoglobinas, CA1, CAT y ALAD como posible hemólisis, y TFF3, TFF1, AHCY, RETN, LCN2 y SELE como compartidas con inflamación intestinal, en la columna `lista_ruido`.
- Dejar el atlas de Malmström como capa de origen tisular para cuando el suplemento esté en el disco. No está aplicado en esta ejecución.
- Tratar el aporte de péptidos no canónicos y la degradómica como trabajo futuro.

## Corrida del 31 de octubre de 2026

La etapa `tesis_cierre` se ejecutó sobre el disco de datos. El suero de PXD055218 trae los tres grupos, 20 corridas cada uno. El umbral de magnitud es el del certificado (log2FC 0,585) y el FDR es 0,05. Alcance declarado: cribado.

Contraste de proteínas con t moderada, unas 700 proteínas con prueba en cada comparación:

| Contraste | Significativas |
|---|---|
| Cáncer frente a sanos | 3: SAA1 (P0DJI8, log2FC 3,26, 40 observaciones), LBP (P18428, log2FC 1,25, 40 observaciones), CHMP4A (Q9BY43, log2FC 5,42, 16 observaciones) |
| Cáncer frente a colonoscopia negativa | 0 |
| Colonoscopia negativa frente a sanos | 1: APOA4 (P06727, log2FC −1,32) |

SAA1 y LBP son reactantes de fase aguda. CHMP4A pasa el umbral con la mitad de las observaciones; no se presenta como marcador hasta revisar el patrón de ausencias. Ninguna candidata de la literatura pasa el FDR en cáncer frente a enfermedad, incluido BAIAP3 (log2FC 2,20, FDR 0,20).

La regla de AUC, más laxa que el FDR, marca 8 candidatas como `especifica`. Solo BAIAP3, TPM4 y FN1 no están en la lista de fase aguda. BAIAP3: AUC 0,817 frente a enfermedad (límite inferior 0,673) y AUC 0,675 de enfermedad frente a sanos. C9, LRG1, SERPINA1, CRP y HP también superan la regla y llevan `lista_ruido = fase_aguda`: el veredicto estadístico no las convierte en específicas de cáncer. CRP no se cita como hallazgo propio.

En tejido CPTAC, 33 de las candidatas se detectan y 7 replican el signo de la literatura con FDR: CEACAM5, CA1, SMPDL3A, HBB, HBD, LCN2 y TOP2A. CA1, HBB y HBD quedan marcadas como hemólisis. LCN2 queda marcada como compartida con enfermedad inflamatoria intestinal. TOP2A sube en tejido y no tiene contraste evaluable en el suero.

Los 19 péptidos de Rai 2025 y los neoepítopos de Xiang 2026 siguen en `peptidos_neoantigenos.tsv` como no evaluables en esta matriz.

## Qué número puede presentarse

Un número de suero de esta etapa solo si `cierre_audit.json` dice `estado_suero = tres_grupos` y el contraste correspondiente dice `calculado`. Un número de la versión 2 solo si sale de la tabla citada en el informe, que ya se genera desde esas tablas. El corte alfabético antiguo de PXD036594 no se presenta. Los 19 neoantígenos no se presentan como resultado propio ni como replicados.
