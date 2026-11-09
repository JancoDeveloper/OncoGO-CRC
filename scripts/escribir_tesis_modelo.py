"""Borrador de tesis con la osamenta de las tesis de maestría de la Universidad de Caldas
y el contenido de OncoGO-CRC. No copia el texto científico de otros trabajos.
"""

from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SALIDA = REPO / "docs" / "borradores"

BLOQUES: list[tuple] = []


def h(nivel: int, texto: str) -> None:
    BLOQUES.append(("h", nivel, texto))


def p(texto: str) -> None:
    BLOQUES.append(("p", texto))


def viñetas(items: list[str]) -> None:
    BLOQUES.append(("ul", items))


def tabla(titulo: str, columnas: list[str], filas: list[list[str]]) -> None:
    BLOQUES.append(("tabla", titulo, columnas, filas))


def portada() -> None:
    p(
        "Nota para quien entrega el documento. Esta organización sigue el modelo de las tesis de maestría "
        "de la Universidad de Caldas: portada, página de título, epígrafe, agradecimientos, resumen, abstract, "
        "capítulos y referencias. El contenido, las cifras y la bibliografía son de OncoGO-CRC, corrida tesis_cierre "
        "del 31 de octubre de 2026. Quite esta nota, complete los corchetes y, si el programa lo pide, pase el texto "
        "a primera persona del plural. En Word, clic derecho sobre la tabla de contenido y elija Actualizar campos."
    )
    h(0, "Qué parte del tumor colorrectal se deja ver en la sangre")
    p("[Nombre completo del estudiante]")
    p("Universidad de Caldas")
    p("Facultad de Ingenierías")
    p("Manizales, Colombia")
    p("2026")
    h(0, "Página de título")
    p("Qué parte del tumor colorrectal se deja ver en la sangre")
    p("[Nombre completo del estudiante]")
    p("Tesis presentada como requisito parcial para optar al título de:")
    p("Magíster en Bioinformática y Biología Computacional")
    p("Director:")
    p("[Nombre del director, título]")
    p("Grupo de investigación:")
    p("[Nombre del grupo]")
    p("Línea de investigación:")
    p("[Nombre de la línea]")
    p("Universidad de Caldas")
    p("Facultad de Ingenierías")
    p("Manizales, Colombia")
    p("2026")
    h(1, "Epígrafe")
    p("[Escriba aquí una frase breve y cite a su autor. En las tesis de referencia esta página es personal y no forma parte del argumento científico.]")
    h(1, "Agradecimientos")
    p(
        "A la Universidad de Caldas, por la formación en el programa de Bioinformática y Biología Computacional. "
        "A [nombre del director], por la dirección de este trabajo. A [familia, profesores y compañeros, en las frases que cada quien quiera dejar]."
    )
    h(1, "Resumen")
    p(
        "El tumor colorrectal ya tiene un proteoma distinto del de la mucosa. Esta tesis pregunta cuánto de ese cambio "
        "sigue visible cuando la muestra es suero, y cuánto de lo que aparece en suero es compartido con otra enfermedad del intestino. "
        "Las matrices no se fusionaron. El tejido de anclaje es el de CPTAC, PDC000116, TMT, tumor y mucosa del mismo estudio. "
        "La réplica tisular es PXD034575: 16 pacientes con tumor primario, metástasis hepática y mucosa normal. "
        "El suero es PXD055218, adquisición independiente de datos, con 20 personas con cáncer, 20 con colonoscopia negativa y 20 sanas. "
        "El contraste es una t moderada contra un cambio mínimo de 1,5 veces, FDR 0,05, sin imputar ausencias."
    )
    p(
        "En el contraste oficial del suero, 205 proteínas, ninguna separó el cáncer de la colonoscopia negativa: "
        "105 son equivalentes y 100 quedan indeterminadas. Frente a personas sanas pasan SAA1 y LBP. "
        "Las dos están secretadas y se leen como respuesta del organismo. CHMP4A no entra en las 312: su q global de grupo es 0,037. "
        "El máximo del área bajo la curva, sin contaminantes ni paneles, es 0,865 en SERPINA3. "
        "En tejido, el cierre del 31 de octubre coincide en signo con CPTAC en 687 de 874 genes, correlación 0,83. "
        "Entre las 688 certificadas, la corrida v3.2 da 89 de 90 para el primario y 69 de 82 para la metástasis. "
        "La cuenta 116 de 118 queda en results/v3/. "
        "De los procesos del tumor, este suero solo mide alrededor de doce con al menos diez genes, y esa concordancia no supera el azar. "
        "Una fila de grupo proteico no confirma un péptido mutado. La tesis responde la pregunta de visibilidad. No propone un examen de sangre."
    )
    p("Palabras clave: cáncer colorrectal, proteómica, suero, especificidad, Gene Ontology, neoantígeno.")
    h(1, "Abstract")
    p(
        "The colorectal tumor already has a proteome distinct from that of the mucosa. This thesis asks how much of that change "
        "remains visible in serum, and how much of the serum signal is shared with another bowel disease. Cohorts were not fused. "
        "The anchor tissue is CPTAC PDC000116. The tissue replicate is PXD034575, sixteen patients with primary tumor, hepatic metastasis "
        "and normal mucosa. The serum is PXD055218, data-independent acquisition, twenty people with cancer, twenty with a negative colonoscopy "
        "and twenty healthy people. The contrast is a moderated t test against a minimum 1.5-fold change, FDR 0.05, with no imputation."
    )
    p(
        "In the official serum contrast, 205 proteins, none separated cancer from the negative-colonoscopy group: "
        "105 are equivalent and 100 remain indeterminate. SAA1 and LBP separate cancer from healthy people. "
        "CHMP4A does not enter the 312: its global group q-value is 0.037. "
        "After contaminants and panels are removed, the maximum area under the curve is 0.865, on SERPINA3. "
        "The 31 October close agrees in sign with CPTAC in 687 of 874 genes, correlation 0.83. "
        "Among the 688 certified proteins, run v3.2 gives 89 of 90 for the primary and 69 of 82 for the metastasis. "
        "About twelve tumor Gene Ontology terms are measured in this serum with at least ten genes, and that concordance does not beat the null. "
        "A protein-group row cannot confirm a mutated peptide. The thesis answers the visibility question. It does not propose a blood test."
    )
    p("Keywords: colorectal cancer, proteomics, serum, specificity, Gene Ontology, neoantigen.")
    h(1, "Contenido")
    p(
        "El orden es el de las tesis de maestría de la Universidad de Caldas: introducción, justificación, política pública y objetivos; "
        "marco teórico; materiales y métodos; análisis de resultados; conclusiones, trabajos futuros e impacto; referencias. "
        "En Word, clic derecho sobre el campo de contenido y elija Actualizar campos. En LaTeX, el contenido se genera al compilar dos veces."
    )
    h(1, "Lista de tablas")
    p("Tabla 1. Cohortes y papel de cada una.")
    p("Tabla 2. Contrastes del suero PXD055218.")
    p("Tabla 3. Contrastes pareados del tejido PXD034575.")
    p("Tabla 4. Concordancia del tejido de Tanaka con CPTAC.")
    p("Tabla 5. Candidatas con AUC por encima de la regla y su lectura.")


def capitulo_1() -> None:
    h(1, "Capítulo 1. Introducción")
    h(2, "Introducción")
    p(
        "El cáncer colorrectal se estudia desde hace años en el tejido. El proteoma del tumor y el de la mucosa del mismo estudio "
        "ya fueron medidos en cohortes públicas, entre ellas la de CPTAC publicada por Vasaikar y colaboradores en 2019 @@vasaikar2019@@. "
        "La pregunta clínica que queda abierta es otra: si una parte de ese cambio se puede leer en sangre, y si lo que se lee en sangre "
        "pertenece al cáncer o también aparece cuando el intestino está enfermo por otra causa."
    )
    p(
        "Comparar cáncer solo contra personas sanas enciende proteínas de respuesta del organismo. Esa comparación describe inflamación. "
        "La comparación que decide una señal propia del cáncer es cáncer frente a otra enfermedad intestinal. Rai y colaboradores, en 2025 @@rai2025@@, "
        "midieron suero con tres grupos: cáncer, colonoscopia negativa y sanos. Esta tesis usa esa matriz para auditar la especificidad, "
        "y usa el tejido para preguntar si el cambio del tumor se repite y si llega a la sangre."
    )
    p(
        "El trabajo es de bioinformática sobre datos públicos. No recluta pacientes. No fusiona intensidades de plataformas distintas. "
        "No presenta un panel para la consulta. El producto es una lectura, con reglas escritas antes del resultado, de lo que el suero de esta matriz deja ver."
    )
    h(2, "Justificación")
    p(
        "Integrar proteomas públicos de cáncer colorrectal, tejido y biopsia líquida, ya está publicado. Robles, Prakash, Vizcaíno y Casal @@robles2024@@ "
        "lo hicieron en 2024. Esta tesis no reclama esa prioridad. El hueco que sí cierra es operativo: con una regla fijada de antemano, "
        "medir cuánto del cambio tisular se parece entre dos tejidos, cuánto de ese cambio está en el suero, y si las candidatas publicadas "
        "separan el cáncer de una colonoscopia negativa."
    )
    p(
        "Ese hueco importa porque varias listas de proteínas «de cáncer» se arman contra controles sanos, se solapan con fase aguda o con hemólisis, "
        "y luego se leen como marcadores. Una tesis que deje el cero por escrito, con el control de enfermedad y sin relajar el FDR, evita esa lectura."
    )
    h(2, "Política pública a la que se da respuesta")
    p(
        "En Colombia, la Ley 1384 de 2010 @@ley1384@@ ordena acciones para la atención integral del cáncer, incluida la detección temprana. "
        "Un examen de sangre sería, en ese marco, una herramienta deseable. Esta tesis no diseña ese examen. Lo que aporta a la discusión pública "
        "es el estándar de evidencia: una proteína que solo separa el cáncer de una persona sana no basta para una prueba de tamización. "
        "Hace falta que separe el cáncer de otra enfermedad del intestino, en una cohorte más grande que veinte personas por grupo, "
        "y a ser posible junto al antígeno carcinoembrionario medido en la misma persona. El resultado de este suero dice que esa barra, hoy, no se cruza."
    )
    h(2, "Objetivos")
    h(3, "Objetivo general")
    p(
        "Cuantificar cuánto del cambio proteico del tumor colorrectal se repite en otro tejido y cuánto es visible y específico en suero, "
        "sin fusionar matrices de plataformas distintas y sin presentar como propio un resultado ya publicado."
    )
    h(3, "Objetivos específicos")
    viñetas([
        "Medir la concordancia de signo y de magnitud entre el tumor primario de PXD034575 y el tejido certificado de CPTAC, con la leyenda real de las columnas.",
        "Contrastar el suero de tres grupos y evaluar si las candidatas publicadas separan el cáncer de la colonoscopia negativa, además de separarlo de los sanos.",
        "Estimar si los procesos de Gene Ontology alterados en el tejido están medidos en el suero y si su dirección concuerda por encima del azar.",
        "Delimitar qué parte de los péptidos y neoantígenos publicados en 2025 y 2026 puede afirmarse con una matriz de grupos proteicos.",
    ])


def capitulo_2() -> None:
    h(1, "Capítulo 2. Marco teórico")
    h(2, "El tumor y la sangre son dos preguntas")
    p(
        "Una proteína alterada en el tumor sería útil en sangre si el mismo cambio se pudiera medir en el suero y si ese cambio no lo produjera "
        "también otra enfermedad. Malmström y colaboradores, en 2025 @@malmstrom2025@@, describen el plasma como una mezcla cuyo origen es sobre todo el organismo, "
        "no la célula tumoral. Seo y colaboradores, en 2026 @@seo2026@@, encuentran en vesículas proteínas del huésped. Esa literatura prepara el resultado "
        "de esta tesis: lo que separa el cáncer de una persona sana puede ser respuesta del huésped."
    )
    h(2, "Qué mide cada técnica usada aquí")
    p(
        "El tejido de CPTAC llega marcado con TMT. Varias muestras se comparan dentro del mismo experimento y el archivo ya está en log2. "
        "El tejido de Tanaka es LFQ: la intensidad de cada proteína se estima en cada muestra por separado, con MaxQuant. "
        "El suero de Rai es DIA: se intentan ver muchas proteínas a la vez, sin elegirlas de antemano, y el informe las entrega en grupos. "
        "El PRM de Rai, cuatro péptidos, es una medición dirigida que esta tesis no reanaliza. Lesur y colaboradores, en 2023 @@lesur2023@@, miden un panel dirigido de plasma; "
        "esas columnas no entran porque un orden alfabético no es un diagnóstico."
    )
    h(2, "Cómo se lee un cambio de abundancia")
    p(
        "En MaxQuant, el cero significa que la proteína no se identificó en esa muestra. Cox y Mann, en 2008 @@cox2008@@, describen ese programa. "
        "Tratar el cero como cantidad cero hunde la mediana y hace creer que la matriz ya está en logaritmo. Por eso el cero se vuelve un hueco "
        "y solo la intensidad positiva pasa a log2. Un log2FC de 1 es el doble. Un log2FC de −1 es la mitad. El umbral 0,585 es el log2 de 1,5."
    )
    p(
        "Smyth, en 2004 @@smyth2004@@, modera la varianza entre proteínas para que una medición escasa no se declare distinta solo porque su desviación salió pequeña. "
        "McCarthy y Smyth, en 2009 @@mccarthy2009@@, llaman treat a la prueba contra un cambio mínimo, no contra cero. Benjamini y Hochberg, en 1995 @@benjamini1995@@, controlan "
        "la fracción esperada de falsos entre las proteínas que se declaran distintas. El corte usado es 0,05. El AUC ordena personas: 0,50 es el azar. "
        "El AUC no dice si la proteína sube o baja y no sustituye al FDR."
    )
    p(
        "En Tanaka, tumor y mucosa son del mismo paciente: el contraste resta dentro de la persona. En el suero, los grupos son personas distintas: "
        "el contraste compara grupos. Pedir pareo donde no hay pares inventaría una relación."
    )
    h(2, "Gene Ontology")
    p(
        "Ashburner y colaboradores, en 2000 @@ashburner2000@@, proponen Gene Ontology como vocabulario de procesos. Subramanian y colaboradores, en 2005 @@subramanian2005@@, "
        "proponen mirar si los genes de un proceso se acumulan al principio de un ranking, no solo si un gen pasa un corte. "
        "Un suero que mide pocos genes de cada proceso no puede declarar que el proceso del tumor está en la sangre. "
        "Uhlén y colaboradores, en 2015 @@uhlen2015@@, publican concentraciones de referencia en sangre. Esas concentraciones no son las de estos pacientes. "
        "Jumper y colaboradores, en 2021 @@jumper2021@@, predicen un pliegue con AlphaFold. Un pliegue no demuestra que la cantidad haya cambiado."
    )
    h(2, "Neoantígenos")
    p(
        "Un neoantígeno es un péptido cuya secuencia cambió por una mutación y que el tumor puede presentar. "
        "Rai y colaboradores describen 19 péptidos en suero y validan cuatro por PRM, uno de ellos neoantigénico. "
        "Xiang y colaboradores, en 2026 @@xiang2026@@, identifican 96 neoepítopos mutados de MHC-I en tejido de diez pacientes; la mayoría proviene de regiones no codificantes. "
        "Hembrough y colaboradores, en un resumen de congreso de 2026 @@hembrough2026@@, recuperan cinco péptidos de una persona con cáncer colorrectal avanzado, desde vesículas de plasma. "
        "Ninguno de esos trabajos se convierte en una fila de esta tesis. La matriz guarda grupos de proteínas, no la letra mutada."
    )
    h(2, "Estado del arte que cambia la lectura")
    p(
        "Vasaikar y colaboradores, Cell 2019, son el ancla tisular. Tanaka y colaboradores, Cell Reports 2024 @@tanaka2024@@, describen tripletes de primario, "
        "metástasis hepática y mucosa; la columna P de esa matriz es el primario, no un pólipo. "
        "Ögren y colaboradores, en 2026 @@ogren2026@@, publican un panel de plasma para cáncer colorrectal formado por CEACAM5, FLT1, IL19 y ferritina. "
        "Verathamjamras, Chantaraamporn y colaboradores, en 2023 @@verathamjamras2023@@, encuentran la fibronectina más baja en plasma de cáncer que en sanos. "
        "Yu y colaboradores, en 2026 @@yu2026@@, proponen candidatas tempranas que esta tesis no descubre: las audita. "
        "Dhami y colaboradores, en 2026 @@dhami2026@@, muestran proteínas elevadas a la vez en cáncer y en enfermedad inflamatoria intestinal. "
        "García-Foncillas, en 2026 @@garciafoncillas2026@@, revisa la biología molecular del cáncer colorrectal; es una revisión narrativa, no una evidencia de subtipos en proteína. "
        "Zhang y colaboradores, Nature 2014 @@zhang2014@@, son otro proteoma de colon y recto, distinto de PDC000116. "
        "Li y colaboradores, en 2024 @@li2024@@, perfilan el plasma durante el tratamiento con cetuximab: es respuesta a un fármaco, no un estudio de tamización."
    )


def capitulo_3() -> None:
    h(1, "Capítulo 3. Materiales y métodos")
    h(2, "Cohortes")
    p(
        "Cada cohorte se contrasta dentro de sí. No hay corrección de lote entre estudios. No hay personas compartidas que autoricen una fusión. "
        "La configuración que declara PXD015163 descarga, en la práctica, el informe DIA-NN de PXD055218. En la tesis se nombra el archivo analizado."
    )
    tabla(
        "Tabla 1. Cohortes y papel de cada una.",
        ["Accesión", "Qué es", "Papel", "Fuente"],
        [
            ["PDC000116", "Tejido, TMT, ya en log2", "Ancla tumoral", "Vasaikar et al., 2019"],
            ["PXD034575", "16 pacientes; P primario, M metástasis, N mucosa", "Réplica tisular", "Tanaka et al., 2024"],
            ["PXD055218", "Suero DIA, 20+20+20", "Especificidad", "Rai et al., 2025"],
            ["PXD065161", "PRM de cuatro péptidos", "No reanalizado", "Rai et al., 2025"],
            ["PXD036594", "Plasma dirigido", "Excluido: sin etiqueta clínica usable", "Lesur et al., 2023"],
        ],
    )
    h(2, "Reglas fijadas antes de leer el suero")
    viñetas([
        "FDR de Benjamini–Hochberg 0,05 y |log2FC| mayor que 0,585. Los dos cortes se exigen juntos.",
        "Sin imputar ausencias.",
        "En Tanaka, al menos ocho pares. El cero de MaxQuant es ausencia. La intensidad positiva se transforma a log2.",
        "Especificidad: AUC de cáncer frente a colonoscopia negativa mayor que 0,70, límite inferior del intervalo bootstrap mayor que 0,50, y ese AUC mayor que el de enfermedad frente a sanos. Dos mil remuestreos, semilla 20261031.",
        "La marca de fase aguda, hemólisis o enfermedad inflamatoria se escribió en el catálogo antes del contraste. No baja el umbral del AUC. Cambia la lectura.",
        "Vías: términos biológicos con FDR de hasta 0,25 en el tejido de CPTAC, al menos diez genes medidos en el suero, y 200 permutaciones del ranking.",
    ])
    h(2, "Contraste y fichas")
    p(
        "El suero usa los prefijos de corrida CP, CN y Normal. Si esas columnas existen, no se usa otra etiqueta. "
        "El tejido de Tanaka se lee por el sufijo de la columna, no por la etiqueta «pólipo» que una versión anterior del pipeline había asignado a P. "
        "Esa versión certificada no se reescribe: sus números quedan como otra estadística."
    )
    p(
        "La concordancia con CPTAC empareja genes y calcula el mismo signo y la correlación de Pearson de los log2FC. "
        "Las fichas reúnen proteínas que pasan un FDR o que tienen una lectura de AUC evaluable, y les añaden términos GO de 15 a 300 genes, "
        "la anotación de secreción del Human Protein Atlas y los enlaces a UniProt y AlphaFold. El origen de cada fila es literatura o resultado de esta corrida."
    )
    h(2, "Lo que el método no puede hacer")
    p(
        "No busca espectros de péptidos mutados. No descarga el suplemento de origen tisular de Malmström como tabla aplicada paciente a paciente. "
        "No parte PXD036594 en cáncer y control. No estima un examen clínico. Con veinte personas por grupo el alcance es de cribado."
    )


def capitulo_4() -> None:
    h(1, "Capítulo 4. Análisis de resultados")
    h(2, "El suero no separa el cáncer de la otra enfermedad")
    p(
        "PXD055218 tiene los tres grupos, veinte corridas cada uno. El contraste oficial de v3.2 tiene 205 proteínas. "
        "Cáncer frente a colonoscopia negativa: cero proteínas dentro del FDR y del cambio mínimo. "
        "De esas 205, 105 tienen un cambio menor que 1,5 veces demostrado y 100 quedan indeterminadas. "
        "La ausencia solo se afirma en las 105. Cáncer frente a sanos: SAA1 y LBP. "
        "El 687 de 874 y la correlación 0,83 son el cierre del 31 de octubre. Entre las 688 certificadas, "
        "el primario coincide en 89 de 90 y la metástasis en 69 de 82. La cuenta 116 de 118 queda en results/v3/."
    )
    tabla(
        "Tabla 2. Proteínas de la matriz sin filtrar del 31 de octubre que pasan FDR 0,05 y |log2FC| mayor que 0,585.",
        ["Proteína", "Acceso", "Contraste", "log2FC", "FDR", "Observaciones"],
        [
            ["SAA1", "P0DJI8", "Cáncer frente a sano", "3,26", "0,0002", "40"],
            ["LBP", "P18428", "Cáncer frente a sano", "1,25", "0,007", "40"],
            ["CHMP4A", "Q9BY43", "Cáncer frente a sano", "5,42", "3,4e-7", "16 de 40"],
            ["APOA4", "P06727", "Enfermedad frente a sano", "−1,32", "0,022", "40"],
        ],
    )
    p(
        "SAA1 y LBP están anotadas como secretadas a la sangre. En la ficha de LBP aparece el término GO de respuesta de fase aguda. "
        "Frente a la colonoscopia negativa no pasan el FDR. CHMP4A tiene un cambio grande y la mitad de los sueros sin medición. "
        "Un hueco no es un cero biológico. APOA4 baja en la enfermedad intestinal frente a los sanos y no separa el cáncer de esa enfermedad."
    )
    h(2, "El AUC es más laxo que el FDR")
    p(
        "Ocho candidatas de la literatura cumplen la regla de AUC. Cinco ya estaban en la lista de fase aguda: C9, LRG1, SERPINA1, CRP y HP. "
        "Tres no tienen esa marca y tampoco pasan el FDR: BAIAP3, FN1 y TPM4."
    )
    tabla(
        "Tabla 3. Candidatas con la regla de AUC cumplida.",
        ["Gen", "AUC cáncer vs enfermedad", "Lectura", "FDR del mismo contraste"],
        [
            ["BAIAP3", "0,82", "Sin marca de ruido", "0,20"],
            ["FN1", "0,70", "Sin marca de ruido", "No pasa"],
            ["TPM4", "0,80", "Sin marca de ruido", "No pasa"],
            ["C9", "0,81", "Fase aguda", "No pasa"],
            ["LRG1", "0,78", "Fase aguda", "No pasa"],
            ["SERPINA1", "0,73", "Fase aguda", "No pasa"],
            ["CRP", "0,73", "Fase aguda", "No pasa"],
            ["HP", "0,71", "Fase aguda", "No pasa"],
        ],
    )
    p(
        "BAIAP3, acceso O94812, es la proteína que Rai señaló en suero. El log2FC frente a la colonoscopia negativa es 2,20. "
        "El AUC pasa y el FDR no. Queda como candidata de otro estudio, con más personas. No es un hallazgo de esta tesis. "
        "FN1 en este suero está un poco más baja en cáncer que en sanos, log2FC −0,27, en la misma dirección que el plasma publicado en 2023. "
        "CEACAM5 no está medida en esta matriz DIA, así que su panel sérico no se puede auditar aquí."
    )
    h(2, "El tejido, leído con la leyenda del artículo")
    p(
        "PXD034575 tiene 16 pacientes y tres columnas por paciente. P es tumor primario, M es metástasis hepática y N es mucosa. "
        "Después de tratar el cero como ausencia y pasar a log2, la escala queda como logaritmo de intensidad lineal. "
        "CEACAM5 sube en el primario: log2FC 3,87, FDR 0,002, 13 pares. Esa dirección confirma que el contraste quedó orientado."
    )
    tabla(
        "Tabla 4. Contrastes pareados de PXD034575.",
        ["Contraste", "Significativas", "Probadas"],
        [
            ["Primario frente a mucosa", "58", "2490"],
            ["Metástasis frente a mucosa", "160", "2430"],
            ["Metástasis frente a primario", "16", "2710"],
        ],
    )
    tabla(
        "Tabla 5. Concordancia con el tejido certificado de CPTAC.",
        ["Comparación", "Mismo signo", "Genes", "Correlación"],
        [
            ["Primario frente a mucosa", "687 (78,6 %)", "874", "0,83"],
            ["Metástasis frente a mucosa", "597 (70,9 %)", "842", "0,59"],
        ],
    )
    p(
        "Entre primario y metástasis el signo coincide en 708 de 876 genes, el 80,8 %. Entre los que coinciden, solo el 46,2 % tiene magnitud menor en el primario. "
        "No queda un gradiente en el que el primario sea una versión pequeña de la metástasis. "
        "La correlación 0,104 y el 56 % de igual signo, sobre 1 993 proteínas, pertenecen a un cálculo anterior, con la etiqueta vieja de columnas. "
        "Esa cifra se conserva. No se sustituye por 0,83. Ninguna de las dos mide tejido contra sangre."
    )
    h(2, "Las vías del tumor casi no están en este suero")
    p(
        "De los términos GO biológicos con FDR menor o igual a 0,25 en CPTAC, entran al suero los que tienen al menos diez genes medidos. "
        "Son 12 en cáncer frente a sano, 13 en cáncer frente a enfermedad y 12 en enfermedad frente a sano. "
        "La correlación de sus puntuaciones con el tejido es −0,28, −0,31 y −0,04. El nulo, 200 permutaciones, da p 0,34, 0,30 y 0,93. "
        "Este suero no contiene el grueso de los procesos que el tumor altera."
    )
    h(2, "Lo que la versión certificada ya había medido hacia la sangre")
    p(
        "Una corrida anterior, sin rehacerse en este cierre, comparó la dirección del cambio en tejido y en plasma: coincide en 42 de 86 proteínas, 0,488. "
        "Dentro de la ventana de concentración en la que el suero sí puede ver la proteína, coincide en 34 de 63, 0,540, con intervalo 0,409 a 0,666. "
        "Ese intervalo incluye el azar. La falta de transferencia no se explica solo porque el suero no alcance a medir proteínas escasas."
    )
    h(2, "Péptidos")
    p(
        "La tabla de péptidos de esta corrida fija un límite, no una lista descubierta. Los 19 péptidos de Rai, los 96 neoepítopos de Xiang "
        "y los cinco péptidos de Hembrough quedan fuera de las filas de abundancia. En tejido CPTAC, 33 candidatas del catálogo se detectan "
        "y siete repiten el signo de la literatura con FDR: CEACAM5, CA1, SMPDL3A, HBB, HBD, LCN2 y TOP2A. CA1, HBB y HBD quedan marcadas como posible hemólisis. "
        "LCN2 queda marcada como compartida con inflamación intestinal. TOP2A sube en tejido y no tiene contraste evaluable en este suero."
    )
    h(2, "Discusión")
    p(
        "El resultado se lee en dos capas. En tejido, el cambio del primario se parece al de CPTAC. La biología del tumor, medida como abundancia, "
        "se repite cuando las columnas se nombran como el artículo las nombra. En suero, esa biología no aparece como señal propia del cáncer. "
        "Lo que el suero separa de una persona sana es, en dos de los tres casos, una proteína secretada de respuesta del huésped. "
        "Eso encaja con el atlas de Malmström y con la señal de huésped descrita por Seo. No encaja con usar SAA1 o LBP como prueba de tumor."
    )
    p(
        "Rai y colaboradores combinaron proteínas con el antígeno carcinoembrionario y reportaron un AUC de 0,758 en su propio estudio. "
        "Esta matriz no trae ese inmunoensayo, así que no se puede preguntar si alguna proteína añade algo al marcador que ya se usa. "
        "BAIAP3 sobrevive a la regla de AUC sin marca de fase aguda y no sobrevive al FDR. Con veinte por grupo, eso es una candidata, no una validación. "
        "Ögren y colaboradores arman un panel que incluye CEACAM5. Aquí CEACAM5 sube en el primario y no está en el suero no dirigido. "
        "El diálogo con ese panel queda abierto para una medición dirigida, no para esta matriz."
    )
    p(
        "Xiang mide presentación por MHC en el tumor. Hembrough recupera péptidos de una persona. Esta tesis mide grupos de proteínas. "
        "Las tres preguntas son vecinas y no son la misma. Afirmar un neoantígeno a partir de estas filas cruzaría ese límite."
    )


def capitulo_5() -> None:
    h(1, "Capítulo 5. Conclusiones")
    h(2, "Conclusiones")
    p("Cada frase corresponde a un objetivo.")
    viñetas([
        "El tumor primario de Tanaka, leído como primario y no como pólipo, coincide en signo con CPTAC en 687 de 874 genes y correlaciona 0,83. La metástasis correlaciona 0,59. El primario no es una metástasis pequeña: solo el 46,2 % de los genes de igual signo tiene magnitud menor en el primario.",
        "En el suero de tres grupos, ninguna proteína separa el cáncer de la colonoscopia negativa con FDR 0,05 y cambio de al menos 1,5 veces. SAA1, LBP y CHMP4A separan el cáncer de los sanos. Las dos primeras se leen como respuesta del huésped. La tercera tiene demasiadas ausencias para presentarse como marcador.",
        "Los procesos GO del tejido casi no están medidos en este suero. Los doce a trece términos que sí entran no concuerdan por encima de 200 permutaciones.",
        "Los péptidos mutados de Rai, Xiang y Hembrough no se confirman ni se rechazan con esta matriz. El resultado de ese objetivo es el límite del método.",
    ])
    p(
        "En conjunto, el cambio de proteínas del tumor se repite entre tejidos y en este suero no aparece como una señal propia del cáncer. "
        "La tesis sirve para responder esa pregunta. No sirve para proponer un examen de sangre."
    )
    h(2, "Trabajos futuros")
    viñetas([
        "Repetir BAIAP3 en otra cohorte con colonoscopia negativa, más de veinte personas por grupo, y con el antígeno carcinoembrionario de la misma persona.",
        "Copiar del suplemento de Rai las secuencias de los 19 péptidos y rebuscarlas en crudos independientes, con control de falsos por clase de péptido. Eso no se anticipa con un número.",
        "Medir CEACAM5 de forma dirigida en suero, que es la escala del panel de Ögren, no la de esta matriz DIA.",
        "Emparejar la tabla de muestras de Lesur con las columnas de PXD036594 antes de usar ese plasma.",
        "Aplicar el suplemento de origen tisular de Malmström cuando esté en el disco. Hoy la ficha usa secreción y GO.",
        "Incorporar, cuando exista, tejido y sangre del mismo paciente. Estas matrices no los traen juntos.",
    ])
    h(2, "Impacto")
    p(
        "El impacto es de método y de lectura, no de tamización. Queda escrito un control que muchos trabajos omiten: cáncer frente a otra enfermedad del intestino, "
        "con el umbral sin relajar. Queda corregida la leyenda de PXD034575 para cualquier lectura posterior. "
        "Queda una página pública, https://oncogo-crc-tesis.web.app, con las fichas, la teoría, los diccionarios y los enlaces de los estudios, "
        "para estudiar la sustentación y para que otra persona pueda verificar un DOI. "
        "No hay cohorte colombiana en estos datos públicos. La Ley 1384 de 2010 enmarca la detección temprana; este resultado dice qué evidencia faltaría antes de hablar de una prueba."
    )
    h(2, "Limitaciones")
    viñetas([
        "Veinte sueros por grupo. El alcance es de cribado.",
        "TMT, LFQ y DIA no se pueden sumar.",
        "No hay tejido y sangre del mismo paciente.",
        "La colonoscopia negativa agrupa otras enfermedades intestinales; no es una cohorte pura de enfermedad inflamatoria recién diagnosticada.",
        "PXD036594 sigue sin etiqueta clínica usable.",
        "No hay población latinoamericana en estas matrices.",
        "Los neoantígenos no se rebuscaron en los crudos.",
        "La correlación 0,104 de la versión anterior y la 0,83 de este cierre son estadísticas distintas y deben citarse juntas, cada una con su n.",
    ])


def referencias() -> None:
    h(1, "Referencias bibliográficas")
    p(
        "La lista está numerada en el orden de aparición, como en las tesis de la Facultad de Ingenierías. "
        "Cada DOI se confirmó en Crossref el 31 de octubre de 2026 y el título del registro coincide con la cita. "
        "La ficha completa está en referencias.bib."
    )
    p(
        "No entran, porque el identificador no abrió un artículo el 31 de octubre de 2026: doi:10.1016/j.xcrm.2025.101689, "
        "la cita mal formada de ENCODE en Nature, y el DOI 10.64898/2026.08.30.26361622. "
        "Los preprints de Lee, Tan y Ramsbottom pueden ir a un anexo marcados como no revisados por pares."
    )


def latex_escape(texto: str) -> str:
    reemplazos = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    salida = []
    for caracter in texto:
        salida.append(reemplazos.get(caracter, caracter))
    return "".join(salida)


def _titulo_latex(texto: str) -> str:
    import re

    return latex_escape(re.sub(r"^Capítulo \d+\. ", "", texto))


def orden_citas() -> list[str]:
    import re

    vistos: list[str] = []
    marca = re.compile(r"@@([A-Za-z0-9]+)@@")

    def tomar(texto: str) -> None:
        for hallazgo in marca.finditer(texto):
            clave = hallazgo.group(1)
            if clave not in vistos:
                vistos.append(clave)

    for bloque in BLOQUES:
        if bloque[0] == "p":
            tomar(bloque[1])
        elif bloque[0] == "h":
            tomar(bloque[2])
        elif bloque[0] == "ul":
            for item in bloque[1]:
                tomar(item)
        elif bloque[0] == "tabla":
            for fila in bloque[3]:
                for celda in fila:
                    tomar(celda)
    return vistos


def citar_word(texto: str, orden: list[str]) -> str:
    import re

    return re.sub(r"@@([A-Za-z0-9]+)@@", lambda m: f"[{orden.index(m.group(1)) + 1}]", texto)


def leer_bib(ruta) -> dict[str, dict[str, str]]:
    import re

    entradas: dict[str, dict[str, str]] = {}
    clave = None
    campos: dict[str, str] = {}
    for linea in ruta.read_text(encoding="utf-8").splitlines():
        inicio = re.match(r"@\w+\{(\w+),", linea)
        if inicio:
            clave = inicio.group(1)
            campos = {}
            continue
        if linea.strip() == "}" and clave:
            entradas[clave] = campos
            clave = None
            continue
        campo = re.match(r"\s*(\w+)\s*=\s*\{(.*)\},?\s*$", linea)
        if campo and clave:
            valor = campo.group(2).replace(r"\_", "_")
            while valor.startswith("{") and valor.endswith("}"):
                valor = valor[1:-1]
            campos[campo.group(1)] = valor
    return entradas


def _autor_ieee(persona: str) -> str:
    persona = persona.strip()
    if "," not in persona:
        return persona
    apellido, nombres = persona.split(",", 1)
    iniciales = " ".join(f"{parte[0]}." for parte in nombres.split() if parte)
    return f"{iniciales} {apellido.strip()}"


def linea_ieee(numero: int, campos: dict[str, str]) -> str:
    personas = [p for p in campos.get("author", "").split(" and ") if p.strip()]
    if len(personas) == 1:
        autores = _autor_ieee(personas[0])
    elif len(personas) == 2:
        autores = f"{_autor_ieee(personas[0])} and {_autor_ieee(personas[1])}"
    else:
        autores = f"{_autor_ieee(personas[0])} et al."
    revista = campos.get("journal") or campos.get("booktitle") or ""
    trozos = [f"[{numero}] {autores}, “{campos.get('title', '')},”"]
    if revista:
        trozos.append(f" {revista},")
    if campos.get("volume"):
        trozos.append(f" vol. {campos['volume']},")
    if campos.get("number"):
        trozos.append(f" no. {campos['number']},")
    if campos.get("pages"):
        trozos.append(f" pp. {campos['pages']},")
    if campos.get("year"):
        trozos.append(f" {campos['year']},")
    if campos.get("doi"):
        trozos.append(f" doi: {campos['doi']}.")
    texto = "".join(trozos).rstrip(",")
    if not texto.endswith("."):
        texto += "."
    return texto


def a_latex() -> str:
    partes = [
        "% Archivo principal de Overleaf. Compiler: pdfLaTeX. Overleaf ejecuta Biber.",
        "% Estilo IEEE numerado, el de las tesis de la Facultad de Ingenierías.",
        r"\documentclass[12pt,oneside]{report}",
        r"\usepackage[english,spanish]{babel}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{geometry}",
        r"\usepackage{longtable}",
        r"\usepackage{booktabs}",
        r"\usepackage{csquotes}",
        r"\usepackage[backend=biber,style=ieee,sorting=none,doi=true,isbn=false,url=false]{biblatex}",
        r"\addbibresource{referencias.bib}",
        r"\AtBeginBibliography{\selectlanguage{english}}",
        r"\usepackage{hyperref}",
        r"\usepackage{setspace}",
        r"\geometry{margin=2.5cm}",
        r"\onehalfspacing",
        r"\hypersetup{colorlinks=true,linkcolor=black,urlcolor=blue,citecolor=black}",
        r"\title{Qué parte del tumor colorrectal se deja ver en la sangre}",
        r"\author{[Nombre completo del estudiante]\\Universidad de Caldas}",
        r"\date{Manizales, 2026}",
        r"\begin{document}",
        r"\chapter*{Nota para el autor}",
        latex_escape(BLOQUES[0][1]),
        r"\begin{titlepage}",
        r"\centering",
        r"\vspace*{2cm}",
        r"{\Large Universidad de Caldas}\\[0.4cm]",
        r"{\large Facultad de Ingenierías}\\[2cm]",
        r"{\LARGE Qué parte del tumor colorrectal se deja ver en la sangre}\\[2cm]",
        r"{\Large [Nombre completo del estudiante]}\\[2cm]",
        r"Tesis presentada como requisito parcial para optar al título de:\\[0.3cm]",
        r"Magíster en Bioinformática y Biología Computacional\\[2cm]",
        r"Manizales, Colombia\\2026",
        r"\end{titlepage}",
        r"\tableofcontents",
        r"\listoftables",
    ]
    fase = "portada"
    for bloque in BLOQUES:
        if fase == "portada":
            if bloque[0] == "h" and bloque[1] == 1 and bloque[2] == "Epígrafe":
                fase = "cuerpo"
            else:
                continue
        if bloque[0] == "h":
            nivel, texto = bloque[1], _titulo_latex(bloque[2])
            if nivel == 1 and bloque[2] in {
                "Epígrafe",
                "Agradecimientos",
                "Resumen",
                "Abstract",
                "Contenido",
                "Lista de tablas",
            }:
                partes.append(r"\chapter*{" + texto + "}")
                partes.append(r"\addcontentsline{toc}{chapter}{" + texto + "}")
            elif nivel == 1:
                partes.append(r"\chapter{" + texto + "}")
            elif nivel == 2:
                partes.append(r"\section{" + texto + "}")
            else:
                partes.append(r"\subsection{" + texto + "}")
        elif bloque[0] == "p":
            partes.append(latex_escape(bloque[1]) + "\n")
        elif bloque[0] == "ul":
            partes.append(r"\begin{itemize}")
            for item in bloque[1]:
                partes.append(r"\item " + latex_escape(item))
            partes.append(r"\end{itemize}")
        elif bloque[0] == "tabla":
            _, titulo, columnas, filas = bloque
            spec = "p{3.2cm}" * len(columnas)
            # better column spec
            ancho = round(15.5 / len(columnas), 2)
            spec = " ".join([f"p{{{ancho}cm}}" for _ in columnas])
            partes.append(r"\begin{longtable}{" + spec + "}")
            partes.append(r"\caption{" + latex_escape(titulo) + r"}\\")
            partes.append(r"\toprule")
            partes.append(" & ".join(latex_escape(c) for c in columnas) + r"\\")
            partes.append(r"\midrule")
            partes.append(r"\endfirsthead")
            partes.append(r"\toprule")
            partes.append(" & ".join(latex_escape(c) for c in columnas) + r"\\")
            partes.append(r"\midrule")
            partes.append(r"\endhead")
            for fila in filas:
                partes.append(" & ".join(latex_escape(c) for c in fila) + r"\\")
            partes.append(r"\bottomrule")
            partes.append(r"\end{longtable}")
    partes.append(r"\printbibliography[heading=none]")
    partes.append(r"\end{document}")
    return "\n".join(partes)


def a_docx() -> None:
    from docx import Document
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Cm, Pt, RGBColor

    doc = Document()
    for section in doc.sections:
        section.top_margin = Cm(2.5)
        section.bottom_margin = Cm(2.5)
        section.left_margin = Cm(3)
        section.right_margin = Cm(2.5)
    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(12)
    normal.paragraph_format.line_spacing = 1.5
    normal.paragraph_format.space_after = Pt(8)
    for estilo in ("Heading 1", "Heading 2", "Heading 3"):
        fuente = doc.styles[estilo].font
        fuente.name = "Times New Roman"
        fuente.color.rgb = RGBColor(0x1C, 0x24, 0x30)

    def campo_toc(parrafo) -> None:
        run = parrafo.add_run()
        r = run._r
        inicio = OxmlElement("w:fldChar")
        inicio.set(qn("w:fldCharType"), "begin")
        instr = OxmlElement("w:instrText")
        instr.set(qn("xml:space"), "preserve")
        instr.text = r' TOC \o "1-3" \h \z \u '
        separado = OxmlElement("w:fldChar")
        separado.set(qn("w:fldCharType"), "separate")
        fin = OxmlElement("w:fldChar")
        fin.set(qn("w:fldCharType"), "end")
        r.append(inicio)
        r.append(instr)
        r.append(separado)
        parrafo.add_run("Clic derecho y Actualizar campos para ver el contenido.")
        parrafo.add_run()._r.append(fin)

    orden = orden_citas()
    toc_puesto = False
    for bloque in BLOQUES:
        if bloque[0] == "h" and bloque[1] == 0:
            par = doc.add_paragraph(bloque[2])
            par.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for run in par.runs:
                run.bold = True
                run.font.size = Pt(16)
            continue
        if bloque[0] == "h":
            texto = bloque[2]
            if texto == "Contenido" and not toc_puesto:
                doc.add_page_break()
                doc.add_heading(texto, level=1)
                campo_toc(doc.add_paragraph())
                toc_puesto = True
                continue
            if texto.startswith("Capítulo") or texto in {"Epígrafe", "Agradecimientos", "Resumen", "Abstract", "Lista de tablas", "Referencias bibliográficas"}:
                doc.add_page_break()
            doc.add_heading(texto, level=bloque[1])
            continue
        if bloque[0] == "p":
            doc.add_paragraph(citar_word(bloque[1], orden))
            continue
        if bloque[0] == "ul":
            for item in bloque[1]:
                doc.add_paragraph(citar_word(item, orden), style="List Bullet")
            continue
        if bloque[0] == "tabla":
            _, titulo, columnas, filas = bloque
            doc.add_paragraph(titulo).runs[0].bold = True
            tabla_doc = doc.add_table(rows=1, cols=len(columnas))
            tabla_doc.style = "Table Grid"
            for i, col in enumerate(columnas):
                tabla_doc.rows[0].cells[i].text = col
            for fila in filas:
                celdas = tabla_doc.add_row().cells
                for i, valor in enumerate(fila):
                    celdas[i].text = citar_word(valor, orden)
            doc.add_paragraph("")
    fichas = leer_bib(SALIDA / "referencias.bib")
    for indice, clave in enumerate(orden, start=1):
        doc.add_paragraph(linea_ieee(indice, fichas[clave]))
    doc.save(SALIDA / "tesis_modelo_caldas.docx")


def main() -> None:
    portada()
    capitulo_1()
    capitulo_2()
    capitulo_3()
    capitulo_4()
    capitulo_5()
    referencias()
    import re
    import zipfile

    SALIDA.mkdir(parents=True, exist_ok=True)
    tex_texto = re.sub(r"@@([A-Za-z0-9]+)@@", r"\\cite{\1}", a_latex())
    tex = SALIDA / "tesis_modelo_caldas.tex"
    tex.write_text(tex_texto, encoding="utf-8")
    proyecto = REPO / "overleaf"
    proyecto.mkdir(parents=True, exist_ok=True)
    (proyecto / "main.tex").write_text(tex_texto, encoding="utf-8")
    a_docx()
    paquete = SALIDA / "oncogo-crc-overleaf.zip"
    with zipfile.ZipFile(paquete, "w", zipfile.ZIP_DEFLATED) as archivo:
        archivo.write(proyecto / "main.tex", "main.tex")
        archivo.write(proyecto / "referencias.bib", "referencias.bib")
    print(tex)
    print(SALIDA / "tesis_modelo_caldas.docx")
    print(paquete)


if __name__ == "__main__":
    main()
