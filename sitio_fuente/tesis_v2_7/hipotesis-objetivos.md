1.3 Pregunta de investigación

¿Qué parte del cambio proteico del tumor colorrectal puede verse en el suero cuando el contraste es otra enfermedad intestinal y no una persona sana?

1.4 Hipótesis

La hipótesis principal (HP) se desprende de la pregunta. Si el tumor colorrectal deja en la sangre una huella propia, al menos una proteína del suero debería diferir entre las personas con cáncer (CP) y las personas con colonoscopia negativa (CN) en más de 1,5 veces. La hipótesis nula es que ninguna proteína supera ese margen. Se decide con TREAT y TOST con un margen de log2(1,5) = 0,585 y corrección de Benjamini-Hochberg al 5 % (secciones 3.4 y 4.4). En el resto del documento esta prueba aparece también como la «celda preregistrada» del suero.

Alrededor de la HP hay hipótesis que responden partes concretas de la pregunta. Su numeración viene de los preregistros del repositorio (config/preregistro_v2.yaml y config/preregistro_v3.yaml) y se conserva para que el documento, el código y el sitio web usen los mismos nombres:

• H1 (tejido, heredada de v2): entre dos cohortes de tejido, las proteínas individuales coincidirían poco en la dirección del cambio (concordancia de signo menor que 0,70), pero las vías sí (mayor que 0,80).

• H2 (tejido): la estructura de las redes de coexpresión difiere entre cohortes.

• H4 (medibilidad): se da por falsada si al menos el 90 % de las proteínas candidatas ya es medible en sangre.

• H5 (especificidad de órgano): una parte importante del cambio tumoral del colon es común a otros carcinomas.

• H6 (origen celular): la concordancia entre cohortes es distinta en las proteínas epiteliales y en las del estroma o inmunes.

• H7 (visibilidad en sangre): lo que se sabe de una proteína antes de mirar la sangre (localización, secreción, concentración y otras variables) anticipa si su cambio tumoral se verá en sangre en la misma dirección.

• H8 (carga tumoral): lo que sube en el suero en estadio IV frente a estadios I–II está enriquecido en proteínas tumorales epiteliales.

• H9 (dirección causal, exploratoria): una variante genética que modifica el nivel de una proteína en sangre modifica también el riesgo de cáncer colorrectal.

• H10 (especificidad clínica, descriptiva): las proteínas de partida son específicas del cáncer colorrectal frente a otras enfermedades.

La HP, H1 y H5 son las hipótesis confirmatorias; H6, H7 y H8 son secundarias; H9 es exploratoria y H10 es descriptiva. Por objetivo específico (sección 1.7.2), H1, H5 y H6 responden al objetivo 1, de tejido; la HP, H4, H7, H8, H9 y H10 responden al objetivo 2, de suero, porque todas usan el dato de la sangre; y H2, sobre las redes de coexpresión, queda en el objetivo 3, de interpretación funcional. La Tabla 5 (sección 3.9) reúne la regla de decisión de cada una y la Tabla 8 (sección 4.1), el resultado por objetivo; el detalle está en las secciones 4.2, 4.3, 4.4 y 4.8.

Sobre la numeración: H0 y H3 no aparecen en este documento y H4 no tiene veredicto. H0 era una auditoría de ficheros compartidos entre depósitos de PRIDE prevista en el preregistro v2; no tiene resultado en la versión v3.7 y aquí se trata como un supuesto de los datos. H3 se planteó en una etapa temprana del proyecto, pero no entró en ninguno de los dos preregistros. H4 tiene una regla escrita (umbral de 0,90), pero su veredicto no se calculó en la versión v3.7. No se renumeró para no romper la correspondencia con el repositorio y con el sitio web.

En pocas palabras, los resultados fueron estos. La HP no se rechazó: ninguna proteína del suero fue diferente entre CP y CN y 166 fueron equivalentes. H1 se rechazó, pero no porque el tejido coincidiera poco, sino porque coincidió mucho más de lo que H1 predecía (147 de 147 proteínas con el mismo signo). H5 quedó no refutada, porque solo se cumplió una de las dos condiciones de refutación. H7 se refutó en los dos recursos que se pudieron evaluar; el tercero no fue evaluable. H6 no fue informativa, H8 no fue evaluable, H2 no se evaluó, H9 no mostró señal y H10 es descriptiva.

1.5 Justificación

Integrar proteomas públicos de cáncer colorrectal, de tejido y de biopsia líquida, no es nuevo; ya se ha hecho antes (Robles et al., 2024). El vacío que se aborda es operativo y metodológico. Con un umbral de cambio fijado de 1,5 veces (incremento del 50 % o una reducción a dos tercios en la concentración de una proteína respecto al grupo de control) se mide la concordancia entre dos tejidos, se clasifica cada proteína del suero como diferente, equivalente o indeterminada, y se pone a prueba si la separación entre cáncer y colonoscopia negativa sobrevive cuando se tiene en cuenta la calidad de la muestra.

El vacío importa porque varias listas de proteínas candidatas se construyen frente a controles sanos y luego se leen como marcadores. Las alteraciones de la muestra (coagulación, contaminación por eritrocitos o plaquetas, manipulación) pueden producir diferencias aparentes (Geyer et al., 2019), y un efecto de lote puede imitar una diferencia biológica (Leek et al., 2010). Documentar el resultado nulo, con el control de enfermedad y con control de la tasa de descubrimientos falsos, evita una lectura optimista y aporta un procedimiento reproducible para otras cohortes.

El aporte también es formativo para la línea de investigación: el trabajo deja un pipeline, de modo que cualquier cifra del documento se puede rastrear hasta el repositorio (Anexo B).

1.6 Política pública e impacto esperado

En Colombia, la Ley 1384 de 2010 establece acciones para la atención integral del cáncer, incluida la detección temprana (Congreso de la República de Colombia, 2010). En ese marco, una prueba de sangre para el cáncer colorrectal sería una herramienta deseable. Este trabajo no diseña esa prueba. Lo que aporta a la discusión pública es un estándar de evidencia: una proteína que solo separa el cáncer de una persona sana no basta para proponer una prueba de tamización; debe separarlo de otra enfermedad intestinal, en una cohorte más grande que veinte personas por grupo, con el orden de adquisición aleatorizado y validada en una población externa (Pepe et al., 2001).

Las matrices públicas analizadas no incluyen una cohorte colombiana identificada como tal, esa limitación refuerza la necesidad de cohortes locales antes de trasladar cualquier hallazgo a la práctica clínica del país.

1.7 Objetivos

1.7.1 Objetivo general

Cuantificar qué parte del cambio proteico del tumor colorrectal observado en el tejido se repite en una cohorte independiente de tejido y se ve de forma específica en el suero cuando la comparación es con personas con colonoscopia negativa, mediante el reanálisis de datos proteómicos públicos.

1.7.2 Objetivos específicos

1. Caracterizar el cambio proteico del tumor colorrectal frente a la mucosa y medir su concordancia entre dos cohortes independientes de tejido, CPTAC y Tanaka, incluidas su especificidad de órgano y su relación con el origen celular (H1, H5 y H6).

2. Evaluar si ese cambio aparece en el suero al comparar cáncer con colonoscopia negativa, con una regla de diferencia y de equivalencia escrita de antemano, y separar el efecto de la calidad de la muestra y del orden de adquisición con los modelos M1 a M6 (HP, H4, H7, H8, H9 y H10).

3. Interpretar funcionalmente el cambio del tejido con Gene Ontology y documentar un pipeline reproducible, con sello, versión v3.7 y sitio web documental, que permita rastrear cada cifra hasta su fichero (H2).

Cada objetivo corresponde a una línea del diseño (sección 3.1): el tejido, el suero y la interpretación con la reproducibilidad. La Tabla 8 (sección 4.1) reúne, por objetivo, las hipótesis y el resultado que las responde.

1.8 Alcance y estructura del documento

El trabajo es de bioinformática sobre datos públicos. No se reclutaron pacientes, no se fusionaron intensidades de distintas plataformas y no se propone un panel para la consulta. Con veinte sueros por grupo, el alcance es exploratorio y corresponde a la fase 1 del desarrollo de biomarcadores (Pepe et al., 2001). La colonoscopia negativa reúne a personas con síntomas intestinales cuya colonoscopia no mostró cáncer, con otras enfermedades del intestino; no es una cohorte pura de enfermedad inflamatoria.

Tampoco es un estudio de exactitud diagnóstica. El tejido y el suero vienen de personas distintas, no se valida ningún biomarcador, no se estima utilidad clínica y la tesis no asigna las muestras a un segmento del colon. Ningún resultado de esta tesis sirve para decidir sobre un paciente. Las piezas que el preregistro pedía y que no se ejecutaron (vías del suero, cláusula de vías de H1, unión de HC por XIC, deriva del PRTC, redes de coexpresión, catálogo de candidatas y neoantígenos) no se presentan como resultados: se declaran en la sección 5.4 (Tabla 16).

El documento se organiza así. El capítulo 2 presenta los conceptos que se usan para leer los resultados: el proteoma del tejido y del suero, las técnicas de espectrometría de masas, la inferencia de cambio y de equivalencia, el enriquecimiento funcional, el aprendizaje automático y los sesgos preanalíticos. El capítulo 3 describe los datos, las reglas y el entorno computacional. El capítulo 4 presenta los resultados. El capítulo 5 los discute junto con las desviaciones y las limitaciones. El capítulo 6 resume las conclusiones, las contribuciones y el trabajo futuro. Los anexos documentan la correspondencia con la propuesta, la reproducibilidad y las tablas suplementarias.
