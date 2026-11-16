---
titulo: Pregunta e hipótesis
ruta: /investigacion/pregunta-e-hipotesis
seccion: Investigación
orden: 2
resumen: La pregunta, las hipótesis y los objetivos, con el texto de la tesis y el estado de cada una.
---

# Pregunta e hipótesis

{{resumen}}
La pregunta pregunta qué parte del cambio del tumor se ve en el suero cuando el contraste es otra enfermedad intestinal.
{{/resumen}}

## Pregunta

¿Qué parte del cambio proteico del tumor colorrectal puede verse en el suero cuando el contraste es otra enfermedad intestinal y no una persona sana?

## Hipótesis principal

La hipótesis principal (HP) se desprende de la pregunta. Si el tumor colorrectal deja en la sangre una huella propia, al menos una proteína del suero debería diferir entre las personas con cáncer (CP) y las personas con colonoscopia negativa (CN) en más de 1,5 veces. La hipótesis nula es que ninguna proteína supera ese margen. Se decide con TREAT y TOST con un margen de log2(1,5) = 0,585 y corrección de Benjamini-Hochberg al 5 % (secciones 3.4 y 4.4). En el resto del documento esta prueba aparece también como la «celda preregistrada» del suero.

{{estado:HP}}

## Las demás hipótesis

Alrededor de la HP hay hipótesis que responden partes concretas de la pregunta. Su numeración se conserva para que el documento, el código y el sitio web usen los mismos nombres.

- H1 (tejido, heredada de v2): entre dos cohortes de tejido, las proteínas individuales coincidirían poco en la dirección del cambio (concordancia de signo menor que 0,70), pero las vías sí (mayor que 0,80). {{estado:H1}}
- H2 (tejido): la estructura de las redes de coexpresión difiere entre cohortes. {{estado:H2}}
- H3 se planteó en una etapa temprana del proyecto, pero no entró en ninguno de los dos preregistros. {{estado:H3}}
- H4 (medibilidad): se da por falsada si al menos el 90 % de las proteínas candidatas ya es medible en sangre. {{estado:H4}}
- H5 (especificidad de órgano): una parte importante del cambio tumoral del colon es común a otros carcinomas. {{estado:H5}}
- H6 (origen celular): la concordancia entre cohortes es distinta en las proteínas epiteliales y en las del estroma o inmunes. {{estado:H6}}
- H7 (visibilidad en sangre): lo que se sabe de una proteína antes de mirar la sangre (localización, secreción, concentración y otras variables) anticipa si su cambio tumoral se verá en sangre en la misma dirección. {{estado:H7}}
- H8 (carga tumoral): lo que sube en el suero en estadio IV frente a estadios I–II está enriquecido en proteínas tumorales epiteliales. {{estado:H8}}
- H9 (dirección causal, exploratoria): una variante genética que modifica el nivel de una proteína en sangre modifica también el riesgo de cáncer colorrectal. {{estado:H9}}
- H10 (especificidad clínica, descriptiva): las proteínas de partida son específicas del cáncer colorrectal frente a otras enfermedades. {{estado:H10}}

La HP, H1 y H5 son las hipótesis confirmatorias; H6, H7 y H8 son secundarias; H9 es exploratoria y H10 es descriptiva. Por objetivo específico (sección 1.7.2), H1, H5 y H6 responden al objetivo 1, de tejido; la HP, H4, H7, H8, H9 y H10 responden al objetivo 2, de suero, porque todas usan el dato de la sangre; y H2, sobre las redes de coexpresión, queda en el objetivo 3, de interpretación funcional.

Sobre la numeración: H0 y H3 no aparecen en este documento y H4 no tiene veredicto. H0 era una auditoría de ficheros compartidos entre depósitos de PRIDE prevista en el preregistro v2; no tiene resultado en la versión v3.7 y aquí se trata como un supuesto de los datos. H3 se planteó en una etapa temprana del proyecto, pero no entró en ninguno de los dos preregistros. H4 tiene una regla escrita (umbral de 0,90), pero su veredicto no se calculó en la versión v3.7. No se renumeró para no romper la correspondencia con el repositorio y con el sitio web.

{{especialistas}}
La enmienda de cierre registra la priorización así: {{dato:config/preregistro_v3.yaml#enmienda_v3_7_cierre.descriptivos_que_se_ejecutan.priorizacion_secuencial.regla}}
{{/especialistas}}

H0: era una auditoría de ficheros compartidos entre depósitos de PRIDE prevista en el preregistro v2; no tiene resultado en la versión v3.7 y aquí se trata como un supuesto de los datos. {{estado:H0}}

## Objetivos

Cuantificar qué parte del cambio proteico del tumor colorrectal observado en el tejido se repite en una cohorte independiente de tejido y se ve de forma específica en el suero cuando la comparación es con personas con colonoscopia negativa, mediante el reanálisis de datos proteómicos públicos.

- 1. Caracterizar el cambio proteico del tumor colorrectal frente a la mucosa y medir su concordancia entre dos cohortes independientes de tejido, CPTAC y Tanaka, incluidas su especificidad de órgano y su relación con el origen celular (H1, H5 y H6).
- 2. Evaluar si ese cambio aparece en el suero al comparar cáncer con colonoscopia negativa, con una regla de diferencia y de equivalencia escrita de antemano, y separar el efecto de la calidad de la muestra y del orden de adquisición con los modelos M1 a M6 (HP, H4, H7, H8, H9 y H10).
- 3. Interpretar funcionalmente el cambio del tejido con Gene Ontology y documentar un pipeline reproducible, con sello, versión v3.7 y sitio web documental, que permita rastrear cada cifra hasta su fichero (H2).

El estado de cada hipótesis y la figura están en [Hipótesis](/resultados/hipotesis).
