---
titulo: Discusión
ruta: /discusion
seccion: Discusion
orden: 50
resumen: Qué se puede leer de la corrida y qué queda confundido con el lote.
---

# Discusión

## El diseño en bloques

Los grupos se midieron en bloques. El orden de {{termino:CP}} va de {{dato:results/v3_7/preanalitica.json#posicion_por_grupo.CP.orden_min}} a {{dato:results/v3_7/preanalitica.json#posicion_por_grupo.CP.orden_max}}, el de {{termino:CN}} de {{dato:results/v3_7/preanalitica.json#posicion_por_grupo.CN.orden_min}} a {{dato:results/v3_7/preanalitica.json#posicion_por_grupo.CN.orden_max}}, y el de {{termino:HC}} de {{dato:results/v3_7/preanalitica.json#posicion_por_grupo.HC.orden_min}} a {{dato:results/v3_7/preanalitica.json#posicion_por_grupo.HC.orden_max}}. Hay fecha en {{dato:results/v3_7/preanalitica.json#n_con_fecha}} de {{dato:results/v3_7/preanalitica.json#n_crudos_validos}} crudos. El estado de la preanalítica es {{dato:results/v3_7/preanalitica.json#estado}}. [@geyer2019] [@leek2010]

Con estos datos no se separa el efecto del grupo del efecto del lote de adquisición. Eso no crea una proteína del tumor: refuerza que la separación del suero no se atribuye al tumor. [@rai2025]

## Coagulación: tubo o tumor

La coagulación difiere entre cáncer y colonoscopia negativa, con p {{dato:results/v3_7/suero_calidad.json#confusion_preanalitica.0.p_kruskal}}. La piel, con p {{dato:results/v3_7/suero_calidad.json#confusion_preanalitica.1.p_kruskal}}, apunta a manipulación. El cáncer puede activar la coagulación, y el tubo y el tiempo hasta centrifugar también la cambian. Con estos datos no se separan. [@geyer2019]

## El resultado negativo

En la celda, {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.diferente}} de {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.n}} proteínas separan el cáncer de la colonoscopia negativa. {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.equivalente}} son equivalentes y {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.indeterminada}} quedan indeterminadas. Con este tamaño, no detectar una diferencia no equivale a demostrar que no existe. [@rai2025]

## Coherencia y lo que no se alcanza

SAA1, LBP y LRG1 forman parte de la respuesta de fase aguda. Esa lectura es plausible y no se separa del lote, porque HC se midió al final. [@gabay1999] CEACAM5 se usa en el seguimiento, no como cribado. [@locker2006] H8 queda en {{dato:results/v3_7/h8_carga_tumoral.json#estado}}, con {{dato:results/v3_7/h8_carga_tumoral.json#n_tumorales_epiteliales}} proteínas tumorales epiteliales en un universo de {{dato:results/v3_7/h8_carga_tumoral.json#n_universo}}. Un DIA sin depleción no alcanza esas proteínas. [@hermoso2026]

## Otros estudios y limitaciones

Otros estudios preguntan cosas distintas y usan diseños distintos. No se comparan áreas bajo la curva entre diseños. H7 queda en {{dato:results/v3_7/h7_visibilidad_en_sangre.json#estado}} dentro de sus propias reglas. [@rai2025]

El suero es de un solo centro y de grupos pequeños. No hay tejido y sangre de la misma persona. H6 queda en {{dato:results/v3_7/h6_origen_celular.json#estado}}. H9 queda en {{dato:results/v3_7/h9_wald.json#estado}} y no decide la tesis. [@vasaikar2019]

## Qué diseño lo evitaría

Aleatorizar el orden de adquisición, intercalar los grupos y repetir el contraste en otra cohorte de personas con síntomas. Sigue siendo la fase 1 de ese esquema, y no es un examen. [@pepe2001]

La enmienda declara como desviación, no como resultado, el catálogo de candidatas, la cláusula de vías de H1, los neoantígenos, la deriva de PRTC, las redes de coexpresión, la unión de HC por XIC y las vías del suero. Esa casilla está en {{dato:results/v3_7/certificado.json#casillas.6.estado}}. [@rai2025]
