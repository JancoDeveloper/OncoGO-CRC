---
titulo: Hipótesis
ruta: /resultados/hipotesis
seccion: Resultados
orden: 35
resumen: Doce tarjetas con el enunciado de la tesis y el estado en lenguaje llano.
---

# Hipótesis

{{resumen}}
Cada hipótesis conserva el enunciado de la tesis. El chip dice el estado, no un umbral.
{{/resumen}}

{{figura:n11}}

## Tablero

- HP. Si el tumor deja una huella propia, al menos una proteína del suero debería diferir entre cáncer y colonoscopia negativa en más de 1,5 veces. {{estado:HP}} [La celda](/resultados/suero).
- H0. Auditoría de ficheros compartidos entre depósitos. Sin resultado en v3.7: se trata como un supuesto de los datos. {{estado:H0}}
- H1. Las proteínas individuales coincidirían poco en la dirección del cambio y las vías, mucho. {{estado:H1}} [Tejido](/resultados/tejido).
- H2. La estructura de las redes de coexpresión difiere entre cohortes. {{estado:H2}}
- H3. La priorización no usa pesos. {{estado:H3}} [Priorización](/resultados/priorizacion).
- H4. Se daría por falsada si al menos el 90 % de las candidatas ya fuera medible en sangre. El número es un umbral, no un estado. {{estado:H4}}
- H5. Una parte importante del cambio tumoral del colon es común a otros carcinomas. {{estado:H5}}
- H6. La concordancia entre cohortes es distinta en proteínas epiteliales y en las del estroma o inmunes. {{estado:H6}}
- H7. Lo que se sabe de una proteína antes de mirar la sangre anticipa si su cambio se verá en la misma dirección. {{estado:H7}}
- H8. Lo que sube en estadio IV está enriquecido en proteínas tumorales epiteliales. {{estado:H8}} [Por qué no se evaluó](/datos/pxd064655).
- H9. Una variante que modifica el nivel de una proteína en sangre modifica también el riesgo. {{estado:H9}}
- H10. Las proteínas de partida son específicas del cáncer colorrectal frente a otras enfermedades. {{estado:H10}}

{{especialistas}}
H1 lee {{dato:results/v3_7/k1_etiquetas.json#h1_veredicto}}. H2 lee el estado heredado del preregistro. H4 lee el umbral {{dato:config/preregistro_v3.yaml#hipotesis.h4.umbral}} y no un veredicto. H5 a H10 leen el campo estado de su fichero. H8 quedó no evaluable: el depósito nombrado no tiene matriz en el manifiesto.
{{/especialistas}}
