---
titulo: Aprendizaje automático
ruta: /metodos/aprendizaje-automatico
seccion: Métodos
orden: 24
resumen: El modelo de visibilidad se juzga por su AUC, y el veredicto ya está escrito.
---

# Aprendizaje automático

## Qué entra

Entran variables conocidas antes de mirar la sangre. La predicción de H7 es la que el preregistro escribe. [@rai2025]

## Qué se hace

Un modelo intenta anticipar qué proteínas del tumor cambian en la sangre en la misma dirección. El veredicto no usa el modelo secundario. [@rai2025]

## Qué sale

H7 cerró en {{dato:results/v3_7/h7_visibilidad_en_sangre.json#estado}}. El AUC principal es {{dato:results/v3_7/h7_visibilidad_en_sangre.json#modelo_secundario.auc_principal}} y el secundario es {{dato:results/v3_7/h7_visibilidad_en_sangre.json#modelo_secundario.auc_secundario}}. [@rai2025]

## Por qué así

El umbral de sostén está en el preregistro: el límite inferior del AUC por encima de {{dato:config/preregistro_v3.yaml#hipotesis.h7.umbral_limite_inferior}}, en al menos {{dato:config/preregistro_v3.yaml#hipotesis.h7.recursos_minimos}} recursos. [@vasaikar2019]

## Qué tiene que ver con la bioinformática

El modelo comprime variables de anotación. No crea una medición nueva de la sangre. [@tanaka2024]

## Qué no demuestra

Un AUC cercano a un medio no es una prueba diagnóstica. El propio preregistro lo marca: {{dato:config/preregistro_v3.yaml#hipotesis.h7.no_es_prueba_diagnostica}}. [@rai2025]
