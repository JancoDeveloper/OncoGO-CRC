---
titulo: Hipótesis
ruta: /metodos/hipotesis
seccion: Métodos
orden: 27
resumen: Cada hipótesis tiene su regla y su estado. No se promedian.
---

# Hipótesis

## Qué entra

Entra la regla escrita en el preregistro y la tabla de la corrida que le corresponde. [@vasaikar2019]

## Qué se hace

Cada hipótesis se declara por separado. H5 usa la fuente {{dato:config/preregistro_v3.yaml#hipotesis.h5.fuente}}. H8 usa {{dato:config/preregistro_v3.yaml#hipotesis.h8.fuente}}. [@tanaka2024] [@rai2025]

## Qué sale

Los estados son H5 {{dato:results/v3_7/h5_especificidad_organo.json#estado}}, H6 {{dato:results/v3_7/h6_origen_celular.json#estado}}, H7 {{dato:results/v3_7/h7_visibilidad_en_sangre.json#estado}}, H8 {{dato:results/v3_7/h8_carga_tumoral.json#estado}}, H9 {{dato:results/v3_7/h9_wald.json#estado}} y H10 {{dato:results/v3_7/h10_especificidad.json#estado}}. [@zhang2014]

## Por qué así

Una hipótesis que no se calculó no recibe un número aproximado. H2 sigue en {{dato:config/preregistro_v3.yaml#hipotesis.h2.estado_heredado}}. [@vasaikar2019]

## Qué tiene que ver con la bioinformática

La hipótesis dice qué comparación se acepta como respuesta. El programa solo aplica esa comparación. [@rai2025]

## Qué no demuestra

Ningún estado de esta lista sirve para decidir sobre un paciente. [@rai2025]
