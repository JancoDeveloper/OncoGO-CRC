---
titulo: Pregunta e hipótesis
ruta: /investigacion/pregunta-e-hipotesis
seccion: Investigación
orden: 2
resumen: Qué se preguntó antes de mirar, y el estado con el que cerró cada hipótesis.
---

# Pregunta e hipótesis

La pregunta es qué parte del cambio del tumor se deja ver en la sangre cuando el contraste es otra enfermedad intestinal. [@rai2025] [@garciafoncillas2026]

## Lo que se fijó antes

H1 conserva el umbral {{dato:config/preregistro_v3.yaml#hipotesis.h1.umbral_concordancia_proteina}}. La nota del preregistro dice: {{dato:config/preregistro_v3.yaml#hipotesis.h1.nota}} [@vasaikar2019] [@tanaka2024]

H2 quedó con el estado heredado {{dato:config/preregistro_v3.yaml#hipotesis.h2.estado_heredado}}. [@zhang2014]

## Cómo cerró cada una

H5, especificidad de órgano, cerró en {{dato:results/v3_7/h5_especificidad_organo.json#estado}}. [@vasaikar2019]

H6, origen celular, cerró en {{dato:results/v3_7/h6_origen_celular.json#estado}}. [@tanaka2024]

H7, visibilidad en sangre, cerró en {{dato:results/v3_7/h7_visibilidad_en_sangre.json#estado}}. El AUC principal de ese modelo es {{dato:results/v3_7/h7_visibilidad_en_sangre.json#modelo_secundario.auc_principal}}. [@rai2025]

H8, carga tumoral, cerró en {{dato:results/v3_7/h8_carga_tumoral.json#estado}}. [@rai2025]

H9 cerró en {{dato:results/v3_7/h9_wald.json#estado}}. H10 cerró en {{dato:results/v3_7/h10_especificidad.json#estado}}. [@zhang2014]

Las hipótesis no se promedian entre sí. [@vasaikar2019]
