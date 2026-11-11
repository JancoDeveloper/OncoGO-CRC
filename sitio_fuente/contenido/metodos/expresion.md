---
titulo: Expresión diferencial
ruta: /metodos/expresion-diferencial
seccion: Métodos
orden: 22
resumen: Diferente, equivalente e indeterminada son tres clases, no un solo umbral.
---

# Expresión diferencial

## Qué entra

Entra la {{termino:celda}} del {{termino:suero}} y los contrastes de tejido ya definidos. [@rai2025] [@vasaikar2019]

## Qué se hace

Cada proteína recibe una clase. Diferente y equivalente no se colapsan en «no significativo». El margen TOST declarado, en log2, va de {{dato:config/preregistro_v3.yaml#k4_suero.margenes_tost_log2.0}} a {{dato:config/preregistro_v3.yaml#k4_suero.margenes_tost_log2.1}}. [@rai2025]

## Qué sale

Frente a {{termino:CN}} hay {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.diferente}} diferentes, {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.equivalente}} equivalentes y {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.indeterminada}} indeterminadas. [@rai2025]

## Por qué así

Un valor no significativo no demuestra que el efecto sea pequeño. La clase equivalente exige el margen escrito antes. [@zhang2014]

## Qué tiene que ver con la bioinformática

La prueba se calcula sobre la matriz, con el error estándar y los grados de libertad de cada proteína. No se elige la clase a mano. [@tanaka2024]

## Qué no demuestra

La clase equivalente no dice que la proteína esté ausente, ni que el resultado sirva para un paciente. [@rai2025]
