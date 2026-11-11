---
titulo: Inicio
ruta: /
seccion: Inicio
orden: 0
resumen: La pregunta de la tesis y el resultado de la celda del suero.
---

# ¿Qué parte del cambio del tumor se deja ver en la sangre?

{{llamadas}}

El tejido publicado ya describe el tumor y la mucosa. [@vasaikar2019] [@tanaka2024] Esta tesis pregunta cuánto de ese cambio se deja ver en la sangre cuando el contraste es otra enfermedad intestinal, y no fusiona las matrices. [@rai2025]

## ¿Qué problema se midió?

El cáncer colorrectal se describe, en el tejido, por el proteoma del tumor y el de la mucosa. [@vasaikar2019] [@zhang2014] El {{termino:suero}} es otra matriz y otras personas. El depósito que lee esta tesis es {{dato:config/preregistro_v3.yaml#k4_suero.accesion}}. [@rai2025]

## El colon, como contexto

El dibujo nombra la anatomía. No es un resultado y no es un mapa de pacientes. La tesis no dice en qué segmento estaban las muestras.

{{anatomia:colon}}

## ¿Cuál es la pregunta?

La pregunta es qué parte del cambio del tumor se deja ver en la sangre cuando el contraste es otra enfermedad intestinal. [@rai2025] [@garciafoncillas2026]

H7, visibilidad en sangre, cerró en {{dato:results/v3_7/h7_visibilidad_en_sangre.json#estado}}. [@rai2025] Las demás hipótesis y su estado están en la página de la pregunta. [@vasaikar2019]

## ¿Qué muestra cada matriz?

{{matrices}}

El cambio del tumor primario coincide en signo entre Tanaka y CPTAC en {{dato:results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.numerador}} de {{dato:results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.denominador}}. [@tanaka2024] [@vasaikar2019]

La clase diferente del suero es la de la clasificación de más abajo. Qué no demuestra: no es un examen y no detecta una proteína del tumor. [@rai2025]

## ¿Cómo se recorrió el análisis?

{{diagrama:D1}}

El embudo de priorización sale de `priorizacion_embudo.json`. El ancho de cada barra es el cociente entre el n de ese paso y el n mayor de la lista. Nadie pasa el paso 6: {{dato:results/v3_7/priorizacion_embudo.json#nadie_pasa_el_paso_6}}. [@rai2025] [@tanaka2024]

{{embudo:priorizacion}}

## ¿Cómo se clasificaron las proteínas frente a CN?

{{comparacion}}

{{termino:CN}} son pacientes sintomáticos con otra enfermedad intestinal. {{termino:HC}} es otro grupo, medido aparte, y no es el contraste que decide la pregunta. [@rai2025]

## ¿Dónde se consulta una proteína?

Una ficha se busca por identificador o por gen en el campo de arriba, o en el listado. La función es el texto de UniProt. Los contrastes de suero y de tejido son resultados calculados. La etiqueta de candidata, cuando sale, es una hipótesis para cohortes futuras.

{{explorar}}

## ¿Qué no demuestra este resultado?

El cero de la clase diferente no dice que el tumor no cambie en el tejido. Las {{dato:results/v3_7/priorizacion_embudo.json#pasos.4.n}} del paso 5 quedan como hipótesis para cohortes futuras. No hay una cohorte nueva medida. Las redes no se calcularon. [@vasaikar2019] [@tanaka2024] [@zhang2014]
