---
titulo: Calidad de la muestra
ruta: /aprende/calidad-de-la-muestra
seccion: Aprende
orden: 6
resumen: El tubo, el orden y el lote pueden imitar una diferencia biológica.
---

# Calidad de la muestra

{{resumen}}
En este suero, el grupo y el orden de medición van juntos. No se puede decir que la separación sea del tumor.
{{/resumen}}

La {{termino:preanalítica}} es todo lo que pasa antes de que el espectrómetro mida: cómo se sacó la sangre, cuánto tardó en centrifugar, si hemolizó, si se contaminó de piel. [@geyer2019]

La coagulación cambia el suero. Los eritrocitos, si se rompen, sueltan proteínas que no son del tumor. Las plaquetas también. Esas señales pueden diferir entre grupos sin que el tumor las haya secretado.

El {{termino:efecto de lote}} es una diferencia técnica que se alinea con el grupo. Si todas las muestras de un grupo se midieron juntas, el modelo no sabe si aprendió la biología o el día. [@leek2010]

El {{termino:orden de adquisición}} de esta corrida pone a CP, luego a CN y al final a HC. El estado de la preanalítica es {{dato:results/v3_7/preanalitica.json#estado}}. Hay fecha en {{dato:results/v3_7/preanalitica.json#n_con_fecha}} de {{dato:results/v3_7/preanalitica.json#n_crudos_validos}} crudos. [@rai2025]

{{analogia}}
Si fotografías todas las manzanas por la mañana y todas las naranjas por la tarde, un modelo puede «distinguirlas» por la luz, no por la fruta.
{{/analogia}}

{{figura:n10}}

## Lo que aprendiste

- Coagulación, eritrocitos y plaquetas cambian el suero.
- Un lote alineado con el grupo no se puede separar del grupo.
- Aquí el orden y el grupo van juntos.

## Términos clave

{{termino:preanalítica}} {{termino:efecto de lote}} {{termino:orden de adquisición}} {{termino:proteína de fase aguda}}

## Compruébalo

<details>
<summary>¿Se puede atribuir la separación del suero al tumor?</summary>
<p>No con estos datos. El orden de adquisición y el grupo están confundidos.</p>
</details>
