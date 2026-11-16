---
titulo: Comparar grupos sin engañarse
ruta: /aprende/comparar-grupos
seccion: Aprende
orden: 4
resumen: Un cambio pequeño, un p pequeño y una equivalencia no son la misma frase.
---

# Comparar grupos sin engañarse

{{resumen}}
Esta tesis no pregunta solo si el cambio es distinto de cero. Pregunta si supera 1,5 veces, si cabe dentro de ese margen, o si no se puede decidir.
{{/resumen}}

El {{termino:log2FC}} dice cuántas veces cambia una proteína, en escala de logaritmo en base 2. Un margen de 1,5 veces es un aumento de la mitad o una bajada a dos tercios. En log2, ese margen es 0,585.

El p-valor dice qué tan raro sería el dato si no hubiera cambio. Con miles de proteínas, un p pequeño aparece por azar. El {{termino:FDR}} y el procedimiento de {{termino:Benjamini-Hochberg}} ponen un techo a esa tasa. El {{termino:q-valor}} es el número ya ajustado.

{{termino:limma}} estima el cambio con una {{termino:t moderada}}, que toma prestada información de las otras proteínas para no fiarse de una varianza inestable. [@smyth2004] [@ritchie2015]

{{termino:TREAT}} exige que el cambio supere el margen, no solo que sea distinto de cero. {{termino:TOST}} declara {{termino:equivalencia}} cuando el intervalo cabe dentro del margen. Si no ocurre ninguna de las dos, la proteína queda indeterminada: no ver un cambio, con esta potencia, no demuestra que no exista. [@mccarthy2009] [@lakens2017]

{{juguete:clases}}

{{analogia}}
Un p pequeño es oír un ruido. TREAT pregunta si el ruido pasa de un volumen acordado. TOST pregunta si estamos seguros de que el volumen se quedó por debajo. Si el oído no alcanza, la respuesta honesta es «no se sabe», no «hay silencio».
{{/analogia}}

En la celda, frente a CN, hay {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.diferente}} diferentes, {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.equivalente}} equivalentes y {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.indeterminada}} indeterminadas, de {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.n}}. [@rai2025]

{{comparacion}}

## Lo que aprendiste

- 1,5 veces no es «cualquier cambio».
- FDR corrige el hecho de mirar miles de proteínas.
- Indeterminada no es equivalente.

## Términos clave

{{termino:log2FC}} {{termino:FDR}} {{termino:Benjamini-Hochberg}} {{termino:limma}} {{termino:TREAT}} {{termino:TOST}} {{termino:equivalencia}}

## Compruébalo

<details>
<summary>Si el intervalo cruza el margen, ¿la clase de juguete es equivalente?</summary>
<p>No. Equivalente exige que todo el intervalo quepa dentro del margen. Si lo cruza, queda indeterminada.</p>
</details>
