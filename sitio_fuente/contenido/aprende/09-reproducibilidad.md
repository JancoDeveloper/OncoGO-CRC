---
titulo: Reproducibilidad
ruta: /aprende/reproducibilidad
seccion: Aprende
orden: 9
resumen: Lo que se escribió antes de mirar, y el sello de lo que sí se corrió.
---

# Reproducibilidad

{{resumen}}
Un preregistro anota la regla antes del resultado. Una enmienda anota el cambio y por qué. El sello dice qué código produjo los números.
{{/resumen}}

El {{termino:preregistro}} separa lo que se iba a confirmar de lo que se explora después. La {{termino:enmienda}} no esconde el cambio: lo fecha y dice qué no se movió. [@nosek2018]

El {{termino:multiverso}} repite el análisis con otras especificaciones razonables, para ver si la conclusión depende de una sola elección. [@steegen2016]

El {{termino:sello}} guarda el resumen criptográfico de la corrida. El {{termino:certificado}} dice si esa corrida se puede publicar. Su estado es {{dato:results/v3_7/certificado.json#estado}}. [@zhang2014]

{{analogia}}
El preregistro es la receta escrita antes de probar la sopa. La enmienda es la nota al margen cuando cambiaste un ingrediente. El sello es la foto de la cocina el día que serviste.
{{/analogia}}

{{especialistas}}
La repetición D/E comprueba que dos corridas con el mismo código dejan el mismo sello. El commit que el certificado guarda es {{dato:results/v3_7/certificado.json#sello.commit}}. El árbol limpio es {{dato:results/v3_7/certificado.json#sello.arbol_limpio}}. La semilla es {{dato:results/v3_7/certificado.json#sello.semilla}}. Lo que la enmienda declara no ejecutado —vías del suero, redes de coexpresión, neoantígenos— es una desviación, no un resultado a favor.
{{/especialistas}}

## Lo que aprendiste

- La regla se escribe antes de ver el número.
- Una enmienda se declara, no se disimula.
- El certificado dice el estado de esta corrida, y cada cifra del sitio apunta a su fichero.

## Términos clave

{{termino:preregistro}} {{termino:enmienda}} {{termino:multiverso}} {{termino:sello}} {{termino:certificado}}

## Compruébalo

<details>
<summary>¿Una cifra sin fichero puede publicarse en este sitio?</summary>
<p>No. Si el atajo no encuentra la fuente, la construcción se detiene.</p>
</details>

<details>
<summary>¿Las redes de coexpresión quedaron calculadas?</summary>
<p>No. La enmienda de cierre las declara no ejecutadas. H2 queda no evaluada.</p>
</details>
