# Qué cambia de v1 a v2, número por número

Esta tabla sale de `certification/changelog_numeros.tsv`. Cada fila pone el valor de v1 al lado del de v2, con el motivo y el código de la corrección que lo movió. Si una fila dice que no cambia, el número es el mismo y lo que cambia es lo que se publica junto a él.

> **Ningún resultado de esta tesis sirve para decidir sobre un paciente.**

| Número | v1 | v2 | ¿Cambia? | Corrección | Motivo |
| --- | --- | --- | --- | --- | --- |
| `R1_concordancia_de_signo` | 0,8739 | 0,8739 | no | C15 | Mismo cálculo. Lo que cambia es que ahora se publica junto a su cobertura: el denominador es solo el 17 % de los hallazgos. |
| `R1_cobertura` | no disponible | 0,1730 | no aplica | C15 | v1 no publicaba este número. Sin él, R1 parece una concordancia general y es una concordancia condicionada a que la réplica mida la proteína. |
| `R2_detectabilidad` | 0,3285 | 0,3285 | sí | C3 | Mismo valor. Lo que se corrige es que v1 tenía tres denominadores distintos circulando (688, 702 y 244) para cosas que parecían la misma. |
| `R5_p_empirico` | 0,0099 | 0,000100 | sí | C1 | 100 permutaciones pasan a 10 000. En las dos el nulo nunca alcanza lo observado, así que el número es el mínimo alcanzable y se declara como tope. |
| `R4_composicion` | 0,4462 | 0,4477 | sí | C7 | v1 usaba las 282 proteínas de la matriz extracelular como índice y a la vez las contrastaba. Ahora cada proteína se excluye de su propio índice. |
| `concordancia_de_via` | 1,0000 | 0,9216 | sí | C4 | v1 daba 1,0 con cuatro términos. Ahora el mínimo es diez términos, la prueba es competitiva y se repite en tres colecciones. |
| `n_terminos_de_via` | 4,0000 | 51,0000 | sí | C4 | El número de términos es parte del resultado, no un detalle: con cuatro, el intervalo iba de 0,398 a 1,0. |
| `H1_veredicto` | rechazada_por_clausula_proteina | rechaza_h1 | sí | C4 | Mismo veredicto, pero ahora sale de una tabla de verdad explícita con el estado no_evaluable incluido. |
| `polipo_magnitud_intermedia` | 0,9195 | 0,5625 | sí | C5 | v1 seleccionaba las proteínas en el brazo tumor y medía su magnitud en ese mismo brazo. Con selección independiente en CPTAC el efecto desaparece. |
| `polipo_esperado_por_regresion_a_la_media` | no disponible | 0,7516 | no aplica | C5 | v1 no tenía este nulo. Dice cuánta de la asimetría la produce la propia selección. |
| `concordancia_tejido_plasma` | 0,4884 | 0,5397 | sí | C8 | v1 medía sobre todas las proteínas comunes. Ahora se mide solo entre las que el método alcanza a detectar en suero. |

De los 11 números del changelog, 8 cambian.

## Lo que no cambia

v1 no se toca. Su manifiesto de hashes está en `results/v1_manifest.sha256` y `make manifest-verify` falla si algún fichero de v1 o alguna matriz de entrada se ha movido. Las matrices crudas se leen en modo solo lectura y su manifiesto está en `results/v2/manifest_inputs.sha256`.

## Enmiendas

El preregistro no se edita. Cuando un umbral o un método cambió después de ver los datos de v1, queda escrito como enmienda en `docs/enmiendas_v2.md`, con la fecha y la justificación.