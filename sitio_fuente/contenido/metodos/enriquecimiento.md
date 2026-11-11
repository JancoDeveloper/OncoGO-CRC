---
titulo: Enriquecimiento
ruta: /metodos/enriquecimiento
seccion: Métodos
orden: 23
resumen: Los términos salen de la tabla local de Gene Ontology, no de una lista escrita a posteriori.
---

# Enriquecimiento

## Qué entra

Entra la tabla `metadata/gene_ontology/go_gene_long.tsv`, que está en el manifiesto, y el ranking del contraste de tejido. [@vasaikar2019]

## Qué se hace

Se pregunta qué términos están más presentes en el cambio que en el universo de esa cohorte. El tejido y el {{termino:suero}} no se mezclan en un solo ranking. [@zhang2014]

## Qué sale

El resultado de tejido queda en los ficheros de la corrida. Esta página no copia una lista de términos. [@tanaka2024]

## Por qué así

Una lista de términos elegida después de ver los nombres no es una prueba. El universo tiene que ser el de la matriz. [@vasaikar2019]

## Qué tiene que ver con la bioinformática

El enriquecimiento resume muchas proteínas en términos. El resumen no sustituye a la proteína ni al contraste. [@rai2025]

## Qué no demuestra

Un término enriquecido en tejido no demuestra que ese proceso se vea en la sangre. [@rai2025]
