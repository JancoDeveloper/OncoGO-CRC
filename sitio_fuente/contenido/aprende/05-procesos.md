---
titulo: Procesos, no solo proteínas
ruta: /aprende/procesos
seccion: Aprende
orden: 5
resumen: Las proteínas pueden no coincidir y, aun así, apuntar al mismo proceso.
---

# Procesos, no solo proteínas

{{resumen}}
H1 esperaba poca concordancia entre proteínas y más concordancia entre vías. En el tejido, el signo de las proteínas coincidió más de lo previsto.
{{/resumen}}

{{termino:GO}} es un vocabulario: cada proteína puede llevar términos que dicen en qué proceso participa. [@ashburner2000] [@goconsortium2023]

{{termino:ORA}} pregunta si un proceso está más presente de lo esperado entre las proteínas que cambian. {{termino:GSEA}} usa toda la lista ordenada, no solo las que pasan un corte. El {{termino:NES}} dice si el proceso se enriquece hacia arriba o hacia abajo, y cuánto. [@subramanian2005]

La concordancia de signo del NES entre CPTAC y Tanaka es {{dato:results/v3_7/go_tejido.json#concordancia_nes.numerador}} de {{dato:results/v3_7/go_tejido.json#concordancia_nes.denominador}}. Eso no emite el veredicto de H1. [@vasaikar2019]

{{termino:concordancia de signo}} en proteínas, para el tumor primario frente a la mucosa, es {{dato:results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.numerador}} de {{dato:results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.denominador}}. H1 se planteó como «coincidirían poco». Coincidieron mucho, y por eso H1 no se sostiene. {{estado:H1}}

{{analogia}}
Dos cocinas pueden no comprar las mismas marcas y, aun así, estar haciendo sopa. GO mira el plato, no solo el frasco.
{{/analogia}}

{{figura:r04}}

{{especialistas}}
La red de coexpresión de H2 no se calculó: exigía al menos mil redes reconectadas y la enmienda de cierre la declara no ejecutada. {{estado:H2}}
{{/especialistas}}

## Lo que aprendiste

- GO nombra procesos, no diagnósticos.
- ORA mira un corte. GSEA mira el ranking entero.
- En este tejido, las proteínas coincidieron de signo más de lo que H1 preveía.

## Términos clave

{{termino:GO}} {{termino:ORA}} {{termino:GSEA}} {{termino:NES}} {{termino:concordancia de signo}}

## Compruébalo

<details>
<summary>¿Un NES concordante decide H1?</summary>
<p>No. El veredicto de H1 es el de las proteínas. GO se informa aparte.</p>
</details>
