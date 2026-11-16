---
titulo: Cómo se mide una proteína
ruta: /aprende/como-se-mide
seccion: Aprende
orden: 3
resumen: De la digestión al número, y por qué TMT, LFQ y DIA no se pegan.
---

# Cómo se mide una proteína

{{resumen}}
La espectrometría de masas no ve la proteína entera: ve péptidos, les asigna una identidad y devuelve una intensidad.
{{/resumen}}

{{paso 1}} La digestión con tripsina corta la proteína en {{termino:péptido}}s.

{{paso 2}} La cromatografía líquida los separa en el tiempo.

{{paso 3}} El espectrómetro mide la masa y fragmenta. Eso es {{termino:LC-MS/MS}}.

{{paso 4}} Un programa identifica el péptido y, con varios péptidos, el {{termino:grupo de proteínas}}.

{{paso 5}} La cuantificación convierte esa identificación en un número comparable dentro de la misma corrida.

## Tres maneras de cuantificar

{{termino:TMT}} pone una etiqueta química. Varias muestras viajan juntas y se distinguen por la etiqueta. CPTAC midió el tejido así. [@vasaikar2019]

{{termino:LFQ}} no pone etiqueta. Compara la intensidad del péptido entre corridas. Tanaka midió el tejido así, y directLFQ es una forma de hacer esa comparación. [@cox2014] [@ammar2023] [@tanaka2024]

{{termino:DIA}} fragmenta ventanas de precursores en lugar de elegir un ion por su intensidad. {{termino:DIA-NN}} lee ese informe. El suero de esta tesis sale de ahí. [@demichev2020] [@rai2025]

{{analogia}}
TMT es marcar las maletas con una cinta de color y pesarlas en el mismo viaje. LFQ es pesar cada maleta en viajes distintos y alinear las básculas después. DIA es abrir todas las maletas de una ventana, no solo la más pesada. Las tres básculas no dan un número que se pueda sumar.
{{/analogia}}

{{diagrama:D5}}

{{especialistas}}
No se aplica ComBat para pegar matrices: un efecto de lote entre plataformas no es un lote dentro de un mismo experimento. ThermoRawFileParser convierte crudos cuando hace falta leerlos; la matriz principal del suero parte del informe ya cuantificado. [@hulstaert2020] [@leek2010]
{{/especialistas}}

## Lo que aprendiste

- Se miden péptidos, no la proteína fotografiada.
- TMT, LFQ y DIA responden preguntas distintas.
- Por eso las matrices de esta tesis no se fusionan.

## Términos clave

{{termino:péptido}} {{termino:LC-MS/MS}} {{termino:TMT}} {{termino:LFQ}} {{termino:DIA}} {{termino:DIA-NN}} {{termino:grupo de proteínas}}

## Compruébalo

<details>
<summary>¿Por qué no se suma el tejido TMT con el suero DIA?</summary>
<p>Porque la etiqueta, la escala y las personas no son las mismas. Fusionarlas crearía un número sin significado.</p>
</details>
