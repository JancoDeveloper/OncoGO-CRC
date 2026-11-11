---
titulo: Proteómica
ruta: /metodos/proteomica
seccion: Métodos
orden: 10
resumen: Qué se midió en el tejido y en el suero, y qué hace con eso el cálculo.
---

# Proteómica

## Qué entra

Entra el {{termino:suero}} medido por {{termino:DIA}} en tres grupos, {{termino:CP}}, {{termino:CN}} y {{termino:HC}}. [@rai2025] Entra también el {{termino:proteoma}} de tejido de descubrimiento y el de la réplica. [@vasaikar2019] [@tanaka2024]

## Qué se hace

La digestión, el cromatógrafo y el espectrómetro son los del artículo de cada cohorte. Esta tesis no repitió el experimento: lee las matrices publicadas y les aplica la regla escrita antes de mirar. [@rai2025] [@zhang2014]

## Qué sale

La {{termino:celda}} conserva {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.n}} proteínas. Frente a {{termino:CN}} hay {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.diferente}} diferentes. [@rai2025]

## Por qué así

Las intensidades de una plataforma no se suman a las de otra. El tejido de descubrimiento y el suero no son la misma medición ni las mismas personas. [@zhang2014] [@rai2025]

## Qué tiene que ver con la bioinformática

El cálculo empieza cuando la matriz ya existe: identificar las proteínas, separar los grupos y contar con la regla fijada. No diseña el ensayo de laboratorio. [@vasaikar2019]

## Qué no demuestra

El resultado del suero no demuestra que el cambio del tumor esté ausente en el tejido, ni que el cero sirva para descartar una enfermedad en una persona. [@tanaka2024] [@rai2025]
