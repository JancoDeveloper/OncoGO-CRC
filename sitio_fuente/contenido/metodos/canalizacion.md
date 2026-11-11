---
titulo: Canalización
ruta: /metodos/canalizacion
seccion: Métodos
orden: 20
resumen: De la matriz publicada a la tabla que entra en el contraste.
---

# Canalización

## Qué entra

Entra el informe ya cuantificado del {{termino:suero}} y las matrices de tejido publicadas. Esta tesis no vuelve a adquirir los espectros. [@rai2025] [@vasaikar2019]

## Qué se hace

Se leen las proteínas, se aplican los filtros escritos en el preregistro y se dejan tres contrastes de suero separados. El diseño declarado del suero es {{dato:config/preregistro_v3.yaml#k4_suero.diseno}}. [@rai2025]

## Qué sale

La {{termino:celda}} queda con {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.n}} proteínas. [@rai2025]

## Por qué así

Un filtro escrito después de ver el resultado cambiaría el universo. El preregistro fija el universo antes. [@zhang2014]

## Qué tiene que ver con la bioinformática

La canalización de esta tesis es el paso que convierte un fichero publicado en una tabla con grupos y con una regla de ausencia. No es el programa del espectrómetro. [@tanaka2024]

## Qué no demuestra

Que una proteína entre en la celda no dice que sirva para decidir sobre un paciente. [@rai2025]

{{diagrama:D6}}
