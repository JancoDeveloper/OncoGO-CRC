---
titulo: Cohortes y diseño
ruta: /metodos/cohortes-y-diseno
seccion: Métodos
orden: 21
resumen: Tres matrices, tres diseños, sin sumar intensidades.
---

# Cohortes y diseño

## Qué entra

El descubrimiento de tejido, la réplica {{dato:config/preregistro_v3.yaml#k1_etiquetas.accesion}} y el {{termino:suero}} {{dato:config/preregistro_v3.yaml#k4_suero.accesion}}. [@vasaikar2019] [@tanaka2024] [@rai2025]

## Qué se hace

En tejido el contraste principal es {{dato:config/preregistro_v3.yaml#k1_etiquetas.contraste_principal.nombre}}, {{dato:config/preregistro_v3.yaml#k1_etiquetas.contraste_principal.diseno}}. En suero el diseño es {{dato:config/preregistro_v3.yaml#k4_suero.diseno}}. [@tanaka2024] [@rai2025]

## Qué sale

K2 registra la cohorte de tejido como {{dato:results/v3_7/k2_independencia.json#cohorte_del_contraste_de_tejido}}. [@vasaikar2019]

## Por qué así

Si no se puede descartar que dos depósitos comparten personas, no se cuentan como dos cohortes. [@zhang2014]

## Qué tiene que ver con la bioinformática

El diseño decide qué filas se pueden restar. Emparejar por participante no es lo mismo que comparar grupos independientes. [@tanaka2024]

## Qué no demuestra

Un diseño correcto no convierte el cero del suero en una prueba de que el tumor no cambia en el tejido. [@rai2025]

{{diagrama:D7}}
