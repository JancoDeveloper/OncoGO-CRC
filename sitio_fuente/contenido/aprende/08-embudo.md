---
titulo: El embudo
ruta: /aprende/el-embudo
seccion: Aprende
orden: 8
resumen: Seis filtros sin pesos, y nadie pasa el sexto.
---

# El embudo

{{resumen}}
La priorización no suma puntos. Cada proteína se queda en el último paso que supera, y el paso 6 se queda en cero.
{{/resumen}}

{{estado:H3}}

{{paso 1}} Diferencial certificada en el tejido de CPTAC. Quedan {{dato:results/v3_7/priorizacion_embudo.json#pasos.0.n}}.

{{paso 2}} Medida en Tanaka, tumor primario frente a mucosa, con el mismo signo. Quedan {{dato:results/v3_7/priorizacion_embudo.json#pasos.1.n}}.

{{paso 3}} La clase de especificidad de órgano se informa y no elimina. Siguen {{dato:results/v3_7/priorizacion_embudo.json#pasos.2.n}}.

{{paso 4}} El atlas de proteínas la marca como secretada a sangre o con concentración conocida. Quedan {{dato:results/v3_7/priorizacion_embudo.json#pasos.3.n}}.

{{paso 5}} Está en la matriz principal del suero. Quedan {{dato:results/v3_7/priorizacion_embudo.json#pasos.4.n}}. Son hipótesis para otra cohorte, no un panel para la consulta.

{{paso 6}} En el suero, frente a CN, es diferente y con el mismo signo que el tejido. Quedan {{dato:results/v3_7/priorizacion_embudo.json#pasos.5.n}}.

{{paso 7}} El atlas de sangre la marca específica de este cáncer. También {{dato:results/v3_7/priorizacion_embudo.json#pasos.6.n}}, porque el paso anterior ya estaba vacío.

Nadie pasa el paso 6: {{dato:results/v3_7/priorizacion_embudo.json#nadie_pasa_el_paso_6}}. [@rai2025]

{{embudo:priorizacion}}

{{analogia}}
Es una aduana de puertas en fila, no una nota media. Quien no pasa una puerta no compensa con las anteriores.
{{/analogia}}

## Lo que aprendiste

- No hay pesos ni un puntaje único.
- Veinte proteínas llegan al paso 5.
- Ninguna pasa el paso 6.

## Términos clave

{{termino:HPA}} {{termino:CN}} {{termino:concordancia de signo}} {{termino:preregistro}}

## Compruébalo

<details>
<summary>¿Las 20 del paso 5 son un panel listo para usar?</summary>
<p>No. Son hipótesis para cohortes futuras. El paso que pedía diferencia en el suero, con el mismo signo que el tejido, quedó en cero.</p>
</details>
