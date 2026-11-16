---
titulo: Inicio
ruta: /
seccion: Inicio
orden: 0
resumen: En 284 proteínas del suero, ninguna distinguió el cáncer de otra enfermedad intestinal.
---

# ¿Qué parte del cambio proteico del tumor colorrectal puede verse en el suero cuando el contraste es otra enfermedad intestinal y no una persona sana?

{{resumen}}
En {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.n}} proteínas del suero, {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.diferente}} distinguieron el cáncer de otra enfermedad intestinal.
{{/resumen}}

<p class="cta"><a class="boton" href="/aprende/cancer-colorrectal">Empieza el recorrido</a> <a class="boton boton-secundario" href="/resultados">Ver los resultados</a></p>

## Para qué sirve y para qué no

<div class="etica">
<div>
<p class="resumen-etiqueta">Sí</p>
<p>Sirve para generar hipótesis que otra cohorte, con síntomas y con el orden de medición mezclado, pueda volver a mirar.</p>
</div>
<div>
<p class="resumen-etiqueta">No</p>
<p>Ningún resultado de esta tesis sirve para decidir sobre un paciente.</p>
</div>
</div>

## Tres puertas de entrada

<div class="tarjetas">
<a class="tarjeta" href="/aprende/cancer-colorrectal"><h2>Biología</h2><p>Qué es el cáncer colorrectal, qué es una proteína y por qué el contraste no es una persona sana.</p></a>
<a class="tarjeta" href="/aprende/como-se-mide"><h2>Bioinformática</h2><p>Cómo se mide una proteína, cómo se comparan grupos y por qué las matrices no se fusionan.</p></a>
<a class="tarjeta" href="/aprende/que-hace-la-ia"><h2>Inteligencia artificial</h2><p>Qué hacen los modelos M1 a M6 y por qué un área bajo la curva interna no es un examen.</p></a>
</div>

## Cuatro cifras

<div class="cifras">
{{cifra:results/v3_7/priorizacion_embudo.json#pasos.0.n|Certificadas en tejido}}
{{cifra:results/v3_7/priorizacion_embudo.json#pasos.4.n|Llegan al paso 5}}
{{cifra:results/v3_7/priorizacion_embudo.json#pasos.5.n|Pasan el paso 6}}
{{cifra:results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.numerador|Mismo signo, tumor primario}}
{{cifra:results/v3_7/ml_suero_permutacion.json#versiones.M2_proteinas_principal.auc_oof_media|AUC de M2, proteínas}}
{{cifra:results/v3_7/ml_suero_permutacion.json#versiones.M3_solo_qc.auc_oof_media|AUC de M3, solo calidad}}
{{cifra:results/v3_7/suero_embudo.json#celda_preregistrada.n|Proteínas de la celda}}
{{cifra:results/v3_7/suero_embudo.json#celda_preregistrada.diferente|Diferentes frente a CN}}
{{cifra:results/v3_7/suero_embudo.json#celda_preregistrada.equivalente|Equivalentes}}
{{cifra:results/v3_7/suero_embudo.json#celda_preregistrada.indeterminada|Indeterminadas}}
</div>

## El mapa

{{diagrama:D1}}

## ¿Cómo se clasificaron las proteínas frente a CN?

{{termino:CN}} son personas con síntomas y otra enfermedad intestinal. {{termino:HC}} es otro grupo y no decide la pregunta. [@rai2025]

{{comparacion}}

## El embudo

Nadie pasa el paso 6: {{dato:results/v3_7/priorizacion_embudo.json#nadie_pasa_el_paso_6}}. [@rai2025]

{{embudo:priorizacion}}

## Qué más hay en el tejido

El signo del tumor primario coincide en {{dato:results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.numerador}} de {{dato:results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.denominador}}. [@vasaikar2019] [@tanaka2024]

{{matrices}}

{{explorar}}
