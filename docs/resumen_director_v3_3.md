# Resumen para el director, corrida v3.3

Las cifras salen de `results/v3_3/certificado_v3_3.json` y de `suero_resumen.json`.

- No se detecta diferencia; 166 proteínas son equivalentes La matriz principal tiene 284 proteínas: 0 diferentes, 166 equivalentes y 118 indeterminadas.
- La separación del suero entre cáncer y colonoscopia negativa se explica por la calidad de muestra. Clave: `separacion_explicada_por_calidad`.
- H1: `rechaza_h1`. Concordancia 147/147. Frágil a la normalización: False.
- H5: `refuta`. H7: `refuta`. H10: `descriptiva`, 5 específicas.
- CA1 pasa. AQP8, SLC26A3 y GUCA2A fallan. Las etiquetas de tumor y mucosa no quedan confirmadas.
- ESM-2, TMbed y Papier 2024 están en `variables_faltantes` de H7.
- Ningún resultado de esta tesis sirve para decidir sobre un paciente.

No está publicado. El checklist de la sección 16 no está en verde.
