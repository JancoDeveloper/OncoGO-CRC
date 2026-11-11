---
titulo: Conclusiones
ruta: /conclusiones
seccion: Conclusiones
orden: 60
resumen: Una conclusión por objetivo, y qué no demuestra cada una.
---

# Conclusiones

El cambio del tumor primario coincide en signo entre Tanaka y CPTAC en {{dato:results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.numerador}} de {{dato:results/v3_7/k1_etiquetas.json#contrastes.primario_frente_a_mucosa.denominador}}. H5 queda en {{dato:results/v3_7/h5_especificidad_organo.json#estado}}, en un universo de {{dato:results/v3_7/h5_especificidad_organo.json#universo_principal_n}} genes. Qué no demuestra: no confirma todas las etiquetas de tumor y mucosa. [@tanaka2024] [@vasaikar2019]

En la celda del suero, {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.diferente}} de {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.n}} proteínas separan el cáncer de la colonoscopia negativa. {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.equivalente}} son equivalentes y {{dato:results/v3_7/suero_embudo.json#celda_preregistrada.indeterminada}} indeterminadas. El modelo con las proteínas de la celda da {{dato:results/v3_7/ml_suero_permutacion.json#versiones.M2_proteinas_principal.auc_oof_media}}, la calidad sola da {{dato:results/v3_7/ml_suero_permutacion.json#versiones.M3_solo_qc.auc_oof_media}} y sin la calidad el área cae a {{dato:results/v3_7/ml_suero_permutacion.json#versiones.M5_principal_residualizada_qc.auc_oof_media}} (p {{dato:results/v3_7/ml_suero_permutacion.json#versiones.M5_principal_residualizada_qc.p_permutacion}}). Qué no demuestra: no es un examen y no detecta una proteína del tumor. [@rai2025]

La calidad de la muestra, el orden de adquisición y el grupo van juntos. Hay fecha en {{dato:results/v3_7/preanalitica.json#n_con_fecha}} de {{dato:results/v3_7/preanalitica.json#n_crudos_validos}} crudos y el estado es {{dato:results/v3_7/preanalitica.json#estado}}. Qué no demuestra: no separa el efecto del grupo del efecto del lote. [@leek2010]

H6 queda en {{dato:results/v3_7/h6_origen_celular.json#estado}} y no informa. H7 queda en {{dato:results/v3_7/h7_visibilidad_en_sangre.json#estado}}: las variables del tejido no anticipan el cambio en sangre. H8 queda en {{dato:results/v3_7/h8_carga_tumoral.json#estado}}, con {{dato:results/v3_7/h8_carga_tumoral.json#n_tumorales_epiteliales}} proteínas tumorales epiteliales en el universo. Qué no demuestra: no dice que el estadio IV no lleve proteína del tumor a la sangre. [@tanaka2024]

H9 queda en {{dato:results/v3_7/h9_wald.json#estado}} y H10 queda en {{dato:results/v3_7/h10_especificidad.json#estado}}, con {{dato:results/v3_7/h10_especificidad.json#n_especificas_crc}} proteínas. Qué no demuestra: no es un hallazgo y no decide la tesis. [@rai2025]

El certificado queda en {{dato:results/v3_7/certificado.json#estado}}. Hasta una medición que aleatorice el orden y repita el contraste en otra cohorte con síntomas, no se propone un examen de sangre. [@pepe2001]
