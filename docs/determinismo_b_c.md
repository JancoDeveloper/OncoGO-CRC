# Determinismo B frente a C

B: `/Volumes/JANCODATA/OncoGO-CRC/results/v3_7`. C: `/Volumes/JANCODATA/OncoGO-CRC/results/archivo/v3_7_final_c`.

Las dos salidas salen del commit `b7f8b8e844f6c484b72e19fecf5fb7d07f53c2ab`. En el sello coinciden el commit, `arbol_limpio`, el preregistro, las enmiendas, el manifiesto, los paquetes, la semilla y la etapa. Solo cambia `fecha_utc`.

Diferencias admitidas: certificado.json, corrida.json, gsea_tejido.tsv, multiverso.json, multiverso_v3_4.json, procedencia.json, tests.txt.

El secundario de H7 es 0.4773843370379944 en B y 0.4882398552597571 en C. El veredicto sigue en refuta y el AUC principal no cambia. En `h7_tabla_proteinas.tsv` los valores presentes de ESM-2 coinciden; cambia en qué filas falta el valor. Ese secundario no entra en el veredicto, y esta comparación no lo admite como diferencia menor.

Diferencias que impiden la etiqueta: h7_tabla_proteinas.tsv, h7_visibilidad_en_sangre.json.

No se crea la etiqueta `v3.7-final`.
