# Banco de preguntas, versión del 6 de noviembre de 2026

Este banco describe la corrida v3.2. La corrida vigente es v3.3 y está en `docs/sustentacion_v3_3.md`.

Cada respuesta usa una tabla. Si la hipótesis no está calculada, el estado es `no_evaluable`. No hay número aproximado. Ningún resultado de aquí es una prueba diagnóstica. Como mucho está en la fase 1 de Pepe et al. 2001: exploración en muestras ya recogidas.

## Lo que ya se puede defender

1. **¿Qué aporta esta tesis que no se supiera?**
La corrección de la etiqueta, el filtro del suero y la regla de equivalencia. El primario frente a la mucosa coincide en signo en 89 de 90, tras excluir 237 filas y pedir 8 pares. La metástasis frente a la mucosa coincide en 69 de 82. En el suero oficial, de 205 proteínas, 0 son diferentes, 105 equivalentes y 100 indeterminadas entre cáncer y colonoscopia negativa. H5 quedó `no_disponible`. H10 cuenta 0 específicas de 1 165 y no tiene veredicto. Fuente: `results/v3_2/k1_etiquetas.json`, `results/v3_2/suero_resumen.json`, `results/v3_2/hipotesis_v3_2.tsv`.

2. **¿Por qué no fusionaste las cohortes?**
Cada contraste ocurre dentro de su cohorte. TMT, LFQ y DIA no comparten escala ni personas. Unirlas contaría dos veces al mismo laboratorio si los depósitos se solapan, y mezclaría un hígado metastásico con una mucosa. El preregistro lo prohíbe. Fuente: `config/preregistro_v3.yaml`, sección de inferencia.

3. **Con 20 sueros por grupo, ¿cómo dices que no hay señal?**
Solo para las 105 proteínas de clase equivalente en el contraste oficial, cáncer frente a colonoscopia negativa, sobre 205 proteínas. En esas, las dos pruebas unilaterales demuestran un cambio menor que 1,5 veces. Las 100 indeterminadas no demuestran ausencia. Las 0 diferentes son el corte de treat y BH, que es otra frase. La lectura anterior, de 876 proteínas, era 145 equivalentes y 731 indeterminadas, y queda en `results/v3/`. Fuente: `results/v3_2/suero_resumen.json`.

4. **¿Por qué la réplica de tejido no es «cualquier tumor contra cualquier normal»?**
Esa pregunta es H5. En v3.2 quedó `no_disponible`: el carcinoma de endometrio no tiene 8 pares con el mismo identificador de muestra. El mínimo no se bajó. La réplica que sí está hecha compara el primario de colon con la mucosa del mismo paciente, y aparte nombra la metástasis hepática. No se presenta como específica de colon. Fuente: `results/v3_2/h5_especificidad_organo.json`, `results/v3_2/hipotesis_v3_2.tsv`.

5. **¿Qué pasó con el pólipo?**
La columna P de PXD034575 es el tumor primario. Una versión anterior la guardó como pólipo y usó la metástasis (M) como tumor. El fichero viejo no se editó. La leyenda nueva está en `metadata/sample_metadata/pride_pxd034575_v3.tsv`. El contraste principal es primario frente a mucosa. Fuente: `docs/auditoria_v3.md`, sección 2.1; `results/v3/k1_comparacion_contrastes.tsv`.

6. **¿PXD019504 y PXD031705 son independientes de Tanaka?**
No está recerrado en v3. El informe de v2 ya dijo que PXD031705 comparte ficheros con PXD019504, así que solo uno entró. El solapamiento de personas con PXD034575 no se ha vuelto a medir. Por eso R3, el 1 − I² de tres cohortes, queda `no_evaluable`. El 0,094 de v2 se conserva con su nombre y no se reutiliza. Fuente: `results/v3/k1_etiquetas.json`, campo `r3`; `docs/reportes/informe_v2.md`, sección 7.

7. **¿Por qué no usaste R ni limma?**
No hay R en el equipo. treat, el meta-análisis REML, camera y Reitsma están reimplementados en Python, con tests contra casos de resultado conocido. No se garantiza el último dígito frente a R. Fuente: `docs/enmiendas_v2.md`, enmienda E2.

8. **¿Qué es treat y por qué 0,585?**
treat pregunta si el cambio supera un umbral, no si es distinto de cero. El umbral es log2(1,5) = 0,585: un cambio de al menos una vez y media. Va con la tasa de descubrimientos falsos de Benjamini–Hochberg en 0,05. No se movió después de ver los datos. Fuente: `config/preregistro_v3.yaml`; McCarthy y Smyth 2009 en `tesis-ucaldas/latex/references.bib`.

9. **¿Cómo controlas el sesgo de abundancia en Gene Ontology?**
En esta versión, no. El plan es un fondo emparejado por abundancia y una prueba competitiva con el universo de proteínas detectadas en la cohorte, no el genoma. Esa sensibilidad no está corrida. Las vías del suero que sí hay son las del cierre del 31 de octubre: 12 a 13 términos y 200 permutaciones. No se repitieron con 10 000. Fuente: `results/v3/hipotesis_v3.tsv`; `tesis-ucaldas/latex/chapters/09-resultados.tex`.

10. **¿SAA1 sirve como marcador? ¿Y BAIAP3?**
No. SAA1 y LBP separan el cáncer de las personas sanas (log2FC 3,26 y 1,25) y no separan el cáncer de la colonoscopia negativa. Se leen como respuesta del huésped. Están en fase exploratoria. BAIAP3, en el cierre del 31 de octubre, tiene un AUC plegado de 0,82 y no pasa el control de falsos (0,20). Es una candidata de otro estudio, no un hallazgo de esta tesis. Fuente: `tesis-ucaldas/latex/chapters/09-resultados.tex`; `results/v3/k4_crc_frente_a_sanos.tsv`.

11. **¿Qué limitaciones tiene la aleatorización mendeliana aquí?**
H9 no está calculada. El diseño, ya escrito, usa solo certificadas secretadas con un instrumento cis, y no decide la tesis. Si la fuente pide registro, el estado será `requiere_acceso`. No hay un efecto causal que defender hoy. Fuente: `results/v3/hipotesis_v3.tsv` (H9); `config/preregistro_v3.yaml`.

12. **¿Qué estudio clínico propondrías y con cuántos pacientes?**
El protocolo ya escrito no es un examen proteómico nuevo. Pregunta si la trayectoria del CEA mejora un umbral fijo de 5 ng/ml. El tamaño escrito es 330 pacientes para una diferencia de AUC de 0,10, y 1 340 para una de 0,05, con potencia del 80 % y alfa del 5 %. Un tamaño para la regla tejido → sangre de H7 no está calculado. Fuente: `docs/fuentes_datos.md`, ficha del protocolo de CEA.

13. **¿Cuántas versiones hubo y cómo evitaste elegir la que más convenía?**
v1, v2 y el cierre del 31 de octubre están congelados. v3 no los reescribe. El preregistro se commiteó antes de la corrida (`09b47dc`, y un arreglo de comillas en `d0b8441` que no cambió el texto). Los umbrales no se movieron después de ver los números. El cambio de cada cifra está en el changelog. Fuente: `results/v3/changelog_numeros.tsv`; `docs/enmiendas_v3.md`.

14. **¿Por qué un AUC de 0,99 contra sanos y de 0,80 contra sintomáticos?**
No son AUC de esta matriz. Álvez 2023, contra personas sanas, da 0,994 con las 1 463 proteínas del panel. Wannberg 2025, contra pacientes con síntomas y sin cáncer, da 0,80 en descubrimiento y 0,82 en réplica, con el mismo tipo de panel. La diferencia está en contra quién se compara. Los controles sanos no se tratan como si fueran sintomáticos. Es una comprobación externa de dirección publicada, no una validación de un examen. Fuente: `docs/fuentes_datos.md`, fichas de Álvez 2023 y Wannberg 2025; `external/alvez_2023/estado.json`.

15. **¿Qué parte del cambio tumoral es estroma?**
No está repartido en v3. R4 de v2, 308 de 688, usa un índice sin referencia celular propia de colon. Repetirlo con Pelka 2021 es H6 y A3, y sigue `no_evaluable`. Fuente: `results/v3/hipotesis_v3.tsv` (H6).

16. **¿Cuál es el contraste principal de PXD034575?**
Primario frente a mucosa, pareado, 16 pares. Concordancia de signo 116/118, intervalo 0,940 a 0,998, Spearman 0,741, concordancia de Lin 0,734. Fuente: `results/v3/k1_comparacion_contrastes.tsv`.

17. **¿Y la metástasis?**
Metástasis hepática frente a mucosa: 93/106, Spearman 0,710, Lin 0,440. El signo se parece; la magnitud, menos. No se llama tumor a esa columna. Fuente: `results/v3/k1_etiquetas.json`.

18. **¿El diseño pareado cambió el resultado, o fue la etiqueta?**
El mismo primario, sin parear, da 124/128. La lectura vieja, con la metástasis como tumor, era 104/119. El cambio grande es la etiqueta. Fuente: `results/v3/k1_comparacion_contrastes.tsv`.

19. **¿R2 cambió?**
No. Sigue en 226 de 688. Cuenta péptidos únicos del proteinGroups. La etiqueta de la columna no entra. Fuente: `results/v3/changelog_numeros.tsv`.

20. **¿Cuántas de las 688 entran en la concordancia?**
118 de 688 tienen un log2FC evaluable en el primario frente a la mucosa. La cobertura es 0,172. La concordancia está condicionada a lo que la réplica alcanza a medir. Fuente: `results/v3/k1_etiquetas.json`, bloque `cobertura_principal`.

21. **¿H1 se sostiene?**
No. H1 predecía concordancia de signo por debajo de 0,70. Se sostiene solo si el límite superior del intervalo queda bajo ese umbral. El límite observado es 0,998. Veredicto: `rechaza_h1`. La cláusula de vías no se recalculó con P frente a N y queda `no_evaluable`. Fuente: `results/v3/hipotesis_v3.tsv`.

22. **¿CHMP4A está medida?**
En abundancia, pasa el corte de cáncer frente a sanos, con 16,85 grados de libertad. En presencia: 13 de 20 sueros de cáncer, 4 de 20 de colonoscopia negativa y 3 de 20 de personas sanas. El Fisher ajustado es 0,171. No se imputa y no se presenta como marcador. Fuente: `results/v3/k4_crc_frente_a_sanos.tsv`; `results/v3/k6_presencia.tsv`.

23. **¿El suero está libre de diferencia técnica?**
No. La fracción lineal de albúmina e inmunoglobulinas tiene medianas 0,253, 0,323 y 0,327, con p de Kruskal–Wallis de 1,23 × 10⁻⁵. La coagulación también difiere, p = 0,0063. Hemólisis, profundidad, mediana de intensidad y plaquetas no difieren. Fuente: `results/v3/k7_control_tecnico.json`.

24. **¿Por qué la fracción de albúmina no se sumó en log2?**
Sumar log2 comprime el rango y no es una fracción de masa. La corrida que se reporta recupera la intensidad lineal. Una suma previa en log2 dio p = 0,056 y no se usó en ningún documento. Fuente: `docs/enmiendas_v3.md`; `results/v3/k7_control_tecnico.json`, campo `escala_fraccion_alta_abundancia`.

25. **¿Dónde están los precursores y el orden de corrida?**
`DIA-NN_report.tsv` no está en el disco: la matriz ancha lo sustituyó. El número de precursores queda `no_disponible`. El orden de corrida no se reconstruye. Fuente: `results/v3/k7_control_tecnico.json`.

26. **¿El AUC de 0,94 demuestra un marcador?**
No. En la matriz humana filtrada el máximo plegado, cáncer frente a colonoscopia negativa, es 0,939 (FGB, P02675), con p = 0,00490 en 10 000 permutaciones. Al quitar contaminantes y paneles baja a 0,865 (SERPINA3, P01011), con p = 0,0464. SERPINA3 no pasa treat. La lectura anterior, 0,9375 en KRT5, queda en `results/v3/k5_auc.json`. Fuente: `results/v3_2/k5_auc.json`, `results/v3_2/k5_auc.tsv`.

27. **¿Por qué no usas DeLong para comparar los dos AUC?**
DeLong compara curvas sobre las mismas personas. Cáncer frente a colonoscopia negativa y colonoscopia negativa frente a sanos no comparten los dos grupos. La diferencia va por bootstrap de los tres grupos. Fuente: `config/preregistro_v3.yaml`, sección `k5_auc`; `results/v3/k5_diferencia_auc.tsv`.

28. **¿Por qué no imputas los huecos?**
En MaxQuant el cero es ausencia de identificación. En DIA-NN la ausencia es ausencia. Rellenar CHMP4A fabricaría una abundancia. La presencia se prueba aparte, con Fisher, y la abundancia solo sobre lo observado. Fuente: `config/preregistro_v3.yaml`; `results/v3/k6_presencia.tsv`.

29. **¿Qué proteínas quedan por debajo del 70 % de observaciones?**
358 de las del suero. Once tienen una diferencia de presencia que pasa el control de falsos. CHMP4A no está entre esas once. Fuente: `results/v3/changelog_numeros.tsv`; `results/v3/k6_presencia.tsv`.

30. **¿De dónde sale el 688?**
De CPTAC, PDC000116, ya en log2. treat a 0,585 y BH a 0,05 dejan 688 de 7 394. v3 no repite un segundo log2. La regresión de v1 sigue reproduciendo 688, 140, 28, 104/119, 226 y 307. Fuente: `results/v3/k1_etiquetas.json`, campo `n_certificadas`.

31. **¿687 de 874, 116 de 118 y 89 de 90 son la misma cuenta?**
No. 687 de 874 y la correlación de Pearson 0,83 son el cierre del 31 de octubre, sobre los genes con efecto en las dos matrices. 116 de 118 es la concordancia de `results/v3/`, antes del filtro de filas CON__/REV__ y con otro mínimo de pares. 89 de 90 es la corrida v3.2: 237 filas excluidas, mínimo de 8 pares, intervalo 0,940 a 0,9997. Fuente: `results/v3_2/k1_etiquetas.json`.

32. **¿Qué cita iba mal en PXD036434?**
El DOI 10.1186/s13073-023-01176-5 es Cascianelli et al., Genome Medicine 2023. El depósito es Praus et al., Journal of Translational Medicine 2023;21:41, doi:10.1186/s12967-022-03855-0. Fuente: `docs/auditoria_v3.md`; `results/v3/auditoria_citas.tsv`.

33. **¿Álvez 2023 se llama «proteograph»?**
No. El título en Crossref es «Next generation pan-cancer blood proteome profiling using proximity extension assay», Nature Communications 2023;14:4308. Fuente: `docs/auditoria_v3.md`.

34. **¿Wannberg es de 2026?**
Crossref da el 29 de diciembre de 2025, Nature Communications 17:151. Se cita 2025. El nombre `wannberg_2026.py` se conserva porque v2 ya escribió en esa ruta. Fuente: `docs/fuentes_datos.md`.

35. **¿PXD055218 es PXD015163?**
No. La entrada vieja `plasma_crc_2024` declara PXD015163 y descarga PXD055218. v3 lee `serum_pxd055218`, con la accesión real. La entrada vieja se deja porque v1 y v2 la usan. Fuente: `config/datasets.yaml`; `docs/auditoria_v3.md`, sección 2.5.

36. **¿La tesis usa el pipeline que imputa con knn?**
No. `config/config.yaml` conserva semilla 42 e imputación knn. Ningún número de v3 sale de `results/differential/`, `results/go/`, `results/gsea/` ni `results/ml/`. Un test falla si una salida de v3 apunta ahí. Fuente: `tests/test_v3_correcciones.py`; `docs/RESULTADOS_Y_LIMITACIONES_TESIS.md`.

37. **¿H7, el modelo de visibilidad, ya predice la sangre?**
No. Está preregistrado: variables fijadas antes de mirar la sangre, validación dejando un recurso fuera, y se sostiene solo si el límite inferior del AUC pasa de 0,55 en al menos dos de tres recursos. Hoy el estado es `no_evaluable`. No es una prueba clínica. Fuente: `config/preregistro_v3.yaml`; `results/v3/hipotesis_v3.tsv`.

38. **¿H8 ya dice que el estadio IV se ve en el suero?**
No. PXD064655, si entra, es un eje dentro del cáncer (estadio IV frente a I–II), no cáncer frente a control. La matriz de precursores no se ha agregado. Estado: `no_evaluable`. Fuente: `results/v3/hipotesis_v3.tsv`; `docs/auditoria_v3.md`, sección 2.4.

39. **¿H10 ya dice que SAA1 no es específica de colon?**
La ficha de Wannberg y la de Álvez ya separan «contra sanos» de «contra sintomáticos». La tabla de candidatas, una fila por proteína y por enfermedad, no está hecha. H10 queda `no_evaluable`. Fuente: `results/v3/hipotesis_v3.tsv`; `docs/fuentes_datos.md`.

40. **¿En qué fase de desarrollo está cada hallazgo?**
Fase 1 de Pepe et al. 2001: muestras ya recogidas, sin una cohorte de cribado. No hay biomarcador validado ni prueba diagnóstica. Fuente: la frase está en el encabezado de este archivo y en `tesis-ucaldas/latex/chapters/90-conclusiones.tex`.

41. **¿Por qué un p de permutación no sale cero?**
Aunque ninguna permutación supere el dato, el p es (b + 1) / (B + 1). Con 10 000 permutaciones el suelo es 1/10001. En el nulo del máximo, b = 244, así que p = 245/10001. Fuente: `results/v3/k5_auc.json`; `results/v3/changelog_numeros.tsv`.

42. **¿APOA4 separa el cáncer?**
No de la colonoscopia negativa. Baja en colonoscopia negativa frente a sanos: es la única proteína diferente de ese contraste. Fuente: `results/v3/k4_potencia.json`; `results/v3/k4_colonoscopia_negativa_frente_a_sanos.tsv` (P06727).

43. **¿Qué semillas se usaron?**
La base es 20260914. Los desplazamientos de v3 están en el preregistro, entre 610 y 670: alfa efectivo, efecto mínimo, bootstrap del AUC, permutaciones y la diferencia de AUC. Fuente: `config/preregistro_v3.yaml`.

44. **¿Dónde está cada número de la sustentación?**
Tejido: `results/v3/k1_etiquetas.json` y `results/v3/k1_comparacion_contrastes.tsv`. Suero: `results/v3/k4_potencia.json`. AUC: `results/v3/k5_auc.json`. Presencia: `results/v3/k6_presencia.tsv`. Control técnico: `results/v3/k7_control_tecnico.json`. Estados: `results/v3/hipotesis_v3.tsv`. Cambios: `results/v3/changelog_numeros.tsv`.

45. **¿Qué no debes decir en la defensa?**
Que no hay señal en las 100 indeterminadas. Que la columna P es un pólipo. Que 89/90 demuestra un cambio propio del colon. Que SAA1, la logística o el modelo H7 son una prueba de consulta. Que Álvez 2023 valida especificidad clínica: sus controles son sanos.
