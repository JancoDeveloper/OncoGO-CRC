# Banco de preguntas, corrida v3.3

Cada cifra sale de `results/v3_3/`. Este archivo lo escribe `scripts/v3/escribir_documentacion_v3_3.py`. El banco del 6 de noviembre describe v3.2 y no es la corrida vigente.

1. **¿Cuál es la pregunta de la corrida v3.3?**
Si el suero separa el cáncer de la colonoscopia negativa más allá de la calidad de muestra. Fuente: certificado, campo titular_ml.

2. **¿Cuántas proteínas tiene el contraste principal?**
284. Fuente: suero_resumen.json, n_principal.

3. **¿Cuántas son diferentes, equivalentes e indeterminadas?**
0, 166 y 118. No se detecta diferencia; 166 proteínas son equivalentes

4. **¿Qué dice el titular del modelo?**
La separación del suero entre cáncer y colonoscopia negativa se explica por la calidad de muestra.

5. **¿Por qué el titular no es que el suero separa el cáncer?**
M2 tiene AUC 0,880 y p 0,000999, pero M5, ya sin la calidad, tiene AUC 0,647 y p 0,119. La regla pide las dos p menores que 0,05 y un AUC de M5 de al menos 0,70.

6. **¿Qué pasó con los AUC 0,996, 0,981 y 0,819 de v3.2?**
Quedan retirados. Esas corridas duplicaban cada muestra. La demostración con ruido pasa de 0,475 a 0,931.

7. **¿La media del nulo cae donde debía?**
La de M2 es 0,473. La regla pedía entre 0,45 y 0,55.

8. **¿Cuántas filas entran al modelo?**
40: 20 de cáncer y 20 de colonoscopia negativa. Fuente: ml_suero_permutacion.json, n_por_clase.

9. **¿H1 se sostiene?**
No. El estado es rechaza_h1. La concordancia es 147/147, y la regla solo sostiene H1 si el límite alto del intervalo queda bajo 0,70.

10. **¿Esa concordancia es frágil al centrado?**
No. h1_fragil_a_normalizacion es False.

11. **¿Las etiquetas P y N quedaron confirmadas?**
No. CA1 pasa. AQP8, SLC26A3 y GUCA2A fallan. Fuente: controles.tsv.

12. **¿Por qué fallan AQP8 y GUCA2A?**
No están en la salida de MaxQuant. No se inventa una fila.

13. **¿Por qué falla SLC26A3?**
Tiene una fila y el log2FC no es finito. El control lo marca no evaluable y la etapa lo cuenta como falla.

14. **¿Qué es la columna P?**
El tumor primario de PXD034575. No es un pólipo.

15. **¿Cuántos pares tiene el tejido?**
16 de primario frente a mucosa y 16 de metástasis frente a mucosa. El mínimo sigue en 8.

16. **¿Cuántas filas se excluyen del tejido?**
237, por CON__, REV__ y las marcas de MaxQuant.

17. **¿Cuál es la concordancia de la metástasis?**
116/143, con Lin 0,434.

18. **¿Cuántos grupos trae el informe de suero?**
897. Tras el FDR de proteína y la especie humana quedan 312.

19. **¿Por qué la matriz principal no es 312?**
312 es la matriz humana con FDR. El contraste principal quita contaminantes de piel y paneles y se queda en 284.

20. **¿Qué proteínas separan el cáncer de las personas sanas?**
LRG1 (P02750), SAA1 (P0DJI8), LBP (P18428). Fuente: k4_crc_frente_a_sanos.tsv, clase diferente, cruzada con la traza.

21. **¿Esa separación frente a sanos es la pregunta clínica?**
No. La colonoscopia negativa reúne otras enfermedades intestinales. Comparar solo con sanos infla la señal inflamatoria.

22. **¿SAA1 o LBP son un examen?**
No. Se separan de las personas sanas y el contraste principal, frente a colonoscopia negativa, no tiene proteínas diferentes.

23. **¿H5 qué dice?**
Estado refuta. La fracción pan-carcinoma es 0,119.

24. **¿Por qué UCEC no entra en H5?**
Sus columnas normales se llaman NX y no comparten identificador con el tumor. El mínimo de ocho pares no se bajó.

25. **¿H7 qué dice?**
Estado refuta. AUC 0,503 en Álvez 2023 y 0,502 en el atlas del HPA. Los dos intervalos incluyen 0,50.

26. **¿H7 es una prueba diagnóstica?**
No. Es un modelo de fase 1, fuera de recurso y fuera de proteína. El certificado lo dice.

27. **¿Qué le falta a H7?**
ESM-2, TMbed, Papier 2024 y otras variables listadas en variables_faltantes. No se sustituyen.

28. **¿H10 qué dice?**
Estado descriptiva. Cuenta 5 proteínas específicas. No hay umbral de sostiene o refuta.

29. **¿Por qué H10 no se declara refutada?**
El preregistro no fijó ese umbral. Fijarlo después de ver el atlas sería elegirlo por el resultado.

30. **¿El multiverso ya está completo?**
No. Lo calculado es la rejilla K4, 48 celdas. K5, el modelo y la normalización del tejido no se repitieron en cada celda.

31. **¿Se puede publicar ya?**
No. El checklist de la sección 16 no está en verde. No hay push ni Firebase.

32. **¿Se imputó el suero para el contraste?**
No. El cero de MaxQuant y la ausencia de DIA-NN se quedan como ausencia. La mediana del modelo solo se calcula dentro del pliegue de entrenamiento.

33. **¿Se ejecutó DIA-NN otra vez?**
No. Se usa el report ya hecho, de 213 313 241 bytes.

34. **¿Cuántas permutaciones tiene K5?**
10 000, en tres matrices. Fuente: k5_auc.json.

35. **¿Cuántas permutaciones tiene el modelo?**
1 000 por modelo. La p no puede ser menor que 1/1001.

36. **¿Qué semilla se usó?**
La base 20260914. Los desplazamientos de v3.3 están en la enmienda, no se eligieron después de ver el resultado.

37. **¿La albúmina humana quedó marcada como contaminante?**
No. El cruce es por acceso y especie. La prueba de plasma clásico pasa.

38. **¿Cuántas entradas tiene la biblioteca de contaminantes?**
381. Fuente: suero_controles.json, libreria_n_entradas.

39. **¿Quién falla el control de sexo?**
CRC15 y CRC16. Fuente: qc_sexo_por_participante.tsv.

40. **¿Quién falla el hígado?**
CRC13, CRC16 y CRC17. Fuente: qc_etiquetas_higado.tsv.

41. **¿La mucosa de los 16 participantes es evaluable?**
Sí, y los 16 pasan. Fuente: qc_etiquetas_mucosa.tsv.

42. **¿A quién saca la sensibilidad de etiquetas?**
A CRC16, el único que falla dos de tres controles. Fuente: k1_etiquetas.json.

43. **¿v3 y v3.2 se reescribieron?**
No. Siguen congelados por sus manifiestos. v3.3 escribe solo en results/v3_3.

44. **¿De dónde sale cada cifra de la página?**
De scripts/v3/exportar_sitio.py, que lee el certificado y el changelog.

45. **¿De dónde sale cada cifra de la sección nueva de la tesis?**
De cifras.tex, escrito por scripts/v3/exportar_cifras_tesis.py.

46. **¿Qué fase de Pepe es esta?**
La fase 1. Exploración en muestras ya recogidas. No hay una cohorte de cribado.

47. **¿El intervalo de los pliegues es un intervalo de confianza?**
No. La dispersión entre pliegues se informa aparte. El intervalo de LeDell es el que lleva ese nombre en el JSON.

48. **¿CHMP4A volvió a entrar?**
Está en la traza, acceso Q9BY43;Q9BY43-2, con Global.PG.Q 0,0374 y fracción de péptidos proteotípicos 0. No está entre las proteínas diferentes del contraste principal.

49. **¿Qué no se debe decir en la sustentación?**
Biomarcador validado, prueba diagnóstica, y pólipo para la columna P.

50. **¿Qué queda abierto a propósito?**
Los tres controles sin log2FC, el multiverso fuera de la rejilla K4, y los modelos de secuencia que no se ejecutaron.

51. **¿El resumen de v3.2 sigue en la tesis?**
Sí, en el capítulo de resultados, marcado como corrida del 6 de noviembre. La sección vigente es la de v3.3.

52. **¿Cuántas diapositivas hay?**
18, en tesis-ucaldas/sustentacion/diapositivas.tex, escritas por este mismo guion.
