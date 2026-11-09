# Qué se ve del tumor en la sangre

La página con las gráficas, las fichas, la teoría y los diccionarios está en [divulgacion.html](divulgacion.html). Lo que hay que estudiar está en [sustentacion.md](sustentacion.md). El mapa de la teoría y de los enlaces está en [diccionario.md](diccionario.md). Esta nota deja los números por escrito.

El suero de PXD055218 tiene tres grupos de 20. En el contraste oficial de v3.2, 205 proteínas, con FDR 0,05 y |log2FC| mayor que 0,585, ninguna separa el cáncer de la colonoscopia negativa. 105 son equivalentes y 100 quedan indeterminadas. Frente a sanos pasan SAA1 y LBP. CHMP4A no entra en las 312: su q global de grupo es 0,037.

En el catálogo del 31 de octubre, ocho candidatas de la literatura superan la regla de AUC. Cinco están marcadas como fase aguda. BAIAP3, FN1 y TPM4 no llevan esa marca y tampoco pasan el FDR. Esa lista no se recalculó en v3.2.

PXD034575, Tanaka et al. 2024, es tumor primario (P), metástasis hepática (M) y mucosa (N). No es un pólipo. Tras tratar el cero de MaxQuant como ausencia y pasar a log2, el primario coincide en signo con CPTAC en 687 de 874 genes, correlación 0,83. La metástasis correlaciona 0,59. Esa cuenta es el cierre del 31 de octubre, sobre todos los genes con efecto en las dos matrices. La cifra 0,104 de la versión 2 es otro cálculo y sigue en los resultados de reproducibilidad.

La corrida v3.2 del 6 de noviembre está en `results/v3_2/`. Entre las 688 proteínas certificadas, el primario frente a la mucosa coincide en signo en 89 de 90 (intervalo 0,940 a 0,9997; Lin 0,704). La metástasis frente a la mucosa coincide en 69 de 82 (Lin 0,390). El contraste oficial del suero tiene 205 proteínas. El corte de 1,5 veces deja cero diferentes entre cáncer y colonoscopia negativa. De esas 205, 105 son equivalentes y 100 quedan indeterminadas. La ausencia solo se afirma en las 105. Fuente: `results/v3_2/k1_etiquetas.json` y `results/v3_2/suero_resumen.json`.

De los términos GO del tejido, el suero solo comparte 12 con al menos 10 genes en cáncer frente a sano. Su correlación de NES no supera el nulo (p 0,34).

En el cierre del 31 de octubre, primario frente a mucosa: 58 proteínas significativas de 2490. El signo entre primario y metástasis coincide en el 80,8 % de 876 genes.
