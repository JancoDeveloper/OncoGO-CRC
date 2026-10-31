# Qué se ve del tumor en la sangre

La página con las gráficas, las fichas, la teoría y los diccionarios está en [divulgacion.html](divulgacion.html). Lo que hay que estudiar está en [sustentacion.md](sustentacion.md). El mapa de la teoría y de los enlaces está en [diccionario.md](diccionario.md). Esta nota deja los números por escrito.

El suero de PXD055218 tiene tres grupos de 20. Con FDR 0,05 y |log2FC| mayor que 0,585, ninguna proteína separa el cáncer de la colonoscopia negativa. Frente a sanos pasan SAA1, LBP y CHMP4A. CHMP4A tiene 16 de 40 observaciones. APOA4 baja en la enfermedad intestinal frente a los sanos.

Ocho candidatas de la literatura superan la regla de AUC. Cinco están marcadas como fase aguda. BAIAP3, FN1 y TPM4 no llevan esa marca y tampoco pasan el FDR.

PXD034575, Tanaka et al. 2024, es tumor primario (P), metástasis hepática (M) y mucosa (N). No es un pólipo. Tras tratar el cero de MaxQuant como ausencia y pasar a log2, el primario coincide en signo con CPTAC en 687 de 874 genes, correlación 0,83. La metástasis correlaciona 0,59. La cifra 0,104 de la versión 2 es otro cálculo y sigue en los resultados de reproducibilidad.

De los términos GO del tejido, el suero solo comparte 12 con al menos 10 genes en cáncer frente a sano. Su correlación de NES no supera el nulo (p 0,34).

Primario frente a mucosa: 58 proteínas significativas de 2490. El signo entre primario y metástasis coincide en el 80,8 % de 876 genes.
