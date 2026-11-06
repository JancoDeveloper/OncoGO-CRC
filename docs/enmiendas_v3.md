# Enmiendas v3

Fecha de todas: 6 de noviembre de 2026. Ningún umbral de treat, BH, TOST, permutaciones o AUC se movió después de ver los datos. Los valores nuevos están en `results/v3/changelog_numeros.tsv`.

## Cambios de forma, no de umbral

La nota de H1 en `config/preregistro_v3.yaml` tenía dos puntos sin comillas. YAML no leía el fichero. El texto es el mismo. Commit `d0b8441`, anterior a la corrida de K1 a K7.

Los códigos nuevos de R1 no entran en `config/registro_subindices.yaml`. Ese fichero tiene el contrato de cinco códigos de v2 y un test lo fija. Los de v3 están en `config/registro_subindices_v3.yaml`.

La rejilla del efecto mínimo detectable usa 20 pasos entre el margen de treat y 4. Ese número no está en el preregistro. Se fijó en el código antes de la primera corrida y no se cambia para perseguir un resultado.

## Fracción de alta abundancia

La primera corrida sumó la fracción de albúmina e inmunoglobulinas sobre valores ya en log2. El Kruskal–Wallis de esa suma dio p = 0,05596. Ese número no se usó en ningún documento. La corrida que queda recupera la intensidad lineal con exp2 y suma ahí. El p que se reporta es 1,23 × 10⁻⁵ (`results/v3/k7_control_tecnico.json`, campo `escala_fraccion_alta_abundancia`).

## Correcciones

| Código | Valor anterior | Valor v3 | Motivo |
| --- | --- | --- | --- |
| K1, R1 | 104/119, M frente a N | 116/118, P frente a N, pareado | La columna P es el primario. M frente a N queda nombrada: 93/106. |
| K1, R2 | 226/688 | 226/688 | Los péptidos únicos no dependen de la etiqueta. |
| K1, R3 | 0,094 con tres cohortes | `no_evaluable` | PXD019504 y PXD031705 comparten ficheros. El solapamiento con PXD034575 no está recerrado en v3. |
| K4, cáncer frente a colonoscopia negativa | 0 proteínas, sin potencia | 0 diferentes, 145 equivalentes, 731 indeterminadas, de 876 | La ausencia solo se afirma en las equivalentes. |
| K4, cáncer frente a sanos | sin TOST | 3 diferentes, 153 equivalentes, 720 indeterminadas, de 876 | SAA1, LBP y CHMP4A. CHMP4A con 16,8 grados de libertad. |
| K4, colonoscopia negativa frente a sanos | sin TOST | 1 diferente, 217 equivalentes, 658 indeterminadas, de 876 | APOA4. |
| K5 | AUC plegado sin nulo | máximo plegado 0,9375; media nula 0,843; p = 245/10001 | 10 000 permutaciones. El AUC de cada proteína va con signo. |
| K6 | faltantes analizados solo sobre observados | 358 proteínas bajo el 70 %; Fisher de presencia con BH | CHMP4A: 13/20, 4/20 y 3/20. Presencia no significativa tras BH (0,171). No se imputa. |
| K7 | sin pantalla de hemolisis ni de orden | confusión en fracción de alta abundancia y en coagulación; precursores `no_disponible`; orden como limitación | Hemólisis, profundidad y plaquetas no difieren. |
| K9 | el informe podía leer etapas con imputación knn | v3 no lee `results/differential/`, `results/go/`, `results/gsea/` ni `results/ml/` | Test en `tests/test_v3_correcciones.py`. |
| K10 | fracciones sin unidad explícita | numerador, denominador y unidad `proteinas_unicas` | `results/v3/indice_R_v3.tsv`. |
