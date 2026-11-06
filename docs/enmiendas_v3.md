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

## Enmienda v3.2 — 6 de noviembre de 2026

No borra los números de v3. Los deja con su nombre. Lo nuevo se calcula en `results/v3_2/`.

El report de PXD055218 se vuelve a bajar a `data/raw/plasma/DIA-NN_report.tsv`. El SHA-256 esperado es `A8A5AADD74391AFEF68CF4D6C3A7D3860023101E5C92590E2E87C152A39E59E8` y el tamaño es 213 313 241 bytes. El fichero queda solo lectura. La matriz ancha ya no lo sustituye.

Filtro principal del suero: `Q.Value` ≤ 0,01, `PG.Q.Value` ≤ 0,01 en la corrida, `Global.PG.Q.Value` ≤ 0,01, al menos un miembro `_HUMAN`, cuantificación `PG.MaxLFQ` por `Protein.Group`. Tres sensibilidades: los filtros por defecto de DIA-NN sin el q por corrida; `Genes.MaxLFQ` con `GG.Q.Value`; y el filtro principal sin el corte de especie.

Contaminantes: FASTA universal de Frankenfield et al. 2022, commit `bddf6da9334f4300caed75f7a87c096dd96b6bca` (3 de julio de 2026), SHA-256 `c0782146be6dad6ecf53e23e5a21a7ed4fea9b4802cb01e8d00091404ec66926`, 201 365 bytes, 381 entradas. Queratinas de piel y DCD, lista cerrada: KRT1, KRT2, KRT5, KRT9, KRT10, KRT14, KRT16, KRT6A, KRT6B, KRT6C, KRT17, KRT73, KRT78 y DCD. El contraste principal las quita. La librería también trae KRT8, KRT18, KRT19 y KRT20. Salen del contraste principal porque la librería se fija entera, y vuelven en la sensibilidad. No se retiran a mano.

Paneles de Geyer et al. 2019. La tabla suplementaria no se obtuvo: la URL probada devolvió 404. Se usa la lista del texto, sin genes añadidos. Eritrocitos: HBA1, HBB, HBD, CA1, CA2, PRDX2, CAT, SLC4A1. Plaquetas: PF4, PPBP, TLN1, MYH9, ITGA2B, y los ya nombrados en el cuerpo del artículo en el preregistro de v3 (GP1BA, THBS1, GP5, PF4V1). Coagulación: FGA, FGB, FGG, F13A1, SERPINC1. La puntuación de cada muestra es la mediana del z de log2. Esas proteínas salen del contraste principal y se reportan aparte. La puntuación de coagulación entra como covariable en una sensibilidad.

Tejido. Se excluyen identificadores `CON__` y `REV__`, y las filas con Reverse, Potential contaminant u Only identified by site igual a `+`. `min_pares_observados` vuelve a 8, el valor del cierre del 31 de octubre. v3 lo había bajado a 3 sin enmienda. En 16 pares, un efecto visto en menos de la mitad no entra en el contraste principal. No se elige por el número de proteínas que salen.

Semillas de v3.2: base 20260914, desplazamientos 710 a 770. Están en el bloque `enmienda_v3_2` del preregistro.

Declaración de transparencia. El 6-nov-2026, una revisión externa analizó el report original de PXD055218 con filtros de FDR y especie, y vio un ranking exploratorio (Mann-Whitney, sin umbral). Esta enmienda se motiva en las guías de DIA-NN y del HPP y en la presencia de especies no humanas, no en ese ranking. El análisis oficial es el preregistrado.

Los 60 crudos (65,6 GB) no se descargan. La búsqueda con AlphaDIA queda como extensión documentada. La web no se despliega.
