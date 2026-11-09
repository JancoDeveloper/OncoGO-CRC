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

Los 60 crudos (65,6 GB) no se descargan. La búsqueda con AlphaDIA queda como extensión documentada. La web no se despliega en esta corrida.

## Resultado de esa corrida

Fuente: `results/v3_2/suero_resumen.json` y `results/v3_2/k1_etiquetas.json`. Sello del suero: commit `293e6ee`.

Suero, cáncer frente a colonoscopia negativa, matriz principal de 205 proteínas: 0 diferentes, 105 equivalentes, 100 indeterminadas. La frase de la tabla es «No se detecta diferencia; 105 proteínas son equivalentes». Frente a sanos, en esa matriz, son diferentes SAA1 y LBP. El informe tiene 897 grupos. Con FDR de proteína y especie humana quedan 312.

Tejido, después de excluir 237 filas y con mínimo de 8 pares: primario frente a mucosa 89/90, intervalo 0,940 a 0,9997, Lin 0,704. Metástasis frente a mucosa 69/82. H1 queda en `rechaza_h1`.

El máximo del AUC plegado baja de 0,939 (FGB, matriz humana filtrada, p = 0,00490) a 0,865 (SERPINA3, sin contaminantes ni paneles, p = 0,0464). Piel y coagulación difieren entre grupos. La logística y la coagulación sola rechazan el nulo de permutación. No es una prueba diagnóstica.

H5 quedó `no_disponible`: el carcinoma de endometrio no tiene 8 pares con el mismo identificador. H10 cuenta 0 específicas de 1 165, en 60 enfermedades, y el veredicto queda `no_evaluable` porque el preregistro no fija un umbral. H6, H7, H8 y H9 quedan sin número, con el motivo de `hipotesis_v3_2.tsv`.

La página pública se armó después, a pedido del autor: https://oncogo-crc-tesis.web.app. La frase de arriba sobre no desplegar la web es la regla de la corrida.

## Enmienda v3.3 — 9 de noviembre de 2026

No borra v3.2. Esos números quedan en `results/v3_2/` y en la etiqueta `v3.2-final`. Lo nuevo se calculará en `results/v3_3/`. El autor ya había visto los resultados de v3.2 cuando se escribió esta enmienda. Ningún umbral de inferencia, TOST, AUC, FDR, pares mínimos ni regla de H1, H5 o H7 se mueve.

Contaminantes. El valor anterior cruzaba por símbolo de gen, sin mirar la especie. La albúmina bovina de la librería marcaba la albúmina humana, y lo mismo pasaba con el fibrinógeno, la transferrina, el complemento y las apolipoproteínas: 98 de 312 grupos humanos. El valor nuevo cruza por acceso, sin `Cont_`, sin `sp|…|` y sin el sufijo de isoforma. Un grupo humano es contaminante solo si uno de sus accesos está en la librería con `OX=9606`. Un grupo sin miembro humano es contaminante si alguno de sus accesos está en la librería. El símbolo solo se usa para entradas humanas y solo si el grupo no trae `Protein.Ids`. Los grupos mixtos se quedan en el principal con la marca `comparte_con_contaminante` y salen en la sensibilidad `sin_grupos_mixtos`. KRT7, KRT8, KRT18, KRT19 y KRT20 siguen fuera del principal, porque la librería se fija entera, y vuelven en `con_queratinas_epiteliales`. ALB, APOA1, C3, TF, SERPINA1, A2M, APOE, FGA, FGB, FGG y LRG1 humanas no pueden salir como contaminantes.

ML del suero. El valor anterior añadía cada muestra dos veces, una por cada matriz, y la validación cruzada veía la misma persona en el entrenamiento y en la prueba. Por eso el nulo no era interpretable. El valor nuevo usa columnas únicas, 20 CP y 20 CN, entra en log2 si hace falta, y calcula seis modelos fijados aquí: proteínas humanas, proteínas del principal, solo calidad de muestra, solo coagulación, principal residualizado por la calidad dentro del pliegue, y cinco componentes principales dentro del pliegue. El titular sale de una regla escrita antes de correr. La dispersión entre pliegues no se llama intervalo. El intervalo es el de LeDell et al. 2015.

Tejido. El valor anterior buscaba el gen por acceso exacto y, si no lo encontraba, escribía `no_evaluable` sin detener la etapa. CA1, AQP8, SLC26A3 y GUCA2A no se evaluaban. El valor nuevo usa un solo mapa: Gene names, el `GN=` de las cabeceras FASTA y todos los Protein IDs pasados por el mapa de UniProt. Un control positivo no evaluable cuenta como fallo. El análisis principal no normaliza. Se añade la sensibilidad de centrado por mediana. Si la concordancia cruza 0,70 solo con esa sensibilidad, H1 se declara frágil a la normalización.

H5. El valor anterior devolvía `no_disponible` en cuanto un carcinoma, UCEC, no tenía ocho pares. El valor nuevo excluye ese carcinoma con su diagnóstico y evalúa la hipótesis si quedan al menos seis. Pan-carcinoma sigue siendo mismo signo y treat en al menos cuatro de los disponibles. Con siete carcinomas, la regla es la de v3.

H10. El valor anterior tomaba la primera fila de colorrectal e ignoraba la columna Control. El valor nuevo lee esa columna y llama específica a la proteína que es significativa, con el mismo signo, frente a sanos, frente a la misma clase y frente a todas las demás. El informe empieza por SAA1, LBP, CHMP4A, BAIAP3, FN1 y TPM4. H10 queda descriptiva: el preregistro no fijó umbral, y fijarlo después de ver 0/1165 sería elegir por resultado.

H7. El valor anterior se escribía a mano como `no_evaluable`. El valor nuevo es un modelo por proteína, fuera de recurso y fuera de proteína, con las variables que ya estaban en el preregistro. ESM-2 y TMbed entran solo como modelo secundario. No se sustituye una variable que falte.

El autor pidió, el mismo día, desplegar la web solo si el checklist queda en verde, guardar ProtT5-XL en el disco externo y descargar los 60 crudos al disco externo. Esos crudos no entran en el cálculo y DIA-NN no se ejecuta.

Estas correcciones cambian cifras que ya estaban en la tesis y en la web. Las cifras v3.2 no se borran: pasan al anexo de historial con su motivo.
