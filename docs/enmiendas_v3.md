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

## Enmienda v3.4 — 9 de noviembre de 2026

No borra v3.3. Esos números quedan en `results/v3_3/` y en la etiqueta `v3.3-final`. Lo nuevo se calculará en `results/v3_4/`. Los resultados v3.3 ya se habían visto cuando se escribió esta enmienda. La revisión externa del 9-nov-2026 leyó, además, Sample_Metadata.txt de PXD064655 (100 muestras; estadio 4 en 50, 3 en 21 y 1-2 en 29). No leyó la matriz proteómica de PXD064655. Leyó en el artículo de Álvez 2023 que el análisis diferencial es un cáncer frente a los demás cánceres. Vio en las tablas v3.3 que CEACAM5 sube en Tanaka y es específica de CRC en el HPA; por eso la triangulación es exploratoria y post hoc. Ningún umbral de inferencia, de K4, de K5, del filtro principal, de los pares mínimos, del ML, de H1, de H5 ni de H7 se mueve.

H5, universo. El valor anterior leía la concordancia de colon en las 147 proteínas certificadas de K1 y cruzaba los otros carcinomas con todos los genes de Tanaka. El valor nuevo usa un solo universo: genes con treat significativo en Tanaka, gen único y log2FC finito en COAD y en todos los carcinomas disponibles. La selección se hace solo con Tanaka. La cifra 147/147 se conserva como referencia y no entra en la regla. La sensibilidad `G_amplio` usa todos los genes de Tanaka con log2FC finito medidos en todos.

H5, submuestreo. El valor anterior pedía, a cada carcinoma, los 96 pares de COAD. CCRCC tiene 80, el bucle salía con `break` y se escribía 0,0 como si fuera un resultado. El valor nuevo submuestrea todos los carcinomas, COAD incluido, a `n_sub` igual al mínimo de pares comunes. Si `n_sub` es menor que 8, el estado es `no_calculado`. Hay 200 repeticiones con semilla fija. Se informa la fracción al mismo lado de 0,20, la mediana de la fracción pan-carcinoma y los percentiles 2,5 y 97,5. La sensibilidad sin OV usa los otros cinco carcinomas y la regla de 4 de 5. No decide el veredicto. El JSON incluye, por carcinoma, pares y proteínas significativas en la corrida completa y en la submuestreada.

Controles del tejido. El valor anterior contaba como falla un control positivo no evaluable. AQP8 y GUCA2A no tienen filas en el contraste, y SLC26A3 tiene menos de 8 pares. El valor nuevo cuenta primero las filas físicas en `Cohort1_proteinGroups.txt`, antes y después de los filtros de MaxQuant. Los estados son `no_detectado`, `filtrado`, `no_evaluable_por_pares`, `pasa` y `falla`. El bloque pasa si hay al menos un control evaluable, ningún evaluable falla y el QC de mucosa pasa en todos los participantes del contraste principal. Esta regla se escribe después de ver los resultados de v3.3. Se informan los dos veredictos. CA1 es marcador de colonocito y de eritrocito, así que su descenso no basta solo.

Certificado. El valor anterior listaba solo los JSON ausentes, y podía decir que no faltaba nada con casillas en rojo. El valor nuevo lista cada casilla del checklist con estado verde, rojo o limitación declarada. El estado es `publicable` solo si ninguna casilla está en rojo. Las limitaciones aceptadas son H9 con `requiere_acceso`, ESM-2 o TMbed no ejecutados si el autor dijo que no, un recurso de H7 `no_disponible` con su motivo, y H6 `no_evaluable` si ni Pelka ni el HPA tienen leyenda verificable.

H7. El valor anterior describía a Álvez 2023 como una validación frente a sanos. El contraste del artículo es colorrectal frente a los demás cánceres. El HPA sigue siendo colorrectal frente a Healthy. Papier 2024, si el fichero trae todas las proteínas, es riesgo prospectivo: la etiqueta es 1 si p < 0,05/639 y el signo de log(HR) coincide con CPTAC. Si solo hay asociaciones significativas, el recurso queda `no_disponible`. La regla de sostiene y refuta no cambia. Se comprueba si los casos de 2023 están en el atlas de 2025. Se escriben la tabla 2×2, el kappa y `modelo_completo`.

H6. El valor anterior quedaba `no_evaluable`. El valor nuevo clasifica cada gen certificado como epitelial, no epitelial o mixto, con Pelka 2021 si la leyenda se verifica y, si no, con el tipo celular único del HPA. La prueba es un Fisher bilateral. El sentido no se predice.

H8. No se había calculado. La fuente es PXD064655, agregada con directLFQ. El contraste principal es estadio IV frente a I–II. El estadio III queda fuera del principal. La prueba es un Fisher unilateral después de quitar las listas de ruido. Si el mapa de corridas no es uno a uno, el estado es `no_evaluable`.

H10 y triangulación. Wannberg entra solo como descripción, para las seis proteínas de partida y las catorce de la triangulación. La triangulación es exploratoria y post hoc. No hay pruebas nuevas. El estado de H10 no cambia.

Multiverso. El valor anterior calculó 48 celdas. El valor nuevo calcula las 288 del suero, el ML en las 72 combinaciones de cuantificación, contaminantes, paneles y presencia, y R1 en las dos normalizaciones. La celda preregistrada conserva las 1000 permutaciones de la corrida principal.

H9. No hay un GWAS de CRC descargable sin cuenta y verificado en el GWAS Catalog. Queda `requiere_acceso`.

El mismo día, después de esa frase, se verificó `GCST90255676`: europeos, sin UK Biobank, GRCh37, md5 `919638df0721391b4414e8873a2c87da`. El condicional ya escrito —Wald con el cis-pQTL líder de Sun 2023— queda habilitado. SAA1 no entra si no tiene fila cis en la tabla de descubrimiento. El cálculo es exploratorio y no decide la tesis. La frase anterior no se borra.

Las preguntas iniciales se cierran con el valor por defecto: Firebase solo con el checklist en verde; el PDF como zip de Overleaf; ESM-2 y TMbed sí, porque el disco externo tiene más de 20 GB libres; los metadatos de adquisición sí cuando estén los 60 crudos.

Las cifras v3.3 no se borran: pasan al anexo de historial con su motivo.

## Enmienda v3.5 — 10 de noviembre de 2026

No borra v3.4. Esos números quedan en `results/v3_4/` y en la etiqueta `v3.4-final`. Lo nuevo se calculará en `results/v3_5/`. Los resultados v3.4 ya se habían visto cuando se escribió esta enmienda. La revisión del 10-nov-2026 leyó el código de `bff0c5e`, los JSON de v3.4 y el encabezado de `report.pr_matrix.tsv`. No calculó el contraste de PXD064655 ni leyó los JSON de ThermoRawFileParser. Ningún umbral de inferencia, de K4, de K5, del filtro principal, de los pares mínimos, del ML, de H1, de H5, de H6, de H7 ni de H8 se mueve.

H8, mapa. El valor anterior usaba `Path(columna).name` en macOS. Esa llamada no corta una ruta de Windows, el patrón no encontraba `SH_` y el mapa unía 0 de 100. El valor nuevo toma `PureWindowsPath` si hay barra invertida y `PurePosixPath` si no, y exige `SH_NN.mzML` completo. Las columnas `ST_` son estándares de lote: se excluyen y se cuentan. La columna de estadio es la que contiene `stage`, sin distinguir mayúsculas, sin retorno de carro ni espacios.

H8, agregación. El valor anterior comprobaba que el paquete se pudiera importar y escribía que la agregación no se había ejecutado. El valor nuevo llama a directLFQ 0.3.3 con el formato genérico `protein` / `ion` y los parámetros por defecto. El Fisher, los filtros y las listas de ruido siguen siendo los de la enmienda v3.4. Si después de cablearlo el estado sigue siendo `no_evaluable`, queda como limitación declarada. No vuelve rojo el certificado. Lo decidió el autor el 10-nov-2026.

Orden de adquisición. El valor anterior contaba `*.raw` (120, porque hay 60 ficheros AppleDouble) y, si encontraba el binario, decía que el analizador no estaba cableado. El valor nuevo cuenta solo `(CP|CN|HC)_NN_DIA.raw`. El campo de fecha, leído en un JSON del disco antes de esta frase, es `Content Creation Date`. Se calculan el AUC y los Spearman con los crudos que tienen fecha. Con menos de 60 el estado es `no_evaluable`. No se inventan tiempos. No cambia el titular y no vuelve rojo el certificado. Los seis crudos rechazados no se vuelven a descargar hoy. Si el autor los descarga otro día, una v3.6 cierra la prueba de 60.

Certificado. El valor anterior aceptaba un `no_evaluable` si la lista contenía el texto `H6`, así que H8 pasaba con la regla de H6. La casilla de H9 apuntaba a `hipotesis_v3_4.tsv`, que no existe. El valor nuevo compara el identificador de la pieza y su estado con una lista estructurada. Añade casillas de preanalítica, triangulación, reproducción del ML, regresión v1, manifiestos y la fuente real de H9. Publicar exige que cada limitación declarada esté en esa lista y que exista cada fuente citada.

H7. El valor anterior metía a Papier 2024, con 1 positivo y 125 negativos, en el entrenamiento de Álvez y del HPA. El valor nuevo saca de los recursos entrenables a quien tenga menos de 10 observaciones en alguna clase y lo deja `no_evaluable`, con su conteo. La regla de sostiene y de refuta no cambia. Se informan las dos versiones.

H6. El valor anterior devolvía `no_sostiene` con una tabla sin signos distintos. El estado no se cambia después de ver el resultado: lo sigue dando la regla escrita. El valor nuevo añade `tabla_degenerada` cuando una columna de la tabla 2×2 es cero, y el motivo de que la prueba no puede detectar una diferencia. La lectura no puede decir que el tipo celular no influye.

Aprendizaje. El valor anterior de la celda preregistrada en v3.4 usaba la semilla del multiverso y daba un AUC distinto de 0,88025. El valor nuevo, solo para esa celda, usa las semillas de v3.3 y el mismo orden de muestras. El resto de celdas conserva la semilla de v3.4. Si no reproduce 0,88025, la casilla queda en rojo.

Controles. El valor anterior copiaba el estado `calculado` de v3.3. El veredicto de v3.3 era falla cuando `n_fallos` era mayor que cero. El valor nuevo lo deriva de `n_fallos`.

Triangulación. El valor anterior dejaba solo la columna `gen`. El valor nuevo rellena las columnas de la enmienda con los ficheros que ya están en disco. Lo que falta se escribe `no_medida`. No hay pruebas nuevas.

H9. La enmienda v3.4 decía que H9 no se intentaba. El mismo día, el commit `f5a26fb` habilitó el Wald exploratorio de `GCST90255676`. Se declara como desviación de procedimiento. H9 no decide la tesis. La frase anterior no se borra.

Las preguntas del 10-nov-2026 quedan cerradas: Firebase solo si el certificado v3.5 es publicable; el PDF como zip de Overleaf, sin instalar MacTeX.

Las cifras v3.4 no se borran: pasan al anexo de historial con su motivo.

## Enmienda v3.6 (definitiva), 10 de noviembre de 2026

Corrida limpia. La v3.5 copió H5, H9, H10, el multiverso y la secuencia desde `results/v3_4/`, y la web citaba ficheros de `results/v3_3/`. La v3.6 recalcula cada pieza en `results/v3_6/` y no copia ninguna. Un vigilante registra cada apertura y falla si se lee una corrida anterior, salvo `results/certified/` y las excepciones de `config/entradas_v3_6.yaml`. No cambia ningún umbral ni ninguna regla de decisión.

Presentación del orden de adquisición. El intervalo de LeDell queda en nulo cuando el AUC es 0 o 1, porque esa separación completa degenera el intervalo. Se añade la p exacta de Mann-Whitney y una cota descriptiva por si los crudos sin fecha se leyeran. Se describe la posición por grupo y si los bloques no se solapan. El estado lo sigue dando la regla escrita. La lectura obligatoria es que calidad, orden y grupo van juntos. La clave del titular sigue siendo `separacion_explicada_por_calidad`.

Los seis crudos. Se hace un intento con la ThermoRawFileParser ya instalada, versión 2.0.0.0, primero metadatos en JSON y, si falla, en texto. No se descargan y no se imputa ninguna fecha. Si se leen los sesenta, la prueba preregistrada se aplica tal cual. Si no, el estado sigue `no_evaluable` y la cota queda escrita.

H8. Se añade un embudo descriptivo del universo. El universo sigue siendo toda la tabla de Tanaka del contraste primario frente a mucosa, no un subconjunto de proteínas concordantes. Si el código ya cumple esa regla, el motivo explica el solapamiento. La regla de Fisher no se mueve.

Formato. Tesis, web y figuras usan un solo módulo de formato. Los JSON guardan el valor completo.

Las corridas v3.2 a v3.5 no se borran. Quedan con su etiqueta y su manifiesto.

## Enmienda v3.7 (definitiva), 10 de noviembre de 2026

Qué cambia. La v3.7 recalcula todas las piezas en `results/v3_7/`. Corrige fallos de implementación que ya se habían visto en v3.6: la unión de los sanos con la tabla de calidad, el mapa de genes de Tanaka en H8 y en la triangulación, los embudos para que cada paso sea un subconjunto del anterior, y el sello, que pasa a llevar el SHA-256 del manifiesto de entradas y la etapa `v3_7`.

Por qué. Esos fallos cambian qué filas entran en una prueba sin cambiar la regla. La v3.6 dejó el Spearman de los sanos en `no_evaluable` porque los crudos se llaman `HC_` y la tabla de calidad `Normal_`. H8 asignaba el gen con un mapa de UniProt revisado y dejaba fuera proteínas que el mapa de MaxQuant sí nombra. El embudo contaba cada paso sobre el mismo conjunto, no sobre el anterior.

Qué no cambia. No se mueve la inferencia, el contraste del suero, el AUC de K5, las hipótesis, las vías del suero tal como se preregistraron, ni las enmiendas v3.3, v3.4 y v3.5. El titular y el formato numérico de la v3.6 se conservan.

Crudos. Los 60 se comparan con el SHA-1 de `checksum.txt` de PRIDE. El que no coincida se vuelve a descargar hasta tres veces. El que coincida y siga ilegible se lee con ThermoRawFileParser 2.0.0 y con la 1.4.5. Un lector de cabecera solo se acepta si la fecha coincide en las 54 corridas que el lector oficial ya leyó. Si después de eso falta la fecha, el fichero se declara dañado en origen y el estado sigue `no_evaluable`, con la cota.

Piezas que el preregistro ya pedía y no se habían ejecutado. Las vías del suero, con la semilla 670 y 10 000 permutaciones, y la cláusula de vías de H1, con la semilla 101 y Tanaka como réplica.

Descriptivos que no deciden. La correlación dentro de cada grupo con su q de Benjamini-Hochberg, la deriva del estándar interno, la sobrerrepresentación de GO en el tejido, la auditoría del catálogo de candidatas y los neoantígenos como frontera.

Las corridas v3.2 a v3.6 quedan congeladas con su etiqueta. Solo se leen después de la corrida, para el changelog.
