# Auditoría v3

Fecha: 6 de noviembre de 2026. Rama local `v3`. Los números nuevos salen de `results/v3/` en el workspace de datos. Lo que no se calculó queda `no_evaluable`.

## 2.1 Etiquetas de PXD034575

El fichero viejo sigue igual. `metadata/sample_metadata/pride_merged_metadata.tsv`, línea 4: `LFQ intensity CRC02P` tiene `condition=polyp`. El recuento es 16 `tumor`, 16 `normal` y 16 `polyp`.

`docs/fuentes_datos.md`, líneas 40 a 47, ya decía lo contrario: en Tanaka et al. 2024, P es el primario, M la metástasis hepática y N la mucosa. La versión 2 calculó la réplica con M frente a N.

v3 no edita ese fichero. Escribe `metadata/sample_metadata/pride_pxd034575_v3.tsv`: 16 `primary`, 16 `liver_metastasis`, 16 `mucosa`. El contraste principal es primario frente a mucosa.

Resultado en `results/v3/k1_etiquetas.json` y `results/v3/k1_comparacion_contrastes.tsv`. Universo: 688 proteínas certificadas. Dieciséis pares en cada contraste.

| Contraste | Mismo signo | Spearman | Lin |
| --- | --- | --- | --- |
| Primario frente a mucosa, pareado | 116/118 (IC 0,940 a 0,998) | 0,741 | 0,734 |
| Metástasis hepática frente a mucosa, pareado | 93/106 (IC 0,799 a 0,933) | 0,710 | 0,440 |
| Primario frente a mucosa, no pareado | 124/128 (IC 0,922 a 0,991) | 0,737 | 0,746 |

El valor congelado de v1 y v2 era 104/119, calculado con M como si fuera el tumor. El contraste no pareado con P aísla el cambio de etiqueta del cambio de diseño.

R2 sigue en 226/688. Cuenta péptidos únicos y la etiqueta no entra.

R3 queda `no_evaluable`. Ver 2.2.

H1 predice concordancia de signo menor que 0,70 (`config/preregistro_v2.yaml`, cláusula de proteína). La cláusula se sostiene solo si el límite superior del intervalo queda por debajo de 0,70. El límite observado es 0,998. Veredicto: `rechaza_h1`. La cláusula de vía no se recalculó con P frente a N y queda `no_evaluable`.

## 2.2 Independencia de PXD034575, PXD019504 y PXD031705

K2 no se ha vuelto a correr en v3.

`docs/reportes/informe_v2.md`, líneas 159 a 163, ya dejó escrito un resultado de v2: la correlación muestra a muestra no llegó a 0,98 (máximo 0,877 en ese informe) y PXD031705 comparte ficheros con PXD019504, así que v2 metió solo a PXD019504 como tercera cohorte.

Eso no cierra la pregunta de si PXD019504 comparte personas con PXD034575. Los tres depósitos son del mismo grupo y los nombres crudos comparten el patrón `CRC<n>`. Hasta repetir los dos niveles (checksums de PRIDE y correlación > 0,98 con al menos 200 proteínas) dentro de v3, el 1 − I² de tres cohortes no se recalcula. En `results/v3/indice_R_v3.tsv`, R3 está `no_evaluable`. El 0,094 de v2 se conserva con su nombre y no se reutiliza.

## 2.3 Citas

`scripts/v3/verificar_citas.py` consultó Crossref el 6 de noviembre de 2026. La caché está en `results/v3/crossref_cache/`. La auditoría final, `results/v3/auditoria_citas.tsv`, cerró con cero errores.

Tres errores eran reales. El resto de la primera pasada eran fallos del extractor y se corrigieron en el mismo script.

PXD036434 se citaba con doi:10.1186/s13073-023-01176-5. Crossref devuelve a Cascianelli, Silvia, 2023, *Genome Medicine* 15:37, título «Multi-label transcriptional classification of colorectal cancer reflects tumor cell population heterogeneity». El depósito corresponde a Praus, Friederike, 2023, *Journal of Translational Medicine* 21:41, doi:10.1186/s12967-022-03855-0, título «Panomics reveals patient individuality as the major driver of colorectal cancer progression».

Álvez 2023, doi:10.1038/s41467-023-39765-y. Crossref: Álvez, María Bueno, 2023, *Nature Communications* 14:4308, título «Next generation pan-cancer blood proteome profiling using proximity extension assay».

Wannberg, doi:10.1038/s41467-025-67688-3. Crossref: Wannberg, Fredrika, fecha de publicación 29 de diciembre de 2025, *Nature Communications* 17:151, título «Plasma protein profiling predicts cancer in patients with non-specific symptoms». El año que se cita es 2025.

El nombre del script de descarga se conserva, con el año que ya usaba v2 en la ruta, porque esa carpeta ya existe.

La corrección de texto está en el commit `101c8f6`.

## 2.4 PXD064655

Una búsqueda en el repositorio el 6 de noviembre de 2026 no encontró la frase «100 sueros daneses, cáncer frente a controles sintomáticos». `Sample_Metadata.txt` y `report.pr_matrix.tsv` no están en el workspace y no se descargaron en esta sesión.

El preregistro ya asigna PXD064655 a un contraste dentro del cáncer (estadio IV frente a I–II), no a especificidad cáncer frente a control (`config/preregistro_v3.yaml`, H8). H8 queda `no_evaluable` hasta tener la matriz y el estado de la fuente.

## 2.5 Configuración

`config/config.yaml`, línea 6, `project.seed: 42`. Líneas 105 y 107: `imputation_method: knn` y `normalize_method: zscore_protein`. Esas claves pertenecen al pipeline antiguo. v3 no las lee. `tests/test_v3_correcciones.py` falla si una salida de v3 apunta a `results/differential/`, `results/go/`, `results/gsea/` o `results/ml/`.

`config/datasets.yaml`, entrada `plasma_crc_2024`, líneas 105 y 106: la accesión declarada es PXD015163 y la matriz descargada es PXD055218. Esa entrada se deja, porque v1 y v2 la usan. v3 lee `serum_pxd055218` (línea 134), accesión PXD055218, la misma matriz, sin nueva descarga.

## 2.6 Debilidades que el jurado va a usar

Suero, PXD055218, 876 proteínas, 20 muestras por grupo. Escala: log2 de intensidad lineal (`results/v3/suero_escala.json`).

Cáncer frente a colonoscopia negativa, treat 0,585 y BH 0,05, con TOST en los mismos márgenes (`results/v3/k4_potencia.json`): 0 diferentes, 145 equivalentes, 731 indeterminadas. La frase de ausencia solo vale para las 145 equivalentes. Las 731 no demuestran ausencia. El alfa efectivo de BH, simulado con las varianzas observadas, es 7,14 × 10⁻⁵.

Cáncer frente a sanos: 3 diferentes, 153 equivalentes, 720 indeterminadas. Las tres, unidas al mapa `metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv`, son SAA1 (P0DJI8, log2FC 3,26), LBP (P18428, log2FC 1,25) y CHMP4A (Q9BY43, log2FC 5,42, 16,8 grados de libertad). Colonoscopia negativa frente a sanos: 1 diferente, APOA4 (P06727, log2FC −1,32).

AUC con signo en `results/v3/k5_auc.tsv`. El nulo del máximo plegado, 10 000 permutaciones de etiquetas, semilla del preregistro (`results/v3/k5_auc.json`): máximo observado 0,9375 en 613 proteínas (KRT5, P13647, y ese AUC ya tiene signo positivo), media nula 0,843, percentil 95 nulo 0,917, b = 244, p = 245/10001 = 0,0245. Con 20 contra 20 el máximo del catálogo nace alto. Ese p controla la búsqueda entre muchas proteínas. No es una prueba diagnóstica.

La diferencia entre el AUC de cáncer frente a colonoscopia negativa y el de colonoscopia negativa frente a sanos está en `results/v3/k5_diferencia_auc.tsv` (759 proteínas con los tres grupos). En 163 el intervalo bootstrap no incluye 0. Esos intervalos no llevan corrección por el número de proteínas. La corrección de la búsqueda es el nulo del máximo.

Faltantes (`results/v3/k6_presencia.tsv`): 358 proteínas por debajo del 70 % de observaciones. La presencia se prueba con Fisher y BH, aparte de la abundancia. Once quedan significativas en presencia. CHMP4A tiene 13/20 presentes en cáncer, 4/20 en colonoscopia negativa y 3/20 en sanos. Su p de Fisher ajustado es 0,171: la abundancia sobre los observados y la presencia no dicen lo mismo, y no se imputa.

Control técnico (`results/v3/k7_control_tecnico.json`). Kruskal–Wallis exploratorio a 0,05. Profundidad, mediana de intensidad, hemólisis y plaquetas no difieren. La fracción de albúmina e inmunoglobulinas, sumada en intensidad lineal, sí difiere: medianas 0,253 (cáncer), 0,323 (colonoscopia negativa) y 0,327 (sanos), p = 1,23 × 10⁻⁵. Entran 8 filas de ALB, IGHG2, IGHG3, IGHG4, IGKC e IGLC2. No están IGHA1, IGHA2, IGHG1, IGHM ni IGLC1. La coagulación también difiere (p = 0,0063): la mediana en log2 es 19,28 en cáncer, 18,91 en colonoscopia negativa y 19,12 en sanos. Se declara confusión técnica en esos dos indicadores. Es una pantalla, no un descubrimiento.

`DIA-NN_report.tsv` no está en `data/raw/plasma`. El número de precursores por corrida queda `no_disponible`. El orden de corrida no está en la matriz ancha y queda como limitación declarada.

## Lo que esta auditoría no cierra

H5 a H10 siguen `no_evaluable` en `results/v3/hipotesis_v3.tsv`. No hay número aproximado. LinkedOmicsKB, PXD031705, PXD064655 y el resto de la sección 15 del plan no se han abierto como datos de v3.

## v3.3, arranque

El encargo dice que `35fad60`, `50caef8` y `c2be8e4` no llevaron prefijo. El mensaje de `35fad60` es `feat: la corrida v3.2 queda en el repositorio`. Gana el commit. `50caef8` (`resultados`) y `c2be8e4` (`update doc`) sí están sin prefijo. No se reescriben. Quedan en `docs/cronologia_v3.md`.

El 9 de noviembre de 2026 el autor pidió desplegar la web si el checklist queda en verde, guardar ProtT5-XL en el disco externo y descargar los 60 crudos al disco externo. El disco del sistema tiene 14 GB libres, por debajo de los 20 GB del valor por defecto. El externo tiene 464 GB.

## v3.2

El 6 de noviembre de 2026 el report de DIA-NN coincidió con el SHA-256 `A8A5AADD74391AFEF68CF4D6C3A7D3860023101E5C92590E2E87C152A39E59E8` y la corrida terminó. Las tablas están en `results/v3_2/`. La frase del suero, en `suero_resumen.json`, es «No se detecta diferencia; 105 proteínas son equivalentes», sobre 205 proteínas. El tejido, en `k1_etiquetas.json`, es 89/90. El bloque de abajo es lo que se comprobó antes y durante esa corrida.

Comprobado en esta pasada:

- El parser de `src/oncogo_crc/utils/dia_nn_plasma.py` pivotaba `Genes.MaxLFQ` por `Protein.Group` y `postprocess_dia_nn_plasma_download` escribía encima del propio fichero (línea 60).
- `results/v3/k7_control_tecnico.json` dice que el report no está en el workspace. Al empezar v3.2 seguía sin estar.
- La librería universal de Frankenfield et al. 2022 se bajó del commit `bddf6da9334f4300caed75f7a87c096dd96b6bca`. SHA-256 `c0782146be6dad6ecf53e23e5a21a7ed4fea9b4802cb01e8d00091404ec66926`. 381 entradas. Incluye KRT8, KRT18, KRT19 y KRT20 además de las de piel. La enmienda no las saca a mano.
- La tabla suplementaria de Geyer et al. 2019 no se obtuvo. `https://www.embopress.org/doi/suppl/10.15252/emmm.201910427/suppl_file/emmm201910427-sup-0001-figs.pdf` devolvió 404 el 6 de noviembre de 2026. Los paneles usan la lista del texto, escrita en la enmienda.
- El intérprete del entorno es Python 3.14.3. El prompt dice 3.11. Gana el fichero real: el sello escribe 3.14.3.
- `make manifest-verify` coincidió (11 entradas, 31 de v1). `make regresion-v1` pasó: 688, 140, 28, 104/119, 226 y 307.
- Los 60 crudos no se descargan. En esta pasada la web no se despliega. Después, a pedido del autor, la página pública quedó en https://oncogo-crc-tesis.web.app. El preregistro no se reescribe.
- Lista de ficheros de PRIDE, 6 de noviembre de 2026, `external/pride/PXD055218/files.json`: 61 ficheros. 60 son `.raw` de Thermo (`CP_`, `CN_`, `HC_`) y el otro es `DIA-NN_report.tsv` de 213 313 241 bytes. No hay `.fasta` ni `report.log.txt`. En el report parcial, las corridas se llaman `CP_`, `CN_` y `Normal_`, no `HC_`. El nombre del crudo y el de la corrida no coinciden. La base de búsqueda sigue sin fichero.
- UniProt humano revisado con isoformas, descarga del 6 de noviembre de 2026: `external/uniprot/human_reviewed_isoforms.fasta.gz`, 42 560 entradas, SHA-256 `96a0fc42f6dbba02386eaae1d0241407eef335876bf0ddc6b1aa91144cb961ae`. Sirve para el ortólogo compartido. No es la base con la que los autores buscaron PXD055218.
- El atlas de enfermedades del HPA (`blood_pea_diseases.tsv.zip`, 6 de noviembre de 2026) trae 60 enfermedades, no 59. Gana el fichero. El colorrectal está en el catálogo con la clase Cancer.
- El 6 de noviembre de 2026, con el report completo y el SHA-256 coincidente, `humanos_pg_y_global` es 312, no 307. Los otros controles coinciden: 283 488 filas, 60 corridas, 329 con `Global.PG.Q.Value` ≤ 0,01 y 432 que no pasan `PG.Q.Value` en ninguna corrida. Los 5 de diferencia son P0DOX2, P0DOX3, P0DOX5, P0DOX7 y P0DOX8: humanos, pasan los dos q-valores y la columna `Genes` está vacía. Con gen, el recuento es 307. El análisis principal indexa por `Protein.Group`, así que esos cinco grupos se quedan.

## v3.3, mapa de genes de PXD034575

`Cohort1_proteinGroups.txt` tiene 420 columnas. Las que usa el mapa son `Protein IDs`, `Gene names` y `Fasta headers`. Los cinco primeros `protein_id` de `pride_crc_matrix.tsv` son D2Y6Y7, X2C0Y9, Q4R1L7, A0A141ZSL4 y A0A1B1PEM6. Ninguno tiene forma de Swiss-Prot revisado. Por eso el cruce viejo, acceso exacto contra `gene_uniprot_human_reviewed_bulk.tsv`, dejaba fuera a la mayoría.

La matriz tiene 6824 filas. Con el cruce viejo, 2415 tienen el primer acceso en el mapa revisado. Con el mapa nuevo, 3941 tienen un solo gen, 2688 tienen más de uno y 195 no tienen ninguno. Una fila con varios genes no entra en R1. La columna de pares es `n_obs`.

CA1 está en el grupo V9HWE3, con genes HEL-S-11 y CA1, y el acceso P00915 va dentro del grupo. SLC26A3 está en P40879, junto al símbolo viejo DKFZp686P10213. AQP8 (O94778) y GUCA2A (Q02747) no aparecen en `Protein IDs`, ni en `Gene names`, ni en `GN=` de las cabeceras, ni en la matriz. No es un fallo del mapa: esas dos proteínas no están en la salida de MaxQuant de esta cohorte. El control las escribe como falla. La corrida no las va a dar por buenas.

## v3.3, H5 y la columna Control del HPA

UCEC tiene 95 columnas de tumor (`C3N-00858`, `C3L-01252`, …) y 18 de normal (`NX1`, `NX2`, …). La intersección es 0. No es el mismo participante con otro sufijo, así que no se normaliza el identificador. UCEC queda excluido. COAD tiene 96 pares. CCRCC 80, HNSCC 62, LSCC 99, LUAD 101, OV 10 y PDAC 44. Con el mínimo de 8 pares hay 6 carcinomas además del colon.

La columna `Control` del atlas trae exactamente `Healthy`, `Class` y `All other diseases`. `Class` contiene la subcadena `all`, así que un mapeo por subcadena la mandaría al control equivocado. El mapeo explícito es Healthy → sanos, Class → misma_clase, All other diseases → todas_las_demas. Un texto distinto detiene la etapa.

## v3.3, recursos de H7

Álvez 2023 está en `external/alvez_2023/alvez_2023_colorrectal_por_proteina.tsv`. Columnas: proteina, olink_id, diferencia_npx, p_ajustado, significativa, direccion, etiqueta_del_articulo. El contraste es colorrectal en el panel de Olink. La etiqueta de H7 es 1 si p_ajustado < 0,05 y el signo coincide con el de CPTAC.

El atlas del HPA usa el contraste de cáncer colorrectal frente a `Healthy`. La concentración sale de `external/hpa/hpa_concentracion_sangre.tsv`, columna `conc_sangre_ms_pg_l`, en log10. `Secreted to blood` es la anotación de secreción. No hay una columna de especificidad tisular: esa variable queda faltante y no se sustituye.

Papier 2024 no está en el disco. H7 se evalúa con dos recursos, y la regla de al menos dos exige los dos. El solapamiento de casos entre Álvez y el HPA no está comprobado.

ESM-2 y TMbed no entran en el modelo principal. La rejilla de regularización tiene 10 valores desde 0,1 hasta 10 000. El entero `Cs=10` de scikit-learn empieza en 10⁻⁴ y, cuando todas las C empatan en AUC, se queda la más pequeña: las probabilidades quedan en 0,5 y el AUC agregado entre pliegues se desarma.

## v3.4, H5

`Fraccion` no tiene el atributo `ci_low`. El intervalo sale de `intervalo()`, que es Clopper–Pearson. El resto de `concordancias_mismo_universo` y de `submuestreo` sigue el encargo.

## v3.5, comprobaciones del 10 de noviembre de 2026

Los manifiestos de `results/v3_2/`, `results/v3_3/` y `results/v3_4/` coinciden con los ficheros: 0 hashes distintos y 0 extras. `make regresion-v1` pasa y conserva 688, 140, 28, 104/119, 226 y 307. El SHA-256 de `data/raw/plasma/DIA-NN_report.tsv` es `a8a5aadd74391afef68cf4d6c3a7d3860023101e5c92590e2e87c152a39e59e8`.

`report.pr_matrix.tsv` de PXD064655 pesa 6300706 bytes, SHA-256 `b51a062c23a99cf93296ada8ba80061823289cdf0c5c872706e0e41bd005b025`. `Sample_Metadata.txt` pesa 3961 bytes, SHA-256 `c53e5d279118fff2b0a4891ed9873cafb12f4baec503ee2769522fdf36274862`.

Los JSON de ThermoRawFileParser en el disco se llaman `CN_02_DIA.json`, no `<crudo>-metadata.json`. El campo de fecha es `FileProperties` → `Content Creation Date`, por ejemplo `10/02/2023 04:14:10`. El instrumento leído en `CN_02_DIA.json` es Orbitrap Fusion, serie `FSN10392`.

Hay 60 crudos cuyo nombre cumple `(CP|CN|HC)_NN_DIA.raw` y 60 ficheros AppleDouble `._*`. Un glob `*.raw` cuenta 120. Los seis rechazados por el lector tienen el mismo tamaño que la API de PRIDE: CN_01 1049764670, CN_08 1000697970, CN_16 1205043346, CP_18 1226017242, HC_17 1065175806, HC_20 1088171988. El tamaño coincide, así que la descarga no quedó corta. El lector de Thermo, con `-l 1` sobre CN_01, escribió `RAW file is likely corrupted` y `OutOfMemoryException` al cargar MsDevice en el offset 973984150. Hoy no se vuelven a descargar.

Wannberg está en `external/wannberg_2026/wannberg_2026_tabla_s5.tsv`.

## v3.6, comprobaciones del 10 de noviembre de 2026

ThermoRawFileParser disponible: `external/bin/trfp/ThermoRawFileParser`, versión `2.0.0.0`. No aparece otra instalación. La v3.6 hará un intento con esta versión sobre los seis crudos sin JSON, primero `-m 0` y, si falla, `-m 1`. No se descargan.

Inventario de funciones que leen una corrida anterior, para el refactor de la etapa limpia:

- `stages/v3_5.py`: copia `results/v3_4`, lee `results/v3_3` para el ML y para R1, y abre los manifiestos de v3.2 a v3.4.
- `v3/h8.py`: `results/v3_3/k1_primario_frente_a_mucosa.tsv` y, si falta en la salida, `results/v3_4/h6_origen_celular.tsv`.
- `v3/triangulacion.py`: K1 y K4 de v3.3, y H5 y H10 de v3.4.
- `stages/v3_5_preanalitica.py`: `results/v3_3/suero_calidad_muestra.tsv`.
- `v3/salida.py` `entrada_v3`: si el fichero no está en la salida, abre `results/v3_3`. Lo usan H5, H6 y H7, y el multiverso de v3.4.
- `v3/secuencia_features.py`: escribe `tmbed_fuente` y `esm2_fuente` como ruta absoluta.

`k1_primario_frente_a_mucosa.tsv` tiene 6587 filas. El código de H8 usa esa tabla entera, no el subconjunto de 147 concordantes. El embudo de v3.6 contará por qué el universo del suero es más pequeño.

matplotlib 3.10.9 y Jinja2 3.1.6 ya están en el entorno. No se instala nada.
