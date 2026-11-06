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

## v3.2

El 6 de noviembre de 2026 el prompt de v3.2 describe un report de DIA-NN sin filtro de proteína, especies no humanas y queratinas de piel. Esas cifras se vuelven a medir cuando `data/raw/plasma/DIA-NN_report.tsv` esté en disco y su SHA-256 coincida con `A8A5AADD74391AFEF68CF4D6C3A7D3860023101E5C92590E2E87C152A39E59E8`. Hasta entonces no se copian como resultado de esta corrida.

Comprobado en esta pasada:

- El parser de `src/oncogo_crc/utils/dia_nn_plasma.py` pivotaba `Genes.MaxLFQ` por `Protein.Group` y `postprocess_dia_nn_plasma_download` escribía encima del propio fichero (línea 60).
- `results/v3/k7_control_tecnico.json` dice que el report no está en el workspace. Al empezar v3.2 seguía sin estar.
- La librería universal de Frankenfield et al. 2022 se bajó del commit `bddf6da9334f4300caed75f7a87c096dd96b6bca`. SHA-256 `c0782146be6dad6ecf53e23e5a21a7ed4fea9b4802cb01e8d00091404ec66926`. 381 entradas. Incluye KRT8, KRT18, KRT19 y KRT20 además de las de piel. La enmienda no las saca a mano.
- La tabla suplementaria de Geyer et al. 2019 no se obtuvo. `https://www.embopress.org/doi/suppl/10.15252/emmm.201910427/suppl_file/emmm201910427-sup-0001-figs.pdf` devolvió 404 el 6 de noviembre de 2026. Los paneles usan la lista del texto, escrita en la enmienda.
- El intérprete del entorno es Python 3.14.3. El prompt dice 3.11. Gana el fichero real: el sello escribe 3.14.3.
- `make manifest-verify` coincidió (11 entradas, 31 de v1). `make regresion-v1` pasó: 688, 140, 28, 104/119, 226 y 307.
- Los 60 crudos no se descargan. La web no se despliega.
