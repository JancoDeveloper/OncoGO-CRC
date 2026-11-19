# Fuentes de datos

Cuatro cohortes, analizadas por separado. Los números de muestras son los
del capítulo 3 de la tesis (tabla 7) y coinciden con lo que el código lee
de los ficheros.

| Accesión | Matriz | Muestras | Papel | Fuente |
| --- | --- | --- | --- | --- |
| PDC000116 | Tejido, TMT-10, razón log2 | 95 tumor, 100 normal, 3 sin tipo | Ancla del cambio tumoral | Vasaikar et al. 2019; Proteomic Data Commons |
| PXD034575 | Tejido, LFQ, MaxQuant | 16 pacientes con primario, metástasis hepática y mucosa | Segunda cohorte de tejido | Tanaka et al. 2024 |
| PXD055218 | Suero, DIA, DIA-NN, Orbitrap Fusion | 20 CP, 20 CN, 20 HC; 60 crudos | Contraste frente a colonoscopia negativa | Rai et al. 2025 |
| PXD064655 | Suero o plasma, DIA-NN 1.8.1, matriz de precursores | 100 muestras: 50 estadio IV, 21 estadio III, 29 estadios I–II | Validación por estadio | Hermoso-Durán et al. 2026; McKay 2026 |

El metadato de PXD064655 dice plasma y el protocolo depositado en PRIDE dice
suero. El código no resuelve esa discrepancia.

## Dónde está cada fichero

| Fichero | URL | Huella sellada |
| --- | --- | --- |
| `Cohort1_proteinGroups.txt` | `https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/01/PXD034575/Cohort1_proteinGroups.txt` | `da9ec08e…` en `data/checksums.sha256`. 69 272 674 bytes |
| `pride_crc_matrix.tsv` | Se reconstruye con `scripts/descargar_datos.py` desde el proteinGroups | `a5baca5a…`. Coincide byte a byte con pandas 3.0.2 |
| `DIA-NN_report.tsv` | `https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/06/PXD055218/DIA-NN_report.tsv` | `a8a5aadd…`. 213 313 241 bytes |
| `report.pr_matrix.tsv` y `Sample_Metadata.txt` | `https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/01/PXD064655/` | No están en el manifiesto sellado. Las huellas calculadas el 19 de noviembre de 2026 están en `data/checksums.sha256` |
| `cptac_crc_matrix.tsv` | Estudio PDC000116 en `https://pdc.cancer.gov/pdc/study/PDC000116` | `a1272baf…`. No hay un único URL de archivo en este repositorio |
| Tres ficheros de `results/certified` y `results/certification` | Regeneración con `make certified`, o el depósito futuro `https://zenodo.org/deposit/PENDIENTE` | En `manifest_entradas.sha256` |

Los crudos de PXD055218 son 60 archivos `.raw`, 65,56 GB en total (entre
687,6 MB y 1,3 GB cada uno, consulta a la API de PRIDE el 19 de noviembre
de 2026). `make crudos-v3-7` los compara con el SHA-1 de `checksum.txt`.

## Anotaciones

| Recurso | Fichero | Huella |
| --- | --- | --- |
| Gene Ontology | `metadata/gene_ontology/go_gene_long.tsv` | `560b33fa…`, en el repositorio |
| UniProt revisado | `metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv` | `2d1619fa…`, en el repositorio |
| Human Protein Atlas, concentración en sangre | `external/hpa/hpa_concentracion_sangre.tsv` | `9e9e3239…`. La genera `scripts/download/human_protein_atlas.py` desde `proteinatlas.tsv.zip` |
| Biblioteca de contaminantes de Frankenfield et al. 2022 | `external/frankenfield_2022/0602_Universal_Contaminants.fasta` | Huella calculada en esta preparación, no sellada |
| Embeddings ESM-2 | `external/modelos/h7_esm2_embeddings.tsv` | `2b1e48e3…`. No está en git. El modelo principal no lo usa |

LinkedOmics no lo lee la etapa v3.7. No forma parte de las entradas
declaradas.

Las licencias de los depósitos son las de PRIDE, ProteomeXchange, PDC, GO,
UniProt y el Human Protein Atlas (CC BY-SA 4.0 para el atlas). Las citas
están en `docs/references.bib`.
