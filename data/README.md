# Datos

Las matrices públicas no se versionan. `data/checksums.sha256` guarda la huella
SHA-256 que la corrida sellada dejó escrita. `scripts/descargar_datos.py`
descarga lo que tiene URL pública y comprueba esa huella.

| Accesión | Papel | Dónde queda |
| --- | --- | --- |
| PDC000116 | Tejido TMT, ancla del cambio | `data/raw/cptac/cptac_crc_matrix.tsv` y, tras `make certified`, `results/certified/` |
| PXD034575 | Tejido LFQ, segunda cohorte | `data/raw/pride/Cohort1_proteinGroups.txt` y `pride_crc_matrix.tsv` |
| PXD055218 | Suero DIA, contraste CP frente a CN | `data/raw/plasma/DIA-NN_report.tsv` |
| PXD064655 | Suero por estadio | `external/pxd064655/report.pr_matrix.tsv` |

La corrida v3.7 no recalcula el contraste de CPTAC: lo lee de
`results/certified/cptac_unshared_treat.tsv`,
`results/certified/inference_audit.json` y
`results/certification/certificado_resumen.json`. Esos tres ficheros no están en
git. `make descargar` los trae, junto con `data/raw/cptac/cptac_crc_matrix.tsv`
y `external/modelos/h7_esm2_embeddings.tsv`, del depósito
<https://doi.org/10.5281/zenodo.24398807> (registro
<https://zenodo.org/records/24398807>, licencia CC-BY-4.0) y comprueba el
SHA-256. La matriz CPTAC se deriva del estudio de acceso abierto PDC000116
(Vasaikar et al. 2019) y se redistribuye con atribución. `make certified`
puede regenerar los tres ficheros certificados cuando esa matriz ya está en
el workspace y coincide con la huella.

Los sesenta crudos `.raw` de PXD055218 suman 65,56 GB en PRIDE. `make crudos-v3-7`
los comprueba por SHA-1. La tabla sellada, con las fechas de adquisición, está en
`results/v3_7/crudos_integridad.tsv`.
