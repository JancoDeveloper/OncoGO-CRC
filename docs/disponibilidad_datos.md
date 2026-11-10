# Disponibilidad de los datos

El código y los resultados de la corrida v3.5 están en este repositorio, en `results/v3_5/`. El manifiesto es `results/v3_5_manifest.sha256`. La licencia del código es MIT.

Los crudos no están en el repositorio. Viven en el disco de trabajo y en los repositorios públicos de origen. Cada cohorte se lee sola.

| Accesión | Uso en esta tesis | Dónde está |
| --- | --- | --- |
| PXD055218 | Suero principal | PRIDE y `data/raw/plasma/` en el disco de trabajo |
| PXD034575 | Tejido de Tanaka | PRIDE |
| PDC000116 | Tejido de CPTAC | PDC y LinkedOmicsKB |
| PXD064655 | Estadio, hipótesis H8 | `external/pxd064655/` en el disco de trabajo |
| PXD033236 | Vesículas | El fichero usado está en el disco de trabajo |
| HPA, Álvez 2023, Papier 2024 | Atlas de sangre | Las tablas usadas están en el disco de trabajo |
| GCST90255676 | Genética, hipótesis H9 | GWAS Catalog |

Los seis crudos de PXD055218 que ThermoRawFileParser no leyó siguen en el disco. Su tamaño coincide con el de PRIDE. No se volvieron a descargar en esta corrida. Una corrida posterior puede cerrar esa prueba si se descargan de nuevo.

Las cifras de la página y de la tesis salen de `results/v3_5/` y, donde la corrida anterior sigue vigente, de `results/v3_3/`. No se editan a mano.
