> Documento de la v1. No forma parte de los resultados certificados de la tesis; ver results/v3_7.

# Diccionario y teoría

La explicación completa, con las tablas y los enlaces, está en [divulgacion.html](divulgacion.html), en las secciones Teoría, Datos, Procesos y Estudios.

Cada DOI de esa página se confirmó en Crossref el 31 de octubre de 2026: el título registrado coincide con la cita. El DOI 10.1016/j.xcrm.2025.101689 devolvió 404 y no se enlaza.

## Cómo está armada la página

1. El resultado, en una frase.
2. Las tres cohortes, cada una con el artículo y el depósito: [Vasaikar 2019](https://doi.org/10.1016/j.cell.2019.03.030) y [PDC000116](https://pdc.cancer.gov/pdc/study/PDC000116); [Tanaka 2024](https://doi.org/10.1016/j.celrep.2024.113810) y [PXD034575](https://www.ebi.ac.uk/pride/archive/projects/PXD034575); [Rai 2025](https://doi.org/10.1021/acs.jproteome.5c00285) y [PXD055218](https://www.ebi.ac.uk/pride/archive/projects/PXD055218).
3. Los números de tejido y de suero.
4. Las fichas, con [QuickGO](https://www.ebi.ac.uk/QuickGO/), el [Human Protein Atlas](https://www.proteinatlas.org/) y [AlphaFold](https://alphafold.ebi.ac.uk/entry/O94812).
5. La teoría de cada procedimiento: por qué no se suman las cohortes, qué es el cero de MaxQuant, el log2, el pareo, la t moderada, el FDR, el AUC, Gene Ontology y el límite de los neoantígenos.
6. El diccionario de datos y el de procesos.
7. Las preguntas para practicar la sustentación.

## Procedimientos, en corto

| Paso | Idea |
|---|---|
| Suero | Tres grupos de Rai. t moderada, umbral 0,585, FDR 0,05, sin imputar. |
| Especificidad | AUC mayor que 0,70, intervalo por encima de 0,50 y mayor que el AUC de enfermedad frente a sanos. La marca de fase aguda no baja el umbral. |
| Tanaka | P primario, M metástasis, N mucosa. El cero es un hueco. Luego log2 y contraste pareado. |
| Concordancia | Mismo signo y correlación contra la tabla certificada de CPTAC. 0,83 y 0,104 son estadísticas distintas. |
| Vías | Unos 12 términos GO del tejido tienen al menos 10 genes en el suero. No superan 200 permutaciones. |
| Fichas | GO, atlas de referencia, UniProt y AlphaFold. La estructura no valida la abundancia. |

Los métodos enlazados son [Smyth 2004](https://doi.org/10.2202/1544-6115.1027), [McCarthy y Smyth 2009](https://doi.org/10.1093/bioinformatics/btp053), [Benjamini y Hochberg 1995](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x), [Cox y Mann 2008](https://doi.org/10.1038/nbt.1511), [Ashburner et al. 2000](https://doi.org/10.1038/75556) y [Subramanian et al. 2005](https://doi.org/10.1073/pnas.0506580102).
