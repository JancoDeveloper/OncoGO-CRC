# OncoGO-CRC — Resultados actuales y limitaciones explícitas para la tesis

Este documento resume **lo que el pipeline OncoGO-CRC ya puede afirmar con evidencia reproducible** y **lo que queda explícitamente como trabajo en progreso o limitación metodológica**. Está pensado como guía honesta para redactar el capítulo de resultados / discusión de la tesis sin sobre-interpretar.

> **Encuadre:** la propuesta enviada al Dr. Bonilla (oncólogo investigador) y a García-Foncillas et al. (2026) enfatiza un trabajo **estrictamente investigativo y computacional**. Los resultados se interpretan como **hipótesis generadoras**, no como biomarcadores listos para uso clínico — alineado con las recomendaciones de prudencia traslacional en estudios ómicos (Nature 2013, 12564) y con el meta-análisis de reuso de proteómica pública de CCR (Robles et al., 2024, *PLOS Comp Biol*).

---

## 1. Lo que sí se puede afirmar (con evidencia en `results/`)

### 1.1 Descubrimiento diferencial tisular robusto en CPTAC
- `results/differential/cptac_tumor_vs_normal_diff.tsv` reporta proteínas diferenciales con `padj < 0.05` y separación clara en PCA (`results/qc/cptac_colorectal_pca_log2norm.tsv`).
- El contraste `tumor vs normal` en tejido CPTAC se considera el **núcleo de descubrimiento**.

### 1.2 Reproducibilidad parcial CPTAC ∩ PRIDE (multi-cohorte tisular)
- Tras la armonización con el mapa efectivo (`protein_mapping_master.tsv` + bulk UniProt humano revisado en `gene_uniprot_human_reviewed_bulk.tsv`), el cruce por puente UniProt arroja:
  - **`n_cptac_with_join_uni` ≈ 6 000**, **`n_pride_with_join_uni` ≈ 6 700**, **intersección ≈ 1 967 genes puente**, **≈ 430 doble-significativos `padj < 0.05`**, **`same_sign_fraction` ≈ 0.56**.
- Lectura prudente: existe una **firma compartida parcial** entre cohortes tisulares; la **concordancia de magnitud es modesta** (Pearson ≈ 0.10, Spearman ≈ 0.12), lo cual es consistente con la heterogeneidad esperada entre TMT (CPTAC) y LFQ/DIA (PRIDE) entre cohortes de origen distinto.
- Esto está alineado con Robles et al. (2024) sobre reuso multi-cohorte: la **curación de IDs y metadata es el paso que más cambia la posibilidad de integrar PRIDE/CPTAC**.

### 1.3 Capa plasmática como señal exploratoria
- `plasma_crc_vs_healthy_diff.tsv` mantiene una señal más corta pero plausible.
- Lectura honesta: se reporta **estrictamente como hallazgo exploratorio**, alineado con la literatura plasmática reciente (Nat Commun 2024, 44911-1; Cell Rep Med 2025, 101689) que recomienda interpretación prudente de proteómica en sangre.

### 1.4 Priorización interna ML (CPTAC)
- `results/ml/cptac_colorectal_nested_cv_*.tsv` reporta AUROC anidada (nested CV con 5×outer / 3×inner), bootstrap de estabilidad L1 y ranking SHAP/RF (`combined_score`).
- `results/clinical_translation/cptac_colorectal_translation_ranked.tsv` integra capas YAML curadas (`crc_context_layers.yaml`) y produce un **score de priorización traslacional (TPS)** — explícitamente como ranking investigativo de hipótesis.

### 1.5 Trazabilidad reproducible
- `results/run_manifest.json` registra `run_id`, timestamp, etapas ejecutadas/omitidas y hashes SHA-256 de archivos clave.
- `publish_results` y `thesis_web_digest` **leen este manifest** y rechazan publicar resultados de etapas omitidas en la corrida actual (gating estricto).

---

## 2. Lo que NO se debe afirmar todavía como cerrado

### 2.1 GO/GSEA — módulos funcionales del CCR
- Tras descargar la anotación GOA humana real (`scripts/fetch_go_human_annotation.py` → ~334 000 anotaciones, ~38 800 genes, ~17 900 términos), la **infraestructura GO/GSEA es funcional**.
- **Pendiente:** validar en la próxima corrida no-demo que al menos `cptac_tumor_vs_normal` produzca **conjuntos GO significativos** (`NES`/`FDR < 0.25`). Los audits (`results/gsea/<comp>_gsea_input_*.json` y `results/go/<comp>_go_input_audit.json`) ahora reportan `n_ranked_genes_total`, `n_ranked_genes_with_go`, `n_terms_pre_filter`, `n_terms_post_filter`, `top_missing_genes_from_go` para diagnóstico cuantitativo.
- **Cómo redactar:** "los módulos funcionales se reportarán tras validación de la cobertura GO real en la próxima corrida no-demo; los audits cuantitativos están disponibles para cualquier término que no pase el filtro de tamaño".

### 2.2 Transferencia ML CPTAC → PRIDE
- El bug histórico (`common=0`) está resuelto: ahora `harmonize_protein_table` hace doble pasada (merge directo + lookup UniProt→gene) y `normalized_matrix_to_samples_x_genes` detecta cuando un `gene_symbol` post-harmonize aún parece UniProt y lo traduce. La intersección de features debe ser del **orden de cientos a algunos miles** en datos reales.
- **Pendiente:** correr el pipeline real (no-demo) y verificar que `results/ml/cptac_colorectal_pride_crc_meta_feature_intersection.tsv` muestre `common > 0` (idealmente > 500) y `*_train_*_test_transfer.json` reporte AUROC plausible.
- **Cómo redactar:** "la validación externa por transferencia ML quedó habilitada metodológicamente tras alinear el espacio canónico `gene_symbol`; las métricas finales corresponden a la última corrida reproducible y se reportarán con `run_id` explícito".

### 2.3 Evidencia externa transcriptómica / supervivencia
- Mientras `external_evidence.tcga_rnaseq_log2fc_tsv` y `survival_gene_cox_tsv` estén `skipped_no_path`, el reporte web **declara explícitamente que no hay evidencia externa cargada** y borra los TSV previos para no arrastrar plantillas.
- **Pendiente:** cargar tablas reales (GDC/Xena para TCGA COAD/READ; GEPIA/cBioPortal para Cox por gen) bajo `metadata/tcga_auxiliary/` y referenciarlas en `config.yaml`.
- **Cómo redactar:** "la concordancia con RNA-seq externa y supervivencia se reporta como **trabajo en progreso**: el pipeline tiene la etapa habilitada pero requiere TSV curados externos, cuya carga se documenta en `docs/proceso_datos_y_pipeline.md`".

---

## 3. Limitaciones metodológicas explícitas

1. **Heterogeneidad tecnológica entre cohortes.** CPTAC (TMT) y PRIDE (LFQ/DIA) miden abundancia con escalas y errores distintos; por eso la concordancia de magnitud (Pearson/Spearman de log2FC) es modesta. La interpretación se centra en **signo y reproducibilidad cualitativa**, no en magnitudes absolutas.
2. **Cobertura clínica fragmentada.** Solo CPTAC trae metadata clínica consolidada (edad, estadio, MSI, CMS opcional via supplement). PRIDE y plasma se manejan con anotación mínima (tumor/normal o crc/healthy) — limitando contrastes opcionales como CMS4 vs resto o respondedor vs no.
3. **Plasma como hipótesis exploratoria.** La señal proteómica en sangre es intrínsecamente más débil que en tejido (Nat Commun 2024, 44911-1). Se interpreta como **generadora de hipótesis**, no como biomarcador validado.
4. **Sin validación experimental.** El pipeline produce rankings computacionales de proteínas y procesos; **no sustituye** validación de laboratorio (Western blot, IHC, ELISA, etc.) ni guías terapéuticas.
5. **GOA / supplementary maps son snapshots.** Los archivos `go_gene_long.tsv` y `gene_uniprot_human_reviewed_bulk.tsv` son descargas estáticas; deben re-generarse periódicamente (`scripts/fetch_go_human_annotation.py`, `scripts/fetch_uniprot_human_reviewed_gene_map.py`) y su versión queda registrada en `run_manifest.json`.

---

## 4. Plan de validación final antes de cerrar el capítulo de resultados

Orden recomendado (alineado con la Fase A-D del informe técnico):

1. **Fase A — outputs estables.** ✅ Hecho. `run_manifest` enforcement + gating de `external_evidence` + auditorías GO/GSEA con campos cuantitativos.
2. **Fase B — interpretación funcional.** Correr `oncogo-crc run --stage gsea` y `--stage go` con la anotación GOA real y revisar `results/gsea/*_gsea_input_*.json` y `results/go/*_go_input_audit.json`. Si `n_ranked_genes_with_go > 1000` y `n_terms_post_filter > 100`, se espera detectar módulos significativos. Si no, documentar cuantitativamente el porqué (típicamente: símbolos sin alias HGNC actualizado o filtros `min_term_size` demasiado estrictos).
3. **Fase C — validación externa por ML.** Verificar `results/ml/*_feature_intersection.tsv` y `*_train_*_test_transfer.json` para AUROC plausible.
4. **Fase D — evidencia externa real.** Cargar tablas TCGA y Cox curadas; re-correr `external_evidence`.

---

## 5. Cómo citar este estado del pipeline en la tesis

- **Para resultados que ya soportan afirmaciones:** cita el `run_id` específico que aparece en `results/run_manifest.json`, junto con los hashes SHA-256 de los archivos clave (también en el manifest).
- **Para limitaciones:** cita este documento (`docs/RESULTADOS_Y_LIMITACIONES_TESIS.md`) y el manifest correspondiente.
- **Para reproducibilidad:** referencia `GUIA_CONFIGURACION.md` y `docs/proceso_datos_y_pipeline.md`.

---

## 6. Referencias clave usadas en la tesis

- Guinney et al. (2015). *Consensus molecular subtypes of colorectal cancer.* Nat Med.
- Zhang et al. (2014). *Proteogenomic characterization of human colon and rectal cancer.* Nature.
- Robles et al. (2024). *Meta-analysis of public proteomics datasets for colorectal cancer.* PLOS Comp Biol. doi:10.1371/journal.pcbi.1011828
- Ferrari et al. Sensitizing cancer cells to immune checkpoint inhibitors by microbiota-mediated upregulation of HLA class I. *Cancer Cell* 2023. doi:10.1016/j.ccell.2023.08.014. Esta línea se había anotado como proteogenómica pan-cáncer de CPTAC. El registro no dice eso, así que no se usa como tal.
- Li et al. Comprehensive proteogenomic profiling reveals the molecular characteristics of colorectal cancer at distinct stages of progression. *Cancer Research* 2024. doi:10.1158/0008-5472.CAN-23-1878.
- Li Y. et al. (2024). Longitudinal plasma proteome profiling of colorectal cancer during cetuximab. *Nat Commun*. doi:10.1038/s41467-024-44911-1. Es respuesta a tratamiento, no un estudio de cribado.
- García-Foncillas J. (2026). doi:10.51643/22562915.853. Revisión narrativa. No es evidencia de un término GO.
- Tanaka A. et al. (2024). Proteogenomic characterization of primary colorectal cancer and metastatic progression. *Cell Reports* 43:113810. doi:10.1016/j.celrep.2024.113810. PXD034575, cohorte 1: 16 tripletes de tumor primario, metástasis hepática y mucosa normal.
