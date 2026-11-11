> Documento de la v1. No forma parte de los resultados certificados de la tesis; ver results/v3_7.

# Síntesis para tesis

> **Workspace de datos:** `data/`, `results/`, `reports/`, `metadata/` y `docs/web/` apuntan a `/Volumes/JANCODATA/OncoGO-CRC` (el clon conserva `config/`, `src/` y el resto de `docs/`). Trazabilidad de corrida: `results/run_manifest.json`.

## Qué hace este pipeline (en una frase)

Integra matrices proteómicas públicas (CPTAC, validación PRIDE, plasma), normaliza e imputa, contrasta grupos (tumor vs normal y comparaciones opcionales MSI/EOCRC si hay datos), interpreta con **ORA + GSEA pre-ranked** sobre Gene Ontology offline (mismo `go_gene_long.tsv`), cuantifica **firma núcleo CPTAC ∩ PRIDE**, resume redes, prioriza con ML explicable y exporta un **bundle** listo para capas contextuales (STRING, UniProt, AlphaFold, ProteinsLanguage).

## Arquitectura conceptual (3 niveles)

1. **Nivel 1 — OncoGO-CRC:** pipeline estadístico reproducible (diferencial, ORA, GSEA, KEGG opcional, ML, redes, estabilidad).
2. **Nivel 2 — ProteinsLanguage / recursos externos:** GO profundo, estructura, embeddings, grafos, visualización (consume `results/export/protein_context_bundle.json` + tablas en `results/`).
3. **Nivel 3 — Traslación clínica:** CMS / MSI / ECM / plasma / EOCRC vía metadata y `crc_context_layers.yaml` (interpretación siempre acotada a lo medido o inferido con criterio explícito).

## Enfoque metodológico: ORA, GSEA y señal distribuida

En proteómica, **listas muy grandes o heterogéneas** pueden hacer que la **ORA (sobrerrepresentación binaria)** sea **poco informativa** (pocos o ningún término GO con FDR estricto), **sin que ello invalide** el contraste diferencial. Por ello el proyecto **prioriza también una estrategia basada en ranking continuo (GSEA pre-ranked)** y en **reducción semántica** de términos GO (`*_go_reduced.tsv`), para capturar **cambios distribuidos a nivel de proceso**. Se reportan **ORA y GSEA juntos**: la ORA describe conjuntos discretos bajo reglas explícitas; el GSEA explora si la ranking completa enriquece vías cuando la ORA colapsa. En la memoria se puede enmarcar un **0 términos ORA** como **decisión metodológica defendible** (lista no específica / universo próximo al foreground), no como “fallo”.

## Preguntas de investigación: módulos funcionales, no solo “proteína top”

Además de identificar proteínas punteras, el marco sugiere preguntar **qué procesos o módulos son consistentes** entre cohortes o con la literatura (p. ej. EMT/ECM, metabolismo lipídico/glicolítico, MAPK/EGFR, inmunidad/MSI-like, reparación del ADN, angiogénesis/hipoxia). Las etiquetas `functional_module_hints` y `cms_pathway_hints` en `results/clinical_translation/*.tsv` orientan esa lectura; la síntesis clínica debe citar **García-Foncillas *et al.*, 2026** y tu revisión narrativa donde corresponda.

## Validación externa y transferencia entre cohortes

Diseño explícito: **descubrimiento en CPTAC**, **validación de dirección en PRIDE**, **apoyo traslacional en plasma**. El resumen `results/reproducibility/cptac_pride_differential_overlap_summary.tsv` incluye **correlación de `log2fc`** (Pearson/Spearman) y **fracción de mismo signo** entre cohortes unidas por **UniProt canónico**. Con `ml.external_transfer.enabled`, el pipeline escribe **`results/ml/*_train_*_test_transfer.json`**: AUC/accuracy en PRIDE tras alinear **genes** (post-`harmonize`) y entrenar en CPTAC. Un **rendimiento ML alto dentro de CPTAC** pero **transferencia modesta a PRIDE** es defendible científicamente: subraya batch/tecnología y la necesidad de **priorizar señales biológicas robustas** frente a clasificadores sobreajustados.

## ML: ranking estable y explicable (no “el mejor AUC”)

Los modelos **L1**, **Elastic Net (SGD)** y **RF** con métricas de hold-out sirven para **ordenar genes con evidencia multivariante**, no para vender un clasificador cerrado. Con `ml.nested_cv.enabled`, se generan **`results/ml/<cohorte>_nested_cv_per_fold.tsv`** (AUC por pliegue externo; `GridSearchCV` interno para `alpha` o hiperparámetros de RF). Opcional: **`ml.feature_stability_bootstrap_n`** para frecuencias bootstrap del top‑k en L1. La narrativa honesta favorece **nested CV + transferencia** frente a un único split optimista.

## Translational Priority Score (TPS)

La tabla de traslación incluye **`translational_priority_score`**: combinación explícita de magnitud/FDR, red, ML, contexto CMS/MAPK, módulos funcionales y capas literarias (`crc_context_layers.yaml`). Sirve como **marco cuantitativo** para priorizar hipótesis, no como decisión clínica.

## Druggability e hipótesis de reposicionamiento (in silico)

Sin experimentos de laboratorio, el **bundle** JSON incluye plantillas de consulta (p. ej. Open Targets / ChEMBL) para **hipótesis de diana o reposicionamiento**. La redacción debe evitar afirmaciones tipo “este fármaco sirve en CRC”; el formulación correcta es **“priorización de diana / hipótesis generadora”**, alineada con revisiones de farmacoproteómica.

## Validaciones in silico de bajo costo (opcionales)

La etapa **`external_evidence`** (tras `clinical`) escribe **`results/external_evidence/evidence_matrix.tsv`** (columnas `hit_*` desde `crc_context_layers.yaml` + fusión con traslación). Si configuras `external_evidence.tcga_rnaseq_log2fc_tsv` o `survival_gene_cox_tsv` (plantillas en `metadata/tcga_auxiliary/`), se generan **`tcga_rnaseq_concordance.tsv`** y **`survival_gene_cox_merged.tsv`** con Spearman/mismo signo o Cox fusionado. Rutas adicionales típicas: soporte en **TCGA/GEO**, **supervivencia** pública, **meta-análisis direccional** (≥2 cohortes).

## Qué está **demostrado** aquí (evidencia cuantitativa)

- **Tumor vs tejido adyacente normal (CPTAC)** (`cptac_tumor_vs_normal`): **5031** proteínas con FDR por debajo de **0.05** (de **6237** filas en el contraste).
  - GO (ORA, universo = proteínas detectadas en la matriz): **34** términos (FDR término < **0.05**).
  - GSEA (GO offline, ranking continuo; métrica principal `t_stat`): **446** términos con FDR < **0.25**.
  - GSEA (métrica `signed_neglog10_p`): **353** términos con FDR < **0.25** (`results/gsea/cptac_tumor_vs_normal_gsea_go_signed_neglog10_p.tsv`).
- **Plasma CRC vs control sano** (`plasma_crc_vs_healthy`): **119** proteínas con FDR por debajo de **0.05** (de **517** filas en el contraste).
  - GO (ORA, universo = proteínas detectadas en la matriz): **11** términos (FDR término < **0.05**).
  - GSEA (GO offline, ranking continuo; métrica principal `t_stat`): **16** términos con FDR < **0.25**.
  - GSEA (métrica `signed_neglog10_p`): **10** términos con FDR < **0.25** (`results/gsea/plasma_crc_vs_healthy_gsea_go_signed_neglog10_p.tsv`).
- **PRIDE validación externa: tumor vs normal** (`pride_validation_tumor_vs_normal`): **1951** proteínas con FDR por debajo de **0.05** (de **6782** filas en el contraste).
  - GO (ORA, universo = proteínas detectadas en la matriz): **0** términos (FDR término < **0.05**).
  - GSEA (GO offline, ranking continuo; métrica principal `t_stat`): **31** términos con FDR < **0.25**.
  - GSEA (métrica `signed_neglog10_p`): **19** términos con FDR < **0.25** (`results/gsea/pride_validation_tumor_vs_normal_gsea_go_signed_neglog10_p.tsv`).

**Comparaciones opcionales** (solo si `differential.auto_run_optional_comparisons` y datos):

- **Estadio avanzado vs temprano (si stage está disponible)** (`cptac_stage_advanced_vs_early`): **24** proteínas con FDR por debajo de **0.05** (de **6237** filas en el contraste).
  - GO (ORA, universo = proteínas detectadas en la matriz): **0** términos (FDR término < **0.05**).
  - GSEA (GO offline, ranking continuo; métrica principal `t_stat`): **7** términos con FDR < **0.25**.
  - GSEA (métrica `signed_neglog10_p`): **7** términos con FDR < **0.25** (`results/gsea/cptac_stage_advanced_vs_early_gsea_go_signed_neglog10_p.tsv`).
- **Hispano/latino vs no hispano/latino (etnicidad NIH; contexto US — ver clinical_demographics.md)** (`cptac_hispanic_latino_vs_non`): **33** proteínas con FDR por debajo de **0.05** (de **6237** filas en el contraste).
  - GO (ORA, universo = proteínas detectadas en la matriz): **0** términos (FDR término < **0.05**).
  - GSEA (GO offline, ranking continuo; métrica principal `t_stat`): **32** términos con FDR < **0.25**.
  - GSEA (métrica `signed_neglog10_p`): **37** términos con FDR < **0.25** (`results/gsea/cptac_hispanic_latino_vs_non_gsea_go_signed_neglog10_p.tsv`).

**Firma núcleo CPTAC ∩ PRIDE** (clave común **UniProt** vía mapa `gene_symbol`↔`uniprot_accession`; filas con UniProt en ambas cohortes: **1993.0**; doble significativo `padj`: **432.0**). Cruce ingenuo por el mismo texto en `protein_id` suele ser **1993.0** (CPTAC=genes, PRIDE=UniProt). Correlación Pearson log2FC CPTAC–PRIDE (UniProt común): **0.104**; mismo signo: **56%**. Tablas: `results/reproducibility/cptac_pride_differential_overlap.tsv` y `..._summary.tsv`.

---

## Contexto clínico‑epidemiológico (Colombia) y subtipos CMS

La revisión **García-Foncillas *et al.*, 2026** (*Rev. colomb. hematol. oncol.*, DOI [10.51643/22562915.853](https://doi.org/10.51643/22562915.853)) subraya el aumento de CCR en población joven y el papel de clasificaciones moleculares (**CMS**, MSI, MAPK/EGFR, biopsia líquida). OncoGO-CRC incorpora:

- **Edad / EOCRC:** si el Excel clínico PDC trae edad, la metadata CPTAC incluye `age`, `age_bin` (`young`/`older` según `clinical_filters.age_bins`) y la comparación opcional `cptac_eocrc_vs_lorc` (se ejecuta sola si `differential.auto_run_optional_comparisons` y hay estratos con N mínimo).
- **MSI vs MSS:** comparación opcional `cptac_msi_vs_mss`; con `auto_run_optional_comparisons` corre si `msi_status` tiene suficientes MSI y MSS.
- **CMS4 vs resto:** columna `cms_subtype_bin` si aportas un TSV suplementario (`harmonize_public_metadata.cptac.cms_subtype_supplement_tsv`) y la comparación opcional `cptac_cms4_vs_rest`.
- **Traslación:** `metadata/clinical_layers/crc_context_layers.yaml` etiqueta genes priorizados (biomarcadores, pistas CMS, eje MAPK/EGFR para discusión de plasma).

Interpretación sugerida: hubs en matriz extracelular / EMT pueden discutirse como **coherente con fenotipo CMS4** (TGF‑β, angiogénesis) si el contraste y la literatura lo sustentan — citar siempre el subtipo molecular **medido o inferido con criterio explícito** en tesis.

---

## Qué **no** sustituye este documento

- No interpreta mecanismos ni discute contraste con la bibliografía clínica: eso es texto de tesis.
- Los avisos de SciPy («catastrophic cancellation») afectan filas con varianza casi nula; no invalidan todo el contraste.

## Cómo cerrar el capítulo de resultados

1. Lee [Resultados base](results_base.md) y valida que las tablas coinciden con tus preguntas de investigación.
2. Abre los TSV completos en `results/` (misma jerarquía bajo tu `data_workspace` si lo usas).
3. Usa el prompt en [`docs/PROMPT_SINTESIS_TESIS.md`](../PROMPT_SINTESIS_TESIS.md) con un LLM pegando extractos de tablas y el contexto de la cohorte.

## Archivos clave para citar en la memoria

- `results/differential/*_diff.tsv` — estadísticos por proteína (todas las comparaciones corridas).
- `results/gsea/*_gsea_go.tsv` y `results/gsea/*_gsea_go_<métrica>.tsv` — GSEA pre-ranked (p. ej. `t_stat` y `signed_neglog10_p`).
- `results/go/*_go.tsv` — ORA por término (lista significativos + universo detectado).
- `results/reproducibility/cptac_pride_differential_overlap*.tsv` — firma núcleo multi-cohorte.
- `results/export/protein_context_bundle.json` — puente a ProteinsLanguage / STRING / UniProt / AlphaFold.
- `results/ml/<cohorte>_feature_ranking.tsv` — priorización multivariante.
- `results/ml/<cohorte>_nested_cv_*.tsv` — nested CV (AUC por pliegue externo; hiperparámetros por grid interno).
- `results/ml/*_train_*_test_transfer.json` — evaluación **CPTAC → PRIDE** (genes alineados).
- `results/ml/<cohorte>_feature_stability_bootstrap.tsv` — estabilidad de ranking L1 (si `ml.feature_stability_bootstrap_n>0`).
- `results/external_evidence/evidence_matrix.tsv` — matriz de evidencia bibliográfica / módulos.
- `results/external_evidence/tcga_rnaseq_concordance.tsv` — concordancia direccional proteómica vs RNA-seq (opcional).
- `results/validation/input_validation_report.json` — trazabilidad de calidad.
