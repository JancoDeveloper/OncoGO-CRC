# Evidencia externa e in silico

Esta etapa (`external_evidence`) **no sustituye** experimentos de laboratorio: integra el diferencial núcleo con **capas YAML** (`crc_context_layers.yaml`), tablas **opcionales** de referencia transcriptómica (TCGA u otras) y **Cox por gen** si aportas un TSV (p. ej. export GEPIA / cBioPortal).

## Qué se genera

1. **`evidence_matrix.tsv`**: cada fila = feature del contraste `external_evidence.primary_comparison`; columnas `hit_*` = pertenencia a listas de la revisión (tiers, CMS, módulos funcionales, eje MAPK/plasma). Si existe traslación clínica, se fusionan columnas de scoring (`translational_priority_score`, etc.).
2. **`tcga_rnaseq_concordance.tsv`** (si config apunta a un TSV con `gene_symbol` y `log2fc_rnaseq` o alias): unión por gen con el `log2fc` proteómico; columnas de concordancia de signo.
3. **`survival_gene_cox_merged.tsv`** (si config apunta a `gene_symbol`, `cox_hr`, `cox_p`): unión con el diferencial.
4. **`external_evidence_summary.json`** + **`external_evidence_summary.tsv`**: métricas agregadas (Spearman signo, n matched, …).

## Resumen de la última corrida

```json
{
  "primary_comparison": "cptac_tumor_vs_normal",
  "n_rows_evidence_matrix": 120,
  "n_annotation_columns": 13,
  "tcga_rnaseq": {
    "status": "skipped_no_path"
  },
  "survival_cox": {
    "status": "skipped_no_path"
  }
}
```

## Extracto — matriz de evidencia (capas YAML internas)

| protein_id | log2fc | t_stat | p_value | stat_method | padj | significant | comparison_id | matrix_input | hit_tier_tier1_validated_clinical | hit_tier_tier2_emerging | hit_tier_tier3_pathway_anchor | hit_cms_CMS1_MSI_immune | hit_cms_CMS2_canonical_WNT | hit_cms_CMS3_metabolic | hit_cms_CMS4_mesenchymal_ECM | hit_module_EMT_ECM_invasion | hit_module_metabolism_lipid_glycolytic | hit_module_DNA_repair_damage_response | hit_module_angiogenesis_hypoxia | hit_module_antigen_presentation_immune | hit_hallmark_MAPK_EGFR_axis | evidence_score | combined_score | evidence_tier | translational_priority_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GENE001 | 0.9235921277208908 | 6.304887075924769 | 1.125820612595871e-07 | welch_ttest | 4.5032824503834825e-06 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.9677083333333332 | 0.895 | tier_A_cuantitativo_fuerte | 0.6761666666666667 |
| GENE002 | 0.909495664070206 | 6.385451446560606 | 9.086670483674025e-08 | welch_ttest | 4.5032824503834825e-06 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.8910416666666665 | 0.9041666666666668 | tier_A_cuantitativo_fuerte | 0.6155 |
| GENE003 | 0.3871333247671575 | 1.9854245051103492 | 0.0537739688714026 | welch_ttest | 0.2581150505827326 | False | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.6818749999999999 | 0.3370833333333333 | tier_B_cuantitativo_moderado | 0.4558333333333333 |
| GENE004 | 0.890719220493756 | 4.961971017399509 | 9.987660410227048e-06 | welch_ttest | 0.0001089562953842 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.8627083333333332 | 0.92125 | tier_A_cuantitativo_fuerte | 0.5988333333333333 |
| GENE005 | 0.4586332723064352 | 2.516241222670177 | 0.0155239036789905 | welch_ttest | 0.0887080210228031 | False | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.7877083333333333 | 0.8279166666666666 | tier_B_cuantitativo_moderado | 0.5475000000000001 |
| GENE006 | 0.4346464936536058 | 2.1899342191594164 | 0.0341724543034385 | welch_ttest | 0.1782910659309834 | False | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.7735416666666666 | 0.8104166666666667 | tier_B_cuantitativo_moderado | 0.5378333333333334 |
| GENE007 | 0.9860413371216624 | 5.0900551932527 | 6.570769374679962e-06 | welch_ttest | 7.884923249615955e-05 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.868125 | 0.8270833333333334 | tier_A_cuantitativo_fuerte | 0.5976666666666667 |
| GENE008 | 0.7488073431992144 | 3.789292696818511 | 0.0004393480854779 | welch_ttest | 0.0029289872365197 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.8035416666666666 | 0.7479166666666667 | tier_A_cuantitativo_fuerte | 0.5551666666666666 |
| GENE009 | 1.202307290707827 | 6.57097891681167 | 5.93461308610679e-08 | welch_ttest | 4.5032824503834825e-06 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.9047916666666664 | 0.8783333333333334 | tier_A_cuantitativo_fuerte | 0.623 |
| GENE010 | 1.041265054134849 | 5.408041305382885 | 2.2880197419963825e-06 | welch_ttest | 4.576039483992765e-05 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.8835416666666666 | 0.8616666666666667 | tier_A_cuantitativo_fuerte | 0.6085 |
| GENE011 | 0.7983035393362492 | 3.975448977646191 | 0.0002532874450346 | welch_ttest | 0.0018996558377595 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.819375 | 0.77 | tier_A_cuantitativo_fuerte | 0.5655 |
| GENE012 | 0.7071579084285715 | 4.285613904820984 | 9.416389630266764e-05 | welch_ttest | 0.0007533111704213 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.813125 | 0.7887500000000001 | tier_A_cuantitativo_fuerte | 0.5623333333333334 |

## Concordancia RNA-seq externa (TCGA u otro)

_No cargada en esta corrida._ Para activarla, fija `external_evidence.tcga_rnaseq_log2fc_tsv` (formato: columnas `gene_symbol`, `log2fc_rnaseq` o alias). Mientras tanto, **no se publican extractos**: cualquier `tcga_rnaseq_concordance.tsv` previo se borra al inicio de la etapa para no arrastrar plantillas.


## Supervivencia externa (Cox por gen)

_No cargada en esta corrida._ Fija `external_evidence.survival_gene_cox_tsv` con columnas `gene_symbol`, `cox_hr`, `cox_p` (export GEPIA / cBioPortal). Sin TSV real, no se renderiza tabla.


---

**Cómo preparar las tablas de referencia:** coloca TSV bajo `metadata/tcga_auxiliary/` (o ruta absoluta) y fija `external_evidence.tcga_rnaseq_log2fc_tsv` / `survival_gene_cox_tsv`. Los archivos `*_template.tsv` en el repo son **plantillas de columnas**; sustituye valores por datos curados de tu fuente (GDC, Xena, GEPIA, etc.).
