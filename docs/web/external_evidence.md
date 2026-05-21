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
  "n_rows_evidence_matrix": 6237,
  "n_annotation_columns": 13,
  "tcga_rnaseq": {
    "status": "file_not_found",
    "path": "/Volumes/JANCODATA/OncoGO-CRC/metadata/tcga_auxiliary/tcga_coadread_log2fc.tsv"
  },
  "survival_cox": {
    "status": "file_not_found",
    "path": "/Volumes/JANCODATA/OncoGO-CRC/metadata/tcga_auxiliary/tcga_coadread_cox_os.tsv"
  }
}
```

## Extracto — matriz de evidencia (capas YAML internas)

| protein_id | log2fc | t_stat | p_value | stat_method | padj | significant | comparison_id | matrix_input | hit_tier_tier1_validated_clinical | hit_tier_tier2_emerging | hit_tier_tier3_pathway_anchor | hit_cms_CMS1_MSI_immune | hit_cms_CMS2_canonical_WNT | hit_cms_CMS3_metabolic | hit_cms_CMS4_mesenchymal_ECM | hit_module_EMT_ECM_invasion | hit_module_metabolism_lipid_glycolytic | hit_module_DNA_repair_damage_response | hit_module_angiogenesis_hypoxia | hit_module_antigen_presentation_immune | hit_hallmark_MAPK_EGFR_axis | evidence_score | combined_score | evidence_tier | translational_priority_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Mean | -2.685546859119723 | -20.38680996904215 | 1.8910823210506067e-55 | welch_ttest | 3.574145586785647e-53 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.9847242263908929 | 0.9180535513868848 | tier_A_cuantitativo_fuerte | 0.6879910213243546 |
| Median | -2.1726396281523104 | -18.305103742637815 | 3.773875692718863e-54 | welch_ttest | 5.9023527582814135e-52 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.9887085137085138 | 0.8418310084976752 | tier_A_cuantitativo_fuerte | 0.6908882475549142 |
| StdDev | 0.1092110039519997 | 3.946070899616374 | 9.46413766994114e-05 | welch_ttest | 0.0001512370654558 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.335357543690877 | 0.4293730960397627 | tier_C_hipotesis | 0.2488151354818021 |
| A1BG | -0.8650758249918371 | -7.0373093370808215 | 9.325087489367572e-12 | welch_ttest | 2.675279239704947e-11 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.7091229757896425 | 0.7302068302068302 | tier_A_cuantitativo_fuerte | 0.5068686868686869 |
| A1CF | 0.5595804388109276 | 4.330104535999345 | 1.90203757670288e-05 | welch_ttest | 3.2053521658729695e-05 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.4673400673400673 | 0.5544332210998878 | tier_C_hipotesis | 0.3409042809042809 |
| A2M | -1.1620525014399423 | -9.218894890256031 | 3.282464381400479e-18 | welch_ttest | 1.676718292120785e-17 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.8099005932339265 | 0.4990700657367324 | tier_A_cuantitativo_fuerte | 0.5666602533269199 |
| AAAS | 1.26633305403571 | 10.168841021744711 | 1.1278088191473348e-21 | welch_ttest | 7.51511068912599e-21 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.7687349687349687 | 0.5005611672278338 | tier_A_cuantitativo_fuerte | 0.5284688151354817 |
| AACS | 0.7441865940547165 | 6.066351257584723 | 3.540316877590428e-09 | welch_ttest | 8.208533964881598e-09 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.5076879910213242 | 0.4293730960397627 | tier_B_cuantitativo_moderado | 0.3465640532307198 |
| AAGAB | 1.6029800731958983 | 11.126773061240607 | 3.798134765463995e-25 | welch_ttest | 3.403587145430882e-24 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.8340388007054673 | 0.7729357062690396 | tier_A_cuantitativo_fuerte | 0.5752413019079686 |
| AAK1 | -0.3776032229259676 | -2.6962472345819366 | 0.0073183106172676 | welch_ttest | 0.0097053589878585 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.2959595959595959 | 0.4293730960397627 | tier_D_exploratorio | 0.2134776334776335 |
| AAMDC | -1.1882764997492985 | -7.804741104855018 | 7.248711393758534e-14 | welch_ttest | 2.497801821153147e-13 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.7309042809042808 | 0.4293730960397627 | tier_A_cuantitativo_fuerte | 0.5007118807118808 |
| AAMP | 1.3393668830747605 | 10.321775072424243 | 3.565043151670334e-22 | welch_ttest | 2.4542134809015314e-21 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | False | False | False | False | False | 0.7208193041526374 | 0.4293730960397627 | tier_A_cuantitativo_fuerte | 0.4805323071989739 |

## Concordancia RNA-seq externa (TCGA u otro)

_No cargada en esta corrida._ Para activarla, fija `external_evidence.tcga_rnaseq_log2fc_tsv` (formato: columnas `gene_symbol`, `log2fc_rnaseq` o alias). Mientras tanto, **no se publican extractos**: cualquier `tcga_rnaseq_concordance.tsv` previo se borra al inicio de la etapa para no arrastrar plantillas.


## Supervivencia externa (Cox por gen)

_No cargada en esta corrida._ Fija `external_evidence.survival_gene_cox_tsv` con columnas `gene_symbol`, `cox_hr`, `cox_p` (export GEPIA / cBioPortal). Sin TSV real, no se renderiza tabla.


---

**Cómo preparar las tablas de referencia:** coloca TSV bajo `metadata/tcga_auxiliary/` (o ruta absoluta) y fija `external_evidence.tcga_rnaseq_log2fc_tsv` / `survival_gene_cox_tsv`. Los archivos `*_template.tsv` en el repo son **plantillas de columnas**; sustituye valores por datos curados de tu fuente (GDC, Xena, GEPIA, etc.).
