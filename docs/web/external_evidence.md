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
  "n_rows_evidence_matrix": 6239,
  "n_annotation_columns": 8,
  "tcga_rnaseq": {
    "status": "skipped_no_path"
  },
  "survival_cox": {
    "status": "skipped_no_path"
  }
}
```

## Extracto — matriz de evidencia

| protein_id | log2fc | t_stat | p_value | stat_method | padj | significant | comparison_id | matrix_input | hit_tier_tier1_validated_clinical | hit_tier_tier2_emerging | hit_tier_tier3_pathway_anchor | hit_cms_CMS1_MSI_immune | hit_cms_CMS2_canonical_WNT | hit_cms_CMS3_metabolic | hit_cms_CMS4_mesenchymal_ECM | hit_hallmark_MAPK_EGFR_axis | evidence_score | combined_score | evidence_tier | translational_priority_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Mean | -2.685546859119723 | -20.38680996904215 | 1.8910823210506067e-55 | welch_ttest | 3.575291697283253e-53 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.984596890527328 | 0.9184444622535662 | tier_A_cuantitativo_fuerte | 0.6878762622215099 |
| Median | -2.1726396281523104 | -18.305103742637815 | 3.773875692718797e-54 | welch_ttest | 5.9042454479586336e-52 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.9817278410001602 | 0.7933883635197948 | tier_A_cuantitativo_fuerte | 0.6843676871293476 |
| StdDev | 0.1092110039519997 | 3.946070899616374 | 9.46413766994114e-05 | welch_ttest | 0.000151324333477 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.3356026606828017 | 0.4290671581984292 | tier_C_hipotesis | 0.2490270876743067 |
| A1BG | -0.8650758249918371 | -7.0373093370808215 | 9.325087489367572e-12 | welch_ttest | 2.6724492809446155e-11 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.6396738259336432 | 0.4290671581984292 | tier_B_cuantitativo_moderado | 0.4420740503285782 |
| A1CF | 0.5595804388109276 | 4.330104535999345 | 1.90203757670288e-05 | welch_ttest | 3.207246605688991e-05 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.4685246033018111 | 0.5574050328578297 | tier_C_hipotesis | 0.3419362077255971 |
| A2M | -1.1620525014399423 | -9.218894890256031 | 3.282464381400479e-18 | welch_ttest | 1.6786307602916054e-17 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.8107268793075812 | 0.4997876262221509 | tier_A_cuantitativo_fuerte | 0.5673986215739703 |
| AAAS | 1.26633305403571 | 10.168841021744711 | 1.1278088191473209e-21 | welch_ttest | 7.525560665946671e-21 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.7711051450552973 | 0.5037305657958006 | tier_A_cuantitativo_fuerte | 0.5306202917134157 |
| AACS | 0.7441865940547165 | 6.066351257584723 | 3.540316877590428e-09 | welch_ttest | 8.21421978404116e-09 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.5082024362878665 | 0.4290671581984292 | tier_B_cuantitativo_moderado | 0.3469466260618689 |
| AAGAB | 1.6029800731958983 | 11.126773061240607 | 3.798134765463995e-25 | welch_ttest | 3.3997938022567966e-24 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.7569442218304214 | 0.4290671581984292 | tier_A_cuantitativo_fuerte | 0.5032986055457606 |
| AAK1 | -0.3776032229259676 | -2.6962472345819366 | 0.0073183106172676 | welch_ttest | 0.0097022821804362 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.2962053213656034 | 0.4290671581984292 | tier_D_exploratorio | 0.2136912966821606 |
| AAMDC | -1.1882764997492985 | -7.804741104855018 | 7.248711393758534e-14 | welch_ttest | 2.497223102466012e-13 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.7308623176791152 | 0.4290671581984292 | tier_A_cuantitativo_fuerte | 0.5007340920019233 |
| AAMP | 1.3428178843383811 | 10.354361264707157 | 2.7303828384900956e-22 | welch_ttest | 1.899092366704538e-21 | True | cptac_tumor_vs_normal | data/processed/cptac_colorectal_log2norm.tsv | False | False | False | False | False | False | False | False | 0.7219386119570443 | 0.4290671581984292 | tier_A_cuantitativo_fuerte | 0.4812950793396378 |

## Extracto — concordancia proteómica vs RNA-seq de referencia

| protein_id | log2fc | t_stat | p_value | padj | significant | log2fc_rnaseq_ref | direction_concordant |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GENE001 | 1.3617148998166964 | 6.30488707628926 | 1.1258207826320949e-07 | 4.5032831305283785e-06 | True | 0.18 | True |
| GENE002 | 1.3709664342585903 | 6.3854515433405785 | 9.086670688840952e-08 | 4.5032831305283785e-06 | True | -0.12 | False |
| GENE003 | 0.5618892953656645 | 1.9854243100898143 | 0.0537739910531516 | 0.258115157055128 | False | 0.09 | True |
| GENE004 | 1.1809100041797955 | 4.961970932467973 | 9.987663322871891e-06 | 0.0001089563271586 | True | -0.21 | False |
| GENE005 | 0.695666010625086 | 2.5162411097404855 | 0.0155239079268501 | 0.0887080452962867 | False | 0.14 | True |
| GENE006 | 0.6145356930178294 | 2.1899342190268447 | 0.0341724547616757 | 0.1782910683217864 | False | -0.08 | False |
| GENE007 | 1.2004992385624895 | 5.090055289959083 | 6.570767546214481e-06 | 7.884921055457378e-05 | True | 0.11 | True |
| GENE008 | 0.9754793459720978 | 3.789292804952487 | 0.0004393479338694 | 0.0029289862257961 | True | -0.15 | False |
| GENE009 | 1.3916556610766888 | 6.570979222860722 | 5.9346042633666355e-08 | 4.5032831305283785e-06 | True | 0.07 | True |
| GENE010 | 1.246882863365986 | 5.408041263984643 | 2.288019822140102e-06 | 4.576039644280204e-05 | True | -0.19 | False |
| GENE011 | 1.011363360108103 | 3.975449123064156 | 0.000253287312066 | 0.0018996548404956 | True | 0.22 | True |
| GENE012 | 1.068348980284073 | 4.285613954067783 | 9.416388455484864e-05 | 0.0007533110764387 | True | -0.06 | False |

## Extracto — Cox externo fusionado

| protein_id | log2fc | t_stat | p_value | padj | significant | cox_hr_external | cox_p_external |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GENE001 | 1.3617148998166964 | 6.30488707628926 | 1.1258207826320949e-07 | 4.5032831305283785e-06 | True | 1.12 | 0.42 |
| GENE002 | 1.3709664342585903 | 6.3854515433405785 | 9.086670688840952e-08 | 4.5032831305283785e-06 | True | 0.88 | 0.51 |
| GENE003 | 0.5618892953656645 | 1.9854243100898143 | 0.0537739910531516 | 0.258115157055128 | False | 1.05 | 0.63 |
| GENE004 | 1.1809100041797955 | 4.961970932467973 | 9.987663322871891e-06 | 0.0001089563271586 | True | 0.92 | 0.38 |
| GENE005 | 0.695666010625086 | 2.5162411097404855 | 0.0155239079268501 | 0.0887080452962867 | False | 1.18 | 0.29 |

---

**Cómo preparar las tablas de referencia:** coloca TSV bajo `metadata/tcga_auxiliary/` (o ruta absoluta) y fija `external_evidence.tcga_rnaseq_log2fc_tsv` / `survival_gene_cox_tsv`. Los archivos `*_template.tsv` en el repo son **plantillas de columnas**; sustituye valores por datos curados de tu fuente (GDC, Xena, GEPIA, etc.).
