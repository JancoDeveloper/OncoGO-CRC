# Resultados base (generado automáticamente)

> **Workspace de datos:** `data/`, `results/`, `reports/`, `metadata/` y `docs/web/` apuntan a `/Volumes/JANCODATA/OncoGO-CRC` (el clon conserva `config/`, `src/` y el resto de `docs/`). Trazabilidad de corrida: `results/run_manifest.json`.

Resumen **cuantitativo** del último `oncogo-crc run`: proteínas por `padj`, **ORA** (lista binaria), **GSEA pre-ranked** (ranking continuo, mejor tolerancia cuando casi todo el proteoma cambia), firma **CPTAC ∩ PRIDE** y extractos ML/red/traslación. La narrativa interpretativa va en [**Síntesis para tesis**](thesis_synthesis.md).

**Umbrales en esta corrida:** FDR diferencial **0.05**; columna `significant` incluye |log2fc| ≥ **0** si `differential.min_abs_log2fc` > 0; ORA puede exigir además |log2fc| ≥ **0** (`go.ora_min_abs_log2fc`); GSEA usa métrica(s) de ranking **t_stat, signed_neglog10_p** y FDR término **0.25** (exploratorio). Cada métrica adicional genera `*_gsea_go_<métrica>.tsv`; el primero se refleja también en `*_gsea_go.tsv`. **Universo ORA/GSEA:** genes mapeados desde la matriz normalizada ∩ `metadata/gene_ontology/go_gene_long.tsv`.

## Tumor vs tejido adyacente normal (CPTAC)

`cptac_tumor_vs_normal`

**Proteínas significativas** (FDR **0.05**; ver `significant` y filtros en config): **5031** / **6239**

### Top filas diferenciales (por `padj`)

| protein_id | log2fc | t_stat | p_value | padj | significant |
| --- | --- | --- | --- | --- | --- |
| EPB41L3 | -1.8613918231373312 | -26.87052799125075 | 4.134102127300767e-88 | 2.579266317222949e-84 | True |
| ITGA3 | -2.106697767337489 | -22.840903044193876 | 5.5520568690343056e-73 | 1.7319641402952518e-69 | True |
| ARRB1 | -2.1594262825166552 | -22.803097384555667 | 1.4731679090589213e-71 | 3.0636981948728703e-68 | True |
| MPI | -2.1075478363097195 | -22.042600932060516 | 4.914633074090869e-70 | 7.665598937313233e-67 | True |
| NUP205 | -2.2509836007801036 | -22.64001468906791 | 1.6724669221618915e-69 | 2.0869042254736085e-66 | True |
| PRKAR1A | -2.331044288938557 | -21.763635457226513 | 4.124900715419702e-69 | 4.2892092605839215e-66 | True |
| GAP43 | -2.0003465337746693 | -21.53271083260022 | 3.5078433001955365e-66 | 3.1264906214171355e-63 | True |
| LPAR1 | -1.572830114329287 | -21.935679636398643 | 3.0256995377216557e-65 | 2.359667426980676e-62 | True |
| AGL | -2.176319620522187 | -20.94307650390909 | 9.476029669195104e-65 | 6.56899434512314e-62 | True |
| GNA11 | -2.157946086762687 | -21.019892326210424 | 1.11465927825504e-63 | 6.954359237033194e-61 | True |
| SCGN | -1.7710195175534358 | -20.443348563741 | 1.9417951645989544e-63 | 1.0095716693277398e-60 | True |
| MATN2 | -1.8022054832646208 | -21.021769572361464 | 1.8510390959488733e-63 | 1.0095716693277398e-60 | True |
| MECP2 | -2.0480745659741344 | -20.597173480859183 | 2.796481667577185e-62 | 1.3420960864626197e-59 | True |
| PLP1 | -1.8239672819187784 | -20.400207066351797 | 1.3144849954960609e-61 | 5.857908490642804e-59 | True |
| NIT1 | -1.9052015155076207 | -19.894332248958698 | 3.267261129816389e-61 | 1.3589628125949634e-58 | True |

### GO — ORA offline (FDR término **0.05**; **0** términos pasan)

_Sin filas._

_Reducción semántica:_ `results/go/cptac_tumor_vs_normal_go_reduced.tsv`


_Sin GSEA (etapa `gsea`, requiere `gseapy` en `pip install -e ".[go-remote]"`)._

## Plasma CRC vs control sano

`plasma_crc_vs_healthy`

**Proteínas significativas** (FDR **0.05**; ver `significant` y filtros en config): **120** / **518**

### Top filas diferenciales (por `padj`)

| protein_id | log2fc | t_stat | p_value | padj | significant |
| --- | --- | --- | --- | --- | --- |
| P18428 | 1.2543458337810949 | 8.539477097584285 | 9.18521065610093e-10 | 4.757939119860282e-07 | True |
| P06396 | -0.7530256002908366 | -7.122878739713295 | 3.889622040419932e-08 | 1.0074121084687624e-05 | True |
| P0DJI8 | 3.2565411464710934 | 6.74675151507289 | 6.610912092487636e-08 | 1.1414841546361982e-05 | True |
| P02750 | 1.23288985345199 | 6.6088412889283505 | 1.5672394302055288e-07 | 2.02957506211616e-05 | True |
| P20305-2 | -0.6368224726005849 | -6.073943392676239 | 4.806443565570877e-07 | 4.979475533931429e-05 | True |
| P01011 | 0.7518343786754791 | 6.150226909393234 | 6.856064657322412e-07 | 5.326421877279825e-05 | True |
| P02748 | 0.7440439040564399 | 6.22887306775164 | 8.249711901898953e-07 | 5.326421877279825e-05 | True |
| Q5VU43;Q5VU43-10;Q5VU43-11;Q5VU43-12;Q5VU43-3;Q5VU43-4;Q5VU43-6;Q5VU43-7;Q5VU43-8;Q5VU43-9;Q8K389;Q96SN8;Q96SN8-2;Q96SN8-3;Q96SN8-4;Q9JLH5 | 1.2276585393464394 | 5.849152207479296 | 9.254400945080776e-07 | 5.326421877279825e-05 | True |
| P06727 | -1.1790429765372783 | -6.1557288915122275 | 8.726264140518917e-07 | 5.326421877279825e-05 | True |
| P13647 | 1.3028418464955251 | 5.561569807531528 | 2.287048847888162e-06 | 0.0001184691303206 | True |
| P41981 | 1.3803026589305532 | 5.37381187863695 | 5.423148218600252e-06 | 0.0002553809797486 | True |
| P51884 | -0.589602593541219 | -5.244095157341267 | 6.580828324707079e-06 | 0.0002721870423071 | True |
| P24266 | 1.208052469847999 | 5.238647737641073 | 6.830948938209503e-06 | 0.0002721870423071 | True |
| P49913;Q1KLX1;Q1KLY3 | 1.049102858501854 | 5.172800296324417 | 8.48636329584187e-06 | 0.0003139954419461 | True |
| P35527 | 1.3125976157681445 | 5.067647000302324 | 1.1038077055333566e-05 | 0.0003573577446664 | True |

### GO — ORA offline (FDR término **0.05**; **0** términos pasan)

_Sin filas._

_Reducción semántica:_ `results/go/plasma_crc_vs_healthy_go_reduced.tsv`


_Sin GSEA (etapa `gsea`, requiere `gseapy` en `pip install -e ".[go-remote]"`)._

## PRIDE validación externa: tumor vs normal

`pride_validation_tumor_vs_normal`

**Proteínas significativas** (FDR **0.05**; ver `significant` y filtros en config): **1933** / **6824**

### Top filas diferenciales (por `padj`)

| protein_id | log2fc | t_stat | p_value | padj | significant |
| --- | --- | --- | --- | --- | --- |
| A8K7I4 | -9.374870970463215 | -16.31773684002654 | 1.8785002581162063e-16 | 1.2818885761384993e-12 | True |
| Q59EG0 | -4.531526049606704 | -11.851229147921847 | 8.54047892886888e-13 | 2.1713404767548854e-09 | True |
| K7ES52 | -4.527587596686441 | -11.856034128871062 | 9.545752389016202e-13 | 2.1713404767548854e-09 | True |
| Q15124 | -6.774102607883439 | -11.369835746882387 | 2.128720662793474e-12 | 3.631597450725667e-09 | True |
| B3KWQ6 | -5.574200220491463 | -12.624336913299048 | 9.839478277197546e-12 | 1.3428919952719211e-08 | True |
| H0UI04 | 3.6390853286570537 | 11.07889391437977 | 1.6455099001006225e-11 | 1.4036199447858313e-08 | True |
| A5A3E0 | -2.466397941820901 | -11.820128512020831 | 1.6445865904938168e-11 | 1.4036199447858313e-08 | True |
| P06312 | -4.508443041018882 | -11.048505737574969 | 1.3256135374311349e-11 | 1.4036199447858313e-08 | True |
| Q9NXS2 | -2.1451568214075607 | -11.178540446101218 | 2.3374577965581e-11 | 1.6205030139561094e-08 | True |
| M0R0N8 | -5.003719774394643 | -10.462054277165638 | 2.3747113334644043e-11 | 1.6205030139561094e-08 | True |
| Q9HC39 | -5.0087109323726935 | -10.427387406562444 | 2.6413252599510654e-11 | 1.6385821430823702e-08 | True |
| Q96FN9 | -5.017043738990921 | -10.368134825839386 | 3.171761279189124e-11 | 1.8036749140988817e-08 | True |
| V9HWH9 | 1.7368660206291224 | 10.330818451343132 | 3.812507881219514e-11 | 1.8913900504736072e-08 | True |
| P23141-3 | -4.439025302017971 | -10.563649005238153 | 3.880343010936475e-11 | 1.8913900504736072e-08 | True |
| P17947 | -4.923157865375089 | -10.190938948996145 | 4.6454573422454136e-11 | 2.1133733935655136e-08 | True |

### GO — ORA offline (FDR término **0.05**; **0** términos pasan)

_Sin filas._

_Reducción semántica:_ `results/go/pride_validation_tumor_vs_normal_go_reduced.tsv`


_Sin GSEA (etapa `gsea`, requiere `gseapy` en `pip install -e ".[go-remote]"`)._

## Comparaciones opcionales ejecutadas

### Estadio avanzado vs temprano (si stage está disponible)

`cptac_stage_advanced_vs_early`

**Proteínas significativas** (FDR **0.05**; ver `significant` y filtros en config): **22** / **6239**

### Top filas diferenciales (por `padj`)

| protein_id | log2fc | t_stat | p_value | padj | significant |
| --- | --- | --- | --- | --- | --- |
| MED17 | 1.3778376359094904 | 5.5199751797897 | 1.070113026003464e-06 | 0.0066764351692356 | True |
| HEATR5A | 1.2649719823516756 | 5.160129192617317 | 3.770980657329293e-06 | 0.0117635741605387 | True |
| SLC30A9 | 1.370474717768844 | 5.087781085907412 | 6.983736715994562e-06 | 0.01452384445703 | True |
| FRZB | 1.4318963848296484 | 4.92505575432773 | 1.3427119807896968e-05 | 0.0167543600962938 | True |
| GATB | 1.140655503387193 | 4.845247203296126 | 1.1637485750676762e-05 | 0.0167543600962938 | True |
| PBK | 1.0617117605739714 | 4.708032411542605 | 1.794318440031229e-05 | 0.0186579212455914 | True |
| YTHDC2 | 1.634305495262996 | 4.650732653578827 | 3.815796465080896e-05 | 0.0340096487794853 | True |
| SYAP1 | 1.4327494341179905 | 4.3976363074958975 | 5.321208914623405e-05 | 0.0410848818543092 | True |
| ELMO1 | 1.0695274352990272 | 4.3608893701204625 | 5.926653897880797e-05 | 0.0410848818543092 | True |
| DAB2 | 1.0666068708092058 | 4.245490645775058 | 9.355339984896598e-05 | 0.044696460222188 | True |
| TUBGCP3 | 1.0923464799983384 | 4.347877672339777 | 7.990944994655354e-05 | 0.044696460222188 | True |
| HPS6 | 0.9403025010576992 | 4.23626354195643 | 0.0001002965929012 | 0.044696460222188 | True |
| GIT2 | 1.1334764393061345 | 4.35076909318187 | 7.21714302652643e-05 | 0.044696460222188 | True |
| SPEN | 1.0520220355641805 | 4.247193684565189 | 8.902074101015138e-05 | 0.044696460222188 | True |
| PFKP | 1.1261945208931476 | 4.055869631506288 | 0.0001647153727384 | 0.0481329939358396 | True |

### GO — ORA offline (FDR término **0.05**; **0** términos pasan)

_Sin filas._

_Reducción semántica:_ `results/go/cptac_stage_advanced_vs_early_go_reduced.tsv`


_Sin GSEA (etapa `gsea`, requiere `gseapy` en `pip install -e ".[go-remote]"`)._

## Firma núcleo CPTAC ∩ PRIDE (reproducibilidad multi-cohorte)

Cruce **por UniProt canónico** (primera accesión del grupo MaxQuant; isoformas `-n` colapsadas) usando el mapa maestro `metadata/protein_id_maps/protein_mapping_master.tsv` para traducir símbolos CPTAC → UniProt. Sin mapa amplio, `n_cptac_with_join_uni` será bajo y la firma núcleo se subestima. Métricas `n_intersection_raw_protein_id` vs `n_intersection_uniprot_bridge` diagnostican el problema histórico «0 intersección» por IDs incompatibles.

| metric | value | scope |
| --- | --- | --- |
| n_proteins_cptac | 6239.0 | strict |
| n_proteins_pride | 6824.0 | strict |
| n_intersection_raw_protein_id | 0.0 | strict |
| n_cptac_with_join_uni | 6018.0 | strict |
| n_pride_with_join_uni | 6690.0 | strict |
| n_intersection_uniprot_bridge | 1967.0 | strict |
| n_both_significant_padj_lt_0.05 | 430.0 | strict |
| n_cptac_sig_only | 1177.0 | strict |
| n_pride_sig_only | 91.0 | strict |
| n_pairs_log2fc_cptac_pride | 1967.0 | strict |
| log2fc_pearson_r_cptac_pride | 0.0991049168092914 | strict |
| log2fc_spearman_r_cptac_pride | 0.1180704007724229 | strict |
| log2fc_same_sign_fraction | 0.5592272496187087 | strict |

_Detalle ancho:_ `results/reproducibility/cptac_pride_differential_overlap.tsv` (columna `_join_uni`, etapa `cohort_overlap`).

## ML — priorización (cohorte primaria)

| protein_id | combined_score | rf_importance | l1_abs_coef | elasticnet_abs_coef |
| --- | --- | --- | --- | --- |
| NIT1 | 0.9926991505048886 | 0.0047874629055723 | 1.3009708322666909 | 0.0826624431488116 |
| HIST2H2AA3 | 0.9911043436448148 | 0.005382617279725 | 0.4571256090467452 | 0.0606264972509511 |
| ABI3BP | 0.9893091841641288 | 0.0044391674269996 | 0.0789241088372907 | 0.0565640340737168 |
| LYVE1 | 0.9880990543356308 | 0.0037838862364111 | 2.018609090952526 | 0.0673369539492552 |
| PRKACB | 0.9847651867286424 | 0.0034525730589938 | 0.1169090082364016 | 0.0458256111422921 |
| AHCYL2 | 0.9731046642090078 | 0.0008683151950578 | 0.942912486568058 | 0.0330467252504904 |
| SOD2 | 0.9711412085270076 | 0.0006417547581914 | 1.5783517909403293 | 0.1199307366944412 |
| ILKAP | 0.9674467062029172 | 0.0006463186099697 | 1.435203508411173 | 0.0698393838492519 |
| MRPS12 | 0.9661724635358232 | 0.0006158969269426 | 4.458248135441045 | 0.0691023055896343 |
| DHRS9 | 0.9652989261099536 | 0.0005121799784656 | 0.2923918094612863 | 0.0774118811156348 |
| MTIF3 | 0.9637602179836512 | 0.0006025021639246 | 1.3986057130160858 | 0.0687322259613703 |
| GSTT1 | 0.9618288187209488 | 0.0004653536387714 | 0.0996878442997384 | 0.0865620405524969 |
| OAS2 | 0.943716941817599 | 0.0001401887559813 | 1.0378663064138944 | 0.0727193463607293 |
| MSI2 | 0.9266468985414328 | 0.0001278226651943 | 5.361237041934597 | 0.1088600663290099 |
| AGL | 0.9237337714377304 | 0.0184107470729565 | 0.0 | 0.0871633217481486 |

### Métricas test (hold-out)

| model | accuracy_test | roc_auc_test | roc_auc_mean | roc_auc_std | n_folds | n_outer_folds |
| --- | --- | --- | --- | --- | --- | --- |
| sparse_linear_sgd | 0.9489795918367347 | 0.9804166666666666 |  |  |  |  |
| elastic_net_sgd | 1.0 | 1.0 |  |  |  |  |
| random_forest_cv |  |  | 1.0 | 0.0 | 5.0 |  |
| random_forest | 1.0 | 1.0 |  |  |  |  |
| random_forest_nested |  |  | 1.0 | 0.0 |  | 5.0 |

### Nested cross-validation (sesgo reducido frente a un solo hold-out)

Pliegues **externos** estiman AUC ROC; en cada entrenamiento externo, **GridSearchCV** con pliegues **internos** elige hiperparámetros. Con `ml.nested_cv.only_rf: true` solo aparece **random_forest_nested** (más rápido); con `only_rf: false` se incluyen también L1 y Elastic Net. Los pliegues se **recortan** automáticamente si la clase minoritaria es pequeña (`split_adjustment` en el JSON de resumen). Ver `ml.nested_cv` en `config.yaml`.

| outer_fold | model | roc_auc | best_params_json |
| --- | --- | --- | --- |
| 1 | random_forest_nested | 1.0 | {"max_depth": null, "min_samples_leaf": 1} |
| 2 | random_forest_nested | 1.0 | {"max_depth": null, "min_samples_leaf": 1} |
| 3 | random_forest_nested | 1.0 | {"max_depth": null, "min_samples_leaf": 1} |
| 4 | random_forest_nested | 1.0 | {"max_depth": null, "min_samples_leaf": 1} |
| 5 | random_forest_nested | 1.0 | {"max_depth": null, "min_samples_leaf": 1} |

_Resumen numérico:_ `results/ml/cptac_colorectal_nested_cv_summary.tsv`.


### Estabilidad bootstrap de features (L1)

| feature_id | bootstrap_topk_frequency |
| --- | --- |
| NUP93 | 0.34 |
| GNG2 | 0.32 |
| NOC2L | 0.32 |
| RRP9 | 0.3 |
| PRPF19 | 0.3 |
| EDF1 | 0.28 |
| SCG2 | 0.28 |
| CROT | 0.28 |
| UTP14A | 0.28 |
| NCBP2 | 0.28 |
| ABCA8 | 0.28 |
| RANGAP1 | 0.28 |
| PLSCR1 | 0.26 |
| RPRD2 | 0.26 |
| COA7 | 0.26 |
| GPD1L | 0.26 |
| PGRMC2 | 0.26 |
| ATP2B1 | 0.26 |
| RBM28 | 0.26 |
| ATP5L | 0.26 |

_Frecuencia en la que cada gen entra en el top‑k por |coef| en réplicas bootstrap._ Config: `ml.feature_stability_bootstrap_n`.


## Traslación clínica (extracto)

| feature_id | translational_priority_score | evidence_tier | functional_module_hints | cms_pathway_hints | log2fc | padj | combined_score | evidence_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| KRAS | 0.915952877063632 | tier_A_cuantitativo_fuerte |  | CMS3_metabolic | -1.656143164061718 | 4.551912147824223e-40 | 0.8918937329700272 | 0.95232008334669 |
| AHNAK2 | 0.7666180477640647 | tier_A_cuantitativo_fuerte |  | CMS4_mesenchymal_ECM | -2.3104276904342207 | 1.4692704973789453e-52 | 0.905517711171662 | 0.995596249398942 |
| LAMC1 | 0.7608254527969227 | tier_A_cuantitativo_fuerte |  | CMS4_mesenchymal_ECM | -2.028127435625014 | 1.5074172295950068e-48 | 0.8416973874018272 | 0.987874659400545 |
| MSH2 | 0.7589485494470267 | tier_A_cuantitativo_fuerte |  | CMS1_MSI_immune | -1.3863334782655832 | 6.28270497772663e-17 | 0.4290671581984292 | 0.7902187850617086 |
| COL5A3 | 0.7575364641769515 | tier_A_cuantitativo_fuerte |  | CMS4_mesenchymal_ECM | -2.1547768038290167 | 1.0423726113414053e-57 | 0.8295800609071966 | 0.9849615322968424 |
| POLD1 | 0.7454351658919699 | tier_A_cuantitativo_fuerte |  |  | 1.2242057898103251 | 2.1157445043864924e-21 | 0.769746754287546 | 0.8486175669177752 |
| MET | 0.7104712293636801 | tier_A_cuantitativo_fuerte |  |  | 1.3734738646898148 | 2.5348877220703807e-25 | 0.7516028209648983 | 0.8126502644654592 |
| SMAD4 | 0.7092851418496554 | tier_A_cuantitativo_fuerte |  |  | -1.4392615916658769 | 1.2120586326146172e-26 | 0.4290671581984292 | 0.8257493188010899 |
| NUP205 | 0.6976358390767752 | tier_A_cuantitativo_fuerte |  |  | -2.2509836007801036 | 2.0869042254736085e-66 | 0.9229724314794036 | 0.9969786824811668 |
| PRKAR1A | 0.6973297002724796 | tier_A_cuantitativo_fuerte |  |  | -2.331044288938557 | 4.2892092605839215e-66 | 0.9130509697066838 | 0.9967743228081424 |
| PFKM | 0.6947988459689054 | tier_A_cuantitativo_fuerte |  |  | -2.218544469169721 | 5.107957705360924e-51 | 0.8934083987818561 | 0.9932841801570764 |
| AGL | 0.6935758935726879 | tier_A_cuantitativo_fuerte |  |  | -2.176319620522187 | 6.56899434512314e-62 | 0.9237337714377304 | 0.991677352139766 |
| EML1 | 0.6932168616765507 | tier_A_cuantitativo_fuerte |  |  | -2.079798421239502 | 7.161256431419385e-46 | 0.9035141849655394 | 0.9902227921141208 |
| AHNAK | 0.6929668216060265 | tier_A_cuantitativo_fuerte |  |  | -2.1628062179446865 | 5.947284298917679e-50 | 0.9148140727680718 | 0.9904552011540312 |
| GCN1 | 0.6924955922423466 | tier_A_cuantitativo_fuerte |  |  | -2.3811856981260258 | 1.357371813219395e-50 | 0.916873697707966 | 0.9901546722231126 |

## Red — hubs (extracto)

| protein_id | degree | is_hub |
| --- | --- | --- |
| EML1 | 566 | True |
| ATP2B4 | 561 | True |
| MCAM | 544 | True |
| NID1 | 543 | True |
| GNAO1 | 531 | True |
| DMD | 531 | True |
| INA | 529 | True |
| LAMC1 | 525 | True |
| AHNAK2 | 520 | True |
| PARVA | 520 | True |
| DIP2C | 519 | True |
| LAMA5 | 519 | True |
| EFTUD2 | 516 | True |
| CRYAB | 516 | True |
| RPS7 | 502 | True |

## Export contextual (ProteinsLanguage / STRING / AlphaFold)

Bundle JSON generado por la etapa `export_bundle`: `results/export/protein_context_bundle.json`. Incluye `external_links` (STRING, UniProt, AlphaFold), vista previa de *leading edges* GSEA si existen, métricas de reproducibilidad CPTAC∩PRIDE y, si activas `export_bundle.uniprot_top_n`, resúmenes UniProt cacheados.

