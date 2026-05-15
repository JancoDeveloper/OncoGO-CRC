# Limitaciones (auto)

- **Modelo estadístico:** t-test por proteína + FDR Benjamini–Hochberg; asume razonablemente muestreo independiente entre grupos; **no** modela emparejamiento tumor-normal por sujeto (mejora futura: t pareado o modelos mixtos si hay `participant_id` confiable).
- **GO ORA:** sensible a listas enormes de significativos; usa `go.ora_min_abs_log2fc` y `ora_foreground_max_genes`; universo = proteínas detectadas; FDR término a 0.05.
- **GSEA:** depende de la calidad del ranking y del archivo `go_gene_long.tsv`; FDR término habitualmente más laxo (0.25 en config); permutaciones finitas (200).
- **ML:** hold-out simple (`test_size=0.25`) **además** de `results/ml/*_nested_cv_*.tsv` si `ml.nested_cv.enabled=true` (estimación menos optimista que un único split). **Transferencia CPTAC→PRIDE** (`ml.external_transfer.enabled`) reporta AUC en cohorte externa con genes alineados; opcional **`calibration_method`** (`sigmoid`/`isotonic`) recalibra probabilidades con CV **solo en CPTAC** antes de PRIDE. AUC modesta es interpretable (cambio de dominio). Bootstrap de features: `ml.feature_stability_bootstrap_n`.
- **Datos públicos:** versiones de matrices y metadata dependen de PDC/PRIDE en la fecha de descarga; documenta `data/versions/` y manifiestos.
- **α diferencial:** 0.05; efecto mínimo opcional `min_abs_log2fc` en la escala **log2** de la matriz **log2norm**.
- **Tejido vs plasma:** el contraste tumor/normal en tejido suele separar mejor que plasma circulante; AUC más altos en CPTAC frente a capas líquidas débiles es **esperable** y defendible frente a revisiones que citan AUC modestos en cohortes sanguíneas.
- **CMS / KEGG:** CMS suele provenir de transcriptómica; aquí solo se usa si aportas suplemento explícito. KEGG (`kegg_enrichment`, gseapy) es opcional y depende de conectividad Enrichr.
- **Druggability / reposicionamiento:** solo enlaces in silico en el bundle; no equivalen a evidencia de respuesta clínica.
- **UniProt REST:** el enriquecimiento en el bundle requiere red y respeta límites de frecuencia; resultados cacheados en `data/interim/`.
