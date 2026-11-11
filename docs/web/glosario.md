> Documento de la v1. No forma parte de los resultados certificados de la tesis; ver results/v3_7.

# Glosario — OncoGO-CRC

Glosario de los términos estadísticos, bioinformáticos y traslacionales más usados en este sitio. Está pensado para acompañar la lectura de las páginas **Resultados base**, **Síntesis para tesis** e **Interpretación biológica** sin necesidad de saltar a textbooks externos.

_Las definiciones priorizan precisión operativa sobre exhaustividad; las referencias completas están en `docs/RESULTADOS_Y_LIMITACIONES_TESIS.md`._

## log2FC (log<sub>2</sub> fold change)

Logaritmo en base 2 del cociente de medias entre grupos. log2FC = +1 → grupo A duplica al B; log2FC = −1 → grupo A es la mitad. En matrices `log2norm`, las abundancias ya están en escala log2, así que el log2FC es directamente la diferencia de medias log2.

## p-value vs. padj (FDR Benjamini–Hochberg)

El **p-value** mide la probabilidad de observar la diferencia bajo la hipótesis nula. El **padj** (q-value) corrige por las miles de pruebas simultáneas (una por proteína) usando la tasa de descubrimiento falso (FDR-BH). En proteómica de alta dimensión el padj es el umbral correcto, no el p crudo. Convención habitual: padj < 0.05 = 'significativo'.

## ORA (Over-Representation Analysis, Fisher)

Pregunta: ¿hay **más genes de la lista significativa de los esperados por azar** en cada término GO? Usa lista binaria (sig vs no-sig). Sensible al tamaño de la lista: si casi todo el proteoma es significativo (p. ej. tumor vs normal en CPTAC), la ORA pierde poder porque el 'foreground' se parece al universo. Por eso el pipeline siempre reporta también GSEA.

## GSEA (Gene Set Enrichment Analysis, pre-ranked)

Pregunta: ¿los genes de un término GO se **concentran en los extremos** del ranking (más arriba/abajo que por azar)? Usa el ranking continuo (t-stat o `signed_neglog10_p`), no listas binarias. Mucho más robusto que ORA cuando la señal es distribuida (caso típico de tumor vs normal proteómico). Subramanian et al., PNAS 2005.

## NES (Normalized Enrichment Score)

Versión normalizada del enrichment score de GSEA. **NES > 0** → término enriquecido en la cabeza del ranking (genes UP en el primer grupo). **NES < 0** → enriquecido en la cola (DOWN). |NES| > 1.5 con FDR q-val < 0.25 es la convención exploratoria; se reportan valores más permisivos que el padj diferencial porque el universo es de términos GO, no de genes.

## Leading edge

Subconjunto de genes de un término GO que aporta el grueso del enrichment score (genes más cerca del extremo del ranking dentro de ese término). Es el conjunto **más interpretable**: si ves MRPL/MRPS en el leading edge, ese término refleja realmente traducción mitocondrial.

## Reducción semántica GO (`*_go_reduced.tsv`)

Muchos términos GO se solapan (p. ej. 'mitochondrial translation' incluye 'mitoribosome'). El pipeline agrupa términos con Jaccard ≥ 0.45 sobre sus genes y conserva el representante con menor p — evita inflar listas con redundancia.

## ROC AUC

Área bajo la curva ROC. 0.5 = clasificador al azar; 1.0 = perfecto. En proteómica con N moderado y muchas features, AUC > 0.9 en hold-out puede ser optimista por azar — por eso el pipeline reporta nested CV y transferencia externa.

## Nested cross-validation

Doble bucle de CV: el **externo** estima el AUC de generalización; el **interno** elige hiperparámetros con GridSearchCV. Evita la fuga de datos típica de elegir hiperparámetros y evaluar en el mismo split (que infla el AUC). Es la métrica honesta para el capítulo de resultados ML.

## Transferencia externa (CPTAC → PRIDE)

Entrena en una cohorte (CPTAC), evalúa en otra independiente (PRIDE) **sin reajustar nada**. Es el experimento que pone a prueba la **transportabilidad** del modelo a dominio distinto (tecnología, lab, batch). Un AUC modesto aquí (p. ej. 0.65-0.80 frente a 0.95 dentro de CPTAC) es esperable y honesto: refleja el batch/tech shift, no fracaso biológico.

## Calibración (sigmoid / isotonic)

Las probabilidades de scikit-learn no siempre son fieles: un RF puede dar 0.9 cuando la frecuencia real es 0.7. La calibración (Platt sigmoid o isotonic) ajusta esa relación usando CV en train. En transfer, calibrar a veces baja el AUC pero mejora la utilidad clínica de umbrales — útil cuando se quiere decir 'probabilidad estimada de ser tumor', no solo ranking.

## Bootstrap estabilidad de features (L1)

Re-muestrea N veces el train, entrena un modelo L1 (sparse linear) y cuenta cuántas veces cada gen entra en el top-k por |coef|. Una **frecuencia alta** (p. ej. > 0.7) indica un gen robusto al re-muestreo; baja indica que su selección depende del split concreto.

## TPS (Translational Priority Score)

Score compuesto que combina: magnitud diferencial, FDR, centralidad en red, importancia ML (RF, L1, EN, SHAP), pertenencia a capas YAML (CMS, MAPK, módulos funcionales, tiers de evidencia literaria). **No es una decisión clínica**: es un **ranking explícito** de hipótesis investigativas, defendible ante un comité de tesis porque cada componente es trazable.

## Firma núcleo CPTAC ∩ PRIDE

Genes/proteínas con `padj < 0.05` y mismo signo de log2FC en **ambas cohortes** tras alinear por UniProt canónico. La fracción con mismo signo (`same_sign_fraction`) cuantifica reproducibilidad direccional; correlación Pearson/Spearman de log2FC cuantifica concordancia de magnitud (típicamente modesta entre TMT y LFQ, ~0.1–0.3).

## TMT vs LFQ (¿por qué bajan las correlaciones entre cohortes?)

**TMT** (Tandem Mass Tag, CPTAC) cuantifica multiplexado en un mismo MS; tiene buena precisión relativa pero compresión de ratios. **LFQ** (label-free, PRIDE/MaxQuant) cuantifica por intensidad de péptidos independientes; mayor varianza pero menos sesgo. Comparar log2FC entre ambos siempre baja la correlación absoluta — lo importante es la **concordancia de dirección** (signo) y la reproducibilidad de los **términos GO/GSEA**, no la magnitud bruta.

## padj < 0.25 en GSEA (¿no es muy laxo?)

Sí, intencionalmente. La convención original de Subramanian et al. (PNAS 2005) usa FDR q-val < 0.25 como umbral **exploratorio** porque la corrección se hace sobre miles de términos GO y GSEA es conservador por diseño (permutaciones). Para módulos individuales se discute |NES| alto + leading edge biológicamente coherente, no el q-val aislado.

## Catastrophic cancellation (warning SciPy)

Aparece cuando dos números casi iguales se restan y se pierde precisión. En el pipeline ocurre para proteínas con varianza ~0 dentro de un grupo (todas las muestras con el mismo valor tras imputación). No afecta al contraste global; sí marca esas filas como poco informativas — el ranking las pone al final por padj=NaN o p~1.

## Modo demo vs corrida real

`--demo` genera matrices sintéticas con genes `GENE001…` y mapa GO mínimo para validar la **lógica** del pipeline sin descargar nada. Una corrida 'real' (sin `--demo`) lee `datasets.yaml` y `download_real` baja de PDC/PRIDE. El audit JSON (`*_gsea_input_*.json`) reporta el flag `contains_demo_gene_ids` para evitar publicar resultados demo como reales.

## GOA (Gene Ontology Annotations)

Anotación oficial gen→término GO. El pipeline usa `goa_human.gaf.gz` del Gene Ontology Consortium (~14 MiB, ~38 800 genes humanos, ~17 900 términos). Regenerar con `python scripts/fetch_go_human_annotation.py`. Sin GOA real (solo demo de 10 genes), ORA y GSEA producen 0 términos: es bug de cobertura, no biológico.

## Mapping efectivo (`protein_mapping_master.tsv` + bulk UniProt)

Tabla canónica gen↔UniProt usada para colapsar isoformas y para puentear CPTAC (genes) con PRIDE (UniProt). El bulk humano reviewed (~20 000 entradas) se descarga con `scripts/fetch_uniprot_human_reviewed_gene_map.py`. Sin él, la intersección CPTAC∩PRIDE se subestima drásticamente.

## Workspace de datos (`paths.data_workspace`)

Carpeta externa (típicamente disco USB/SSD) donde residen `data/`, `results/`, `reports/`, `metadata/` y `docs/web/`. Permite mantener el clon git ligero sin perder reproducibilidad. El pipeline sincroniza metadata desde el repo al workspace en cada `run` (mtime-based). Si el workspace no es accesible, la CLI cae a rutas del repo con un aviso.

## run_manifest.json

Manifest que se reinicializa al inicio de cada `oncogo-crc run`. Contiene `run_id`, timestamp, lista de etapas ejecutadas/omitidas, hashes SHA-256 de outputs clave y razones de skip. El publicador web (`publish_results`) **lee este manifest** y se niega a publicar resultados de etapas no ejecutadas en la corrida actual — evita publicar artefactos viejos como nuevos.

## stage_bin (temprano vs avanzado)

Agrupación binaria del estadio patológico AJCC: **early** (I–II) vs **advanced** (III–IV), derivada del texto de `stage` con tokens configurables en `clinical_filters.stage_mapping`. Alimenta la comparación `cptac_stage_advanced_vs_early`.

## stage_group (I–IV)

Grupo AJCC romano extraído automáticamente del campo `stage` (p. ej. «Stage IIIA» → `III`). Permite tabular estadios sin perder granularidad del Excel clínico PDC.

## hispanic_latino_bin

Estrato derivado de `ethnicity` según categorías NIH: `hispanic_latino`, `non_hispanic_latino` o `unknown`. **No** identifica nacionalidad colombiana ni residencia en Latinoamérica; en CPTAC refleja autodeclaración en cohortes reclutadas principalmente en EE.UU. Ver [Demografía clínica](clinical_demographics.md).

## geographic_region_bin

Región inferida desde `country` y, si falta, desde `enrollment_site`: `latin_america`, `united_states`, `other` o `unknown`. Útil para transparencia; en muchos estudios CPTAC el país no está informado y la variable queda mayoritariamente `unknown`.

