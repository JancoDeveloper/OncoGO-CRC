# Prompt único para SciSpace v2 — marco metodológico OncoGO-CRC 2.0

Búsqueda bibliográfica enfocada en los **nuevos requerimientos metodológicos** (homogeneización multi-ómica,
análisis de red, marco de priorización por evidencia convergente, nivel peptídico, interpretabilidad y
reproducibilidad). Dos versiones **autocontenidas** (español e inglés): pega **un** bloque en
[SciSpace](https://scispace.com/), sin archivos extra, para **optimizar créditos y tokens**. Prioriza literatura
**2025–2026** (2024 si es seminal/metodológica).

> **Encuadre obligatorio para la búsqueda:** OncoGO-CRC es una **metodología/herramienta de análisis
> reproducible** que **prioriza candidatos como HIPÓTESIS con evidencia convergente**. La tesis **NO** promete
> dianas terapéuticas validadas, biomarcadores clínicos, ni cura. Se busca bibliografía que **valide el diseño
> metodológico** y respalde el **marco teórico**, con prudencia traslacional.

---

## 🇪🇸 Prompt maestro (ESPAÑOL — copiar todo)

```
Actúa como asistente experto en bioinformática traslacional, proteómica computacional, biología de sistemas
e integración multi-ómica. Ayúdame a encontrar literatura reciente (prioriza 2025–2026; incluye 2024 solo si
es seminal o metodológica), métodos y bases de datos que VALIDEN Y RESPALDEN el MARCO METODOLÓGICO de mi tesis
de maestría "OncoGO-CRC 2.0". IMPORTANTE: mi tesis es una HERRAMIENTA/METODOLOGÍA de análisis reproducible que
prioriza candidatos como HIPÓTESIS con evidencia convergente; NO promete dianas terapéuticas validadas, ni
biomarcadores clínicos, ni cura. Busco respaldo metodológico riguroso y prudencia traslacional, no afirmaciones
clínicas.

== QUÉ ES ONCOGO-CRC (ya implementado) ==
Pipeline reproducible que integra proteómica pública de cáncer colorrectal (CCR): CPTAC colon (PDC000116, TMT10),
PRIDE tejido (LFQ) y plasma (DIA-NN, exploratorio); armonización de identificadores proteína→gen con UniProt;
diferencial (t-test Welch/pareado + FDR-BH); reproducibilidad multi-cohorte por puente UniProt (firma núcleo con
concordancia de signo, ~1967 genes intersección, ~430 doble-significativos, Pearson de magnitud ~0.10);
enriquecimiento GO (ORA + GSEA pre-ranked, GOA humana real); redes de coexpresión; ML interpretable
(L1/ElasticNet/RandomForest, nested CV, SHAP, bootstrap); trazabilidad completa (run_manifest con hashes SHA-256).

== EVOLUCIÓN PROPUESTA (OncoGO-CRC 2.0) — 3 CAPAS PARA LAS QUE NECESITO BIBLIOGRAFÍA ==
Capa A. Homogeneización robusta e integración multi-ómica: corrección de efecto batch a nivel proteína
  (ComBat/HarmonizR), integración proteoma+transcriptoma con MOFA+/VAE, y (avanzado) domain adaptation
  adversarial para representación invariante a plataforma; evaluación con métricas antes/después (kBET, iLISI/LISI,
  gPCA, silueta, preservación de señal biológica).
Capa B. Análisis de red (módulos, no listas planas): redes clásicas (WGCNA/detección de comunidades) y Graph
  Neural Networks (GAT para grafos pequeños, GCN/GraphSAGE para grandes, Graph Transformer) sobre PPI/vías
  (STRING/Reactome/OmniPath), tomando la firma reproducible como señal; interpretabilidad con GNNExplainer y
  pesos de atención; enriquecimiento con MSigDB (Hallmark/Reactome/KEGG).
Capa C. Marco de priorización por EVIDENCIA CONVERGENTE (metodología, NO promesa de dianas): anotar cada
  candidato con evidencia pública objetiva — esencialidad CRISPR (DepMap), tractabilidad/fármacos (Open Targets,
  DGIdb, DrugBank, ChEMBL), especificidad tumoral/seguridad (Human Protein Atlas, GTEx), contexto genómico
  (cBioPortal) — con procedencia y versión fijadas; salida = tabla de HIPÓTESIS priorizadas exploratorias.
Nivel peptídico: considerar análisis a nivel de péptido usando TABLAS PEPTÍDICAS YA PUBLICADAS (robustez del
  rollup péptido→proteína, proteoformas/isoformas, detectabilidad); el reprocesado desde archivos CRUDOS con
  deep learning (DIA-BERT, AlphaDIA) se plantea como TRABAJO FUTURO, no como parte prometida.
Transversal: interpretabilidad (SHAP, Integrated Gradients, GNNExplainer) y reproducibilidad/FAIR (versionado de
  bases de datos, hashes, depósito en Zenodo/GitHub).

== TEMAS A BUSCAR (agrupa la respuesta EXACTAMENTE por estos 12) ==
1. Integración multi-ómica con análisis de factores (MOFA/MOFA+/MEFISTO) y autoencoders/VAE multi-vista en
   proteogenómica de cáncer (proteoma+transcriptoma), 2024–2026.
2. Corrección de efecto batch a NIVEL PROTEÍNA y armonización entre plataformas TMT/LFQ/DIA (ComBat, HarmonizR,
   RUV, métodos recientes), con benchmarks.
3. MÉTRICAS de evaluación de integración/batch (kBET, iLISI/cLISI/LISI, gPCA, silueta) y de preservación de
   señal biológica; buenas prácticas para no eliminar biología real.
4. Domain adaptation / aprendizaje adversarial / contrastive learning para representaciones invariantes a
   plataforma o a lote en datos ómicos.
5. Análisis a NIVEL DE PÉPTIDO: métodos de rollup péptido→proteína, proteoformas/isoformas, detectabilidad de
   péptidos, e impacto de imputación/rollup en listas diferenciales (usando datos ya procesados).
6. Deep learning para identificación/cuantificación proteómica desde datos crudos (p. ej. DIA-BERT, AlphaDIA,
   modelos tipo transformer) — SOLO como contexto de TRABAJO FUTURO y sus requisitos de cómputo/datos.
7. Graph Neural Networks (GAT, GCN, GraphSAGE, Graph Transformer) sobre PPI/vías (STRING/Reactome/OmniPath) para
   descubrimiento de módulos/subredes y clasificación en cáncer; comparación GAT vs GCN según tamaño de grafo.
8. Interpretabilidad de GNN y de modelos ómicos (GNNExplainer, atención, Integrated Gradients, SHAP) y su
   validación.
9. Análisis de red/módulos y priorización basada en biología de sistemas (WGCNA, detección de comunidades,
   enriquecimiento con MSigDB Hallmark/Reactome/KEGG) en CCR.
10. MARCOS COMPUTACIONALES DE PRIORIZACIÓN de candidatos/dianas que REUTILIZAN bases públicas (DepMap CRISPR,
    Open Targets, DGIdb, DrugBank, ChEMBL) integrando esencialidad + tractabilidad; enfoque metodológico, no
    afirmaciones clínicas.
11. Evaluación de ESPECIFICIDAD TUMORAL y seguridad de candidatos con GTEx y Human Protein Atlas; contexto
    genómico con cBioPortal en CCR.
12. Estándares de reporte, interpretabilidad, reproducibilidad y prudencia traslacional (TRIPOD+AI, DOME, REMARK,
    principios FAIR) aplicables a pipelines de reuso de datos públicos.

== BASES DE DATOS A VERIFICAR/EXPANDIR (confirma vigencia 2025–2026, cobertura y acceso; sugiere otras) ==
CPTAC/PDC; PRIDE/ProteomeXchange; TCGA COAD/READ (UCSC Xena, GDC); cBioPortal; UniProt; STRING; Reactome;
OmniPath; MSigDB (Hallmark/Reactome/KEGG); Gene Ontology/GOA; DepMap; Open Targets; DGIdb; DrugBank; ChEMBL;
Human Protein Atlas; GTEx. Añade recursos/benchmarks recientes que mejoren cualquiera de estas capas.

== FILTROS DE CALIDAD (aplícalos al seleccionar) ==
- Recencia: prioriza 2025–2026; 2024 solo si seminal/metodológico.
- Preferencia METODOLÓGICA: benchmarks, comparaciones de métodos, herramientas con código abierto y datos
  públicos, por encima de estudios de una sola cohorte sin replicación.
- Rigor: validación externa, análisis de sensibilidad, control de sobreajuste con N pequeño.
- Honestidad: incluye referencias que señalen LIMITACIONES (eliminación de señal por batch correction,
  sobreajuste de deep learning con N chico, sesgo de líneas celulares en DepMap, cobertura de bases de datos,
  generalización clínica y a poblaciones latinoamericanas) para mi capítulo de limitaciones.

== FORMATO DE RESPUESTA (obligatorio) ==
- Agrupa por los 12 temas.
- Por cada referencia: título; autores; año; revista/venue; DOI/enlace; 1–2 frases de por qué es relevante para
  MI metodología (qué capa A/B/C soporta, qué métrica o método valida, o qué limitación documenta).
- Etiqueta cada una: [DATOS] recurso/dataset descargable; [MÉTODO] metodológica; [MARCO] marco teórico;
  [VALIDACIÓN] valida diseño/resultados; [LIMITACIÓN] documenta un riesgo a declarar.
- Prioriza acceso abierto y herramientas con código.
- Termina con TRES secciones:
  (a) "Mapa metodológico": para cada capa (A homogeneización, B red, C priorización) y para el nivel peptídico,
      lista los 3–5 métodos/herramientas mejor respaldados y la métrica de evaluación recomendada.
  (b) "Factible vs futuro": qué es realista para una tesis de maestría con datos PÚBLICOS ya procesados vs qué
      requiere datos crudos/GPU/validación húmeda (trabajo futuro).
  (c) "Contradicciones y riesgos": hallazgos que contradigan mis supuestos o limitaciones que deba declarar
      explícitamente para mantener prudencia traslacional.
```

---

## 🇬🇧 Master prompt (ENGLISH — copy all)

```
Act as an expert assistant in translational bioinformatics, computational proteomics, systems biology and
multi-omics integration. Help me find recent literature (prioritize 2025–2026; include 2024 only if seminal or
methodological), methods and databases that VALIDATE AND SUPPORT the METHODOLOGICAL FRAMEWORK of my master's
thesis "OncoGO-CRC 2.0". IMPORTANT: my thesis is a reproducible ANALYSIS TOOL/METHODOLOGY that prioritizes
candidates as HYPOTHESES with convergent evidence; it does NOT promise validated therapeutic targets, clinical
biomarkers, or a cure. I seek rigorous methodological support and translational prudence, not clinical claims.

== WHAT ONCOGO-CRC IS (already implemented) ==
A reproducible pipeline integrating public colorectal cancer (CRC) proteomics: CPTAC colon (PDC000116, TMT10),
PRIDE tissue (LFQ) and plasma (DIA-NN, exploratory); protein→gene identifier harmonization via UniProt;
differential analysis (Welch/paired t-test + BH-FDR); multi-cohort reproducibility via UniProt bridge (core
signature with sign concordance, ~1967 intersecting genes, ~430 double-significant, magnitude Pearson ~0.10);
functional enrichment GO (ORA + pre-ranked GSEA, real human GOA); co-expression networks; interpretable ML
(L1/ElasticNet/RandomForest, nested CV, SHAP, bootstrap); full traceability (run_manifest with SHA-256 hashes).

== PROPOSED EVOLUTION (OncoGO-CRC 2.0) — 3 LAYERS I NEED LITERATURE FOR ==
Layer A. Robust homogenization & multi-omics integration: protein-level batch correction (ComBat/HarmonizR),
  proteome+transcriptome integration with MOFA+/VAE, and (advanced) adversarial domain adaptation for
  platform-invariant representations; evaluation with before/after metrics (kBET, iLISI/LISI, gPCA, silhouette,
  biological signal preservation).
Layer B. Network analysis (modules, not flat lists): classical networks (WGCNA/community detection) and Graph
  Neural Networks (GAT for small graphs, GCN/GraphSAGE for large, Graph Transformer) over PPI/pathways
  (STRING/Reactome/OmniPath), using the reproducible signature as node signal; interpretability via GNNExplainer
  and attention weights; enrichment with MSigDB (Hallmark/Reactome/KEGG).
Layer C. CONVERGENT-EVIDENCE prioritization framework (methodology, NOT a promise of targets): annotate each
  candidate with objective public evidence — CRISPR essentiality (DepMap), tractability/drugs (Open Targets,
  DGIdb, DrugBank, ChEMBL), tumor specificity/safety (Human Protein Atlas, GTEx), genomic context (cBioPortal) —
  with pinned provenance and version; output = table of exploratory prioritized HYPOTHESES.
Peptide level: consider peptide-level analysis using ALREADY-PUBLISHED PEPTIDE TABLES (peptide→protein rollup
  robustness, proteoforms/isoforms, detectability); raw-file reprocessing with deep learning (DIA-BERT, AlphaDIA)
  is framed as FUTURE WORK, not a promised deliverable.
Cross-cutting: interpretability (SHAP, Integrated Gradients, GNNExplainer) and reproducibility/FAIR (database
  version pinning, hashes, Zenodo/GitHub deposition).

== TOPICS TO SEARCH (group the answer EXACTLY by these 12) ==
1. Multi-omics integration with factor analysis (MOFA/MOFA+/MEFISTO) and multi-view autoencoders/VAE in cancer
   proteogenomics (proteome+transcriptome), 2024–2026.
2. PROTEIN-LEVEL batch-effect correction and cross-platform TMT/LFQ/DIA harmonization (ComBat, HarmonizR, RUV,
   recent methods), with benchmarks.
3. Integration/batch evaluation METRICS (kBET, iLISI/cLISI/LISI, gPCA, silhouette) and biological signal
   preservation; best practices to avoid removing real biology.
4. Domain adaptation / adversarial learning / contrastive learning for platform- or batch-invariant
   representations in omics.
5. PEPTIDE-LEVEL analysis: peptide→protein rollup methods, proteoforms/isoforms, peptide detectability, and the
   impact of imputation/rollup on differential lists (using already-processed data).
6. Deep learning for proteomic identification/quantification from raw data (e.g., DIA-BERT, AlphaDIA,
   transformer-type models) — ONLY as FUTURE-WORK context and its compute/data requirements.
7. Graph Neural Networks (GAT, GCN, GraphSAGE, Graph Transformer) over PPI/pathways (STRING/Reactome/OmniPath)
   for module/subnetwork discovery and cancer classification; GAT vs GCN by graph size.
8. Interpretability of GNNs and omics models (GNNExplainer, attention, Integrated Gradients, SHAP) and its
   validation.
9. Network/module analysis and systems-biology prioritization (WGCNA, community detection, MSigDB
   Hallmark/Reactome/KEGG enrichment) in CRC.
10. COMPUTATIONAL PRIORITIZATION FRAMEWORKS for candidates/targets REUSING public databases (DepMap CRISPR,
    Open Targets, DGIdb, DrugBank, ChEMBL) integrating essentiality + tractability; methodological focus, no
    clinical claims.
11. TUMOR SPECIFICITY and safety assessment of candidates with GTEx and Human Protein Atlas; genomic context
    with cBioPortal in CRC.
12. Reporting standards, interpretability, reproducibility and translational prudence (TRIPOD+AI, DOME, REMARK,
    FAIR principles) for pipelines reusing public data.

== DATABASES TO VERIFY/EXPAND (confirm 2025–2026 currency, coverage, access; suggest others) ==
CPTAC/PDC; PRIDE/ProteomeXchange; TCGA COAD/READ (UCSC Xena, GDC); cBioPortal; UniProt; STRING; Reactome;
OmniPath; MSigDB (Hallmark/Reactome/KEGG); Gene Ontology/GOA; DepMap; Open Targets; DGIdb; DrugBank; ChEMBL;
Human Protein Atlas; GTEx. Add recent resources/benchmarks that improve any of these layers.

== QUALITY FILTERS (apply when selecting) ==
- Recency: prioritize 2025–2026; 2024 only if seminal/methodological.
- METHODOLOGICAL preference: benchmarks, method comparisons, open-source tools with public data, over
  single-cohort studies without replication.
- Rigor: external validation, sensitivity analysis, overfitting control with small N.
- Honesty: include references highlighting LIMITATIONS (batch-correction removing signal, deep-learning
  overfitting at small N, cell-line bias in DepMap, database coverage, clinical and Latin-American generalization)
  for my limitations chapter.

== OUTPUT FORMAT (mandatory) ==
- Group by the 12 topics.
- For each reference: title; authors; year; journal/venue; DOI/link; 1–2 sentences on why it is relevant to MY
  methodology (which layer A/B/C it supports, which metric/method it validates, or which limitation it documents).
- Tag each: [DATA] downloadable resource/dataset; [METHOD] methodological; [FRAMEWORK] theoretical; [VALIDATION]
  validates design/results; [LIMITATION] documents a risk to declare.
- Prioritize open access and tools with code.
- End with THREE sections:
  (a) "Methodological map": for each layer (A homogenization, B network, C prioritization) and for the peptide
      level, list the 3–5 best-supported methods/tools and the recommended evaluation metric.
  (b) "Feasible vs future": what is realistic for a master's thesis using already-processed PUBLIC data vs what
      requires raw data/GPU/wet-lab validation (future work).
  (c) "Contradictions and risks": findings contradicting my assumptions or limitations I must state explicitly to
      keep translational prudence.
```
