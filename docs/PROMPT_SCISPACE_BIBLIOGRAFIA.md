# Prompt único para SciSpace — bibliografía, bases de datos y métodos (OncoGO-CRC)

Dos versiones **autocontenidas** (español e inglés). Cada una incluye TODO en un solo bloque
(contexto, resultados, temas de búsqueda, bases de datos, filtros de calidad, formato de salida y
mapa de hallazgos) para **optimizar créditos y tokens** de [SciSpace](https://scispace.com/): pega **un**
bloque, sin necesidad de otros archivos. Prioriza literatura **2025–2026** (2024 si es seminal).

---

## 🇪🇸 Prompt maestro (ESPAÑOL — copiar todo)

```
Actúa como asistente experto en oncología molecular, proteómica clínica y bioinformática traslacional.
Ayúdame a encontrar literatura reciente (prioriza 2025–2026; incluye 2024 solo si es seminal), bases de
datos y métodos que SOPORTEN el marco teórico, VALIDEN los resultados y AUMENTEN la utilidad traslacional
de mi tesis de maestría "OncoGO-CRC". Busco generar hipótesis rigurosas y una base metodológica reutilizable,
NO biomarcadores clínicos listos para uso.

== QUÉ ES ONCOGO-CRC (ya implementado) ==
Pipeline reproducible que integra proteómica pública de cáncer colorrectal (CCR):
- Cohortes: CPTAC colon (PDC000116, TMT10, tumor/normal), PRIDE tejido (MaxQuant LFQ) y plasma/suero
  (DIA-NN) como capa exploratoria.
- Armonización de identificadores proteína→gen con UniProt (mapa maestro + bulk humano revisado con alias)
  para puentear CPTAC (símbolos gen) con PRIDE (accesiones UniProt).
- Diferencial (t-test Welch/pareado + FDR-BH); reproducibilidad multi-cohorte por puente UniProt (firma
  núcleo con concordancia de signo de log2FC); enriquecimiento GO (ORA + GSEA pre-ranked con GOA humana real);
  redes de coexpresión; ML interpretable (L1/ElasticNet/RandomForest, nested cross-validation, SHAP, bootstrap
  de estabilidad) y un "score de priorización traslacional" (TPS).
- Metadata clínica CPTAC: estadio AJCC I–IV, sexo, etnicidad NIH (hispano/latino), raza.
- Concordancia externa en desarrollo: TCGA COAD/READ (RNA-seq log2FC y Cox de supervivencia vía UCSC Xena)
  y subtipos moleculares consenso (CMS).
- Trazabilidad reproducible completa (run_manifest con hashes SHA-256).

== RESULTADOS PRELIMINARES (para orientar la búsqueda) ==
- Descubrimiento diferencial robusto en tejido CPTAC tumor vs normal (separación clara en PCA).
- Reproducibilidad parcial CPTAC∩PRIDE: ~6000 y ~6700 proteínas con puente UniProt, ~1967 genes en
  intersección, ~430 doble-significativos (padj<0.05), concordancia de signo ~0.56; concordancia de magnitud
  modesta (Pearson ~0.10) coherente con heterogeneidad TMT vs LFQ/DIA.
- Plasma: señal más corta pero plausible, estrictamente exploratoria.
- ML interno con nested CV + bootstrap + SHAP produce rankings estables; TPS integra capas curadas.

== MARCO CONCEPTUAL BASE (ya usado) ==
Subtipos CMS (Guinney 2015, Nat Med); proteogenómica colon/recto y desacople proteoma–transcriptoma
(Zhang 2014, Nature; CPTAC); reuso riguroso de proteómica pública de CCR con armonización de IDs
(Robles 2024, PLOS Comput Biol); proteogenómica pan-cáncer CPTAC (Cancer Cell 2023); proteómica plasmática
prudente (Nat Commun 2024; Cell Rep Med 2025); contexto latinoamericano/colombiano, CCR de inicio temprano
(EOCRC), MSI y biopsia líquida (García-Foncillas 2026); prudencia traslacional en ómicas (Nature 2013;
estándares TRIPOD/REMARK/DOME).

== TEMAS A BUSCAR (agrupa la respuesta por estos 10) ==
1. Proteómica de CCR (tejido y plasma) 2025–2026 que reutilice/integre datasets públicos (CPTAC, PRIDE,
   ProteomeXchange) y reporte firmas reproducibles entre cohortes.
2. Métodos 2024–2026 de armonización de IDs proteína/gen, integración multi-cohorte, corrección de efecto
   batch, y evaluación de reproducibilidad entre TMT/LFQ/DIA.
3. Subtipos moleculares del CCR a nivel PROTEÓMICO (CMS proteómico, iCMS, estroma/EMT) y su relación con
   pronóstico y microambiente tumoral.
4. Cáncer colorrectal de inicio temprano (EOCRC): biología/proteómica diferencial vs inicio tardío;
   epidemiología en Latinoamérica/Colombia.
5. Concordancia multi-ómica (proteoma vs transcriptoma vs supervivencia); uso de TCGA COAD/READ como
   validación externa; priorización por evidencia convergente.
6. ML interpretable (SHAP y afines) en ómicas oncológicas: validación externa, calibración, transferencia
   entre cohortes, y estándares de reporte (TRIPOD+AI, DOME, REMARK).
7. Proteómica plasmática / biopsia líquida en CCR: reproducibilidad, limitaciones, generación de hipótesis.
8. Bases de datos y recursos públicos actualizados útiles para este proyecto (ver lista abajo) y equivalentes
   más recientes o mejores.
9. Atlas de célula única / espacial de CCR 2025–2026 y cohortes latinoamericanas de CCR.
10. Vacíos de conocimiento (gaps) que una tesis computacional de reuso de proteómica pública de CCR puede
    llenar de forma realista, y cómo maximizar su utilidad e impacto oncológico.

== BASES DE DATOS A VERIFICAR/EXPANDIR (confirma vigencia 2025–2026, cobertura y acceso; sugiere otras) ==
CPTAC/PDC; PRIDE/ProteomeXchange; TCGA COAD/READ (UCSC Xena, GDC); cBioPortal; GEPIA2; Human Protein Atlas;
STRING; Reactome; MSigDB (Hallmark); KEGG; Gene Ontology/GOA; UniProt; DepMap; Open Targets; DGIdb; Pharos;
TCIA (CPTAC-COAD imaging). Añade recursos de single-cell/espacial y cohortes latinoamericanas si existen.

== FILTROS DE CALIDAD (aplícalos al seleccionar) ==
- Recencia: prioriza 2025–2026; 2024 solo si seminal/metodológico.
- Reutilizabilidad: favorece estudios con datos públicos y código disponible.
- Rigor: revisiones sistemáticas, consorcios (CPTAC) y trabajos con validación externa por encima de
  estudios de una sola cohorte sin replicación.
- Pertinencia traslacional: que aporten a EOCRC, CMS, microambiente/estroma o contexto latinoamericano.
- Honestidad: incluye referencias que señalen LIMITACIONES del reuso de proteómica pública (heterogeneidad
  TMT/LFQ/DIA, efecto batch, generalización clínica) para mi capítulo de limitaciones.

== FORMATO DE RESPUESTA (obligatorio) ==
- Agrupa por los 10 temas.
- Por cada referencia: título; autores; año; revista/venue; DOI/enlace; 1–2 frases de por qué es relevante
  para MI proyecto (qué parte del pipeline soporta o qué gap cierra).
- Etiqueta cada una: [DATOS] si ofrece dataset/recurso descargable, [MÉTODO] si es metodológica,
  [MARCO] si es para marco teórico, [VALIDACIÓN] si sirve para validar resultados.
- Prioriza acceso abierto y datos reutilizables.
- Termina con dos secciones:
  (a) "Cómo conectar con OncoGO-CRC": 5–8 acciones concretas (dataset a añadir, método a probar,
      limitación a declarar), y a qué capítulo alimenta cada hallazgo (marco teórico / métodos / discusión /
      limitaciones).
  (b) "Contradicciones y riesgos": cualquier hallazgo que contradiga mis supuestos o sugiera una limitación
      que deba declarar explícitamente.
```

---

## 🇬🇧 Master prompt (ENGLISH — copy all)

```
Act as an expert assistant in molecular oncology, clinical proteomics and translational bioinformatics.
Help me find recent literature (prioritize 2025–2026; include 2024 only if seminal), databases and methods
that SUPPORT the theoretical framework, VALIDATE the results and INCREASE the translational usefulness of my
master's thesis "OncoGO-CRC". I aim to generate rigorous hypotheses and a reusable methodological base,
NOT clinically ready biomarkers.

== WHAT ONCOGO-CRC IS (already implemented) ==
A reproducible pipeline integrating public colorectal cancer (CRC) proteomics:
- Cohorts: CPTAC colon (PDC000116, TMT10, tumor/normal), PRIDE tissue (MaxQuant LFQ) and plasma/serum
  (DIA-NN) as an exploratory layer.
- Protein→gene identifier harmonization via UniProt (master map + reviewed human bulk with aliases) to
  bridge CPTAC (gene symbols) with PRIDE (UniProt accessions).
- Differential analysis (Welch/paired t-test + BH-FDR); multi-cohort reproducibility via UniProt bridge
  (core signature with log2FC sign concordance); functional enrichment GO (ORA + pre-ranked GSEA with real
  human GOA); co-expression networks; interpretable ML (L1/ElasticNet/RandomForest, nested cross-validation,
  SHAP, stability bootstrap) and a "translational prioritization score" (TPS).
- CPTAC clinical metadata: AJCC stage I–IV, sex, NIH ethnicity (Hispanic/Latino), race.
- External concordance in progress: TCGA COAD/READ (RNA-seq log2FC and survival Cox via UCSC Xena) and
  consensus molecular subtypes (CMS).
- Full reproducible traceability (run_manifest with SHA-256 hashes).

== PRELIMINARY RESULTS (to guide the search) ==
- Robust differential discovery in CPTAC tumor vs normal tissue (clear PCA separation).
- Partial CPTAC∩PRIDE reproducibility: ~6000 and ~6700 proteins with UniProt bridge, ~1967 genes in the
  intersection, ~430 double-significant (padj<0.05), sign concordance ~0.56; modest magnitude concordance
  (Pearson ~0.10) consistent with TMT vs LFQ/DIA heterogeneity.
- Plasma: shorter but plausible signal, strictly exploratory.
- Internal ML with nested CV + bootstrap + SHAP yields stable rankings; TPS integrates curated layers.

== BASE CONCEPTUAL FRAMEWORK (already used) ==
CMS subtypes (Guinney 2015, Nat Med); colon/rectal proteogenomics and proteome–transcriptome decoupling
(Zhang 2014, Nature; CPTAC); rigorous reuse of public CRC proteomics with ID harmonization (Robles 2024,
PLOS Comput Biol); CPTAC pan-cancer proteogenomics (Cancer Cell 2023); prudent plasma proteomics (Nat Commun
2024; Cell Rep Med 2025); Latin American/Colombian context, early-onset CRC (EOCRC), MSI and liquid biopsy
(García-Foncillas 2026); translational prudence in omics (Nature 2013; TRIPOD/REMARK/DOME standards).

== TOPICS TO SEARCH (group the answer by these 10) ==
1. CRC proteomics (tissue and plasma) 2025–2026 reusing/integrating public datasets (CPTAC, PRIDE,
   ProteomeXchange) and reporting cross-cohort reproducible signatures.
2. Methods 2024–2026 for protein/gene ID harmonization, multi-cohort integration, batch-effect correction,
   and reproducibility assessment across TMT/LFQ/DIA.
3. Proteomic-level CRC molecular subtypes (proteomic CMS, iCMS, stroma/EMT) and their link to prognosis and
   tumor microenvironment.
4. Early-onset CRC (EOCRC): molecular/proteomic differences vs late-onset; epidemiology in Latin America/
   Colombia.
5. Multi-omic concordance (proteome vs transcriptome vs survival); TCGA COAD/READ as external validation;
   prioritization by convergent evidence.
6. Interpretable ML (SHAP and similar) in cancer omics: external validation, calibration, cross-cohort
   transfer, and reporting standards (TRIPOD+AI, DOME, REMARK).
7. Plasma proteomics / liquid biopsy in CRC: reproducibility, limitations, hypothesis generation.
8. Updated public databases and resources useful for this project (see list below) and better/newer
   equivalents.
9. Single-cell / spatial CRC atlases 2025–2026 and Latin American CRC cohorts.
10. Knowledge gaps that a computational thesis reusing public CRC proteomics can realistically fill, and how
    to maximize its oncological usefulness and impact.

== DATABASES TO VERIFY/EXPAND (confirm 2025–2026 currency, coverage, access; suggest others) ==
CPTAC/PDC; PRIDE/ProteomeXchange; TCGA COAD/READ (UCSC Xena, GDC); cBioPortal; GEPIA2; Human Protein Atlas;
STRING; Reactome; MSigDB (Hallmark); KEGG; Gene Ontology/GOA; UniProt; DepMap; Open Targets; DGIdb; Pharos;
TCIA (CPTAC-COAD imaging). Add single-cell/spatial resources and Latin American cohorts if available.

== QUALITY FILTERS (apply when selecting) ==
- Recency: prioritize 2025–2026; 2024 only if seminal/methodological.
- Reusability: favor studies with public data and available code.
- Rigor: systematic reviews, consortia (CPTAC) and externally validated works over single-cohort studies
  without replication.
- Translational relevance: contributing to EOCRC, CMS, microenvironment/stroma or Latin American context.
- Honesty: include references highlighting LIMITATIONS of public proteomics reuse (TMT/LFQ/DIA heterogeneity,
  batch effect, clinical generalization) for my limitations chapter.

== OUTPUT FORMAT (mandatory) ==
- Group by the 10 topics.
- For each reference: title; authors; year; journal/venue; DOI/link; 1–2 sentences on why it is relevant to
  MY project (which pipeline part it supports or which gap it closes).
- Tag each: [DATA] downloadable dataset/resource, [METHOD] methodological, [FRAMEWORK] theoretical framework,
  [VALIDATION] useful to validate results.
- Prioritize open access and reusable data.
- End with two sections:
  (a) "How to connect with OncoGO-CRC": 5–8 concrete actions (dataset to add, method to test, limitation to
      declare), and which chapter each finding feeds (framework / methods / discussion / limitations).
  (b) "Contradictions and risks": any finding contradicting my assumptions or suggesting a limitation I must
      state explicitly.
```
