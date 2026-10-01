# Diseño técnico OncoGO-CRC 2.0 — herramienta metodológica de análisis, armonización e priorización de hipótesis

> **Lectura vigente.** Este documento describe el diseño de las tres capas (A, B, C) tal como se pensó antes de la capa de certificación. Varias ideas de aquí quedaron **restringidas** en la actualización 3.0: ComBat y Harmony solo para figuras, no para inferencia; el *association score* de Open Targets está prohibido; un GNN no se interpreta si no supera a las líneas base y al nulo de rewiring. La guía que manda ahora es [GUIA_ESTUDIANTE_CERTIFICACION.md](GUIA_ESTUDIANTE_CERTIFICACION.md).

> **Encuadre (leer primero).** OncoGO-CRC es un **pipeline/metodología de investigación computacional, reproducible y reutilizable**. Su producto es un **análisis** y una **tabla de priorización de candidatos como HIPÓTESIS con evidencia convergente**, acompañada de trazabilidad completa.
> **La tesis NO promete** dianas terapéuticas validadas, biomarcadores listos para clínica, ni predicción/tratamiento de pacientes. Toda salida es **exploratoria y generadora de hipótesis**, y requiere validación experimental independiente que queda **fuera del alcance**. Este encuadre sigue la prudencia traslacional estándar en ómicas (p. ej. estándares TRIPOD+AI / DOME / REMARK y principios FAIR).

El objetivo de este documento es especificar, con detalle técnico, las **tres capas** de la evolución 2.0, las **bases de datos** que consume, las **métricas** de evaluación, las **referencias validadas** que respaldan cada decisión, y qué es **factible ahora** vs **trabajo futuro**, maximizando el impacto traslacional *sin comprometer* resultados clínicos.

---

## Principios transversales (aplican a las 3 capas)

1. **Interpretabilidad primero.** Ningún modelo caja negra sin explicación (SHAP / Integrated Gradients / GNNExplainer / pesos de atención) y sin procedencia de evidencia por candidato.
2. **Reproducibilidad y FAIR.** Todo paso se registra en `run_manifest.json` (hashes SHA-256), **incluyendo la versión/fecha de cada base de datos externa** (snapshot pinning). Salidas versionadas y depositables (Zenodo/GitHub).
3. **Encuadre no-promisorio.** Cada artefacto se etiqueta como *hipótesis exploratoria*; los textos y figuras evitan lenguaje clínico afirmativo.
4. **Validación honesta.** Toda capa reporta métricas *antes/después* y análisis de sensibilidad; se declaran limitaciones (heterogeneidad TMT/LFQ/DIA, N pequeño, sesgo de líneas celulares, representación latinoamericana limitada).

---

## Capa A — Homogeneización robusta e integración multi-ómica

**Problema que resuelve:** la concordancia de magnitud modesta (Pearson ~0.10) entre CPTAC (TMT) y PRIDE (LFQ/DIA), y la falta de un espacio común proteoma↔transcriptoma.

### A.1 Entradas
- Matrices proteicas armonizadas (CPTAC TMT10; PRIDE LFQ; plasma DIA-NN, exploratorio).
- Transcriptoma externo: TCGA COAD/READ (UCSC Xena / GDC).
- (Opcional) capa peptídica — ver A.4.

### A.2 Métodos (de base → avanzado)
| Método | Rol | Factibilidad |
|---|---|---|
| **ComBat / limma removeBatchEffect** | Corrección batch a nivel proteína (baseline) | Factible ya |
| **HarmonizR** | Armonización entre datasets independientes con manejo de faltantes | Factible ya |
| **MOFA+ (mofapy2)** | Factores latentes multi-ómicos (proteoma + RNA) interpretables | Factible ya |
| **Autoencoder / VAE multi-vista** | Representación latente compartida (reduce sesgo de plataforma) | Media |
| **Domain adaptation adversarial** | Representación *invariante a plataforma* (capítulo avanzado) | Media-alta |

> **Nota de rigor:** la corrección de batch puede **remover señal biológica** si el diseño está confundido con la plataforma. Por eso se evalúa preservación de señal (separación tumor/normal) *además* de mezcla de lotes.

### A.3 Métricas (antes/después de corrección)
- **Mezcla de lotes:** kBET (rechazo de mezcla), iLISI/LISI, gPCA (varianza atribuible a batch), silueta por plataforma.
- **Preservación biológica:** cLISI, silueta tumor/normal, retención de la separación en PCA/UMAP, correlación de la firma diferencial pre/post.
- **Integración MOFA+:** varianza explicada por factor y por vista; factores asociados a fenotipo (estadio, CMS) con test estadístico.
- **Distribuciones:** CV por proteína, % faltantes, efecto de imputación (benchmark KNN / missForest / low-rank).

### A.4 Nivel peptídico (incorporado, con alcance realista)
- **Factible ahora:** usar **tablas peptídicas ya publicadas** (reportes PSM/peptide de CPTAC/PDC; `peptides.txt` de MaxQuant en PRIDE cuando existan) para:
  - Evaluar robustez del *rollup* péptido→proteína (mediana vs maxLFQ vs iBAQ vs modelo).
  - Detectar **proteoformas/isoformas** (grupos de péptidos del mismo gen con comportamiento divergente).
  - Filtrar péptidos de baja fiabilidad (detectabilidad) antes del diferencial.
- **Trabajo futuro (fuera de alcance de la maestría):** reprocesado desde **archivos crudos** con modelos de deep learning de identificación (**DIA-BERT**, **AlphaDIA**). Son reales y potentes, pero exigen datos crudos (.raw/.d), GPU (A100/V100) y semanas de cómputo → se documentan como **línea abierta / doctorado**, no se prometen.

### A.5 Referencias validadas (Capa A)
- Voß et al. 2022, *Nat Commun* — **HarmonizR**. DOI: 10.1038/s41467-022-31007-x
- Qiaochu et al. 2025, *Nat Commun* — corrección batch a nivel proteína. DOI: 10.1038/s41467-025-64718-y
- Argelaguet et al. 2020, *Genome Biology* — **MOFA+** (verificar DOI: 10.1186/s13059-020-02015-1).
- Chen et al. 2025, *Molecular Biomedicine* — MOFA v2 en CPTAC-2 CRC (factor pronóstico CPLF). DOI: 10.1186/s43556-025-00386-0
- Büttner et al. 2019, *Nat Methods* — **kBET** (verificar DOI: 10.1038/s41592-018-0254-1).
- Korsunsky et al. 2019, *Nat Methods* — **Harmony / LISI** (verificar DOI: 10.1038/s41592-019-0619-0).
- Liu et al. 2025, *Nat Commun* 16:3530 — **DIA-BERT** (futuro). DOI: 10.1038/s41467-025-58866-4
- Wallmann et al. 2025, *Nat Biotechnol* — **AlphaDIA** transfer learning (futuro). DOI: 10.1038/s41587-025-02791-w

---

## Capa B — Análisis de red (módulos, no listas planas)

**Problema que resuelve:** interpretar la firma reproducible en el **contexto de la biología de sistemas** (interacciones, vías) en vez de proteínas aisladas.

### B.1 Entradas
- Firma reproducible (diferencial + concordancia de signo multi-cohorte) como **señal sobre nodos**.
- Grafos de conocimiento: **STRING**, **Reactome**, **OmniPath** (PPI + vías).

### B.2 Métodos
| Método | Uso | Nota |
|---|---|---|
| **WGCNA / detección de comunidades** (baseline robusto) | Módulos de coexpresión / subredes | Interpretable, estable con N pequeño |
| **GAT (Graph Attention)** | Grafos pequeños (<~1000 nodos) | Mejor extracción de features en grafos chicos (evidencia: arXiv 2302.12838) |
| **GCN / GraphSAGE** | Grafos grandes | Escala mejor a redes complejas |
| **Graph Transformer** | Exploratorio | Solo si N/estructura lo justifican |

> **Rigor con N pequeño:** el GNN se usa preferentemente en modo **no supervisado / semi-supervisado informado por la firma** (detección de módulos), no como clasificador con etiquetas escasas, para evitar sobreajuste. El baseline de red clásica (WGCNA/comunidades) se mantiene siempre como control.

### B.3 Salidas
- **Subredes/módulos** priorizados (p. ej. ejes estroma/EMT/CMS4), no listas planas.
- Enriquecimiento de cada módulo con **MSigDB (Hallmark/Reactome/KEGG)** — más estable que solo GO.

### B.4 Métricas
- Estabilidad de módulos (bootstrap de aristas/nodos; índice de reproducibilidad).
- Significancia por permutación (null de red preservando grado).
- Enriquecimiento funcional (FDR) por módulo; solapamiento con firmas conocidas.
- Interpretabilidad: **GNNExplainer**, importancia de nodos/aristas, pesos de atención.

### B.5 Referencias validadas (Capa B)
- Li et al. 2022, *Front Genet* 13:806842 — **MoGCN** (GCN multi-ómico, cáncer). DOI: 10.3389/fgene.2022.806842
- Multimodal GNN para subtipos de cáncer, arXiv:2302.12838 — GAT vs GCN según tamaño de grafo.
- ClassiGraph (repo `collaborativebioinformatics/ClassiGraph`) — GNN sobre CPTAC-COAD + PPI + anotaciones CMS (referencia de implementación, no herramienta madura).
- STRING, Reactome, OmniPath, MSigDB — recursos de red/vías (verificar versiones 2025–2026).

---

## Capa C — Priorización traslacional 2.0 (marco de evidencia convergente)

**Problema que resuelve:** ordenar la firma por **relevancia biológica y traslacional** de forma transparente, **como marco de análisis**, sin afirmar que los candidatos sean dianas válidas.

> **Encuadre explícito de la Capa C:** es un **marco de convergencia de evidencia** que anota cada candidato con banderas objetivas provenientes de bases públicas. El resultado es una **tabla de hipótesis priorizadas**, con procedencia y versión de cada fuente. **No** constituye recomendación clínica ni validación de diana.

### C.1 Capas de evidencia (todas de bases públicas)
| Pregunta | Fuente | Salida |
|---|---|---|
| ¿Es esencial en líneas CRC? | **DepMap** (CRISPR-Cas9) | score de dependencia por línea CRC |
| ¿Es "drogable"/hay fármaco? | **Open Targets** (tractability), **DGIdb**, **DrugBank**, **ChEMBL** | banderas de tractabilidad / interacciones fármaco-gen |
| ¿Es tumor-específica / segura? | **Human Protein Atlas**, **GTEx** | especificidad tumoral, expresión en tejidos normales |
| ¿Contexto genómico? | **cBioPortal** | mutación/CNV en CRC |
| ¿Localización/pronóstico? | **Human Protein Atlas** | localización subcelular, asociación pronóstica |

### C.2 Método
- **Score de convergencia transparente** (pesos explícitos y auditables), NO un único número mágico: cada candidato lleva sus banderas de evidencia y su procedencia.
- **Análisis de sensibilidad** de pesos (¿el ranking es estable ante cambios razonables de pesos?).
- (Avanzado, opcional) meta-modelo que aprende pesos, **siempre** con interpretabilidad y validación de estabilidad.

### C.3 Métricas
- Cobertura de la firma anotada por cada base (¿qué % tiene evidencia?).
- Enriquecimiento vs conjuntos aleatorios de genes (¿la firma está enriquecida en genes esenciales/drogables más que el azar?).
- Estabilidad del ranking (bootstrap) y sensibilidad a pesos.

### C.4 Referencias validadas (Capa C)
- Pacini/Garnett et al. 2024, *Cancer Cell* — mapa de dependencias de 2ª generación y **marco de priorización de dianas** (DepMap + PPI). DOI: 10.1016/j.ccell.2023.11.009 *(verificar)*; nota de prensa Open Targets (11-ene-2024) sobre 370 dianas candidatas.
- Revisión de uso de bases públicas para evaluación de dianas oncológicas, *PMC10346432* (Open Targets, DepMap, DGIdb, tractability).
- Frontiers Immunol 2025, marco integrativo de priorización de genes drogables en CRC (OpenTargets + DrugBank + DGIdb). DOI: 10.3389/fimmu.2025.1604154
- DepMap portal (Broad/Sanger); Open Targets Platform; DGIdb; DrugBank; ChEMBL; Human Protein Atlas; GTEx; cBioPortal (verificar versiones 2025–2026).

---

## Bases de datos incorporadas (todas gratuitas)

| Categoría | Recursos | Uso en el pipeline |
|---|---|---|
| Proteómica | CPTAC/PDC, PRIDE/ProteomeXchange | datos primarios |
| Transcriptoma/clínica | TCGA COAD/READ (UCSC Xena/GDC), cBioPortal | validación externa, contexto genómico |
| Identificadores | UniProt (+ alias) | armonización |
| Vías/red | STRING, Reactome, OmniPath, MSigDB (Hallmark/Reactome/KEGG), GO/GOA | Capa B, enriquecimiento |
| Dianas/fármacos | DepMap, Open Targets, DGIdb, DrugBank, ChEMBL | Capa C |
| Especificidad/seguridad | Human Protein Atlas, GTEx | Capa C |

**Regla de reproducibilidad:** cada base se consume con **versión/fecha fijada** y se registra en `run_manifest.json`.

---

## Factible ahora vs Trabajo futuro (resumen)

| Componente | Factible en la maestría | Trabajo futuro / doctorado |
|---|---|---|
| ComBat/HarmonizR + métricas kBET/iLISI/gPCA | ✅ | — |
| MOFA+ (proteoma + TCGA) | ✅ | MEFISTO / VAE multi-vista más complejo |
| Domain adaptation adversarial | ⚠️ capítulo avanzado | versión robusta multi-dataset |
| Nivel peptídico desde **tablas publicadas** (rollup, proteoformas) | ✅ acotado | — |
| Reprocesado **crudo** DL (DIA-BERT/AlphaDIA) | ❌ | ✅ (GPU + crudos) |
| Red clásica (WGCNA/comunidades) | ✅ | — |
| GNN (GAT/GCN) informado por firma | ✅ exploratorio | GNN supervisado a gran escala |
| Marco de convergencia (DepMap/OpenTargets/DGIdb/HPA/GTEx) | ✅ | meta-modelo aprendido de pesos |
| Validación experimental de candidatos | ❌ (fuera de alcance) | ✅ colaboración húmeda |

---

## Cómo se comunica el impacto (sin comprometer)

- **Contribución de la tesis:** una **metodología reproducible y reutilizable** que integra, homogeneiza y prioriza señales proteómicas públicas de CCR, con evidencia convergente y trazabilidad FAIR.
- **Entregable:** el *pipeline* + una **tabla de hipótesis priorizadas** (candidatos anotados por evidencia) como **base para futuras investigaciones**, no como conjunto de dianas validadas.
- **Frase tipo para la tesis:** «OncoGO-CRC prioriza, de forma exploratoria y reproducible, candidatos proteicos con evidencia convergente multi-fuente, entendidos como **hipótesis** que requieren validación experimental independiente.»
