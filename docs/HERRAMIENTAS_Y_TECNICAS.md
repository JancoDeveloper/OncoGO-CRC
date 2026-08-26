# OncoGO-CRC — herramientas, técnicas y procesos (catálogo exhaustivo)

Este documento cataloga **todas** las herramientas, librerías, técnicas y procesos usados (1.0) y propuestos
(2.0), explicando **para qué** se usan, **por qué** se eligieron y con **qué soporte bibliográfico/documental**.
Es la referencia del marco metodológico de la tesis.

> Convención: **[implementado]** = ya en el pipeline; **[2.0]** = diseñado para la evolución (ver
> `DISENO_TECNICO_ONCOGO_2.0.md` y `ARQUITECTURA_ONCOGO_2.0.md`); **[web]** = capa de visualización.

---

## 1. Lenguaje y entorno

| Herramienta | Para qué | Por qué | Referencia |
|---|---|---|---|
| **Python ≥3.10** [implementado] | Lenguaje base del pipeline y la web | Ecosistema científico maduro, reproducible y multiplataforma | python.org |
| **pip + `pyproject.toml`** [implementado] | Gestión de dependencias con *extras* (`go-remote`, `explain`, `survival`, `docs`, `web`) | Instalación reproducible y modular; solo instalas lo que usas | PEP 621 |
| **entornos virtuales (`.venv`)** [implementado] | Aislamiento de dependencias | Evita contaminación del sistema y fija versiones | PEP 405 |

---

## 2. Manipulación de datos y cómputo numérico

| Herramienta | Para qué | Por qué | Referencia |
|---|---|---|---|
| **pandas ≥2.0** [implementado] | Carga/transformación de matrices y metadata (TSV) | Estándar de facto para datos tabulares; API expresiva | McKinney 2010 |
| **NumPy ≥1.24** [implementado] | Álgebra/vectorización (correlaciones, normalización) | Base numérica eficiente | Harris et al. 2020, *Nature* |
| **SciPy ≥1.11** [implementado] | Pruebas estadísticas (t-test, correlaciones) | Implementaciones validadas y ampliamente citadas | Virtanen et al. 2020, *Nat Methods* |
| **statsmodels ≥0.14** [implementado] | Corrección de comparaciones múltiples (BH-FDR) | `multipletests` robusto y documentado | Seabold & Perktold 2010 |

---

## 3. Técnicas estadísticas (análisis diferencial y reproducibilidad)

| Técnica | Para qué | Por qué | Referencia |
|---|---|---|---|
| **log2 + z-score por proteína** [implementado] | Estabilizar varianza / comparabilidad entre muestras | Práctica estándar en proteómica cuantitativa | Aggarwal & Yadav 2016 |
| **Imputación KNN** [implementado] | Rellenar *missing* (MNAR/MCAR) sin sesgar tanto como media | Preserva estructura local; común en MS-proteómica | Troyanskaya et al. 2001 |
| **t-test de Welch** [implementado] | Contraste tumor vs normal sin asumir varianzas iguales | Robusto a varianzas desiguales y n distintos | Welch 1947 |
| **t-test pareado (mismo participante)** [implementado] | Aprovechar diseño tumor/normal emparejado | Mayor potencia al remover variabilidad inter-paciente | — |
| **Benjamini–Hochberg (FDR)** [implementado] | Controlar falsos positivos en miles de tests | Control de FDR estándar en ómicas | Benjamini & Hochberg 1995 |
| **Puente UniProt + concordancia de signo** [implementado] | Firma reproducible CPTAC∩PRIDE por clave canónica | Evita cruces espurios por identificadores heterogéneos | Robles et al. 2024 |

---

## 4. Armonización de identificadores

| Herramienta/proceso | Para qué | Por qué | Referencia |
|---|---|---|---|
| **UniProt (mapa maestro + bulk humano revisado con alias)** [implementado] | Resolver proteína→gen canónico (símbolo HGNC, accesión base) | Base de referencia de proteínas; alias evitan pérdidas por sinónimos | UniProt Consortium 2023, *NAR* |
| **`uniprot_base` / clave canónica** [implementado] | Colapsar isoformas (`P12345-2`) y grupos (`P1;Q2`) | Cruce 1:1 estable entre cohortes | — |
| **HGNC (símbolos de gen)** [implementado] | Vocabulario de gen unificado | Nomenclatura oficial de genes humanos | Tweedie et al. 2021, *NAR* |
| **diagnóstico de mapeo (`mapping_coverage`)** [implementado] | Auditar cobertura 1:1 vs lenient por dataset | Transparencia sobre cuánta señal se puede cruzar | — |

---

## 5. Enriquecimiento funcional

| Herramienta/técnica | Para qué | Por qué | Referencia |
|---|---|---|---|
| **Gene Ontology + GOA humano (`goa_human.gaf`)** [implementado] | Anotación funcional real (no sintética) | Estándar de anotación funcional; GOA es la fuente oficial humana | Ashburner et al. 2000; GO Consortium 2023 |
| **ORA (over-representation)** [implementado] | Enriquecimiento de la lista significativa | Interpretación simple de conjuntos sobre-representados | Boyle et al. 2004 |
| **GSEA pre-ranked (offline)** [implementado] | Enriquecimiento sobre ranking continuo (sin umbral duro) | Detecta señales coordinadas sutiles; reproducible offline | Subramanian et al. 2005, *PNAS* |
| **gseapy** [implementado, opcional] | Interfaz GSEA/Enrichr (KEGG) | Implementación en Python bien mantenida | Fang et al. 2023 |
| **KEGG (vía Enrichr)** [implementado, opcional] | Rutas metabólicas/de señalización | Complemento de rutas a GO | Kanehisa & Goto 2000 |
| **MSigDB (Hallmark/Reactome/KEGG)** [2.0] | Colecciones de conjuntos para módulos | Referencia curada para enriquecimiento de módulos | Liberzon et al. 2015 |

---

## 6. Redes y biología de sistemas

| Herramienta/técnica | Para qué | Por qué | Referencia |
|---|---|---|---|
| **networkx ≥3.1** [implementado] | Construcción de grafos y métricas (grado, hubs) | Librería de grafos estándar en Python | Hagberg et al. 2008 |
| **Red de coexpresión + hubs** [implementado] | Detectar proteínas centrales (candidatas reguladoras) | Los hubs suelen ser funcionalmente relevantes | Barabási & Oltvai 2004 |
| **STRING / Reactome / OmniPath** [2.0] | Priori de interacción proteína–proteína / rutas | Redes curadas de alta cobertura | Szklarczyk 2023; Gillespie 2022; Türei 2021 |
| **WGCNA (baseline de módulos)** [2.0] | Módulos de coexpresión interpretables | Método consolidado para módulos | Langfelder & Horvath 2008 |
| **Louvain / Leiden (comunidades)** [2.0] | Detección de módulos en el grafo | Leiden corrige defectos de Louvain | Traag et al. 2019 |
| **GNN (GAT/GCN/GraphSAGE)** [2.0] | Aprender representaciones sobre la red y descubrir submódulos | Aprovecha topología + señal ómica; explicable con GNNExplainer | Veličković 2018; Kipf 2017; Hamilton 2017; Ying 2019 |

---

## 7. Machine learning interpretable

| Herramienta/técnica | Para qué | Por qué | Referencia |
|---|---|---|---|
| **scikit-learn ≥1.3** [implementado] | Modelado (L1, ElasticNet, RandomForest), CV, calibración | API consistente, validada y reproducible | Pedregosa et al. 2011, *JMLR* |
| **Lasso / ElasticNet (L1)** [implementado] | Selección dispersa de proteínas | Firmas parsimoniosas e interpretables | Tibshirani 1996; Zou & Hastie 2005 |
| **Random Forest** [implementado] | Modelo no lineal robusto | Buen desempeño con pocas muestras/muchas features | Breiman 2001 |
| **Validación cruzada anidada** [implementado] | Estimación honesta de generalización | Evita optimismo por fuga de selección de hiperparámetros | Cawley & Talbot 2010 |
| **Transferencia + calibración (Platt/isotónica)** [implementado] | Evaluar en cohorte externa (PRIDE) | Mide generalización real entre cohortes | Platt 1999; Niculescu-Mizil 2005 |
| **SHAP** [implementado, opcional] | Explicabilidad de contribuciones por feature | Estándar de interpretabilidad basado en valores de Shapley | Lundberg & Lee 2017, *NeurIPS* |
| **Bootstrap de estabilidad de features** [implementado] | Cuantificar robustez de la firma | La estabilidad es criterio de confianza de biomarcadores | Meinshausen & Bühlmann 2010 |

---

## 8. Integración multi-ómica y homogeneización (Capa A) [2.0]

| Herramienta/técnica | Para qué | Por qué | Referencia |
|---|---|---|---|
| **MOFA+ (mofapy2)** | Factores latentes proteoma + transcriptoma | Integración multi-vista no supervisada e interpretable | Argelaguet et al. 2020, *Genome Biol* |
| **ComBat / HarmonizR** | Corrección de lote (batch) en proteómica | Referencia para batch; HarmonizR maneja *missingness* de MS | Johnson 2007; Voß et al. 2022 |
| **VAE / autoencoders multi-vista** | Homogeneización no lineal (futuro) | Captura estructura no lineal entre plataformas | Kingma & Welling 2014 |
| **kBET / iLISI / silhouette / gPCA** | Métricas de mezcla de lote vs señal biológica | Evaluar corrección **sin borrar biología** | Büttner 2019; Korsunsky 2019 |

---

## 9. Priorización por evidencia convergente (Capa C) [2.0, implementado como scaffold]

| Base de datos | Aporta a la priorización | Por qué | Referencia/acceso |
|---|---|---|---|
| **DepMap (CRISPR, Chronos)** | Esencialidad/dependencia génica en líneas CCR | Prioriza genes cuya pérdida afecta la viabilidad tumoral | Tsherniak 2017; depmap.org |
| **Open Targets** | Tractabilidad y asociación diana–enfermedad | Marco estándar de priorización de dianas | Ochoa 2023, *NAR*; platform.opentargets.org |
| **DGIdb / DrugBank / ChEMBL** | Interacciones gen–fármaco (druggability) | Contexto de accionabilidad farmacológica | Freshour 2021; Wishart 2018; Mendez 2019 |
| **Human Protein Atlas (HPA)** | Especificidad/localización tumoral | Especificidad reduce riesgo de toxicidad off-target | Uhlén 2015, *Science*; proteinatlas.org |
| **GTEx** | Expresión en tejidos normales (seguridad) | Baja expresión normal ↔ mejor ventana terapéutica | GTEx Consortium 2020, *Science* |
| **cBioPortal / TCGA** | Contexto genómico (frecuencia de alteración) | Ancla la hipótesis en genómica del tumor | Cerami 2012; Gao 2013 |

**Score de convergencia** [implementado como scaffold]: media ponderada **transparente** de señales normalizadas
disponibles (diferencial, reproducibilidad, esencialidad, tractabilidad, druggability, especificidad, seguridad,
genómica, hub de red). Diseño no-promisorio: reporta `evidence_count` y `provenance` (versión de cada fuente); la
evidencia ausente **no penaliza** artificialmente. Ver `DISENO_TECNICO_ONCOGO_2.0.md` §Capa C.

---

## 10. Evidencia externa y supervivencia [implementado, opcional]

| Herramienta | Para qué | Por qué | Referencia |
|---|---|---|---|
| **TCGA COAD/READ (UCSC Xena)** | Concordancia transcriptómica de la firma | Confirmación externa de dirección del efecto | Goldman 2020, *Nat Biotechnol* |
| **lifelines (Cox)** | HR de supervivencia por gen | Contexto pronóstico exploratorio | Davidson-Pilon 2019 |

---

## 11. Visualización y web

| Herramienta | Para qué | Por qué | Referencia |
|---|---|---|---|
| **matplotlib / seaborn** [implementado] | Figuras estáticas (volcano, PCA, barras) | Estándar de figuras científicas | Hunter 2007; Waskom 2021 |
| **MkDocs + Material** [implementado] | Sitio estático de documentación y resultados | Reproducible, versionable, fácil de publicar (GitHub Pages) | mkdocs.org; squidfunk Material |
| **Mermaid** [implementado] | Diagramas de arquitectura/flujo en la web | Diagramas como código, versionables | mermaid.js |
| **Streamlit** [web, scaffold] | Dashboard interactivo (tabla, volcano, red) | Prototipado rápido de apps de datos en Python | streamlit.io |
| **Plotly** [web, scaffold] | Gráficas interactivas (zoom, hover, filtro) | Interactividad rica sin JS manual | plotly.com |

**Por qué una web:** hacer los **procesos, análisis y resultados consultables** por cualquier persona
(revisores, clínicos, colaboradores) cumple principios FAIR de accesibilidad y favorece la reutilización.

---

## 12. Reproducibilidad, trazabilidad y estándares

| Herramienta/proceso | Para qué | Por qué | Referencia |
|---|---|---|---|
| **`run_manifest.json` (SHA-256, versiones)** [implementado] | Trazar qué se generó/omitió en cada corrida | Evita mezclar artefactos y garantiza auditoría | — |
| **`data_workspace` externo** [implementado] | Separar datos pesados del repo de código | Portabilidad y limpieza del repositorio | — |
| **Salvaguarda `--demo`** [implementado] | Impedir que datos sintéticos contaminen resultados reales | Integridad de resultados | — |
| **PyYAML (config declarativa)** [implementado] | Parámetros de pipeline versionados | Reproducibilidad y reejecución | yaml.org |
| **FAIR / DOME / TRIPOD+AI / REMARK** [guía] | Marco de reporte de ML/ómica y datos | Estándares aceptados de rigor y reporte | Wilkinson 2016; Walsh 2021; Collins 2024 |

---

## 13. Trazabilidad herramienta → etapa del pipeline

| Etapa (`registry.py`) | Herramientas/técnicas principales |
|---|---|
| `download_real` | requests/FTP, timeouts, manifests |
| `harmonize_public_metadata` | pandas, normalización clínica (age/stage/MSI/etnia/región) |
| `harmonize` (ids) | UniProt bulk + alias, clave canónica |
| `qc` | log2, KNN, z-score, PCA |
| `differential` | Welch/pareado + BH-FDR, volcano |
| `cohort_overlap` | puente UniProt, concordancia de signo |
| `go` / `gsea` / `kegg` | GOA humano, ORA, GSEA pre-ranked, gseapy |
| `network` | networkx, hubs |
| `ml` | scikit-learn, nested CV, SHAP, bootstrap |
| `external_evidence` | TCGA (Xena), Cox (lifelines) |
| **`prioritization` (Capa C)** | DepMap/OpenTargets/DGIdb/HPA/GTEx/cBioPortal + score convergente |
| `publish_results` | MkDocs, Plotly/Streamlit (web) |

---

*Referencias completas y validadas por componente: ver `DISENO_TECNICO_ONCOGO_2.0.md` y
`PROMPT_SCISPACE_BIBLIOGRAFIA_V2.md` (búsqueda bibliográfica dirigida).*
