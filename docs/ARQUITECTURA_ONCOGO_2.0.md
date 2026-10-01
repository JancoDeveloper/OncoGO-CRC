# Arquitectura y diseño de OncoGO-CRC 2.0 (pre-implementación)

> **Actualización 3.0.** La arquitectura de stages de este documento sigue siendo el mapa del código. La capa de certificación se añadió como etapa `certification`, después de `prioritization`. Las reglas nuevas (no fusionar cohortes, no imputar para inferir, índice R, auditoría de independencia) están en [GUIA_ESTUDIANTE_CERTIFICACION.md](GUIA_ESTUDIANTE_CERTIFICACION.md).

> **Estado:** documento de **diseño previo a la implementación**. Define arquitectura, diagramas, flujo de
> datos y la **capa web** para consulta pública. Marca qué **ya existe** en OncoGO-CRC 1.0 y qué es **nuevo**
> (2.0). Los diagramas están en [Mermaid](https://mermaid.js.org/) y se renderizan en la web (MkDocs Material).

> **Encuadre (invariante).** OncoGO-CRC 2.0 es una **metodología/herramienta de análisis reproducible** que
> **prioriza candidatos como HIPÓTESIS con evidencia convergente**. No promete dianas validadas ni uso clínico.
> Toda la arquitectura está diseñada para hacer explícita esa prudencia (trazabilidad, versión de datos,
> etiquetado de salidas como exploratorias).

---

## 1. Visión general

OncoGO-CRC 2.0 se organiza en **tres capas metodológicas** sobre el pipeline base:

- **Capa A — Homogeneización e integración multi-ómica** (corrección batch + MOFA+/VAE + métricas).
- **Capa B — Análisis de red** (módulos/subredes con GNN + WGCNA sobre STRING/Reactome).
- **Capa C — Priorización por evidencia convergente** (DepMap, Open Targets, DGIdb/DrugBank/ChEMBL, HPA, GTEx, cBioPortal).

Todo se ejecuta desde un **orquestador por etapas (stages)**, se registra en `run_manifest.json` (hashes
SHA-256 + versiones de bases de datos) y se publica en una **capa web** consultable por cualquier persona.

### 1.1 Diagrama de contexto (nivel sistema)

```mermaid
flowchart TB
    user["Investigador / Tesista"]
    public["Publico / Jurado / Oncologo / Clinicos"]

    subgraph ext["Fuentes de datos publicas (versionadas)"]
      cptac["CPTAC / PDC (TMT10)"]
      pride["PRIDE / ProteomeXchange (LFQ/DIA)"]
      tcga["TCGA COAD-READ (Xena / GDC)"]
      net["STRING / Reactome / OmniPath / MSigDB"]
      dep["DepMap (CRISPR)"]
      drug["Open Targets / DGIdb / DrugBank / ChEMBL"]
      safe["Human Protein Atlas / GTEx / cBioPortal"]
    end

    sys["OncoGO-CRC 2.0<br/>pipeline reproducible + capa web"]

    user -->|configura y ejecuta| sys
    ext -->|descarga con version fijada| sys
    sys -->|reporte web reproducible + tablas de hipotesis| public
    sys -->|artefactos FAIR + run_manifest| user
```

### 1.2 Diagrama de contenedores (nivel arquitectura)

```mermaid
flowchart LR
    subgraph orch["Orquestador (Python CLI)"]
      conf["config.yaml<br/>comparisons.yaml"]
      engine["Motor de stages"]
      manifest["run_manifest.json<br/>SHA-256 + versiones DB"]
    end

    subgraph store["Almacenamiento (data_workspace)"]
      raw["data/raw"]
      results["results/ (tablas + figuras + JSON)"]
      metad["metadata/ (mapas + snapshots DB)"]
    end

    subgraph webc["Capa Web"]
      site["Sitio estatico<br/>MkDocs Material"]
      dash["Dashboard interactivo<br/>(opcional: Streamlit/Dash)"]
      api["API JSON opcional<br/>(FastAPI)"]
    end

    conf --> engine --> manifest
    engine --> store
    results --> site
    results --> dash
    results --> api
    manifest --> site
```

---

## 2. Arquitectura del pipeline (end-to-end)

### 2.1 Flujo maestro por etapas

```mermaid
flowchart TB
    A0["download_real<br/>(CPTAC/PRIDE/plasma + DB snapshots)"]:::exist
    A1["qc<br/>(missingness, PCA, outliers)"]:::exist
    A2["harmonize_ids<br/>(UniProt master + alias)"]:::exist
    A3["differential<br/>(Welch/pareado + BH-FDR)"]:::exist
    A4["cohort_overlap<br/>(firma reproducible multi-cohorte)"]:::exist

    subgraph LA["Capa A - Homogeneizacion e integracion (NUEVO)"]
      LA1["batch_correct<br/>ComBat / HarmonizR"]:::new
      LA2["mofa_integrate<br/>MOFA+ proteoma + TCGA RNA"]:::new
      LA3["integration_metrics<br/>kBET / iLISI / silueta / gPCA"]:::new
    end

    subgraph LB["Capa B - Analisis de red (NUEVO)"]
      LB1["build_graph<br/>STRING / Reactome / OmniPath"]:::new
      LB2["modules_gnn<br/>GAT/GCN + community detection"]:::new
      LB3["modules_wgcna<br/>coexpresion (baseline)"]:::new
      LB4["module_enrichment<br/>MSigDB Hallmark/Reactome/KEGG"]:::new
    end

    subgraph LC["Capa C - Priorizacion por evidencia convergente (NUEVO)"]
      LC1["evidence_join<br/>DepMap/OpenTargets/DGIdb/HPA/GTEx/cBioPortal"]:::new
      LC2["convergence_score<br/>score transparente + sensibilidad"]:::new
      LC3["hypothesis_table<br/>tabla de hipotesis priorizadas"]:::new
    end

    A5["go_gsea + ml_prioritization"]:::exist
    P["publish_results<br/>(web + figuras + manifest)"]:::exist

    A0 --> A1 --> A2 --> A3 --> A4
    A4 --> LA1 --> LA2 --> LA3
    LA3 --> LB1 --> LB2
    LB1 --> LB3
    LB2 --> LB4
    LB3 --> LB4
    A4 --> A5
    LB4 --> LC1 --> LC2 --> LC3
    A5 --> LC1
    LC3 --> P
    LA3 --> P
    LB4 --> P

    classDef exist fill:#e3f2fd,stroke:#1565c0,color:#0d47a1;
    classDef new fill:#e8f5e9,stroke:#2e7d32,color:#1b5e20;
```

Azul = ya existe en 1.0. Verde = nuevo en 2.0.

### 2.2 Flujo de datos (qué produce cada etapa)

```mermaid
flowchart LR
    raw["Matrices crudas<br/>(procesadas de PDC/PRIDE)"] --> norm["Matrices normalizadas<br/>log2 / zscore"]
    norm --> diff["*_diff.tsv<br/>log2FC + FDR"]
    diff --> sign["core_signature.tsv<br/>concordancia de signo"]
    sign --> corr["matrices corregidas<br/>+ metrics.json"]
    corr --> factors["MOFA factors<br/>+ varianza explicada"]
    sign --> graph["grafo.graphml"]
    graph --> mods["modules.tsv<br/>+ enrichment.tsv"]
    factors --> ev
    mods --> ev["evidence_matrix.tsv"]
    ev --> score["hypothesis_priority.tsv<br/>(candidatos + banderas + score)"]
    score --> web["Reporte web + figuras"]
    corr --> web
    mods --> web
```

---

## 3. Detalle por capa

### 3.1 Capa A — Homogeneización e integración

```mermaid
flowchart TB
    in["Matrices por cohorte<br/>CPTAC / PRIDE / plasma + TCGA RNA"]
    bc{"Estrategia batch"}
    combat["ComBat (protein-level)"]
    harmo["HarmonizR (rescate de faltantes)"]
    adv["(avanzado) VAE / domain adaptation"]
    mofa["MOFA+ : factores latentes"]
    met["Metricas antes/despues"]
    guard{"Preserva senal biologica?<br/>(tumor vs normal)"}
    out["Matriz integrada + factores + reporte de metricas"]

    in --> bc
    bc --> combat --> mofa
    bc --> harmo --> mofa
    bc -.opcional.-> adv --> mofa
    mofa --> met --> guard
    guard -->|si| out
    guard -->|no: sobrecorreccion| bc
```

Métricas: kBET, iLISI/cLISI, ASW batch/bio, gPCA, silueta; **detección de sobrecorrección** (RBET / inspección
de marcadores CRC conocidos) para no borrar biología real.

### 3.2 Capa B — Análisis de red

```mermaid
flowchart TB
    sig["Firma reproducible<br/>(señal sobre nodos)"]
    kb["Grafo de conocimiento<br/>STRING (>=0.4) / Reactome / OmniPath"]
    sz{"Tamano del grafo"}
    gat["GAT (<~1000 nodos)"]
    gcn["GCN / GraphSAGE (grandes)"]
    wg["WGCNA / comunidades (baseline)"]
    comm["Deteccion de comunidades<br/>Louvain / Leiden"]
    enr["Enriquecimiento<br/>MSigDB Hallmark/Reactome/KEGG"]
    exp["Interpretabilidad<br/>GNNExplainer / atencion"]
    out["Modulos / subredes priorizadas<br/>(p.ej. estroma / EMT / CMS4)"]

    sig --> kb --> sz
    sz --> gat --> comm
    sz --> gcn --> comm
    sig --> wg --> enr
    comm --> enr
    gat --> exp
    gcn --> exp
    enr --> out
    exp --> out
```

### 3.3 Capa C — Priorización por evidencia convergente

```mermaid
flowchart TB
    cand["Candidatos<br/>(firma + modulos + ML)"]

    subgraph EV["Anotacion de evidencia (bases publicas versionadas)"]
      e1["DepMap<br/>esencialidad CRISPR (Chronos)"]
      e2["Open Targets<br/>tractabilidad / seguridad"]
      e3["DGIdb / DrugBank / ChEMBL<br/>farmacos"]
      e4["HPA + GTEx<br/>especificidad tumoral / seguridad"]
      e5["cBioPortal<br/>contexto genomico (mut/CNV)"]
    end

    join["evidence_matrix<br/>(candidato x evidencia + procedencia)"]
    score["score de convergencia transparente<br/>+ analisis de sensibilidad de pesos"]
    guard{"Validacion:<br/>dianas CRC conocidas rankean alto?"}
    table["hypothesis_priority.tsv<br/>tabla de HIPOTESIS exploratorias"]

    cand --> e1 --> join
    cand --> e2 --> join
    cand --> e3 --> join
    cand --> e4 --> join
    cand --> e5 --> join
    join --> score --> guard
    guard -->|si| table
    guard -->|no| score
```

> Cada fila lleva **banderas de evidencia + procedencia (base y versión)** y se etiqueta como *hipótesis*.

### 3.4 Ciclo de vida de un candidato (estados)

```mermaid
stateDiagram-v2
    [*] --> Detectado: diferencial + FDR
    Detectado --> Reproducible: concordancia multi-cohorte
    Reproducible --> EnModulo: pertenece a subred/modulo
    EnModulo --> Anotado: evidencia (DepMap/OT/HPA/GTEx)
    Anotado --> Priorizado: score de convergencia
    Priorizado --> HipotesisExploratoria: publicado en web
    HipotesisExploratoria --> [*]: requiere validacion externa (fuera de alcance)
    Reproducible --> Descartado: sin evidencia / no especifico
    Descartado --> [*]
```

---

## 4. Cómo funciona (narrativa de ejecución)

### 4.1 Secuencia de una corrida completa

```mermaid
sequenceDiagram
    autonumber
    actor U as Usuario
    participant CLI as Orquestador CLI
    participant M as run_manifest
    participant DL as Descarga (DB snapshots)
    participant PIPE as Stages A/B/C
    participant FS as Almacenamiento
    participant WEB as Publicador web

    U->>CLI: oncogo-crc run --stage all
    CLI->>M: init run_id + timestamp
    CLI->>DL: descargar datos + fijar versiones
    DL-->>M: registrar versiones + hashes
    CLI->>PIPE: ejecutar etapas en orden
    PIPE->>FS: escribir tablas/figuras/JSON
    PIPE-->>M: registrar entradas/salidas por etapa
    CLI->>WEB: publish_results (solo artefactos del run actual)
    WEB->>FS: leer run_manifest + resultados
    WEB-->>U: sitio web + tabla de hipotesis
    Note over M,WEB: La web solo muestra lo generado en la corrida (gating anti-stale)
```

### 4.2 Modelo de reproducibilidad y trazabilidad

```mermaid
flowchart LR
    subgraph run["run_manifest.json"]
      rid["run_id + timestamp"]
      stg["estado por etapa (ok/omitida)"]
      io["archivos in/out + SHA-256"]
      dbv["version + fecha de cada base de datos"]
      sw["versiones de software"]
    end
    gate{"publish_results<br/>gating"}
    ok["Renderiza artefacto"]
    skip["Marca 'omitido' / 'stale'"]

    run --> gate
    gate -->|generado en este run| ok
    gate -->|no generado / previo| skip
```

---

## 5. Arquitectura de la capa web

Objetivo: que **cualquier persona** (jurado, oncólogo, otros investigadores) pueda **consultar los procesos,
el análisis y los resultados** sin instalar nada. Diseño en **dos niveles** (estático + interactivo), ambos con
opciones de hospedaje **gratuito**.

### 5.1 Vista de la capa web

```mermaid
flowchart TB
    subgraph src["Fuente de verdad (artefactos del pipeline)"]
      man["run_manifest.json"]
      tsv["tablas .tsv/.csv"]
      png["figuras .png"]
      html["figuras interactivas .html (Plotly)"]
      gml["red .graphml / .json"]
    end

    subgraph n1["Nivel 1 - Sitio estatico (MkDocs Material)"]
      pages["Paginas por etapa + metodologia + glosario"]
      embed["Figuras + tablas embebidas"]
      manview["Visor de run_manifest"]
    end

    subgraph n2["Nivel 2 - Dashboard interactivo (opcional)"]
      filt["Explorador de tabla de hipotesis (filtros/orden)"]
      volc["Volcano / scatter interactivo (Plotly)"]
      netv["Explorador de red (Cytoscape.js / pyvis)"]
      det["Ficha por candidato (evidencia + procedencia)"]
    end

    api["API JSON opcional (FastAPI)"]

    src --> n1
    src --> api --> n2
    src --> n2
```

### 5.2 Mapa del sitio (sitemap)

```mermaid
flowchart TB
    home["Inicio / resumen de la corrida"]
    met["Metodologia (3 capas)"]
    data["Datos y cohortes + versiones DB"]
    qc["Control de calidad (PCA, missingness)"]
    dif["Resultados diferenciales"]
    rep["Reproducibilidad multi-cohorte"]
    intg["Capa A: integracion + metricas"]
    redp["Capa B: modulos / red"]
    prio["Capa C: tabla de hipotesis priorizadas"]
    glos["Glosario"]
    lim["Limitaciones y prudencia traslacional"]
    repro["Reproducibilidad (manifest, hashes, FAIR)"]

    home --> met --> data --> qc --> dif --> rep --> intg --> redp --> prio
    prio --> lim
    home --> glos
    home --> repro
```

### 5.3 Interacción del usuario con el dashboard

```mermaid
sequenceDiagram
    autonumber
    actor V as Visitante web
    participant UI as Dashboard (navegador)
    participant API as API JSON / archivos estaticos
    participant FS as Artefactos del run

    V->>UI: abre pagina de priorizacion
    UI->>API: solicita hypothesis_priority + evidencia
    API->>FS: lee tablas/JSON del run actual
    FS-->>API: datos + procedencia (version DB)
    API-->>UI: JSON
    V->>UI: filtra (p.ej. esencial + drogable + tumor-especifico)
    UI-->>V: tabla + volcano + subred interactiva
    V->>UI: clic en un candidato
    UI-->>V: ficha con banderas de evidencia y disclaimer de hipotesis
```

### 5.4 Stack tecnológico propuesto

| Componente | Tecnología | Hospedaje gratuito | Estado |
|---|---|---|---|
| Sitio estático informativo | **MkDocs Material** (ya usado) | GitHub Pages | Existe (ampliar) |
| Figuras estáticas | Matplotlib (ya usado) | — | Existe |
| Figuras interactivas | **Plotly** (HTML embebible) | GitHub Pages | Nuevo |
| Explorador de red | **Cytoscape.js** o **pyvis** | GitHub Pages | Nuevo |
| Dashboard interactivo | **Streamlit** o **Dash** | Streamlit Community Cloud / HF Spaces | Nuevo (opcional) |
| API de datos | **FastAPI** (sirve JSON del run) | Render/Fly.io free / HF Spaces | Nuevo (opcional) |
| Contenedor reproducible | Docker | — | Nuevo (opcional) |

> **Recomendación pragmática:** empezar por **Nivel 1 (estático)** — cubre "ver procesos, análisis y
> resultados" para cualquiera, con costo cero y máxima reproducibilidad — y añadir el **Nivel 2 (dashboard)**
> solo para exploración interactiva de la tabla de hipótesis y las redes.

### 5.5 Contrato de datos web (fuente única de verdad)

La web **no recalcula**: consume artefactos del run. Contrato mínimo:

```mermaid
erDiagram
    RUN_MANIFEST ||--o{ STAGE : contiene
    STAGE ||--o{ ARTIFACT : produce
    ARTIFACT ||--o{ FIGURE : incluye
    RUN_MANIFEST {
      string run_id
      string timestamp
      json db_versions
      json software_versions
    }
    STAGE {
      string name
      string status
      json inputs
      json outputs
    }
    ARTIFACT {
      string path
      string sha256
      string type
    }
    HYPOTHESIS ||--o{ EVIDENCE : anota
    HYPOTHESIS {
      string gene_symbol
      float log2fc
      float fdr
      float reproducibility
      float convergence_score
      int rank
    }
    EVIDENCE {
      string source
      string source_version
      string flag
      float value
    }
```

---

## 6. Vista de despliegue

```mermaid
flowchart TB
    subgraph local["Ejecucion (local o servidor)"]
      py["Python + entorno (pyproject/requirements)"]
      pipe["Pipeline OncoGO-CRC 2.0"]
      ws["data_workspace (volumen externo)"]
    end
    subgraph ci["CI / Publicacion"]
      gh["GitHub Actions (build docs)"]
    end
    subgraph host["Hospedaje"]
      pages["GitHub Pages (sitio estatico)"]
      space["Streamlit Cloud / HF Spaces (dashboard)"]
    end

    py --> pipe --> ws
    ws -->|artefactos + manifest| gh --> pages
    ws -->|artefactos| space
    pages --> users["Consulta publica"]
    space --> users
```

---

## 7. Qué existe vs qué es nuevo (resumen)

| Área | Existe (1.0) | Nuevo (2.0) |
|---|---|---|
| Ingesta + QC + armonización + diferencial + overlap | ✅ | — |
| GO/GSEA + ML interpretable + TPS | ✅ | — |
| run_manifest + gating + FAIR | ✅ | + versión de DB por artefacto |
| Sitio MkDocs con figuras/glosario | ✅ | + figuras interactivas Plotly |
| Capa A (batch + MOFA+ + métricas) | — | ✅ |
| Capa B (GNN/WGCNA + módulos + enrichment) | — | ✅ |
| Capa C (evidencia convergente + tabla hipótesis) | — | ✅ |
| Dashboard interactivo + API | — | ✅ (opcional) |

---

## 8. Secuencia de implementación sugerida (mínimo riesgo)

```mermaid
flowchart LR
    s1["1. Capa C<br/>(mayor impacto, bajo computo)"] --> s2["2. Nivel 1 web<br/>tabla + figuras estaticas"]
    s2 --> s3["3. Capa A<br/>batch + MOFA+ + metricas"]
    s3 --> s4["4. Capa B<br/>WGCNA -> GNN"]
    s4 --> s5["5. Nivel 2 web<br/>dashboard interactivo"]
```

Racional: la **Capa C** entrega el artefacto de mayor valor (tabla de hipótesis con evidencia) con bajo costo
técnico; el **Nivel 1 web** lo hace consultable de inmediato; luego se profundiza integración (A) y red (B), y
por último la interactividad (Nivel 2).

---

## 9. Referencias de diseño

Ver `docs/DISENO_TECNICO_ONCOGO_2.0.md` (métodos, métricas y referencias validadas por capa) y
`docs/PROMPT_SCISPACE_BIBLIOGRAFIA_V2.md` (marco bibliográfico). El soporte metodológico de cada capa está
respaldado por la revisión de literatura 2024–2026 (MOFA+/VAE, corrección batch a nivel proteína, GNN sobre
PPI, marcos de priorización con DepMap/Open Targets, y estándares FAIR/DOME/TRIPOD+AI).
```
