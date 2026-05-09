# Plan de ejecución por fases — OncoGO-CRC y tesis

Este documento organiza el trabajo **después de reconstruir el proyecto**: no como una sola corrida completa, sino como **fases con entregables triples** en cada bloque relevante:

| Entregable | Qué es |
|------------|--------|
| **Técnico** | Archivos, tablas, figuras reproducibles en disco. |
| **Científico** | Respuesta parcial a la pregunta de investigación o validación de supuestos. |
| **Tesis** | Texto, figuras o tablas que pueden ir directo a metodología, resultados o anexos. |

**Comandos del pipeline:** en este repositorio las etapas se ejecutan con la CLI `oncogo-crc` (no `python -m src.01_*`, que no es válido en Python). Equivalencia: `download` ≈ 01, `harmonize` ≈ 02, … Ver [DOCUMENTACION_TESIS.md](DOCUMENTACION_TESIS.md#equivalencia-con-scripts-01_-).

```bash
# Una sola etapa
oncogo-crc run --stage <nombre>

# Orden completo (todas las etapas en secuencia)
oncogo-crc run --stage all
```

Nombres de etapa: `download`, `harmonize`, `qc`, `differential`, `go`, `network`, `ml`, `clinical`, `report`.

---

## Fase 0 — Verificación inicial del proyecto

**Objetivo:** Confirmar que la reconstrucción es utilizable antes de invertir tiempo en datos reales.

### Qué revisar

- Entorno activo (`pip install -e .` o Conda según `environment.yml`).
- Carpetas: `config/`, `data/raw|interim|processed/`, `metadata/`, `src/oncogo_crc/`, `results/`, `reports/`.
- Configuración presente y coherente entre sí.
- Metadatos plantilla en `metadata/sample_metadata/` y mapa `metadata/protein_id_maps/protein_mapping_master.tsv`.
- Que las **claves** en `config/datasets.yaml` (`cptac_colorectal`, `pride_crc_meta`, etc.) coincidan con lo que citas en la tesis y con las rutas que usarás.

### Checklist mínimo

- [ ] `config/config.yaml`
- [ ] `config/datasets.yaml`
- [ ] `config/comparisons.yaml`
- [ ] `metadata/sample_metadata/*`
- [ ] `metadata/protein_id_maps/protein_mapping_master.tsv`
- [ ] `metadata/gene_ontology/go_gene_long.tsv` (para ORA offline)
- [ ] `src/oncogo_crc/` y `src/oncogo_crc/utils/`

### Resultado técnico

Proyecto instalable; `oncogo-crc run --demo --stage all` opcionalmente ejecutable sin datos externos.

### Resultado científico

Ninguno todavía; solo **preparación de reproducibilidad**.

### Material para la tesis

Párrafo breve en metodología: “entorno reproducible” y versión del código (commit / release).

---

## Fase 1 — Definir el MVP de la tesis

**Objetivo:** Acotar el alcance para que sea defendible y terminable.

### Núcleo recomendado (alineado con `Informacion/`)

| Rol | Dataset (clave en `datasets.yaml`) | Función en la narrativa |
|-----|-----------------------------------|-------------------------|
| Backbone tumoral | `cptac_colorectal` | Firma primaria tumor vs normal |
| Validación externa | `pride_crc_meta` | Replicación / consistencia en tejido independiente |
| Capa clínica | `plasma_crc_2024` | Biomarcadores circulantes / traslación |

### Resultado técnico

Lista escrita de datasets activos y cuáles quedan explícitamente **fuera del MVP** (`crcdb_context`, `pharmacoproteomics`, metastásico, etc.).

### Resultado científico

Pregunta e hipótesis acotadas a lo que esos tres bloques pueden sustentar.

### Material para la tesis

Subsección “Alcance y MVP” en introducción o marco teórico; tabla “datasets utilizados”.

---

## Fase 2 — Preparar los datos reales

**Objetivo:** Matrices y metadatos alineados **antes** de presionar el pipeline.

### Paso 2.1 — Descargar y colocar matrices

Los nombres de archivo pueden ser los que prefieras; deben coincidir con **`matrix_path`** en `config/datasets.yaml`. Por defecto el repo espera por ejemplo:

| Dataset | Ruta por defecto en config |
|---------|---------------------------|
| CPTAC | `data/raw/cptac/cptac_crc_matrix.tsv` |
| PRIDE | `data/raw/pride/pride_crc_matrix.tsv` |
| Plasma | `data/raw/plasma/plasma_crc_matrix.tsv` |

Si usas `protein_abundance.tsv`, `pride_merged.tsv`, etc., **actualiza** `datasets.yaml` para que apunte a esos archivos.

### Paso 2.2 — Metadatos por cohorte

Completar al menos:

- `metadata/sample_metadata/cptac_metadata.tsv`
- `metadata/sample_metadata/pride_merged_metadata.tsv`
- `metadata/sample_metadata/plasma_metadata.tsv`

Regla: cada `sample_id` de la matriz debe existir en el TSV correspondiente con las columnas que usarás en `comparisons.yaml` (p. ej. `condition`).

### Paso 2.3 — Identificadores y mapa maestro

Documentar por cohorte: UniProt, símbolo de gen, ID de plataforma; construir o ampliar `protein_mapping_master.tsv` y enriquecer `go_gene_long.tsv` según el universo génico real.

### Resultados de la fase

| Tipo | Entrega |
|------|---------|
| **Técnico** | TSV presentes en `data/raw/…`, metadatos sin muestras huérfanas, mapa ID consistente. |
| **Científico** | Definición explícita de cohorte computacional y criterios de inclusión/exclusión. |
| **Tesis** | Tabla de muestras, fuentes (PDC/PRIDE/publicación), variables clínicas disponibles. |

---

## Fase 3 — Ajustar configuraciones según la tesis

**Objetivo:** Que YAML refleje la pregunta **protein-centric** y el orden CPTAC → PRIDE → plasma.

### Paso 3.1 — `config/config.yaml`

Revisar al menos:

- `project.primary_dataset`: normalmente `cptac_colorectal`
- `preprocessing.*`: missingness, log2, KNN, normalización
- `preprocessing.batch_correction`: hoy `none`; si incorporas COMBAT u otro método, documentarlo como extensión
- `go.semantic_reduction`: `true` para coherencia con el reajuste de la tesis
- `translation.evidence_tiers`: niveles de evidencia traslacional
- `differential.*`, `network.*`, `ml.*`

### Paso 3.2 — `config/datasets.yaml`

- **MVP activo:** rutas y metadatos para `cptac_colorectal`, `pride_crc_meta`, `plasma_crc_2024`
- **Opcional:** comentar o no ejecutar etapas sobre cohortes sin archivo (el harmonize/QC solo procesa datasets con matriz presente)

### Paso 3.3 — `config/comparisons.yaml`

**Mínimo recomendado:**

- `cptac_tumor_vs_normal`
- `plasma_crc_vs_healthy`

**Validación PRIDE** (tumor vs normal en tejido): añade una entrada nueva cuando la metadata lo permita, por ejemplo:

```yaml
pride_validation_tumor_vs_normal:
  dataset: pride_crc_meta
  label: "PRIDE — tumor vs normal (validación externa)"
  condition_column: condition
  group_a: tumor
  group_b: normal
```

Ajusta `group_a` / `group_b` a los valores **reales** de tu TSV. Para comparaciones avanzadas (`stage_bin`, `msi_status`, respondedores), usa `optional: true` hasta tener estratos completos.

### Resultados de la fase

| Tipo | Entrega |
|------|---------|
| **Técnico** | YAML versionados y comentados en Git con mensaje claro. |
| **Científico** | Diseño analítico explícito (qué contrastes responden qué parte de la hipótesis). |
| **Tesis** | Sección “Diseño del estudio computacional” / diagrama de comparaciones. |

---

## Fase 4 — Pipeline núcleo (integración y estadística)

Ejecutar **por etapas**, deteniéndote para QC entre bloques.

| Paso | Comando | Salidas principales (típicas) |
|------|---------|-------------------------------|
| Registro | `oncogo-crc run --stage download` | `data/raw/manifest.tsv` |
| Armonización | `oncogo-crc run --stage harmonize` | `data/interim/<dataset>_harmonized.tsv` |
| QC / normalización | `oncogo-crc run --stage qc` | `data/processed/<dataset>_normalized.tsv` |
| Diferencial | `oncogo-crc run --stage differential` | `results/differential/*_diff.tsv`, `results/figures/*_volcano.png` |

### Por cada sub-fase

| Sub-fase | Resultado técnico | Resultado científico | Material tesis |
|----------|-------------------|----------------------|----------------|
| download | Manifest y trazabilidad | Inventario de fuentes | Tabla de datos / PRISMA-like si aplica |
| harmonize | Matrices unificadas por gen/proteína | Integrabilidad multi-estudio | Párrafo de estandarización de IDs |
| qc | Matrices limpias | Datos aptos para inferencia | Preprocessing detallado + figuras QC si las añades |
| differential | Lista DEP + volcanos | Alteración exprés proteica por contraste | Resultados cuantitativos principales, figuras |

---

## Fase 5 — Gene Ontology como eje central

**Comando:** `oncogo-crc run --stage go`

### Paso 5.1 — Enriquecimiento

- **Técnico:** `results/go/<comparison>_go.tsv`
- **Científico:** Procesos BP/MF/CC asociados al fenotipo contrastado
- **Tesis:** Figuras/tablas de enriquecimiento por comparación

### Paso 5.2 — Reducción semántica

El pipeline escribe `results/go/<comparison>_go_reduced.tsv` cuando `go.semantic_reduction` está activo (proxy por solapamiento de genes). Puedes **complementar** con REVIGO u otra herramienta externa y guardar tablas en `results/go/` o `reports/` citadas en la tesis.

| Resultado técnico | Resultado científico | Material tesis |
|-------------------|----------------------|----------------|
| Lista GO reducida | Ejes funcionales no redundantes | Discusión “procesos representativos” vs listas largas |

---

## Fase 6 — Redes y módulos funcionales

**Comando:** `oncogo-crc run --stage network`

- **Técnico:** `results/networks/<primary>_hub_proteins.tsv`
- **Científico:** Hubs y arquitectura de coexpresión (interpretar como asociación estadística, no PPI curada salvo que integres otra fuente)
- **Tesis:** Figura de red o tabla de hubs ligada a la discusión de vulnerabilidades/módulos

---

## Fase 7 — IA explicable

**Comando:** `oncogo-crc run --stage ml`

- **Técnico:** `results/ml/<primary>_feature_ranking.tsv`, `*_model_metrics.tsv`
- **Científico:** Priorización de características con modelo interpretable + métricas out-of-sample
- **Tesis:** Subsección ML; figura de importancias o ranking combinado

### Explicabilidad avanzada (opcional)

SHAP/LIME no están en el núcleo del MVP; pueden añadirse como extensión (`pip install shap`) y un notebook en `notebooks/`. En la tesis, declara qué implementaste y por qué.

---

## Fase 8 — Traducción clínica

**Comando:** `oncogo-crc run --stage clinical`

- **Técnico:** `results/clinical_translation/<primary>_translation_ranked.tsv`
- **Científico:** Clasificación por niveles de evidencia (`translation.evidence_tiers`) + síntesis traslacional
- **Tesis:** Tabla “validado / exploratorio / validación futura”; fortalece la discusión sin sobreinterpretar

---

## Fase 9 — Reporte técnico consolidado

**Comando:** `oncogo-crc run --stage report`

- **Técnico:** `reports/executive_summary/results_summary.md`
- **Científico:** Resumen transversal del pipeline
- **Tesis:** Base para anexos o material suplementario

---

## Fase 10 — Documentación y despliegue

```bash
pip install mkdocs mkdocs-material
mkdocs serve    # revisión local en http://127.0.0.1:8000
mkdocs build    # sitio estático en site/
```

| Tipo | Entrega |
|------|---------|
| **Técnico** | Sitio `site/` o hosting |
| **Científico** | Paquete reproducible para revisores |
| **Tesis** | Anexo “documentación del proyecto” o enlace a material suplementario |

---

## Fase 11 — Integración narrativa en la tesis

Sin comandos nuevos: redacción.

| Capítulo | Contenido típico basado en el pipeline |
|----------|----------------------------------------|
| Metodología | Fases 0–10: datos, QC, diferencial, GO+reducción, redes, ML, traslación |
| Resultados | DEP, GO reducido, hubs, ML, ranking clínico |
| Discusión | Literatura, limitaciones (batch, cohortes, ausencia validación experimental), validación futura |

---

## Fase 12 — Expansión post-MVP

Cuando el núcleo esté estable:

1. Dataset metastásico / iProx / CRLM  
2. CRCDB / contexto EOCRC  
3. Organoides / farmacoproteómica / drug repurposing como línea de validación futura o discusión

Cada expansión repite la lógica de **Fase 2 (datos) → Fase 3 (config) → subconjunto de etapas 4–9**.

---

## Orden resumido recomendado

1. Fase 0 — Verificación  
2. Fase 1 — MVP  
3. Fase 2 — Datos reales alineados  
4. Fase 3 — YAML acorde a la pregunta  
5. Fases 4 → 9 — Pipeline por etapas (`download` … `report`)  
6. Fase 10 — MkDocs  
7. Fase 11 — Escritura  
8. Fase 12 — Expansiones  

### Qué deberías tener al cerrar el MVP

| Artefacto | Ubicación típica |
|-----------|-------------------|
| Cohorte computacional armonizada | `data/interim/`, `data/processed/` |
| DEP por contraste | `results/differential/` |
| GO + reducción | `results/go/` |
| Hubs | `results/networks/` |
| Priorización ML | `results/ml/` |
| Ranking traslacional | `results/clinical_translation/` |
| Informe ejecutivo | `reports/executive_summary/results_summary.md` |

---

## Recomendación operativa

La mejor forma de avanzar **no** es “correr todo por correr”, sino:

1. Cerrar **Fases 0–3** con datos y configuración estables.  
2. Ejecutar **Fase 4** y revisar QC + diferencial antes de GO.  
3. Ejecutar **Fase 5** y validar interpretación GO antes de redes y ML.  
4. Completar **Fases 6–9** y solo entonces **Fases 10–11**.  
5. Reservar **Fase 12** para cuando la narrativa del MVP esté cerrada.
