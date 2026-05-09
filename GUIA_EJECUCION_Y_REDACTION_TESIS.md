# Guía práctica de ejecución y de redacción (pipeline → tesis)

Este documento hace **dos funciones a la vez**:

1. **Checklist operativo** antes y después de correr el pipeline.  
2. **Mapa de redacción**: qué capítulo de la tesis alimenta cada etapa y qué debes extraer en texto/figuras.

En este repositorio las etapas no son `src.01_*.py` (no son importables en Python con ese nombre). Usa la CLI:

```bash
oncogo-crc run --stage <etapa>
# etapas: download_real | download | validate_inputs | harmonize | qc | differential | go | network | ml | clinical | report | stability | publish_results
```

Equivalencia con la nomenclatura “01–09” del paquete `Informacion/`: ver [DOCUMENTACION_TESIS.md](DOCUMENTACION_TESIS.md#equivalencia-con-scripts-01_-).

---

## Tabla rápida: etapa → módulo histórico

| Comando moderno | Referencia antigua (Informacion/) | Código fuente |
|-----------------|-----------------------------------|---------------|
| `download` | `01_download_data` | `oncogo_crc/stages/download_data.py` |
| `harmonize` | `02_harmonize_ids` | `oncogo_crc/stages/harmonize_ids.py` |
| `qc` | `03_qc_normalize` | `oncogo_crc/stages/qc_normalize.py` |
| `differential` | `04_differential_analysis` | `oncogo_crc/stages/differential_analysis.py` |
| `go` | `05_go_enrichment` | `oncogo_crc/stages/go_enrichment.py` |
| `network` | `06_network_modules` | `oncogo_crc/stages/network_modules.py` |
| `ml` | `07_ml_prioritization` | `oncogo_crc/stages/ml_prioritization.py` |
| `clinical` | `08_clinical_translation` | `oncogo_crc/stages/clinical_translation.py` |
| `report` | `09_generate_report` | `oncogo_crc/stages/generate_report.py` |

---

## 1) Checklist operativo (ultra corto)

### Fase 0 — Antes de correr

**Estructura y entorno**

- [ ] La carpeta del proyecto está reconstruida correctamente.
- [ ] Existen `config/`, `data/`, `metadata/`, `src/`, `results/`, `reports/`.
- [ ] El entorno está creado (Conda con `environment.yml` **o** `python3 -m venv .venv` + `pip install -e .`).
- [ ] El entorno está **activado**.
- [ ] Python y dependencias del proyecto funcionan: al menos `pandas`, `scikit-learn`, `networkx`, `scipy`, `PyYAML` (`pip install -e .` las instala). **GSEApy** no es obligatorio en el núcleo: el enriquecimiento GO por defecto es **ORA offline** con `metadata/gene_ontology/go_gene_long.tsv`; si en tu tesis usas GSEA/Enrichr remoto, instala `gseapy` como extra opcional y documenta el cambio.
- [ ] Existen `config/config.yaml`, `config/datasets.yaml`, `config/comparisons.yaml`.
- [ ] El enfoque del YAML está alineado con el **reajuste** de la tesis (GO central, reducción semántica, niveles de evidencia en `translation.evidence_tiers`).

### Fase 1 — Datasets

- [ ] CPTAC / PDC está descargado y en `data/raw/cptac/` (o la ruta que indique `datasets.yaml` → `cptac_colorectal.matrix_path`).
- [ ] PRIDE tejido está en `data/raw/pride/`.
- [ ] Plasma 2024 está en `data/raw/plasma/`.
- [ ] Los **nombres de archivo** coinciden con `config/datasets.yaml` (si usas otros nombres, actualiza el YAML).
- [ ] Datasets opcionales (CRCDB, farmacoproteómica, metástasis) están marcados como **fase avanzada** y no bloquean el MVP si aún no hay matrices.

### Fase 2 — Metadatos

- [ ] `metadata/sample_metadata/cptac_metadata.tsv` existe y está completo para el MVP.
- [ ] `metadata/sample_metadata/pride_merged_metadata.tsv` existe y está completo (cuando uses PRIDE).
- [ ] `metadata/sample_metadata/plasma_metadata.tsv` existe y está completo (cuando uses plasma).
- [ ] `metadata/protein_id_maps/protein_mapping_master.tsv` existe.
- [ ] Los `sample_id` coinciden **exactamente** entre cada metadata y su matriz (columnas de la matriz).
- [ ] Los valores de la columna usada en comparaciones (`condition`, etc.) coinciden con `group_a` / `group_b` en `comparisons.yaml` (p. ej. `tumor`, `normal`, `crc`, `healthy`).

### Fase 3 — Ejecución núcleo

Marca cuando cada etapa **terminó sin error** y generó las salidas esperadas:

- [ ] `oncogo-crc run --stage download` → `data/raw/manifest.tsv` (y datos en `data/raw/…`).
- [ ] `oncogo-crc run --stage harmonize` → `data/interim/*_harmonized.tsv`.
- [ ] `oncogo-crc run --stage qc` → `data/processed/*_normalized.tsv`.
- [ ] `oncogo-crc run --stage differential` → `results/differential/*_diff.tsv` y `results/figures/*_volcano.png`.
- [ ] `oncogo-crc run --stage go` → `results/go/*_go.tsv` (y `*_go_reduced.tsv` si la reducción semántica está activa).

### Fase 4 — Reajuste funcional

- [ ] Se aplicó **reducción semántica** a GO (automática: `*_go_reduced.tsv`; o complementaria con REVIGO/glue externo, citado en metodología).
- [ ] Se identificaron procesos GO prioritarios y **no redundantes** para narrativa.
- [ ] `oncogo-crc run --stage network` → `results/networks/*_hub_proteins.tsv` (**red de coexpresión proteína–proteína** por correlación en esta implementación).
- [ ] `oncogo-crc run --stage ml` → `results/ml/*_model_metrics.tsv` (y ranking de features).
- [ ] `oncogo-crc run --stage clinical` → `results/clinical_translation/*_translation_ranked.tsv`.
- [ ] Los hallazgos quedaron clasificados según `translation.evidence_tiers` en config (por defecto: `validated`, `hypothesis_generating`, `future_validation`), equivalentes narrativos a: **validados**, **exploratorios / generadores de hipótesis**, **validación futura**.

### Fase 5 — Reportes y tesis

- [ ] `oncogo-crc run --stage report` → `reports/executive_summary/results_summary.md`.
- [ ] La documentación MkDocs puede levantarse: `mkdocs serve` / `mkdocs build`.
- [ ] Sabes **qué tablas/figuras** van a Metodología, Resultados y Discusión (usa la sección 2 más abajo).
- [ ] Existe borrador de narrativa: metodología, resultados, discusión, limitaciones, validación futura.

---

## 2) Ruta exacta: capítulo de tesis ← cada etapa del pipeline

### Mapa general

**Etapa → resultado técnico → capítulo que alimenta**

---

### `download` (antes `01_download_data.py`)

**Qué hace**

- Registra datasets cuya matriz existe en disco (modo real) o genera datos sintéticos (`--demo`).
- Deja trazabilidad vía manifiesto.

**Salidas típicas**

- `data/raw/manifest.tsv`

**Capítulo de tesis**

- **Metodología**: fuentes de datos, organización del proyecto, trazabilidad, criterios de selección inicial.
- **Marco / Estado del arte**: qué repositorios públicos se integraron realmente.

**Cómo redactarlo**

Párrafos sobre origen de datos (PDC/CPTAC, PRIDE, plasma), estructura del repositorio y reproducibilidad del flujo.

---

### `harmonize` (antes `02_harmonize_ids.py`)

**Qué hace**

- Armoniza identificadores con `protein_mapping_master.tsv`, colapsa duplicados por símbolo de gen (mediana entre muestras).

**Salidas típicas**

- `data/interim/<dataset>_harmonized.tsv`

**Capítulo de tesis**

- **Metodología**: armonización de IDs, integración multi-cohorte, preprocesamiento previo al QC.
- **Limitaciones**: discrepancias entre plataformas, pérdidas por mapeo incompleto, genes sin anotación GO.

**Cómo redactarlo**

Cómo pasaste de IDs de plataforma a símbolos comparables y qué fracción del universo proteico quedó mapeada.

---

### `qc` (antes `03_qc_normalize.py`)

**Qué hace**

- Filtra por missingness (umbrales en `config.yaml`), log2, imputación KNN, normalización z-score por proteína.

**Salidas típicas**

- `data/processed/<dataset>_normalized.tsv`

**Capítulo de tesis**

- **Metodología**: QC, imputación, normalización, inclusión/exclusión de proteínas y muestras.
- **Resultados** (opcional): número final de proteínas/muestras, estabilidad post-QC (si añades figuras propias).
- **Limitaciones**: sesgos de imputación, batch effects (si no corriges batch), variación técnica.

**Cómo redactarlo**

Justificación de umbrales y decisión de método de imputación/normalización acorde a `config.yaml`.

---

### `differential` (antes `04_differential_analysis.py`)

**Qué hace**

- Contrasta grupos definidos en `comparisons.yaml` (t-test + FDR por proteína/gén).
- Genera volcanos.

**Salidas típicas**

- `results/differential/<comparison>_diff.tsv`
- `results/figures/<comparison>_volcano.png`

**Capítulo de tesis**

- **Resultados**: núcleo cuantitativo — DEP por comparación (tumor/normal, plasma CRC/sano, etc.).
- **Discusión**: concordancia con literatura, proteínas novedosas vs conocidas, enlaces clinicorracionalistas tentativos.

**Cómo redactarlo**

Subsección por comparación: tabla top up/down + figura volcan; texto citando `padj` y magnitud (`log2fc`).

---

### `go` (antes `05_go_enrichment.py`)

**Qué hace**

- Enriquecimiento tipo ORA (Fisher) sobre genes/proteínas significativas vs universo de la cohorte, usando `go_gene_long.tsv`.
- Reducción semántica interna → `*_go_reduced.tsv` (agrupa términos redundantes por solapamiento de genes).

**Salidas típicas**

- `results/go/<comparison>_go.tsv`
- `results/go/<comparison>_go_reduced.tsv` (si está activo)

**Capítulo de tesis**

- **Resultados**: procesos BP/MF/CC enriquecidos; tras reducción, **ejes funcionales** más claros.
- **Discusión**: interpretación biológica (metabolismo, EMT, inmunidad, etc.) apoyada en términos representativos.
- **Metodología**: definición de universo, umbral FDR, ORA vs otros métodos si amplías con GSEA/GSEApy.

**Reducción semántica (bloque obligatorio en la narrativa de la tesis)**

Aunque también puedas usar REVIGO u otra herramienta externa, en texto debes explicar **cómo** pasaste de una lista larga a procesos representativos y por qué eso refuerza el enfoque protein-centric.

---

### `network` (antes `06_network_modules.py`)

**Qué hace (implementación actual)**

- Red **proteína–proteína** por correlación de expresión entre muestras; identificación de **hubs** por grado.

**Salidas típicas**

- `results/networks/<primary>_hub_proteins.tsv`

**Capítulo de tesis**

- **Resultados**: hubs y lectura de “arquitectura” de coexpresión (no confundir con PPI curado salvo que integres otra fuente).
- **Discusión**: centralidad de ciertas proteínas; posibles módulos ligados a procesos GO si enlazas genes hubs con tu tabla GO.
- **Aplicabilidad / traslación cautelosa**: candidatos más defendibles que una sola DEP aislada.

**Cómo redactarlo**

Dejar claro que la red es de **coexpresión estadística** en tu cohorte computacional; la integración con procesos GO es narrativa y por cruces manuales o tablas suplementarias.

---

### `ml` (antes `07_ml_prioritization.py`)

**Qué hace**

- Modelos interpretables (en el núcleo: clasificación dispersa + Random Forest), métricas held-out, ranking de proteínas.

**Salidas típicas**

- `results/ml/<primary>_feature_ranking.tsv`
- `results/ml/<primary>_model_metrics.tsv`

**Capítulo de tesis**

- **Metodología**: modelo(s), validación, métricas, riesgo de sobreajuste.
- **Resultados**: AUC/accuracy según corresponda, ranking de importancia.
- **Discusión**: coherencia entre biomarcadores ML y DEP/GO.
- **Limitaciones**: tamaño muestral, caracter exploratorio del ML.

**Explicabilidad (SHAP/LIME)**

Si las añades como extensión: **Resultados** (figuras de importancia) y **Discusión** (conexión modelo–proteína–GO). Declara explícitamente qué herramienta usaste.

---

### `clinical` (antes `08_clinical_translation.py`)

**Qué hace**

- Integra evidencias (diferencial, hubs, ML) en un ranking traslacional con **niveles de evidencia** configurables.

**Salidas típicas**

- `results/clinical_translation/<primary>_translation_ranked.tsv`

**Capítulo de tesis**

- **Discusión traslacional**: biomarcadores candidatos, plasma vs tejido si aplica, hipótesis priorizadas.
- **Aplicabilidad**: qué validar en cohorte u organoides.
- **Conclusiones**: síntesis de hallazgos priorizados.

---

### `report` (antes `09_generate_report.py`)

**Qué hace**

- Consolida extractos en `results_summary.md`.

**Salidas típicas**

- `reports/executive_summary/results_summary.md`

**Capítulo de tesis**

- **Anexos / material suplementario**, **resumen ejecutivo**, **presentación de defensa**.

---

## 3) Resumen ultra claro: ¿qué alimenta qué?

| Bloque de la tesis | Se nutre principalmente de |
|--------------------|------------------------------|
| **Metodología** | `download`, `harmonize`, `qc`, configuración YAML, estrategia GO, estrategia ML |
| **Resultados** | `differential`, `go` (+ reducción), `network`, `ml`, `clinical` |
| **Discusión** | DEP, GO reducido, hubs, ranking ML, traslación; confrontación con literatura |
| **Limitaciones** | armonización, missingness, batch, muestra, sobreajuste ML, ausencia de validación experimental |
| **Conclusiones** | priorización de `clinical`, síntesis GO + módulos coexpresivos + ML |

---

## 4) Regla práctica para escribir en paralelo

**Cada etapa del pipeline debe producir dos cosas:**

1. **Salida técnica** (tabla/figura en `data/` o `results/`).  
2. **Salida de redacción**: al menos un párrafo o una subsección que cite explícitamente esa salida.

**Ejemplo**

| Etapa | Salida técnica | Salida de tesis |
|-------|----------------|-----------------|
| `differential` | `cptac_tumor_vs_normal_diff.tsv`, volcanos | Subsección “Proteínas diferencialmente expresadas (tumor vs normal)”; figura volcán; tabla top DEP |

---

## 5) Flujo ideal recomendado

1. Usar el **checklist** antes de cada fase grande.
2. Después de cada etapa, preguntarte: **¿qué capítulo estoy alimentando?** y escribir ese bloque **en paralelo**, no solo al final.
3. Mantener una carpeta de borrador (Word/LaTeX/Notion) con la misma numeración que las etapas para pegar tablas exportadas desde los TSV.

Para el orden por **fases del proyecto** (MVP, datos, configuración), ver también [PLAN_EJECUCION_TESIS.md](PLAN_EJECUCION_TESIS.md).

---

## 6) Plantillas “ficha por etapa” (copiar y pegar)

Duplica el bloque **Ficha universal** cada vez que ejecutes una etapa (`oncogo-crc run --stage <nombre>`). Conserva un archivo por sesión o una sección en tu borrador de tesis.

### Ficha universal (una corrida)

```markdown
## Ficha de ejecución — etapa: _______________

**Fecha:** ____ / ____ / ______  
**Entorno:** Conda `oncogo-crc` / venv `.venv` / otro: _______________  
**Commit o etiqueta del repo:** _______________  
**Comando:** `oncogo-crc run --stage _______________`  
**¿Demo (`--demo`)?** Sí / No  

### Entrada / configuración relevante

- Archivos o YAML tocados hoy: _______________
- Dataset(s) involucrados: _______________

### Salida técnica (qué quedó en disco)

- Rutas principales generadas o actualizadas:
  - _______________
  - _______________
- ¿Sin errores en consola? Sí / No — notas: _______________

### Salida científica (qué permitió decidir)

- Respuesta parcial o validación de supuesto: _______________

### Redacción para la tesis (rellenar el mismo día)

- **Capítulo(es) alimentado(s):** Metodología / Resultados / Discusión / Otro: _______________
- **Subsección o párrafo propuesto (título):** _______________
- **Figura(s) / tabla(s) previstas:** _______________
- **Frase guía (1–2 líneas) para no olvidar el mensaje:** _______________

### Limitaciones observadas en esta etapa

- _______________

### Siguiente paso inmediato

- Etapa o tarea: _______________
```

### Tabla compacta — ciclo completo del MVP

Copia esta tabla y complétala una fila cada vez que termines una etapa del núcleo.

| Etapa | Fecha | Salidas clave (rutas) | Figura/tablas para tesis | Capítulo |
|-------|-------|------------------------|----------------------------|----------|
| `download` | | | | |
| `harmonize` | | | | |
| `qc` | | | | |
| `differential` | | | | |
| `go` | | | | |
| `network` | | | | |
| `ml` | | | | |
| `clinical` | | | | |
| `report` | | | | |

### Plantilla corta — solo recordatorio post-etapa

Si prefieres una línea en el cuaderno de laboratorio:

`[fecha] etapa X → archivos: … → borrador tesis: capítulo Y, subsección Z → pendiente: …`
