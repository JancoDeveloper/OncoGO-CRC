# App web Nivel 2 y Capa C (priorización por evidencia convergente)

Dos entregables complementan la documentación estática (MkDocs) para que **cualquier persona pueda
consultar procesos, análisis y resultados**:

1. **Capa C — priorización por evidencia convergente** (stage `prioritization` del pipeline).
2. **App web Nivel 2** (dashboard **Streamlit** interactivo) que consume las salidas de la Capa C.

> Encuadre no-promisorio: la salida son **hipótesis priorizadas con evidencia convergente**, no dianas
> validadas. La evidencia ausente **no penaliza** el score.

---

## 1. Capa C — stage `prioritization`

Anota la firma diferencial primaria (y la firma reproducible CPTAC∩PRIDE si existe) con evidencia pública
objetiva y calcula un **score de convergencia transparente** (media ponderada de señales normalizadas
disponibles).

### Señales del score

| Señal | Fuente | Favorabilidad |
|---|---|---|
| `diff_effect` | diferencial | mayor \|log2FC\| |
| `diff_significance` | diferencial | menor padj |
| `reproducible` | CPTAC∩PRIDE | reproducible (+ mismo signo) |
| `network_hub` | red de coexpresión | mayor grado |
| `depmap_essentiality` | DepMap (Chronos) | más esencial (Chronos negativo) |
| `opentargets_tractability` | Open Targets | más tratable |
| `dgidb_druggability` | DGIdb | más fármacos que interactúan |
| `hpa_tumor_specificity` | Human Protein Atlas | más específico de tumor |
| `gtex_safety` | GTEx | menor expresión en tejido normal |
| `cbioportal_genomic` | cBioPortal/TCGA | mayor frecuencia de alteración |

Los proveedores externos son **snapshots locales pluggables** en `metadata/evidence/`
(ver el README de esa carpeta). Si falta un snapshot, esa fuente se marca no disponible.

### Salidas (`results/prioritization/`)

| Archivo | Contenido |
|---|---|
| `hypothesis_priority.tsv` | Tabla ancha priorizada (rank, gen, score, tier, evidencias, señales) |
| `convergence_evidence_matrix.tsv` | Formato largo (gen × señal): valor crudo, favorabilidad, peso |
| `prioritization_provenance.json` | Pesos, versión/ruta de cada fuente, parámetros y timestamp |

### Configuración (`config/config.yaml` → `prioritization`)

`enabled`, `primary_comparison`, `candidate_source` (`significant`/`all`/`topfc`), `max_candidates`,
`evidence_dir` y `weights` por señal.

### Ejecución

```bash
oncogo-crc run --stage prioritization   # requiere differential ejecutado
# o dentro del pipeline completo:
oncogo-crc run --stage all
```

---

## 2. App web Nivel 2 (Streamlit)

Dashboard interactivo (código en `webapp/`):

- **Resumen** — métricas globales y distribución por tier.
- **Hipótesis** — tabla filtrable + descarga CSV + perfil de evidencia por gen.
- **Volcano** — Plotly interactivo con umbrales ajustables y resaltado de genes.
- **Red** — explorador de coexpresión (Top-N + umbral de correlación).
- **Procedencia** — versiones de fuentes, pesos y reproducibilidad.

Si no encuentra artefactos reales, arranca en **modo demostración** con datos sintéticos etiquetados.

```bash
pip install -r webapp/requirements.txt
streamlit run webapp/streamlit_app.py
```

Detalles y despliegue (Streamlit Community Cloud / Hugging Face Spaces): ver `webapp/README.md`.

---

## 3. Relación con la arquitectura

Esta capa corresponde a la **Capa C** y a la **capa web interactiva** descritas en
`ARQUITECTURA_ONCOGO_2.0.md` y `DISENO_TECNICO_ONCOGO_2.0.md`. La documentación estática (MkDocs) y el
dashboard (Streamlit) son complementarios: la primera explica **procesos y fundamentos**; el segundo
permite **explorar resultados** de forma interactiva.
