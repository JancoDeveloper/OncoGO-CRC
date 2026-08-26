# OncoGO-CRC — App web Nivel 2 (dashboard interactivo)

Dashboard **Streamlit** para explorar los resultados del pipeline de forma interactiva:

- **🧾 Hipótesis** — tabla priorizada (Capa C) con filtros, descarga CSV y perfil de evidencia por gen.
- **🌋 Volcano** — volcano diferencial interactivo (Plotly): umbrales ajustables y resaltado de genes.
- **🕸️ Red** — explorador de red por correlación de perfiles proteicos (Top-N + umbral ajustable).
- **🔎 Procedencia** — versiones de fuentes, pesos del score y reproducibilidad.

Si no encuentra artefactos reales, arranca en **modo demostración** con datos sintéticos claramente
etiquetados, para que la app sea explorable sin ejecutar el pipeline.

## Ejecutar en local

```bash
pip install -r webapp/requirements.txt
streamlit run webapp/streamlit_app.py
```

Fuente de datos (prioridad): variable `ONCOGO_CRC_RESULTS` → `paths.data_workspace/results` → `./results`.
También puedes cambiar el directorio desde la barra lateral.

```bash
export ONCOGO_CRC_RESULTS=/Volumes/JANCODATA/OncoGO-CRC/results
streamlit run webapp/streamlit_app.py
```

## Artefactos que consume

| Vista | Archivo |
|---|---|
| Hipótesis / Resumen | `results/prioritization/hypothesis_priority.tsv` |
| Perfil de evidencia | `results/prioritization/convergence_evidence_matrix.tsv` |
| Procedencia | `results/prioritization/prioritization_provenance.json` |
| Volcano | `results/differential/*_diff.tsv` |
| Red | `data/processed/<dataset>_zscore.tsv` |

Genera estos artefactos con:

```bash
oncogo-crc run --stage all         # pipeline completo
# o solo la Capa C (requiere differential ya ejecutado):
oncogo-crc run --stage prioritization
```

## Despliegue

- **Streamlit Community Cloud / Hugging Face Spaces:** apunta a `webapp/streamlit_app.py` y sube un
  `results/` de ejemplo (o deja el modo demo). `requirements.txt` está en `webapp/`.
- La documentación estática (procesos y arquitectura) se publica aparte con MkDocs (`mkdocs build`).
