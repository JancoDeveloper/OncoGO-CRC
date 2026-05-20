# Proyecto web — OncoGO-CRC

Sitio generado automáticamente por `oncogo-crc run --stage publish_results` a partir de los TSV bajo `results/`. Las páginas se reescriben en cada corrida; los hashes y `run_id` para trazabilidad están en `results/run_manifest.json`.

## Lectura recomendada (orden para tesis)

1. [**Síntesis para tesis**](thesis_synthesis.md) — qué demuestra el pipeline y qué redactar aparte (narrativa).
2. [**Interpretación biológica**](interpretacion_biologica.md) — módulos funcionales detectados por comparación (OXPHOS, mitorribosoma, proteasoma, ECM/CMS4, fase aguda…) con contexto CRC.
3. [**Figuras**](figuras.md) — volcanos etiquetados, ORA/GSEA top, scatter CPTAC↔PRIDE, PCA, heatmaps, AUC nested CV y transferencia.
4. [**Resultados base**](results_base.md) — tablas cuantitativas: proteínas, ORA, GSEA, firma CPTAC∩PRIDE, ML, red, traslación, export bundle.
5. [**Glosario**](glosario.md) — definiciones operativas: padj/FDR, NES, leading edge, nested CV, transferencia, calibración, TPS, etc.

## Referencia secundaria

- [Datasets](datasets.md) · [Cobertura clínica](clinical_coverage.md) · [Evidencia externa](external_evidence.md) · [Estabilidad](stability.md) · [Inventario](results_inventory.md) · [Limitaciones](limitations.md)
- Documento de limitaciones honestas y plan de validación: [`docs/RESULTADOS_Y_LIMITACIONES_TESIS.md`](../RESULTADOS_Y_LIMITACIONES_TESIS.md)
- Prompt para redacción interpretativa con LLM: [`docs/PROMPT_SINTESIS_TESIS.md`](../PROMPT_SINTESIS_TESIS.md)
