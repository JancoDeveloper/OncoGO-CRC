# Síntesis para tesis

## Qué hace este pipeline (en una frase)

Integra matrices proteómicas públicas (CPTAC, validación PRIDE, plasma), normaliza e imputa, contrasta **tumor vs normal** (y otras comparaciones núcleo), enriquece genes significativos en **Gene Ontology (ORA offline)**, resume redes y prioriza proteínas con modelos interpretables.

## Qué está **demostrado** aquí (evidencia cuantitativa)

- **Tumor vs tejido adyacente normal (CPTAC)** (`cptac_tumor_vs_normal`): **5029** proteínas con FDR por debajo de **0.05** (de **6239** filas en el contraste).
  - GO (ORA): **0** términos enriquecidos (FDR por debajo de **0.05**).
- **Plasma CRC vs control sano** (`plasma_crc_vs_healthy`): **0** proteínas con FDR por debajo de **0.05** (de **739** filas en el contraste).
  - GO (ORA): **0** términos enriquecidos (FDR por debajo de **0.05**).
- **PRIDE validación externa: tumor vs normal** (`pride_validation_tumor_vs_normal`): **788** proteínas con FDR por debajo de **0.05** (de **6824** filas en el contraste).
  - GO (ORA): **0** términos enriquecidos (FDR por debajo de **0.05**).

---

## Qué **no** sustituye este documento

- No interpreta mecanismos ni discute contraste con la bibliografía clínica: eso es texto de tesis.
- Los avisos de SciPy («catastrophic cancellation») afectan filas con varianza casi nula; no invalidan todo el contraste.

## Cómo cerrar el capítulo de resultados

1. Lee [Resultados base](results_base.md) y valida que las tablas coinciden con tus preguntas de investigación.
2. Abre los TSV completos en `results/` (misma jerarquía bajo tu `data_workspace` si lo usas).
3. Usa el prompt en [`docs/PROMPT_SINTESIS_TESIS.md`](../PROMPT_SINTESIS_TESIS.md) con un LLM pegando extractos de tablas y el contexto de la cohorte.

## Archivos clave para citar en la memoria

- `results/differential/*_diff.tsv` — estadísticos por proteína.
- `results/go/*_go.tsv` — enriquecimiento por término.
- `results/ml/<cohorte>_feature_ranking.tsv` — priorización multivariante.
- `results/validation/input_validation_report.json` — trazabilidad de calidad.
