# Prompt para redactar la síntesis de resultados (tesis / artículo)

Copia el bloque siguiente en un asistente (ChatGPT, Claude, Cursor, etc.) y **sustituye** las secciones entre corchetes por tus tablas o rutas reales. El pipeline ya dejó números en `results/` y tablas resumidas en el sitio MkDocs (**Proyecto Web → Resultados base**).

---

## Prompt (copiar desde aquí)

Eres redactor científico en oncología y bioinformática. Tu tarea es producir una **síntesis de resultados** en español para la sección de discusión de una tesis doctoral, usando **solo** la evidencia que te proporciono (tablas TSV / extractos del pipeline OncoGO-CRC). No inventes cohortes ni p-values.

**Contexto del estudio**

- Pregunta de investigación: [PEGA aquí 2–4 frases: hipótesis sobre CRC, proteómica, GO, ML interpretable.]
- Cohortes integradas: [p. ej. CPTAC colorectal PDC000116; PRIDE PXD…; plasma DIA …]
- Comparaciones estadísticas ejecutadas: [lista o pega el bloque “Resumen cuantitativo” de `reports/executive_summary/results_summary.md`.]

**Datos que debes integrar (pegar tablas o resúmenes)**

1. **Diferencial (proteínas):** pega las **15 primeras filas** de `results/differential/cptac_tumor_vs_normal_diff.tsv` (o la comparación principal) ordenadas por `padj`, o el extracto de la web “Resultados base”.
2. **GO (ORA):** pega las **15 primeras filas** significativas de `results/go/cptac_tumor_vs_normal_go.tsv` (columnas `go_id`, `overlap`, `padj`, `odds_ratio`).
3. **Validación externa (si aplica):** mismos extractos para `pride_validation_tumor_vs_normal_*`.
4. **Plasma (si aplica):** extracto `plasma_crc_vs_healthy_*`.
5. **ML / priorización:** primeras **15** filas de `results/ml/cptac_colorectal_feature_ranking.tsv` y las métricas test de `cptac_colorectal_model_metrics.tsv`.
6. **Red:** primeras filas de `results/networks/cptac_colorectal_hub_proteins.tsv`.
7. **Calidad:** un párrafo tomado de `results/validation/input_validation_report.json` (cobertura de variables clínicas, si existe).

**Salida pedida (estructura obligatoria)**

1. **Párrafo guía (150–200 palabras):** qué se comparó, con qué N aproximado de proteínas/muestras, y cuántas proteínas pasaron FDR (cita números del material).
2. **Hallazgos proteicos:** 3–5 frases sobre dirección del efecto (`log2fc`), sin lista interminable; menciona solo **proteínas o familias** que salgan en el top y tengan sentido biológico plausible (marca como *hipótesis generadora* si el apoyo es solo estadístico).
3. **Hallazgos GO:** agrupa términos en **2–4 temas** (p. ej. adhesión, mitosis, respuesta inmune, metabolismo) usando los `go_id` entregados; indica coherencia o tensión con la literatura de CRC **solo** si el usuario pegó referencias; si no, limita a descripción funcional.
4. **Convergencia multi-capa:** un párrafo que relacione CPTAC vs PRIDE vs plasma (concordancia o discrepancia) según los extractos.
5. **ML interpretable:** qué aporta el ranking respecto al univariado (concordancia de `protein_id` en tops); limitaciones (N, batch, hold-out).
6. **Limitaciones metodológicas (4 viñetas):** estadístico simple, ORA/universo, datos públicos congelados en fecha X, generalización clínica.
7. **Frases de cierre “no resultados”:** qué experimento o cohorte faltaría para subir el nivel de evidencia (TRIPOD / REMARK si aplica).

**Restricciones**

- No afirmes “validación clínica prospectiva” si los datos son retrospectivos / públicos.
- Distingue **asociación** de **causalidad**.
- Si falta un bloque de datos, dilo explícitamente y acorta la subsección.

**PEGA AQUÍ TUS TABLAS / JSON**

```

(… tus extractos …)

```

---

## Después de generar el texto

- Cruza cada afirmación con el TSV correspondiente.
- Añade a la tesis la **versión de datos** (`data/versions/…`, manifiestos `download_real`).
- Cita GO por **identificador** (`GO:…`) y, si lo tienes, nombre del término desde una base curada (AmiGO, QuickGO) para la memoria final.
