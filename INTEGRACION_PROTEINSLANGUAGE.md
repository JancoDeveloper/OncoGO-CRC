# Integración OncoGO-CRC ↔ ProteinsLanguage (y publicación tipo “proteína en contexto”)

## Qué aporta cada proyecto

| Aspecto | **OncoGO-CRC** | **ProteinsLanguage** (`ProteinsLanguage/Alpha/`) |
|--------|------------------|---------------------------------------------------|
| Objetivo | Pipeline reproducible multi-cohorte: QC, diferencial, GO **ORA offline**, red, ML, informes | Exploración rica por proteína: Dash, MongoDB, conteos GO por MF/BP/CC, AlphaFold (`dash_bio`), tablas |
| Datos | CPTAC / PRIDE / plasma vía YAML + `data_workspace` | Colección curada en Mongo (`Magnoliopsida`, etc.) |
| Fortaleza | Estadística cohorte → listas priorizadas + trazabilidad | Vista interactiva + estructura + anotación densa |

No son sustitutos: **OncoGO-CRC produce listas y números por contraste**; **ProteinsLanguage enriquece visualmente un universo proteico ya cargado en BD**.

## Mejora reciente en GO (OncoGO-CRC)

Dos causas habituales de **0 términos GO** en ORA:

1. **Genes de foreground demasiado numerosos** (p. ej. ~80 % del proteoma “significativo”): el test de Fisher deja de ver *enriquecimiento* frente al universo.  
   → En `config/config.yaml` se añadió `go.ora_foreground_max_genes` (por defecto **800**): ORA usa las filas significativas **ordenadas por `padj`** hasta ese tope.

2. **Identificadores**: ORA usa anotación en **`gene_symbol`**. El paso `go` ahora mapea explícitamente `protein_id` → `gene_symbol` desde la matriz `*_normalized.tsv`.

Tras cambiar la config, vuelve a ejecutar: `oncogo-crc run --stage go` (y luego `publish_results`).

## Puente de datos (recomendado)

1. Ejecuta el pipeline hasta tener `results/differential/` y `results/ml/`.
2. Genera el bundle JSON:

```bash
python scripts/export_protein_context_bundle.py --comparison cptac_tumor_vs_normal --top 300
```

Salida por defecto: `results/export/protein_context_bundle.json` (respeta `data_workspace` si está activo).

3. **Ingesta en ProteinsLanguage**: adapta un script de carga Mongo (similar a `Alpha/BD/Functions/datos.py`) que lea `gene_symbol` / `uniprot_accession` del JSON y llene o cruce con tu esquema de documentos (`GO_Annotations`, `AlphaFold_Prediction`, etc.).

4. **Seguridad**: en `Alpha/web/app.py` aparece un URI de Mongo con credenciales embebidas; conviene **rotar claves** y usar variables de entorno, no texto plano en el repo.

## Publicación “estilo artículo” unificada (propuesta de flujo)

1. **Cohortes y estadística** → OncoGO-CRC (figuras `results/figures/*_volcano.png`, tablas diferenciales).
2. **Priorización** → `results/ml/*_feature_ranking.tsv` + traslación `results/clinical_translation/`.
3. **Contexto mecanístico por proteína** → ProteinsLanguage o servicios externos (UniProt API, AlphaFold DB, STRING) alimentados con el bundle.
4. **Discusión** → `docs/PROMPT_SINTESIS_TESIS.md` + tablas de **Resultados base** en MkDocs.

## Comparaciones múltiples y “cambiar supuestos”

- En OncoGO-CRC: edita `config/comparisons.yaml` (grupos, columnas) y vuelve a correr `differential` → `go` → `export_protein_context_bundle.py` con otro `--comparison`.
- En ProteinsLanguage: cambia consultas/callbacks de Dash o filtros en Mongo según el subconjunto exportado.

## Próximos desarrollos (si quieres un solo producto)

- Servicio ligero que reciba el JSON del bundle y devuelva paneles reutilizando `dash_cytoscape` / `dash_bio` de ProteinsLanguage.
- GSEA o `fgsea` como etapa opcional post-diferencial (señal continua, mejor con listas largas que ORA).
- Tabla **CPTAC ∩ PRIDE** de proteínas significativas en ambos (script de intersección sobre dos `*_diff.tsv`).
