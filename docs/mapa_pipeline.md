# Mapa del pipeline (Fase 0 de v2)

Documento de reconocimiento. Dice dónde está cada cosa antes de cambiar nada.
Fecha: 3 de octubre de 2026. Commit de referencia: `f88eec3`.

## 1. Lenguaje y punto de entrada

Todo el pipeline es **Python**. No hay R instalado en la máquina y no hay código R
en el repositorio. Esto importa para v2 y se trata en la sección 6.

- Paquete: `src/oncogo_crc/`, instalado como editable (`oncogo-crc` 0.1.0).
- Entrada única: `oncogo-crc run --stage <nombre> [--data-workspace RUTA]`,
  definida en `src/oncogo_crc/cli.py`.
- Registro de etapas: `src/oncogo_crc/pipeline/registry.py`, diccionario `STAGES`.
  Cada valor es una función `run(settings) -> None`.
- Configuración: `config/config.yaml` (umbrales y rutas), `config/datasets.yaml`
  (matrices y escalas), `config/comparisons.yaml` (contrastes).
- Resolución de rutas: `src/oncogo_crc/settings.py`, método `Settings.resolve`.
  Redirige `data/`, `results/`, `reports/`, `metadata/` y `docs/web` al
  `data_workspace` (el disco externo). El código y la configuración siguen en el repo.

### Etapas de la corrida certificada v1

| Etapa (`--stage`) | Módulo | Qué hace |
|---|---|---|
| `certified_inference` | `stages/certified_inference.py` | Contraste primario certificado: CPTAC pareado, PRIDE, cruce, ajuste estromal, permutaciones |
| `translational_axes` | `stages/translational_axes.py` | Pólipos, brecha tejido–plasma, panel externo, selección de modelo |
| `secondary_contrasts` | `stages/secondary_contrasts.py` | Sexo, estadio, plasma e inventario de variables |
| `certification` | `stages/certification.py` | Certificado R, informe de brecha, `indice_R.tsv` |

El resto de etapas (`differential`, `go`, `gsea`, `kegg`, `network`, `ml`,
`prioritization`, `stability`, `report`…) son del pipeline anterior y **no**
alimentan el certificado. `certification` prefiere `results/certified/` sobre
`results/reproducibility/` cuando el primero existe.

## 2. Dónde está la t moderada con `treat`

**No se usa limma ni R.** `treat` está reimplementado en Python en
`src/oncogo_crc/utils/moderated.py`:

| Función | Qué implementa | Referencia |
|---|---|---|
| `fit_f_prior(s2, df)` | Prior F de la varianza por momentos: devuelve `d0` y `s02` | Smyth 2004 |
| `trigamma_inverse(y)` | Inversa de la trigamma por Newton, necesaria para `d0` | Smyth 2004 |
| `posterior_scale(...)` | Varianza posterior `s̃²` y grados de libertad totales | Smyth 2004 |
| `treat_pvalues(coef, se, df_total, lfc)` | p-valor contra `\|log2FC\| > lfc`, mínimo de las dos colas desplazadas | McCarthy y Smyth 2009 |
| `_bh(p, alpha)` | Benjamini–Hochberg | |
| `moderate_table(...)` | Une lo anterior en una tabla con `log2fc`, `se`, `padj`, `significant` |
| `paired_observed_stats(diffs, min_pairs)` | Media de diferencias pareadas sobre valores observados, sin imputar |
| `paired_covariate_stats(diffs, covariate, at)` | `diff = a + b·x`, efecto evaluado en `x = at`; base de R4 |
| `unpaired_observed_stats(...)` | Dos grupos independientes (plasma, estadio) |
| `n_significant_after_signs(...)` | Recuento de significativas para el nulo de permutación |

Consecuencia para v2: cualquier método nuevo que el diseño nombre con un
paquete de R hay que implementarlo en Python o instalar R. Ver sección 6.

## 3. Dónde se calcula cada subíndice y cada eje

| Pieza | Fichero | Función | Salida |
|---|---|---|---|
| R1 concordancia de signo | `certification/contract.py` | `sign_concordance_discovery` | `certificado_resumen.json` |
| R2 detectabilidad | `stages/certified_inference.py` | `_peptide_fraction` | clave `peptide_*` del audit |
| R2 variante por log2FC | `certification/contract.py` | `measurable_fraction_two_cohorts` | `measurable_fraction_log2fc` |
| R3 homogeneidad | — | no implementado | `no_calculable` |
| R4 índice estromal | `stages/certified_inference.py` | `stromal_score`, `purity_adjusted_fraction` | `r4_stroma` del audit |
| R5 nulo de permutación | `stages/certified_inference.py` | `_permute_pair_signs` | `r5_empirical_p` |
| Concordancia de vía | `stages/certified_inference.py` | `_pathway_concordance` | `pathway_*` del audit |
| Cláusulas H1 y H4 | `certification/contract.py` | `h1_protein_clause`, `h1_pathway_clause`, `h1_status`, `h4_proxy_verdict` | |
| Pólipos y secuencia adenoma | `stages/translational_axes.py` | `_paired_arm`, `_adenoma_sequence` | `axes_audit.json` |
| Brecha tejido–plasma | `stages/translational_axes.py` | `_tissue_to_plasma` | `tejido_a_plasma.tsv` |
| Panel externo de plasma | `stages/translational_axes.py` | `_external_panel`, constante `PLASMA_PANEL` | `panel_externo_plasma.tsv` |
| Selección de modelo | `stages/translational_axes.py` | `_model_selection` | `model_selection` del audit |
| Sexo, estadio, plasma | `stages/secondary_contrasts.py` | `_sign_agreement`, `_stage_contrast`, `_plasma_contrast` | `secondary_audit.json` |
| Red (H2) | `stages/network_modules.py` | — | **sin nulo que conserve grados** |

## 4. Dónde se generan las figuras y el PDF

- Figuras del PDF: `scripts/build_thesis_pdf.py`, funciones `fig_architecture`,
  `fig_volcano`, `fig_scatter`, `fig_counts`, `fig_stroma`, `fig_adenoma`,
  `fig_gap`, `fig_models`, `fig_subindices`. Salen a
  `results/secondary/figures/` en PNG a 160 dpi.
- PDF: el mismo script, con `reportlab` y la fuente Arial del sistema.
  Salida: `<workspace>/OncoGO-CRC_tesis_explicada.pdf`.
- **El texto del PDF tiene números escritos a mano.** Las cifras se interpolan
  desde los JSON en la mayoría de los párrafos, pero hay valores literales
  (por ejemplo `−1,14`, `+1,23`, `0,88`) escritos en las cadenas. Eso es lo que
  C14 obliga a corregir.
- Figuras del pipeline anterior: `utils/plots.py` y `utils/web_visualizations.py`,
  a `results/figures/` y `docs/web/figures/`.

## 5. Entradas, manifiestos y entorno

Entradas reales de la corrida certificada, todas en el disco externo:

| Fichero | Tamaño | Papel |
|---|---|---|
| `data/raw/cptac/cptac_crc_matrix.tsv` | 49 MB | PDC000116, log-ratio TMT, ya en log2 |
| `data/raw/pride/pride_crc_matrix.tsv` | 69 MB | PXD034575, LFQ MaxQuant |
| `data/raw/pride/Cohort1_proteinGroups.txt` | — | péptidos únicos para R2 |
| `data/raw/plasma/plasma_crc_matrix.tsv` | 896 KB | PXD055218, suero DIA |
| `metadata/protein_id_maps/*` | — | puente a UniProt |
| `metadata/gene_ontology/go_gene_long.tsv` | — | términos GO, incluido GO:0031012 |

- Manifiestos: `scripts/v2/manifest_hashes.py`. Escribe
  `results/v2/manifest_inputs.sha256` (11 ficheros de entrada) y
  `results/v1_manifest.sha256` (31 ficheros de v1, PDF incluido). Excluye los
  `._*` de AppleDouble que ExFAT crea en el disco externo.
- Entorno: `requirements-lock.txt` con versiones exactas y
  `results/v2/session_info.txt` con Python, plataforma, commit y paquetes.
  Generados por `scripts/v2/session_info.py`.
- Orquestador: `Makefile`, un objetivo por etapa y `all-v2`.
- Regresión: `tests/test_regresion_v1.py`. Verificado el 3 de octubre de 2026:
  re-ejecutar el contraste en un workspace temporal con las entradas enlazadas
  reproduce **688, 140, 28, 104/119, 226 y 307**. v1 no se modifica.

## 6. Discrepancias entre el diseño de v2 y el repositorio

Se documentan aquí y se resuelven en `docs/enmiendas_v2.md`.

1. **R no está instalado.** `limma::treat`, `metafor::rma`, `limma::camera`,
   `mada::reitsma`, `renv.lock` y `testthat` no se pueden usar. `treat` ya está
   reimplementado y validado en Python, así que v2 sigue el mismo camino:
   REML, la prueba competitiva de conjuntos y el meta-análisis diagnóstico se
   implementan en Python con tests propios, y se declara la equivalencia.
2. **R2 circula con tres denominadores distintos.** El valor publicado es
   226/688 = 0,3285 (péptidos únicos en PXD034575). Pero
   `measurable_fraction_log2fc` es 119/702 = 0,1695, y 702 son **filas** del
   cruce, no proteínas: 687 proteínas únicas generan 702 filas porque una
   proteína de CPTAC puede mapear a varios UniProt. Y `n_discovery_significant`
   del certificado es 244, que son las filas significativas restringidas a
   `_merge == "both"`. El prompt de v2 describe la discrepancia como 688 contra
   687; lo real es 688 contra 702 contra 244, y la causa es contar filas
   duplicadas por el puente de identificadores.
3. **R1 se calcula sobre el 17 % de los hallazgos.** 119 proteínas evaluables de
   688 certificadas. El denominador no son las proteínas significativas sino las
   que además tienen log2FC finito y distinto de cero en PRIDE. Una concordancia
   de 0,874 medida solo en lo que la segunda plataforma alcanza a medir está
   condicionada a la detectabilidad. Es la razón de C9.
4. **Las definiciones de R4 y R5 en el JSON están obsoletas.** Confirmado en
   `stages/certification.py`: el bucle actualiza `value`, `status`, `name` y
   `note` del subíndice, pero nunca `definition`. Por eso R4 figura como
   `calculado` y a la vez dice «La pureza tumoral aún no es covariable».
5. **R5 está en el piso del nulo.** `N_PERM = 100` y el p reportado es
   1/101 = 0,0099, el mínimo alcanzable. No es una estimación, es un tope.
6. **H2 no tiene nulo.** `stages/network_modules.py` calcula hubs por grado y no
   tiene ningún modelo nulo, ni de reconexión ni de otro tipo. H2 queda
   `no_evaluada` y no se interpreta ninguna red.
7. **El PDC GraphQL no responde a un GET simple** (tiempo de espera agotado en la
   comprobación de red). GDC, PRIDE, cBioPortal y Human Protein Atlas sí
   responden con código 200. La búsqueda del proteoma retrospectivo de colon
   (D2) depende de esa API, así que puede acabar en `no_disponible`.
8. **Nombres de fichero del diseño que no existen todavía**: no hay
   `config/preregistro_v2.yaml`, ni `docs/enmiendas_v2.md`, ni
   `docs/fuentes_datos.md`, ni `scripts/download/`, ni `results/v2/`. Se crean en
   v2. El diseño también menciona `results/v2_<AAAAMMDD>/` como alternativa; se
   usa `results/v2/` fijo para que el `Makefile` sea reproducible.
9. **El disco externo es ExFAT.** No guarda permisos de POSIX ni enlaces
   simbólicos y crea un `._fichero` por cada fichero. Los manifiestos y los
   recorridos de directorio tienen que excluirlos, y no se puede confiar en el
   bit de solo lectura del sistema de ficheros para proteger v1: la protección
   es el manifiesto de hashes más la separación de carpetas.
