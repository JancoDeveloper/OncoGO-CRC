# Estructura del proyecto OncoGO-CRC

```
OncoGO-CRC/
├── README.md
├── DOCUMENTACION_TESIS.md
├── GUIA_CONFIGURACION.md
├── ESTRUCTURA_DEL_PROYECTO.md   # este archivo
├── pyproject.toml               # empaquetado pip y dependencias
├── requirements.txt
├── environment.yml              # Conda (opcional)
├── Makefile                     # atajos install/demo/docs
├── mkdocs.yml
├── config/                      # YAML de control
├── data/
│   ├── raw/…                    # matrices crudas por cohorte
│   ├── interim/                 # matrices armonizadas (*_harmonized.tsv)
│   └── processed/               # matrices normalizadas (*_normalized.tsv)
├── metadata/
│   ├── sample_metadata/         # fenotipos por muestra
│   ├── protein_id_maps/         # mapa proteína → gen / UniProt
│   └── gene_ontology/           # anotación larga gen → GO (offline)
├── src/oncogo_crc/              # código versionado
│   ├── cli.py                   # punto de entrada `oncogo-crc`
│   ├── settings.py              # carga de YAML y rutas (incl. demo)
│   ├── stages/                  # una etapa = un archivo `run(settings)`
│   ├── utils/                   # IO, QC, estadística, GO, red, ML, informes
│   └── pipeline/registry.py     # orden canónico de etapas
├── docs/                        # fuentes MkDocs
├── notebooks/                   # análisis exploratorio (plantilla)
├── results/                     # salidas tabulares y figuras
├── reports/executive_summary/   # Markdown consolidado
└── scripts/run_pipeline.sh      # envoltorio bash
```

## Paquete `oncogo_crc`

| Módulo | Función |
|--------|---------|
| `settings.load_settings` | Lee `config/*.yaml`, resuelve rutas absolutas y bandera `--demo`. |
| `utils/io.py` | Lectura/escritura TSV y validación de columnas mínimas. |
| `utils/ids.py` | Armonización y colapso por `gene_symbol`. |
| `utils/qc.py` | Filtros de missingness, log2, KNN, z-score por proteína. |
| `utils/stats.py` | t-test por proteína + corrección FDR (Benjamini–Hochberg). |
| `utils/plots.py` | Volcanos para QC visual de resultados diferenciales. |
| `utils/go_utils.py` | ORA hipergeométrica/Fisher y reducción semántica por solapamiento. |
| `utils/network_utils.py` | Correlación proteína–proteína y ranking de hubs. |
| `utils/ml_utils.py` | Modelos baseline interpretables y score combinado. |
| `utils/translation_utils.py` | Fusión de evidencias y niveles traslacionales. |
| `utils/report_utils.py` | Plantilla de resumen ejecutivo Markdown. |

Las etapas en `stages/` solo orquestan: mantienen delgado el código y facilitan la narrativa metodológica en la tesis (cada capítulo puede mapear 1:1 a un archivo).
