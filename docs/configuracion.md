# Configuración

Los archivos YAML en `config/` controlan el comportamiento:

- `config.yaml`: semillas, umbrales de QC, GO, red, ML y niveles de evidencia.
- `datasets.yaml`: rutas de matrices y metadatos por cohorte.
- `comparisons.yaml`: comparaciones A vs B enlazadas a columnas de metadata.

En modo `--demo`, el pipeline escribe metadatos runtime en `data/interim/_demo_<dataset>_samples.tsv` y mapas GO sintéticos en `data/interim/`, sin modificar los TSV versionados bajo `metadata/sample_metadata/`.
