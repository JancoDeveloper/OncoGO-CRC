## Publicación web del proyecto

### Run sugerido

```bash
oncogo-crc run --stage download_real
oncogo-crc run --stage validate_inputs
oncogo-crc run --stage publish_results
mkdocs build
```

### Qué publica la web

- estado de datasets
- cobertura clínica
- resultados base
- resultados estratificados
- estabilidad entre subcohortes
- limitaciones

### Objetivo

Mostrar el avance del proyecto, la reproducibilidad del pipeline y la lógica científica de la tesis.
