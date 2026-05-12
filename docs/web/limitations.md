# Limitaciones (auto)

- **Modelo estadístico:** t-test por proteína + FDR Benjamini–Hochberg; asume razonablemente muestreo independiente entre grupos; no sustituye modelos mixtos si hay batches fuertes no corregidos.
- **GO:** ORA con universo = proteínas detectadas en la matriz; términos muy grandes o redundantes pueden dominar; FDR término a 0.05.
- **ML:** hold-out simple (`test_size={ml.get('test_size', 0.25)}`); métricas en conjunto de prueba, no validación externa salvo que añadas cohortes.
- **Datos públicos:** versiones de matrices y metadata dependen de PDC/PRIDE en la fecha de descarga; documenta `data/versions/` y manifiestos.
- **α diferencial:** 0.05 — ajustar solo con criterio metodológico defendible en tesis.
