# Reportes OncoGO-CRC v2

Documentos listos para leer sin ejecutar el pipeline. Los números salen de las tablas de la corrida v2; no están escritos a mano.

## Para lectores sin formación técnica

| Documento | Contenido |
| --- | --- |
| [OncoGO-CRC_guia_publica_v2.pdf](OncoGO-CRC_guia_publica_v2.pdf) | Qué es el proyecto, qué datos se usaron, qué se encontró y qué **no** significa para un paciente. |
| [index.html](index.html) | Misma información que el informe técnico, en el navegador. |

## Informe técnico completo

| Documento | Contenido |
| --- | --- |
| [OncoGO-CRC_informe_v2.pdf](OncoGO-CRC_informe_v2.pdf) | 17 secciones, 12 figuras, correcciones C1–C14 y fuentes D1–D9. |
| [informe_v2.md](informe_v2.md) | Versión Markdown del informe (generada desde tablas). |
| [comparacion_v1_v2.md](comparacion_v1_v2.md) | Qué cambió respecto a la corrida certificada de octubre de 2026. |

## Tablas de resumen

La carpeta `datos/` contiene el certificado, el índice R y el changelog numérico. La carpeta `figuras/` trae cada gráfico en PNG, SVG y su TSV de origen.

## Aviso

Ningún resultado de esta tesis sirve para decidir sobre un paciente.

## Regenerar

Con el disco de datos montado en `/Volumes/JANCODATA/OncoGO-CRC`:

```bash
make v2-report
.venv/bin/python scripts/build_public_report_v2.py --sync-reportes
```
