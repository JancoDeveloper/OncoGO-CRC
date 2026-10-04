# Entrega OncoGO-CRC v2 — resumen para el repositorio

Corrida integrada en la rama `main` el 3 de octubre de 2026. La corrida certificada v1 queda congelada bajo `results/v1_manifest.sha256`.

## Qué se hizo

1. **Correcciones metodológicas C1–C14** sobre la corrida certificada de v1, sin modificar matrices crudas ni sobrescribir v1.
2. **Nueve fuentes externas D1–D9**: TCGA (CEA, pureza), heterogeneidad intratumoral, deconvolución con referencia unicelular, plasma frente a enfermedad benigna, meta-análisis de marcadores, protocolo de cohorte futura, entre otras.
3. **Informe reproducible** en Markdown, HTML y PDF, con verificador automático de números (C14).
4. **Guía pública** en PDF para lectores sin formación en bioinformática.

## Dónde está cada cosa

| Qué | Dónde |
| --- | --- |
| Guía para público general | `docs/reportes/OncoGO-CRC_guia_publica_v2.pdf` |
| Informe técnico PDF | `docs/reportes/OncoGO-CRC_informe_v2.pdf` |
| Código del pipeline v2 | `src/oncogo_crc/stages/v2_*.py`, `src/oncogo_crc/v2/` |
| Preregistro y enmiendas | `config/preregistro_v2.yaml`, `docs/enmiendas_v2.md` |
| Fichas de fuentes | `docs/fuentes_datos.md` |
| Resultados completos (disco de datos) | `/Volumes/JANCODATA/OncoGO-CRC/results/v2/` |
| Manifiesto de entradas | `results/v2/manifest_inputs.sha256` (en disco de datos) |

## Resultados principales (v2)

- **688** proteínas pasan el contraste en CPTAC (94 pares); la dirección se repite en un **87,4 %** entre cohortes (R1).
- **H1 rechazada** por la cláusula de proteína: la concordancia es mayor de lo esperado, no menor.
- Tras ajustar por composición de tejido, solo el **44,8 %** de los hallazgos sobrevive (R4): buena parte de la señal es mezcla de tejidos.
- El tejido **no predice** la sangre con los datos disponibles; el CEA depende poco de CEACAM5 en el tumor.
- Separar cáncer de **sano** (AUC ≈ 0,99) no es lo mismo que separar cáncer de **enfermo** (AUC ≈ 0,80).

## Cómo reproducir

```bash
make all-v2
make v2-report
.venv/bin/python scripts/build_public_report_v2.py --sync-reportes
```

## Aviso

Ningún resultado de esta tesis sirve para decidir sobre un paciente.
