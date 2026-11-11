# OncoGO-CRC

## Empieza aquí (v3.7)

- [Tesis](../tesis-ucaldas/overleaf_v3_7.zip)
- [Estado](estado_v3_7.md)
- [Fuentes](fuentes_datos.md)
- [Determinismo D frente a E](determinismo_d_e.md) y [B frente a C](determinismo_b_c.md)
- [Desviaciones](../tesis-ucaldas/latex/chapters/11-desviaciones.tex)
- [Compilación](compilacion_tesis_v3_7.md)
- Web: la página pública sigue en la v3.6, [oncogo-crc-tesis.web.app](https://oncogo-crc-tesis.web.app). Lo que falta del sitio v3.7 está en [web_pendiente.md](web_pendiente.md).

## Reportes v2 (PDF)

- [Guía pública para lectores sin formación técnica](reportes/OncoGO-CRC_guia_publica_v2.pdf)
- [Informe técnico completo v2](reportes/OncoGO-CRC_informe_v2.pdf)
- [Índice de reportes y figuras](reportes/README.md)
- [Resumen de entrega v2](ENTREGA_V2.md)

Este sitio resume la arquitectura del proyecto. La **documentación orientada a tesis** (contexto científico, fases, interpretación de resultados y equivalencias con el paquete original) está en la **raíz del repositorio** con estos nombres:

- [GUIA_ESTUDIANTE_CERTIFICACION.md](GUIA_ESTUDIANTE_CERTIFICACION.md) — guía de estudio: conceptos, diagramas, hipótesis y orden de escritura
- `GUIA_EJECUCION_Y_REDACTION_TESIS.md` — checklist + pipeline → capítulos de la tesis
- `PLAN_EJECUCION_TESIS.md` — roadmap por fases
- `DOCUMENTACION_TESIS.md`
- `GUIA_CONFIGURACION.md`
- `ESTRUCTURA_DEL_PROYECTO.md`

## Inicio rápido

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
oncogo-crc run --demo --stage all
```

Los artefactos esperados coinciden con el checklist de `Informacion/Paquete/docs/06_checklist_operativo.txt`.
