# OncoGO-CRC

Este sitio resume la arquitectura del proyecto. La **documentación orientada a tesis** (contexto científico, fases, interpretación de resultados y equivalencias con el paquete original) está en la **raíz del repositorio** con estos nombres:

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
