# Estado de la v3.7

Este archivo lo escribe `scripts/v3/exportar_cifras_v3_7.py` a partir de `results/v3_7/`.

Certificado: Publicable.
Sello: commit `b7f8b8e844f6`, fecha 2026-11-10T21:28:31Z.
`sha256_manifest_entradas` en el sello: 864342f54dd6df16160e65aa85d4b33614cd149f39d131dbb38a7acec2a3db8a.

## Unión de los sanos

En `union_crudos_calidad.tsv`, `verificada` es verdadero en 40 filas y falso en 20.
Las filas falsas son las corridas `Normal_` y la evidencia es «el report no nombra ese crudo». No se inventó `HC_n` = `Normal_n`.

## Ficheros de la especificación que no están

Presentes: `go_ora_tejido.tsv`, `manifest_entradas.sha256`.
Ausentes:
- `vias_suero.json`
- `h1_vias.json`
- `catalogo_auditoria.tsv`
- `neoantigenos.tsv`
- `prtc_por_corrida.tsv`

Esas ausencias están declaradas en `enmienda_v3_7_cierre`. No son un resultado.

La corrida versionada terminó en un proceso. El certificado y el sello son los de `results/v3_7/certificado.json`. La comparación con la corrida archivada está en `docs/determinismo_dos_procesos.md`: los números principales coinciden. Los sellos no, porque la segunda corrida declara `results/certified/inference_audit.json`.
