# Estado de la v3.7

Este archivo lo escribe `scripts/v3/exportar_cifras_v3_7.py` a partir de `results/v3_7/`.

Certificado: Publicable.
Sello: commit `c1c22898e01f`, fecha 2026-11-10T18:31:55Z.
`sha256_manifest_entradas` en el sello: None.

## Unión de los sanos

En `union_crudos_calidad.tsv`, `verificada` es verdadero en 40 filas y falso en 20.
Las filas falsas son las corridas `Normal_` y la evidencia es «el report no nombra ese crudo». No se inventó `HC_n` = `Normal_n`.

## Ficheros de la especificación que no están

Ausentes:
- `vias_suero.json`
- `h1_vias.json`
- `go_ora_tejido.tsv`
- `catalogo_auditoria.tsv`
- `neoantigenos.tsv`
- `prtc_por_corrida.tsv`
- `manifest_entradas.sha256`

Esas piezas quedan pendientes. No hay un resultado que poner en su lugar.

Qué falló en la ejecución: `make v3-7` salió con código 2 al empezar H9, sin traza de Python. H9, H10, la preanalítica, la triangulación y el certificado se calcularon después, con el mismo sello. `lecturas_v3_7.json` solo registra las lecturas de ese segundo proceso. `determinismo.json` solo repite M2.
