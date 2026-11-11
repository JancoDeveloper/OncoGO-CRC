# Estado de la v3.7

Este archivo lo escribe `scripts/v3/exportar_cifras_v3_7.py` a partir de `results/v3_7/`.

Certificado: Publicable.
Sello: commit `71b310a3fcae`, fecha 2026-11-11T00:48:45Z.
Etiqueta: `v3.7-final`, commit `595118ae09ea6006c744c86b28cb0a10d42e0442`.
`sha256_manifest_entradas` en el sello: e03f0c71107e41a3e82f6b105d680ff11db702ed89f8a2f3a146546896062a37.

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

La corrida versionada es E. El certificado es publicable. D y E salen del mismo commit y son idénticas en los números: solo cambian la fecha, el pid, las horas, los tiempos de `tests.txt` y el orden de `gsea_tejido.tsv` (`docs/determinismo_d_e.md`).
B y C, con el sello `b7f8b8e844f6`, no eran idénticas en el secundario de H7 (`docs/determinismo_b_c.md`). La causa era ESM-2 en MPS, que omitía errores sin registro: 181 de 668 proteínas con secuencia tenían vector. La enmienda `enmienda_v3_7_esm2` lo calcula en CPU, fuera de la corrida. El detalle está en `docs/diagnostico_esm2.md` y `docs/cambio_esm2_b_d.md`.
La web pública es la documental, en https://oncogo-crc-tesis.web.app. Hosting `ed65f62acbad33dc`, liberada el 2026-11-11T18:31:00.695Z, sitio `ae7416e`. El despliegue está en `docs/despliegue_sitio_v4.md`. La v3.7 anterior queda en `/archivo/v3_7` y su despliegue, en `docs/despliegue_v3_7.md`.
