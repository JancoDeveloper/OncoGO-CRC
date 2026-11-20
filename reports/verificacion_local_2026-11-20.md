# Verificación local 2026-11-20

Rama nueva: `limpieza-v4.0.0-b`, worktree `/Volumes/JANCODATA/OncoGO-CRC-v4b`.
HEAD: `deca3726f72c186ad8bb46c35020b89f8edeabe1`.
No hay push. `cursor/limpieza-v4.0.0-39fb` sigue en `f617f78`.
`v3.7-final` sigue en `595118ae`.

## Entorno

- macOS 27.2.1, Apple M4, 10 núcleos, 16 GiB
- Python 3.14.3, numpy 2.4.4, pandas 3.0.2, scipy 1.17.1, scikit-learn 1.8.0, openpyxl 3.1.5
- R no se usó
- datos, solo lectura: `/Volumes/JANCODATA/OncoGO-CRC`
- salidas: `/Volumes/JANCODATA/OncoGO-CRC-v4-salidas-v3/results/v3_7`
- metadata leída: la del repositorio, 30 040 bytes, sha256 `c0f1d5d7…`

## Base y commits

`git merge-base origin/main cursor/limpieza-v4.0.0-39fb` = `73e46b9f0d4f534d59c8d6e8559e8ca8417e7d24`.
`v3.2-final` y `v3.3-final` son ancestros de esa base. Entre cada etiqueta y la base, `git rev-list` sobre su carpeta sale vacío.

| commit viejo | commit nuevo | mensaje |
| --- | --- | --- |
| `7d46d97` | `49d4a613d215540cc5e3a37fd7a925aaedaeaa0d` | Release the repository as v4.0.0 around the sealed v3.7 run. Reaplicado sin tocar `results/v3_2` ni `results/v3_3`. |
| `ae2e2ea` | omitido | Ya no hace falta restaurar esas carpetas. |
| `53e77db` | `c632ff1642f4a6b9767ef460d09618ca90267a7f` | Declarar openpyxl 3.1.5 y et-xmlfile 2.0.0. |
| `969c554` | `4c9b8732dcc6472c1322f82e8ee42ae3fe40436a` | Leer los insumos con ONCOGO_DATA_DIR y verificar la corrida en local. |
| `54eb8c6` | `3144cd340ab996c3d416c4afdc08a4deff0c5642` | Conservar la ruta lógica cuando la entrada es un enlace. |
| `f617f78` | `deca3726f72c186ad8bb46c35020b89f8edeabe1` | Documentar qué metadata CPTAC usa la corrida v3.7. |

Ningún commit de la rama nueva toca `results/v3_2` ni `results/v3_3`. `.gitignore` no ignora esas carpetas.

## Paso 4

- `git rev-list v3.2-final..HEAD -- results/v3_2`: vacío
- `git rev-list v3.3-final..HEAD -- results/v3_3`: vacío
- `git diff --stat f617f78 HEAD`: vacío. El árbol final es el de la rama vieja.
- `v3.7-final`: `595118ae09ea6006c744c86b28cb0a10d42e0442`

## Pruebas y corrida

| comando | resultado | tiempo |
| --- | --- | --- |
| ruff | pasa | |
| pytest con cobertura | 134 pasan, 0 fallan, 5 se omiten, cobertura 30% | 49.1 s |
| e2e sintético | 2 pasan | 9.9 s |
| `make v3-7` | código 0, certificado `publicable` | 1068 s |

La pytest de dentro de la etapa: `134 passed, 5 skipped`. Las dos pruebas de cronología pasan. Un primer pytest, antes de borrar 885 ficheros `._*` que el checkout creó en ExFAT, tuvo 3 fallos del sitio. No son fallos del código. Al borrarlos, el pytest completo queda en 134.

## Antes y después

| | `7d46d97` | `f617f78` | `limpieza-v4.0.0-b` |
| --- | --- | --- | --- |
| `rev-list` de v3_2 y v3_3 | no vacío | no vacío | vacío |
| cronología | falla | falla | pasa |
| certificado | `no_publicable` | `no_publicable` | `publicable` |
| casilla tests | falla, código 1 | falla, código 1 | pasa, código 0 |
| 35 métricas | 35, fuera 0 | 35, fuera 0 | 35, fuera 0 |
| suero | 284 = 0 / 166 / 118 | 284 = 0 / 166 / 118 | 284 = 0 / 166 / 118 |
| determinismo M2 | pasa | pasa | pasa |
| idénticos byte a byte | 75 | 76 | 76 |
| sellos | 125 | 125 | 125 |
| distintos | 42 | 41 | 41 |
| faltantes | 0 | 0 | 0 |
| `h7_variables.tsv` | solo el nombre | idéntico | idéntico |
| ruta TMbed | `h7_tmbed.pred` | `external/modelos/h7_tmbed.pred` | `external/modelos/h7_tmbed.pred` |

`gsea_tejido.tsv` tiene las mismas 889 filas en otro orden. `h7_visibilidad_en_sangre.json` y `secuencia.json` solo cambian el sello embebido, la fecha o el pid. `certificado.json` no es idéntico byte a byte porque anota este commit, las horas y el SHA de ficheros que llevan el sello nuevo. Su estado es `publicable`, igual que el de la referencia.

## Subir la rama

No ejecutado:

```bash
cd /Volumes/JANCODATA/OncoGO-CRC-v4b
git -c core.filemode=false push -u origin limpieza-v4.0.0-b
```

Texto para el PR nuevo:

> Esta rama reemplaza al PR #1 (`cursor/limpieza-v4.0.0-39fb`).
>
> El PR #1 deja el certificado de `make v3-7` en `no_publicable`. El commit `7d46d97` borra `results/v3_2` y `results/v3_3` después de `v3.2-final` y `v3.3-final`. Aunque `ae2e2ea` las restaura y el contenido vuelve a coincidir, `git rev-list` sigue viendo esos dos commits y las pruebas de cronología fallan.
>
> `limpieza-v4.0.0-b` parte de la misma base (`73e46b9`) y reaplica el mismo árbol, con una diferencia de historia: ningún commit toca esas carpetas. `ae2e2ea` no se reaplica porque ya no hay nada que restaurar. El `git diff` contra `f617f78` está vacío. `v3.7-final` no se mueve.
>
> En esta rama las dos pruebas de cronología pasan, pytest queda en 134 pasan y 5 se omiten, y `make v3-7` firma el certificado como `publicable`. Las 35 métricas coinciden y la celda de suero sigue en 284 = 0 / 166 / 118.
>
> La rama `cursor/limpieza-v4.0.0-39fb` queda como respaldo y no se reescribe.
