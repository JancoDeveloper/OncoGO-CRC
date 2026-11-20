# OncoGO-CRC 4.0.0

Flujo de trabajo bioinformático reproducible para analizar proteómica pública
de tejido y de suero. El cáncer colorrectal es el caso de estudio. El
repositorio contiene el código de la corrida v3.7, su configuración, las
pruebas, las salidas de referencia y la documentación para repetirla.

La pregunta de la tesis es: ¿cómo desarrollar y evaluar un flujo de trabajo
bioinformático reproducible que caracterice alteraciones proteómicas en tejido
y suero de cáncer colorrectal y que establezca la robustez y las limitaciones
de sus resultados?

Objetivo general: desarrollar y evaluar ese flujo, aplicado al cáncer
colorrectal, para caracterizar las alteraciones y establecer los límites de
su interpretación.

1. Diseñar e implementar el flujo para tejido y suero de plataformas distintas
   (TMT, LFQ y DIA), con selección, control de calidad, análisis independiente
   de cada cohorte y trazabilidad.
2. Aplicarlo con inferencia estadística y aprendizaje automático, y evaluar
   la capacidad discriminativa frente a la calidad de la muestra y al orden
   de medición.
3. Evaluar la robustez y la reproducibilidad: multiverso, repetición de la
   ejecución, permutaciones, remuestreo y concordancia entre cohortes.

La corrida numérica sellada es la etiqueta
[`v3.7-final`](https://github.com/JancoDeveloper/OncoGO-CRC/tree/v3.7-final)
(commit `595118ae`; el sello de la ejecución es `71b310a3`, 11 de noviembre
de 2026, Python 3.14.3). Esta versión 4.0.0 no reescribe esa historia.

## Qué hay aquí

| Ruta | Contenido |
| --- | --- |
| `src/oncogo_crc/` | Etapa `v3_7`, las piezas que esa etapa llama y la regeneración de `results/certified` |
| `config/preregistro_v3.yaml` | Umbrales, semillas y multiverso. El código lo lee con este nombre |
| `config/esperado_v3_7.tsv` | 35 métricas de la corrida de referencia |
| `results/v3_7/` | Salidas selladas. No se editan |
| `docs/` | Diseño, fuentes, configuración, reproducibilidad, pruebas y `references.bib` |
| `web/` | Sitio estático de la tesis. No se despliega desde este repositorio |
| `data/checksums.sha256` | Huellas de las entradas |

El diseño sigue el capítulo 3 de la tesis: requisitos, arquitectura,
componentes, flujo de datos, actividad, despliegue y pruebas. El detalle
está en [`docs/diseno.md`](docs/diseno.md).

## Entorno

La corrida de referencia usó Python 3.14.3, numpy 2.4.4, pandas 3.0.2,
scipy 1.17.1 y scikit-learn 1.8.0. Esas versiones, más statsmodels 0.14.6,
directlfq 0.3.3, gseapy 1.2.1 y shap 0.51.0, están fijadas en
`pyproject.toml`. El archivo `env/requirements-lock.txt` es el `pip freeze`
de un entorno creado con CPython 3.14.3 en Linux. `env/environment.yml` y
`Dockerfile` describen el mismo entorno (`PYTHONHASHSEED=0` y un solo hilo).

```bash
python3.14 -m venv .venv
.venv/bin/pip install -e ".[dev]"
make test
make lint
```

## Repetir la corrida

Las matrices no están en git. Hacen falta las cuatro cohortes, tres ficheros
de la corrida certificada v1 y los embeddings ESM-2. Las huellas están en
`data/checksums.sha256`. `make datos` baja PRIDE y el depósito
[10.5281/zenodo.24398807](https://doi.org/10.5281/zenodo.24398807).

```bash
export WORKSPACE="$PWD/workspace"
make datos          # PRIDE, contaminantes, matriz LFQ y Zenodo 10.5281/zenodo.24398807
make hpa            # tabla de concentración en sangre; la huella sellada puede no coincidir si el HPA cambió
make verificar-datos
make certified      # solo si data/raw/cptac/cptac_crc_matrix.tsv tiene la huella sellada
make all-v3-7       # regresión v1, crudos de PXD055218 y etapa v3_7
make verificar      # 35 métricas y SHA-256 frente a results/v3_7
```

`make all-v3-7` exige un árbol de git limpio: la etapa no firma si hay
cambios sin registrar. Escribe en `$WORKSPACE/results/v3_7` y no modifica
`results/v3_7` del repositorio.

Los sesenta crudos `.raw` de PXD055218 suman 65,56 GB. Sin ellos, la tabla
sellada `results/v3_7/crudos_integridad.tsv` sigue siendo la entrada de
fechas; una corrida nueva no recupera esas fechas.

### Entradas v1 y depósito Zenodo

`results/certified/cptac_unshared_treat.tsv`,
`results/certified/inference_audit.json`,
`results/certification/certificado_resumen.json`,
`data/raw/cptac/cptac_crc_matrix.tsv` y
`external/modelos/h7_esm2_embeddings.tsv` no están en git. Están en el
conjunto de datos
[10.5281/zenodo.24398807](https://doi.org/10.5281/zenodo.24398807)
(registro <https://zenodo.org/records/24398807>, licencia CC-BY-4.0).
`make descargar` los baja a esas rutas y comprueba el SHA-256. Si la huella
no coincide, el script termina con código 1 y deja el fichero para
inspección. Si ya está y la huella coincide, no lo vuelve a bajar.
`--solo` acepta la ruta relativa, igual que el resto de fuentes.

La matriz `cptac_crc_matrix.tsv` (huella
`a1272baf3bb84218b60d48e971f576f3b5aee6f3fdce4c25a46ad85a352117d1`) se deriva
del estudio de acceso abierto PDC000116 (Vasaikar et al. 2019,
<https://doi.org/10.1016/j.cell.2019.03.030>; Proteomic Data Commons,
<https://pdc.cancer.gov/pdc/study/PDC000116>) y se redistribuye con
atribución. Hay que citar el artículo además del DOI del depósito.
`make certified` regenera los tres ficheros certificados a partir de esa
matriz y de la matriz de PXD034575. La matriz de tejido LFQ sí se
reconstruye byte a byte desde `Cohort1_proteinGroups.txt`.

Los embeddings ESM-2 no entran en el modelo principal. El TSV del depósito
es el de la corrida sellada. `scripts/v3/esm2_embeddings.py` puede
regenerarlo, y pide torch.

## Reproducir localmente con disco externo

Hace falta macOS o Linux, Python 3.14, `make` y `curl`. R no entra en esta
corrida: el único sitio que lo nombra es `scripts/v2/session_info.py`, y se
limita a anotar que no está. No hace falta instalarlo.

El disco de datos se lee y no se modifica. Las salidas van a otro directorio.

```bash
export ONCOGO_DATA_DIR=/Volumes/JANCODATA/OncoGO-CRC
export ONCOGO_CRC_DATA_WORKSPACE=/Volumes/JANCODATA/OncoGO-CRC-v4-salidas
export PIP_CACHE_DIR=/Volumes/JANCODATA/pip-cache
export TMPDIR=/tmp
export COPYFILE_DISABLE=1
export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=core.filemode GIT_CONFIG_VALUE_0=false
python3.14 -m venv .venv
.venv/bin/pip install -e ".[dev]"
make verificar-local DATA_DIR="$ONCOGO_DATA_DIR" SALIDAS="$ONCOGO_CRC_DATA_WORKSPACE"
```

`ONCOGO_DATA_DIR` (o `oncogo-crc run --data-dir`) es la raíz de solo lectura.
Ahí se resuelven `data/raw/`, `metadata/`, `external/` y `results/certified/`,
`results/certification/` y `results/v2/`. `ONCOGO_CRC_DATA_WORKSPACE` recibe
`results/v3_7/`, `reports/`, `data/processed/` y
`external/pxd064655/derivados/`, que la etapa reescribe. Los dos directorios
no pueden ser el mismo. No hace falta copiar las matrices ni enlazar `data/`.

Sin el disco, `make descargar` (o `scripts/descargar_datos.py --destino
<carpeta>`) vuelve a bajar los ficheros públicos de PRIDE: PXD034575,
PXD055218 (`DIA-NN_report.tsv`) y PXD064655, y los cinco ficheros del
depósito [10.5281/zenodo.24398807](https://doi.org/10.5281/zenodo.24398807).
Comprueba el SHA-256 cuando está sellado y acepta `--solo` para no repetir
una descarga grande. No baja los crudos `.raw` (65,56 GB). La matriz TMT se
deriva del estudio de acceso abierto PDC000116 (Vasaikar et al. 2019) y se
redistribuye con atribución desde
<https://zenodo.org/records/24398807>. `make certified` puede regenerar los
certificados cuando esa matriz ya está en el destino.

`make verificar-local` comprueba los hashes de los manifiestos, corre ruff,
pytest y la prueba sintética, y, si el árbol de git no tiene cambios de
contenido, repite la regresión v1, el SHA-1 de los crudos en solo lectura y
la etapa v3.7. Escribe `reports/verificacion_local_<fecha>.md`. El árbol sucio
no se firma: es la regla de la corrida, no un fallo del disco.

`pyproject.toml` no declara `openpyxl`. H7, H9 y la triangulación leen
libros xlsx y pandas lo pide. En este Mac hizo falta `openpyxl==3.1.5`
(y `et-xmlfile==2.0.0`); no está en `env/requirements-lock.txt`.

Espacio: el plasma ya descargado ocupa unos 74 GB y `external/` unos 8,6 GB;
no se vuelven a copiar. La carpeta de salidas queda en el orden de
`results/v3_7` (74 MB en la referencia) más los derivados de PXD064655. El
tiempo largo es la etapa v3.7: en este Mac tardó 938 s (15 min 38 s).
El SHA-1 de los 60 crudos, sin descargar, tardó 737 s. `pip install` usa
`TMPDIR=/tmp`: en ExFAT setuptools falla al armar el wheel editable.
`make crudos-v3-7` no se usa contra el disco
de datos: si un SHA-1 no coincide, ese target descarga el crudo y mueve el
viejo a `cuarentena/` dentro de la carpeta de crudos.

## Pruebas

```bash
make test          # pytest y cobertura; omite la prueba e2e
make test-e2e      # la etapa rechaza un report sintético y el certificado tiene 22 casillas
make lint          # ruff 0.15.12, reglas E4, E7, E9 y F
```

No hay integración continua. La verificación de la corrida es local, con
`make verificar-local`. Cómo reportar un secreto está en `SECURITY.md`.

## Cita

El código de esta limpieza no tiene DOI. Para citarlo, `CITATION.cff`
(versión 4.0.0). Para los números de la tesis, la etiqueta `v3.7-final`.
Las entradas que no están en git (matriz CPTAC derivada de PDC000116,
ficheros certificados y embeddings ESM-2) tienen DOI
[10.5281/zenodo.24398807](https://doi.org/10.5281/zenodo.24398807),
licencia CC-BY-4.0.

Licencia MIT. Ver `LICENSE`.

## Sitio

Las páginas están en `web/` y se leen sin compilar, en el navegador o así:

```bash
python3 -m http.server 8000 --directory web
```

Abra `http://127.0.0.1:8000/`. También puede abrir `web/index.html` directo.
El texto está en castellano, para un lector general y para un lector clínico.

Las gráficas salen de `results/v3_7` y de `docs/fuentes_de_datos.md`:

```bash
python3 scripts/web/figuras_sitio.py
```

`firebase.json` publica la carpeta `web/` en el proyecto de Hosting
`oncogo-crc-tesis` (https://oncogo-crc-tesis.web.app). No contiene claves.
Hasta que se despliegue esta versión, lo que haya en esa dirección puede ser
anterior.
