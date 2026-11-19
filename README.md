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

Las matrices no están en git. Hacen falta las cuatro cohortes y tres ficheros
de la corrida certificada v1. Las huellas están en `data/checksums.sha256`.

```bash
export WORKSPACE="$PWD/workspace"
make datos          # PXD034575, PXD055218, PXD064655, contaminantes y matriz LFQ
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

### Entradas v1

`results/certified/cptac_unshared_treat.tsv`,
`results/certified/inference_audit.json` y
`results/certification/certificado_resumen.json` no están en git ni en un
depósito. La etapa v3.7 los lee y comprueba su SHA-256. `make certified`
los regenera con `certified_inference` y `certification` a partir de
`data/raw/cptac/cptac_crc_matrix.tsv` (huella
`a1272baf3bb84218b60d48e971f576f3b5aee6f3fdce4c25a46ad85a352117d1`) y de la
matriz de PXD034575. Esa matriz de tejido LFQ sí se reconstruye byte a byte
desde `Cohort1_proteinGroups.txt`. La matriz TMT de PDC000116 no tiene un
único archivo público identificado en este repositorio: el lugar previsto
para depositarla, junto con los tres ficheros certificados, es
<https://zenodo.org/deposit/PENDIENTE>. Hasta ese depósito, `make all-v3-7`
no puede terminar en un equipo que no tenga ya el workspace de la corrida.

Los embeddings ESM-2 (`external/modelos/h7_esm2_embeddings.tsv`, huella en
el manifiesto) tampoco están en git. El modelo principal no los usa. Se
regeneran con `scripts/v3/esm2_embeddings.py`, que pide torch.

## Pruebas

```bash
make test          # pytest y cobertura; omite la prueba e2e
make test-e2e      # la etapa rechaza un report sintético y el certificado tiene 22 casillas
make lint          # ruff 0.15.12, reglas E4, E7, E9 y F
```

La integración continua está en `.github/workflows/ci.yml` (Ubuntu 24.04,
Python 3.14.3). En noviembre de 2026 la cuenta no asigna máquinas por un
bloqueo de facturación, así que esa corrida hay que mirarla en local.

## Cita

No hay DOI. Para el código de esta limpieza, `CITATION.cff` (versión 4.0.0).
Para los números de la tesis, la etiqueta `v3.7-final`.

Licencia MIT. Ver `LICENSE`.

## Sitio

`web/index.html` se abre en el navegador. Resume la pregunta, los tres
objetivos, el flujo, los métodos, los resultados, la robustez y las
descargas. No incluye fichas de proteínas ni la sección «Aprende»: no
aportan al encargo de la tesis y arrastraban la pregunta anterior. El
glosario es una página con las abreviaturas del capítulo 3.
