PY := .venv/bin/python
ONCOGO := .venv/bin/oncogo-crc
WORKSPACE ?= /Volumes/JANCODATA/OncoGO-CRC
RUN := $(ONCOGO) run --data-workspace $(WORKSPACE) --stage

.PHONY: install demo pipeline docs docs-build lint test
.PHONY: manifest-write manifest-verify session-info regresion-v1
.PHONY: v2-certified v2-axes v2-vias v2-deteccion v2-meta v2-certificado
.PHONY: v2-external v2-report verificar-numeros all-v2 serve-v2 v2-tcga v2-marcadores
.PHONY: v2-cohorte-cea v2-pureza v2-heterogeneidad v2-deconvolucion

install:
	python3 -m venv .venv && . .venv/bin/activate && pip install -e ".[dev]"

demo:
	$(ONCOGO) run --demo --stage all

pipeline:
	$(ONCOGO) run --stage all

docs:
	. .venv/bin/activate && pip install -q mkdocs mkdocs-material && mkdocs serve

docs-build:
	. .venv/bin/activate && pip install -q mkdocs mkdocs-material && mkdocs build

lint:
	.venv/bin/ruff check src scripts

test:
	$(PY) -m pytest -q

# ---------------------------------------------------------------- v1 congelada
# Los manifiestos congelan v1 y las entradas. `manifest-verify` falla si algo cambió.
manifest-write:
	$(PY) scripts/v2/manifest_hashes.py write

manifest-verify:
	$(PY) scripts/v2/manifest_hashes.py verify

session-info:
	$(PY) scripts/v2/session_info.py

# Re-ejecuta el contraste v1 en un workspace temporal. No escribe en results/certified.
regresion-v1:
	ONCOGO_REGRESION_V1=1 $(PY) -m pytest tests/test_regresion_v1.py -q

# ---------------------------------------------------------------------- v2
# Cada objetivo escribe solo bajo results/v2/. v1 no se toca.
# Tiempos medidos en el portátil del autor (M-series, 94 pares CPTAC).

# ~3 min: 10.000 permutaciones en CPTAC y las 65.536 exactas en PXD034575.
v2-certified:
	$(RUN) v2_certified

# ~2 min: pólipos con selección independiente, mitades divididas y nulo de intercambio.
v2-axes:
	$(RUN) v2_axes

# ~1 min: cláusula de vías con tres colecciones, prueba competitiva y nulo de etiquetas.
v2-vias:
	$(RUN) v2_vias

# ~1 min: ventana de detección en plasma con las concentraciones del HPA.
v2-deteccion:
	$(RUN) v2_deteccion

# ~1 min: meta-análisis REML por proteína entre las cohortes de tejido disponibles.
v2-meta:
	$(RUN) v2_meta

# D1: determinantes del CEA sérico en TCGA COAD/READ. Necesita la descarga
# scripts/download/tcga_coadread.py; sin ella la etapa se marca no_disponible.
v2-tcga:
	$(RUN) v2_tcga

# D2: el índice estromal de R4 contra la pureza ABSOLUTE medida. Necesita las
# descargas tcga_cptac_retro.py y tcga_coadread.py; sin las dos no calcula nada.
v2-pureza:
	$(RUN) v2_pureza

# D3: variacion dentro de un mismo tumor. La cohorte no cumple el minimo de
# pacientes del preregistro, asi que no entra en el meta-analisis: lo que aporta
# son ocho regiones del mismo tumor, que ninguna cohorte del certificado tiene.
v2-heterogeneidad:
	$(RUN) v2_heterogeneidad

# D5: deconvolucion con referencia unicelular externa, el metodo de sensibilidad
# que el preregistro llama deconvolucion_referencia_externa. Corre en las
# cohortes LFQ: en una matriz TMT de log-cocientes la deconvolucion no es valida,
# y la etapa lo comprueba con los datos en vez de afirmarlo.
v2-deconvolucion:
	$(RUN) v2_deconvolucion

# D8: tabla 2x2 de la literatura y meta-análisis bivariado de Reitsma.
# La extracción va antes porque la etapa lee su TSV.
v2-marcadores:
	$(PY) scripts/literatura/extraer_2x2.py
	$(RUN) v2_marcadores

# D9: protocolo de cohorte de CEA seriado. El estudio no se ha hecho, así que
# aquí no hay resultados que calcular: lo que se comprueba es que el análisis
# corre y que detecta una diferencia solo cuando existe. La cohorte es sintética
# y se queda en tests/, nunca en results/.
v2-cohorte-cea:
	$(PY) scripts/literatura/generar_sintetico_cohorte_cea.py
	.venv/bin/python -m pytest tests/test_cohorte_cea_d9.py -q

# Descargas idempotentes. Depende de la red; cada fuente queda marcada
# usada / no_disponible / requiere_acceso y nunca se rellena con datos inventados.
v2-external:
	$(PY) scripts/v2/download_all.py

# Certificado v2, índice R, changelog frente a v1 e informe de brecha.
# Va al final: lee las auditorías de las etapas anteriores.
v2-certificado:
	$(RUN) v2_certificado

# Falla si algún número de los informes no sale de una tabla de resumen.
verificar-numeros:
	$(PY) scripts/v2/verificar_numeros.py --workspace $(WORKSPACE)

# ~1 min: figuras a 300 dpi, HTML y PDF del informe de v2.
v2-report:
	$(RUN) v2_report

all-v2: manifest-write session-info regresion-v1 v2-external v2-certified v2-axes \
        v2-vias v2-deteccion v2-meta v2-tcga v2-pureza v2-heterogeneidad v2-deconvolucion v2-marcadores v2-certificado v2-report \
        verificar-numeros manifest-verify
	@echo "v2 completa. Informe en $(WORKSPACE)/results/v2/report/"

serve-v2:
	$(PY) -m http.server --directory $(WORKSPACE)/results/v2/report 8000
