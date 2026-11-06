PY := .venv/bin/python
ONCOGO := .venv/bin/oncogo-crc
WORKSPACE ?= /Volumes/JANCODATA/OncoGO-CRC
RUN := $(ONCOGO) run --data-workspace $(WORKSPACE) --stage

.PHONY: install demo pipeline docs docs-build lint test
.PHONY: manifest-write manifest-verify session-info regresion-v1
.PHONY: v2-certified v2-axes v2-vias v2-deteccion v2-meta v2-certificado
.PHONY: v2-external v2-report verificar-numeros all-v2 serve-v2 v2-tcga v2-marcadores
.PHONY: v2-cohorte-cea v2-pureza v2-heterogeneidad v2-deconvolucion tesis-cierre

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

# PDF explicativo para lectores sin formación técnica + copia de reportes a docs/reportes/.
v2-guia-publica:
	$(PY) scripts/build_public_report_v2.py --workspace $(WORKSPACE) --sync-reportes

all-v2: manifest-write session-info regresion-v1 v2-external v2-certified v2-axes \
        v2-vias v2-deteccion v2-meta v2-tcga v2-pureza v2-heterogeneidad v2-deconvolucion v2-marcadores v2-certificado v2-report \
        verificar-numeros manifest-verify
	@echo "v2 completa. Informe en $(WORKSPACE)/results/v2/report/"

serve-v2:
	$(PY) -m http.server --directory $(WORKSPACE)/results/v2/report 8000

# Tres grupos de suero, auditoría de candidatas y límite de neoantígenos.
# No entra en `all` ni en `all-v2`: no reescribe el certificado.
tesis-cierre:
	$(RUN) tesis_cierre

# ---------------------------------------------------------------------- v3
# Escribe solo bajo results/v3/. v1, v2 y el cierre no se tocan.
.PHONY: v3-manifiestos v3-citas v3-correcciones all-v3

v3-manifiestos:
	$(PY) scripts/v3/congelar_manifiestos.py

v3-citas:
	$(PY) scripts/v3/verificar_citas.py

v3-correcciones:
	$(RUN) v3_correcciones

all-v3: manifest-write v3-manifiestos regresion-v1 v3-citas v3-correcciones
	@echo "v3 parcial. Resultados en $(WORKSPACE)/results/v3/"

# --------------------------------------------------------------------- v3.2
# Escribe solo bajo results/v3_2/. v3 queda congelado por results/v3_manifest.sha256.
.PHONY: v3-2-suero v3-2-tejido v3-2-controles v3-2-multiverso v3-2-hipotesis v3-2-ml v3-2-report v3-2-docs v3-2-tesis all-v3-2

v3-2-suero:
	$(RUN) v3_2_suero

v3-2-tejido:
	$(RUN) v3_2_tejido

v3-2-controles:
	$(RUN) v3_2_controles

v3-2-multiverso:
	$(RUN) v3_2_multiverso

v3-2-hipotesis:
	$(RUN) v3_2_hipotesis

v3-2-ml:
	$(RUN) v3_2_ml

v3-2-report:
	$(RUN) v3_2_report

v3-2-docs:
	@echo "v3.2 docs: los textos leen results/v3_2. No hay una etapa de cálculo."

v3-2-tesis:
	@echo "v3.2 tesis: LaTeX y Word se regeneran desde las tablas. No hay una etapa de cálculo."

all-v3-2: manifest-write v3-manifiestos regresion-v1 v3-2-suero v3-2-tejido v3-2-controles v3-2-hipotesis v3-2-multiverso v3-2-ml v3-2-report
