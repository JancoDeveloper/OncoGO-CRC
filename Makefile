PY := .venv/bin/python
ONCOGO := .venv/bin/oncogo-crc
WORKSPACE ?= /Volumes/JANCODATA/OncoGO-CRC
RUN := $(ONCOGO) run --data-workspace $(WORKSPACE) --stage

.PHONY: install demo pipeline docs docs-build lint test
.PHONY: manifest-write manifest-verify session-info regresion-v1
.PHONY: v2-certified v2-axes v2-vias v2-deteccion v2-meta v2-certificado
.PHONY: v2-external v2-report verificar-numeros all-v2 serve-v2

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
        v2-vias v2-deteccion v2-meta v2-certificado v2-report verificar-numeros manifest-verify
	@echo "v2 completa. Informe en $(WORKSPACE)/results/v2/report/"

serve-v2:
	$(PY) -m http.server --directory $(WORKSPACE)/results/v2/report 8000
