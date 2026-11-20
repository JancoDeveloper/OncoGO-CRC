PY := .venv/bin/python
ONCOGO := .venv/bin/oncogo-crc
# Espacio de escritura. No es el disco de datos.
WORKSPACE ?= $(CURDIR)/workspace
DATA_DIR ?=
SALIDAS ?= /Volumes/JANCODATA/OncoGO-CRC-v4-salidas
export PYTHONHASHSEED := 0
export OMP_NUM_THREADS := 1
export OPENBLAS_NUM_THREADS := 1
export MKL_NUM_THREADS := 1
export ONCOGO_CRC_DATA_WORKSPACE := $(WORKSPACE)
export GIT_CONFIG_COUNT := 1
export GIT_CONFIG_KEY_0 := core.filemode
export GIT_CONFIG_VALUE_0 := false

.PHONY: install lint test test-e2e datos descargar verificar-datos certified regresion-v1 crudos-v3-7 v3-7 all-v3-7 verificar verificar-local hpa

install:
	python3 -m venv .venv && .venv/bin/pip install -U pip && .venv/bin/pip install -e ".[dev]"

lint:
	.venv/bin/ruff check src tests scripts

test:
	$(PY) -m pytest -q --cov=oncogo_crc --cov-report=term-missing --cov-report=xml

test-e2e:
	ONCOGO_E2E=1 $(PY) -m pytest -q -m e2e --override-ini="addopts=" tests/test_e2e_v3_7_sintetico.py

# PRIDE, contaminantes y el depósito Zenodo 10.5281/zenodo.24398807.
datos:
	$(PY) scripts/descargar_datos.py --workspace $(WORKSPACE)

descargar:
	$(PY) scripts/descargar_datos.py --workspace $(WORKSPACE)

verificar-datos:
	$(PY) scripts/verificar_datos.py --workspace $(WORKSPACE)

hpa:
	$(PY) scripts/download/human_protein_atlas.py --workspace $(WORKSPACE)

# Regenera results/certified y results/certification a partir de las matrices del workspace.
certified:
	$(ONCOGO) run --data-workspace $(WORKSPACE) --stage certified_inference
	$(ONCOGO) run --data-workspace $(WORKSPACE) --stage certification

regresion-v1:
	ONCOGO_CRC_DATA_WORKSPACE=$(WORKSPACE) ONCOGO_REGRESION_V1=1 $(PY) -m pytest tests/test_regresion_v1.py -q --override-ini="addopts="

crudos-v3-7:
	$(PY) scripts/download/pxd055218_crudos.py --workspace $(WORKSPACE)

v3-7:
	$(ONCOGO) run --data-workspace $(WORKSPACE) $(if $(DATA_DIR),--data-dir $(DATA_DIR)) --stage v3_7

all-v3-7: regresion-v1 crudos-v3-7 v3-7

verificar:
	$(PY) scripts/comparar_con_referencia.py --nueva $(WORKSPACE)/results/v3_7 --referencia results/v3_7

verificar-local:
	ONCOGO_DATA_DIR="$(DATA_DIR)" ONCOGO_CRC_DATA_WORKSPACE="$(SALIDAS)" TMPDIR="$${TMPDIR:-/Volumes/JANCODATA/tmp}" \
		$(PY) scripts/verificar_local.py --datos "$(DATA_DIR)" --salidas "$(SALIDAS)"
