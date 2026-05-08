.PHONY: install demo pipeline docs lint

install:
	python3 -m venv .venv && . .venv/bin/activate && pip install -e ".[dev]"

demo:
	. .venv/bin/activate && oncogo-crc run --demo --stage all

pipeline:
	. .venv/bin/activate && oncogo-crc run --stage all

docs:
	. .venv/bin/activate && pip install -q mkdocs mkdocs-material && mkdocs serve

docs-build:
	. .venv/bin/activate && pip install -q mkdocs mkdocs-material && mkdocs build

lint:
	. .venv/bin/activate && ruff check src
