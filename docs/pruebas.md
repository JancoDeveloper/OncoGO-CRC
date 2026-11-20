# Pruebas

Niveles, de más pequeño a más grande.

| Nivel | Qué cubre | Cómo se corre |
| --- | --- | --- |
| Unidad | Contrastes, equivalencia, sellos, vigilante, SHAP lineal, directLFQ, multiverso | `make test` |
| Regresión v1 | Los siete números certificados. Se omite si el workspace no tiene `results/certified` | `make regresion-v1` |
| Extremo a extremo sintético | La etapa rechaza un informe que no es el DIA-NN sellado. El emisor del certificado devuelve 22 casillas | `make test-e2e` |
| Corrida completa | `make all-v3-7` y `make verificar` | Pide las entradas de `data/checksums.sha256` |

`make test` informa cobertura con `pytest-cov`. No hay un umbral que falle
la corrida: el informe es para ver qué quedó sin ejercer.

La etapa, al terminar, lanza `pytest` y guarda la salida en `tests.txt`.
Esa segunda pasada no ejecuta la prueba e2e (`ONCOGO_E2E` no está definido).

ruff está fijado en 0.15.12. Las reglas activas son E4, E7, E9 y F.

La integración continua (`.github/workflows/ci.yml`) instala el paquete en
Ubuntu 24.04 con Python 3.14.3 y corre ruff, pytest y la e2e. En noviembre
de 2026 GitHub no asigna máquina a esta cuenta por un bloqueo de
facturación. La validación de esta versión se hizo en local con CPython
3.14.3.
