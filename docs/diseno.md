# Diseño del flujo

Este documento sigue el orden del capítulo 3 de la tesis (secciones 3.1 a
3.1.6): requisitos, arquitectura, desarrollo, ejecución y el lugar de las
pruebas. Los nombres de archivo son los del repositorio. Los diagramas están
en `docs/diagramas/` y en `web/diagramas/`.

## Requisitos

Los requisitos se leen del código, de `config/preregistro_v3.yaml` y de las
pruebas. RF es lo que el flujo hace. RNF es con qué propiedades lo hace.

| Código | Requisito | Dónde se comprueba |
| --- | --- | --- |
| RF1 | Una orden lleva de las entradas declaradas al certificado | `make all-v3-7`; `oncogo-crc run --stage v3_7` |
| RF2 | Cada matriz se lee con las reglas de su plataforma y no se fusionan intensidades | `stages/v3_3_tejido.py`, `stages/v3_3_suero.py`, `v3/h8.py` |
| RF3 | Controles de calidad: valor q, especie, contaminantes, paneles, pares mínimos | `controles.tsv`, `controles_tejido.json`, `suero_calidad_muestra.tsv` |
| RF4 | Contraste de diferencia y de equivalencia, y AUC por proteína | `k4_*.tsv`, `k5_auc.json` |
| RF5 | Concordancia de tejido, integración tejido–suero y ontología genética | `triangulacion.tsv`, `priorizacion.tsv`, `go_ora_tejido.tsv`, `gsea_tejido.tsv` |
| RF6 | Orden de adquisición y modelos M1–M6 con validación cruzada, permutación y SHAP | `preanalitica.json`, `ml_suero_modelos.tsv`, `ml_explicacion.tsv` |
| RF7 | Repetir el análisis con decisiones alternativas | `multiverso_suero.tsv`, `multiverso_ml.tsv` |
| RF8 | Integridad de las entradas antes de calcular | `manifest_entradas.sha256`; SHA-1 de los crudos en `crudos_integridad.tsv` |
| RF9 | Certificado de la ejecución | `certificado.json`, 22 casillas, estado `publicable` o no |
| RNF1 | Mismo commit y mismas entradas, mismos números | `PYTHONHASHSEED=0`, un hilo, 35 métricas, segunda evaluación de M2 |
| RNF2 | Cada resultado se rastrea hasta el código, la configuración y las entradas | Sello por archivo, vigilante, árbol limpio, `procedencia.json` |
| RNF3 | La ejecución cabe en un computador personal | En el equipo del autor, la etapa tardó unos 15 minutos (11 de noviembre de 2026, 00:48 a 01:04 UTC) |
| RNF4 | Python 3.10 a 3.14 | `requires-python` en `pyproject.toml`. La referencia es 3.14.3. El contenedor fija esa versión |
| RNF5 | Umbrales y rutas sin tocar el código | `preregistro_v3.yaml`, `entradas_v3_7.yaml`, `esperado_v3_7.tsv`, variable `WORKSPACE` |
| RNF6 | Una versión no lee ni pisa otra | `results/v3_7` vacío al empezar, salvo el registro de crudos. `tests/test_v3_7_sin_mezcla.py` |
| RNF7 | Código revisable y con pruebas | `make lint`, `make test` |

Los umbrales no se eligen al ejecutar. El margen de equivalencia es ±0,585 en
log2 (1,5 veces). BH usa alfa 0,05. El valor q de corte del suero es 0,01.
La semilla base es 20260914, con un desplazamiento fijo por pieza. El GSEA
usa la semilla 42. Todo eso está en el preregistro y en sus enmiendas.

Hay tres actores. Quien investiga registra el preregistro, prueba los módulos
y ejecuta. Quien reproduce parte de la etiqueta, prepara el entorno, descarga
los datos y compara. Quien evalúa lee el certificado, los sellos y el
manifiesto sin ejecutar. El diagrama de casos de uso es
`docs/diagramas/casos-de-uso.svg`.

## Arquitectura

Cuatro cohortes entran por separado. Una etapa común de preparación escribe
el manifiesto y el sello. Después corren tejido, suero, controles, modelos,
multiverso, suero por estadio, preanalítica, triangulación, ontología
genética, explicación e integración. Las flechas entre cohortes comparan
resultados: no suman intensidades.

`make all-v3-7` encadena tres pasos: la regresión de la corrida certificada
v1, la comprobación de los crudos de PXD055218 y `stages/v3_7.py`.

El orquestador está en `stages/v3_7.py`. Reutiliza etapas de `stages/` y la
biblioteca de `v3/`. Los módulos transversales son el preregistro, la
escritura con sello, el vigilante, el lector de DIA-NN y la t moderada.
Ningún módulo de `v3/` importa una etapa.

| Orden | Módulo | Entradas principales | Salidas en `results/v3_7` |
| --- | --- | --- | --- |
| Preparación | `Makefile`, `stages/v3_7.py` | Entradas declaradas, preregistro, crudos | `manifest_entradas.sha256`, `sello_corrida.json`, `corrida.json` |
| Tejido | `stages/v3_3_tejido.py` | `proteinGroups` de Tanaka, contraste certificado de CPTAC, UniProt | `k1_*.tsv`, `k1_etiquetas.json`, `k3_deming.json` |
| Suero | `stages/v3_3_suero.py`, `v3/suero_v32.py` | Informe DIA-NN, biblioteca de contaminantes, paneles | matrices de suero, `k4_*.tsv`, `k5_auc.json` |
| Controles | `stages/v3_3_controles.py`, `stages/v3_4.py` | Matrices de tejido y de suero | `controles.tsv`, `controles_tejido.json` |
| Modelos | `stages/v3_3_ml.py`, `v3/ml_suero.py` | Matriz de suero y variables de calidad | `ml_suero_modelos.tsv`, `ml_predicciones_oof.tsv` |
| Multiverso | `stages/v3_4.py` | Matrices de suero y de tejido | `multiverso_suero.tsv`, `multiverso_ml.tsv` |
| Suero por estadio | `v3/h8.py` | `report.pr_matrix.tsv` de PXD064655 | `h8_carga_tumoral.json` |
| Preanalítica | `stages/v3_5_preanalitica.py` | Fechas de los 60 crudos y variables de calidad | `preanalitica.json` |
| Triangulación | `v3/triangulacion.py` | Resultados ya calculados | `triangulacion.tsv` |
| Ontología | `v3/cierre_v37.py` | Contrastes de tejido y GO | `go_ora_tejido.tsv`, `gsea_tejido.tsv` |
| Explicación | `v3/cierre_v37.py` | Modelos M2 y M3 | `ml_explicacion.tsv` |
| Integración | `v3/cierre_v37.py` | CPTAC, Tanaka, HPA y matriz de suero | `priorizacion.tsv` |
| Cierre | `stages/v3_7.py`, `v3/cierre_v37.py` | Salidas y `esperado_v3_7.tsv` | `determinismo.json`, `tests.txt`, `certificado.json` |

Los almacenes del flujo de datos son el workspace (matrices y anotaciones),
la configuración versionada y el contraste certificado v1. Los procesos
escriben en `results/v3_7`, cada archivo con su sello.

## Desarrollo y ejecución

El historial llega a `v3.7-final` el 10 de noviembre de 2026. La corrida
certificada v1 es del 1 de octubre de 2026. La versión 4.0.0 no mueve esas
etiquetas.

Para repetir desde cero: clonar, situarse en el código, instalar, descargar
las cuatro matrices y el depósito
[10.5281/zenodo.24398807](https://doi.org/10.5281/zenodo.24398807), comprobar
las huellas, ejecutar con `WORKSPACE` y comparar. Los cinco pasos y las excepciones admitidas (fecha, identificador
de proceso, tiempos de las pruebas y el orden de `gsea_tejido.tsv`) están
en el README y en `scripts/comparar_con_referencia.py`.

La etapa no empieza si `results/v3_7` del workspace tiene restos, y no firma
si el árbol de git está sucio. Después de las piezas vuelve a evaluar M2,
cierra el vigilante, rechaza una salida con ruta local, corre las pruebas y
emite el certificado.

## Despliegue

No hay servicio que desplegar. La ejecución es un proceso local o un
contenedor (`Dockerfile`: Python 3.14.3, `PYTHONHASHSEED=0`, un hilo). El
sitio de `web/` es estático y no se publica desde este repositorio.
