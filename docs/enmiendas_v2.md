# Enmiendas de v2

Toda desviación respecto al diseño de v1 o al preregistro de v2, con fecha y
justificación. Una enmienda escrita después de ver los datos de v1 **no** tiene
el mismo valor que un umbral preregistrado, y así se declara en el reporte.

Preregistro de v2: `config/preregistro_v2.yaml`, commiteado el 2026-10-03 antes
de cualquier análisis de v2.

---

## E1 — El mínimo de términos de la cláusula de vía se fija en 10

- **Fecha**: 2026-10-03
- **Tipo**: enmienda posterior a ver los datos de v1
- **Qué cambia**: v1 declaró `sostiene_clausula_via` con 4 términos evaluables y
  un intervalo de 0,398 a 1,0. El preregistro de v2 fija `n_min_terminos = 10` y
  añade el estado `no_evaluable` cuando no se alcanza.
- **Por qué**: con 4 términos el intervalo cubre casi todo el rango posible. Un
  veredicto de «sostiene» con esa incertidumbre no distingue el resultado del
  azar. El umbral de 10 es el mismo mínimo que el pipeline ya usaba para el
  número de proteínas por término (`MIN_TERM_GENES = 10`), así que no es una cifra
  nueva: es la aplicación del mismo criterio al nivel del conjunto de términos.
- **Consecuencia**: la cláusula de vía de v1 pasa de `sostiene` a `no_evaluable`
  si el recálculo no llega a 10 términos. El veredicto de H1 no cambia, porque la
  cláusula de proteína ya rechaza H1 por sí sola.
- **Lo que se conserva**: el valor de v1 (0,80 como umbral, 1,0 con n = 4) queda
  en `results/v2/changelog_numeros.tsv` junto al de v2.

## E2 — Los paquetes de R se reimplementan en Python

- **Fecha**: 2026-10-03
- **Tipo**: adaptación de implementación, no cambia ninguna regla de decisión
- **Qué cambia**: el diseño de v2 nombra `limma::camera`, `metafor::rma` y
  `mada::reitsma`. No hay R instalado en la máquina ni código R en el repositorio
  (verificado: `which R` y `which Rscript` no lo encuentran).
- **Por qué**: `treat` ya estaba reimplementado en Python en
  `src/oncogo_crc/utils/moderated.py` siguiendo Smyth 2004 y McCarthy y Smyth
  2009, y es lo que produjo los números certificados de v1. Instalar R y
  Bioconductor añadiría unos 2 GB y una dependencia que sin `renv` no se puede
  congelar, para después tener dos implementaciones del mismo estadístico. El
  camino consistente es una sola implementación, en Python, con tests.
- **Cómo se controla**: cada método reimplementado lleva un test que lo compara
  contra un caso con resultado conocido en forma cerrada o publicado:
  - REML y DerSimonian–Laird: `tests/test_meta_reml.py`.
  - Prueba competitiva de conjuntos: `tests/test_set_tests.py`.
- **Lo que no se puede garantizar**: la equivalencia numérica al último dígito
  con las implementaciones de R. Se declara en el reporte y en
  `results/v2/session_info.txt`, que registra `R: no_instalado`.

## E3 — Se añade la corrección C15: el denominador de R1

- **Fecha**: 2026-10-03
- **Tipo**: corrección nueva, no estaba en el diseño de v2
- **Qué cambia**: R1 de v1 es 104/119. Las 119 proteínas evaluables son el 17,3 %
  de las 688 certificadas. El certificado presenta 0,874 sin ese denominador, lo
  que se lee como si la concordancia fuera una propiedad de los 688 hallazgos.
- **Por qué**: el denominador de R1 no son los hallazgos, son los hallazgos que
  la réplica además alcanza a medir con log2FC finito y distinto de cero. Eso
  condiciona la concordancia a la detectabilidad de la segunda plataforma, que es
  precisamente lo que R2 mide por separado. Reportar R1 sin la cobertura mezcla
  las dos cosas.
- **Qué se hace**: el preregistro añade la sección `c15_denominador_r1`. R1 se
  reporta siempre con `n_certificadas`, `n_evaluables`, `n_mismo_signo` y
  `cobertura_evaluable`. Si la cobertura es menor que 0,50, R1 se nombra
  «concordancia condicionada a detectabilidad en la réplica».
- **Relación con C3 y C9**: C3 arregla los denominadores del certificado y C9
  estratifica R1. C15 es lo que obliga a declarar la cobertura en el número
  principal.

## E4 — `results/v2/` fijo en lugar de `results/v2_<AAAAMMDD>/`

- **Fecha**: 2026-10-03
- **Tipo**: decisión de infraestructura
- **Qué cambia**: el diseño ofrece las dos rutas. Se usa `results/v2/`.
- **Por qué**: un nombre con fecha obliga a editar el `Makefile` en cada corrida
  y rompe `make all-v2` como comando único reproducible. La trazabilidad temporal
  la da `results/v2/session_info.txt`, que registra fecha, commit y versiones.

## E5 — La protección de v1 es el manifiesto, no el sistema de ficheros

- **Fecha**: 2026-10-03
- **Tipo**: limitación del entorno
- **Qué cambia**: el disco externo JANCODATA es ExFAT. No guarda permisos POSIX
  ni enlaces simbólicos, y crea un fichero `._nombre` de AppleDouble por cada
  fichero.
- **Consecuencia**: no se puede marcar v1 como solo lectura con `chmod`. La
  protección efectiva son tres cosas: `results/v1_manifest.sha256`, que v2 escribe
  únicamente bajo `results/v2/`, y el test de regresión que falla si los números
  de v1 cambian. Los ficheros `._*` quedan excluidos de los manifiestos.
