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

## E6 — Control negativo con covariables permutadas en C7

- **Fecha**: 2026-10-03
- **Tipo**: control añadido después de ver el primer resultado de v2
- **Qué pasó**: el ajuste conjunto por los cuatro compartimentos (matriz
  extracelular, músculo liso, inmune, epitelial) deja 83 de 688 proteínas, es
  decir el 12,1 %, frente al 44,6 % de v1 con un solo índice. Antes de escribir
  que eso es confusión por composición hay que descartar la explicación trivial.
- **Por qué hacía falta**: ajustar por k covariables gasta k grados de libertad
  en 94 pares e infla el error estándar. Una parte de la caída es pérdida de
  potencia, no confusión. Sin el control, el número no se puede interpretar.
- **Qué se añade**: 200 repeticiones en las que el orden de los pares de cada
  covariable se permuta. Eso conserva la distribución de la covariable y rompe su
  relación con el contraste real del par. La fracción que sobrevive en el control
  es el techo que la pérdida de potencia explica por sí sola; la diferencia con
  el ajuste real es la parte atribuible a composición.
- **Semilla**: desplazamiento 416 sobre la base 20260914, añadido al preregistro.
- **Nota de honestidad**: este control no estaba preregistrado. Se añade porque
  sin él el 12,1 % se leería como un resultado biológico cuando puede ser
  aritmética del modelo. No cambia ningún umbral de decisión.

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

## E7 — El nulo de C5 permuta el brazo por proteína y vuelve a seleccionar

- **Fecha**: 2026-10-03
- **Tipo**: corrección de un nulo mal especificado en el preregistro
- **Qué decía el preregistro**: `c5_polipos.nulo_intercambio` fijaba «permuta las
  etiquetas pólipo y tumor dentro de cada paciente», con unidad paciente.
- **Qué está mal**: ese nulo destruye lo que quiere medir. Al intercambiar las
  etiquetas paciente a paciente, el brazo que hace de tumor queda siendo una
  mezcla de columnas de tumor y de pólipo, así que su media se encoge hacia cero.
  Las proteínas se habían elegido por tener media grande en el tumor, y esa
  ventaja desaparece en el nulo. Resultado: el nulo queda por debajo de lo
  observado incluso cuando los dos brazos son ruido independiente, y el p sale
  significativo sin que haya nada. Se comprobó en
  `tests/test_polipos_c5.py`: con dos matrices de ruido puro, el nulo del
  preregistro daba p = 0,005.
- **Qué se hace en su lugar**: para cada proteína se sortea cuál de los dos
  brazos hace de tumor y después se **vuelve a aplicar la regla de selección**
  sobre ese brazo. El nulo repite todo el procedimiento, selección incluida, que
  es la única forma de que mida el sesgo de selección. Con ruido puro el nulo
  iguala lo observado, como debe ser; con un pólipo de verdad intermedio, no lo
  alcanza.
- **Qué no cambia**: el número de permutaciones (1000), la semilla
  (desplazamiento 206 sobre 20260914) y la regla de decisión de C5.
- **Nota de honestidad**: el preregistro se escribió antes de mirar los datos y
  por eso contiene este error. No se edita el fichero preregistrado; la
  desviación se declara aquí.

## E8 — El umbral que convierte el índice estromal en un sustituto de la pureza

- **Fecha**: 2026-10-03
- **Tipo**: umbral nuevo, fijado después de que la fuente existiera
- **Qué faltaba en el preregistro**: `c7_composicion.metodos_sensibilidad` ya
  incluía `pureza_absolute_tcga`, y el lenguaje obligatorio ya decía que R4 es
  «una cota del efecto de composición, no una pureza medida, salvo cuando la
  covariable es la pureza ABSOLUTE de TCGA». Pero no fijaba **cuánto** tenía que
  parecerse el índice a la pureza para poder ejercer ese «salvo».
- **Qué se fija**: una correlación de rangos de 0,7 en valor absoluto y con signo
  negativo. Por encima, el índice se puede presentar como sustituto de la pureza.
  Por debajo, R4 sigue siendo un análisis de sensibilidad legítimo pero no un
  control de pureza.
- **Por qué 0,7**: es el valor al que el índice y la pureza comparten la mitad de
  la varianza (0,7² ≈ 0,49). Por debajo de ahí, más de la mitad de la variación
  en pureza queda fuera del índice, y llamarlo «control de pureza» describiría
  mal lo que hace.
- **Nota de honestidad**: el umbral se fija con el dato ya calculado delante
  (rho = −0,228), así que no tiene el valor de un umbral preregistrado y se
  declara así en el informe. Con cualquier umbral razonable el veredicto sería el
  mismo: 0,228 queda lejos de 0,7 y también de 0,5 y de 0,4. Lo que el umbral
  decide no es el resultado, es la frase con que se describe.
- **Qué no cambia**: ningún subíndice del certificado. D2 no recalcula R4, porque
  la cohorte que tiene pureza medida no tiene tejido normal pareado.

## E9 — Los umbrales de la deconvolución con referencia externa

- **Fecha**: 2026-10-03
- **Tipo**: reglas y umbrales nuevos, fijados al implementar un método que el
  preregistro listaba sin especificar
- **Qué faltaba en el preregistro**:
  `c7_composicion.metodos_sensibilidad` incluía
  `deconvolucion_referencia_externa`, y `semillas.desplazamientos` reservaba el
  415 para él. Pero no decía de dónde salía la referencia, cómo se elegían los
  genes marcadores ni cuándo había que declarar que el resultado no era
  utilizable. Eran cuatro decisiones sin escribir.
- **Qué se fija**:
  1. **La referencia**: los perfiles por tipo celular del Human Protein Atlas
     (`rna_single_cell_type.tsv`), restringidos a once tipos del colon que se
     listan en `scripts/download/hpa_single_cell.py` antes de mirar ningún
     resultado. La referencia específica de colon (GSE178341, Pelka 2021) queda
     declarada y sin usar.
  2. **La regla de marcadores**: un gen entra si su tipo con más expresión llega
     a 10 nCPM y es al menos 5 veces la media de los demás tipos. Un gen
     expresado igual en todos los tipos no aporta información y sí aporta ruido.
  3. **Los umbrales del diagnóstico**: residuo relativo mediano por encima de
     0,5 significa que la mezcla no explica el perfil; correlación de 0,8 entre
     dos columnas de la referencia significa que esos dos tipos no son
     separables; un tipo ausente (por debajo de 0,005) en más de la mitad de las
     muestras significa que el método lo está descartando, no estimándolo en
     cero.
  4. **Cuándo se usa el eje epitelial como covariable**: solo si baja en el
     tumor en al menos el 75 % de los pares y coincide en sentido con el índice
     estromal. Son dos indicios de que el orden entre muestras lleva señal, que
     es lo que un ajuste necesita.
- **Nota de honestidad**: los umbrales del punto 3 se fijaron **después** de ver
  que el residuo salía alrededor de 0,9 y que varios tipos salían en cero. Es
  decir, se fijaron sabiendo que el diagnóstico iba a fallar. En descargo: la
  primera versión del diagnóstico comprobaba solo el máximo por tipo y declaraba
  las proporciones publicables, y al mirar las tablas se vio que eso era falso.
  El cambio no eligió un resultado: corrigió una comprobación que no detectaba lo
  que decía detectar. El umbral del punto 4 se fijó antes de calcular la
  dirección en cada cohorte, y de hecho una de las dos cohortes no lo pasa.
- **Qué no cambia**: ningún subíndice del certificado. El R4 del certificado
  sigue siendo el de CPTAC con el índice estromal. Esta enmienda solo gobierna un
  análisis de sensibilidad que, además, no se puede correr en CPTAC.
