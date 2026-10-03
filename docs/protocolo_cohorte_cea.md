# Protocolo: ¿sirve de más mirar cómo cambia el CEA que mirar cuánto vale?

**Estado: pendiente de aprobación ética y de datos reales.** Este estudio no se
ha hecho. Lo que hay aquí es el protocolo, el diccionario de datos, el plan de
análisis y el código, probado sobre una cohorte sintética. No hay ni un dato de
ninguna persona en este documento ni en el código.

Ningún resultado de esta tesis sirve para decidir sobre un paciente.

---

## 1. De dónde sale la pregunta

El meta-análisis de D8 dejó un número que llama la atención. Para detectar una
recaída de cáncer colorrectal después de la cirugía:

| Prueba | Cómo se mide | Sensibilidad | Especificidad |
| --- | --- | --- | --- |
| CEA | un valor contra el umbral de 5 ng/ml | 0,605 (0,492–0,708) | 0,879 (0,771–0,939) |
| ctDNA | una sola medición tras la cirugía | 0,51 (0,41–0,61) | 0,91 (0,89–0,93) |
| ctDNA | mediciones repetidas en el seguimiento | 0,84 (0,73–0,91) | 0,95 (0,80–0,99) |

El ctDNA medido una vez no es mejor que el CEA. Medido varias veces, sí: la
sensibilidad pasa del 51 % al 84 % con la misma prueba y la misma población. Lo
que cambia no es la molécula, es que se mira la trayectoria en lugar del valor.

Y al CEA nunca se le ha hecho esa pregunta con el mismo cuidado. En los ocho
estudios de D8 el CEA se evalúa como se usa en la consulta: cada valor contra un
umbral fijo de 5 ng/ml, igual para todo el mundo. Pero el CEA basal de una
persona depende de si fuma, de su función hepática y de su propia biología, así
que un umbral común compara a cada paciente con una referencia que no es la suya.

**Hipótesis.** En una cohorte de seguimiento con CEA seriado, un modelo que use
cómo ha cambiado el CEA de cada paciente respecto a sus propias mediciones
anteriores detecta la recaída mejor que el umbral de 5 ng/ml.

## 2. Diseño

Cohorte prospectiva observacional. No se cambia el manejo de nadie: el CEA
seriado ya forma parte del seguimiento habitual tras la resección de un cáncer
colorrectal. El estudio solo analiza lo que ya se mide.

- **Población**: pacientes con cáncer colorrectal en estadio II o III sometidos a
  resección con intención curativa.
- **Exclusión**: enfermedad metastásica al diagnóstico, segunda neoplasia activa,
  y pacientes con menos de dos determinaciones de CEA en el seguimiento, porque
  con una sola no hay trayectoria que medir.
- **Seguimiento**: 36 meses, con CEA en los meses 3, 6, 9, 12, 18, 24, 30 y 36,
  que es la pauta habitual.
- **Variable de resultado**: recaída confirmada por imagen o histología. La fecha
  es la de la confirmación, no la del primer CEA alterado: usar la segunda haría
  que el CEA predijera su propia definición.

## 3. Diccionario de datos

Una fila por visita. El fichero es `paciente,visita_meses,cea,estadio,fumador,recaida,mes_de_la_recaida,ultimo_seguimiento_meses`.

| Columna | Tipo | Unidad | Obligatoria | Descripción |
| --- | --- | --- | --- | --- |
| `paciente` | texto | — | sí | Identificador seudonimizado. Nunca un número de historia. |
| `visita_meses` | decimal | meses | sí | Meses desde la cirugía hasta esta determinación. |
| `cea` | decimal | ng/ml | sí | CEA sérico de esta visita. Mayor o igual que cero. |
| `estadio` | texto | — | no | II o III, del informe patológico. |
| `fumador` | entero | 0 o 1 | no | Fumador activo en el momento de la cirugía. Entra porque eleva el CEA basal sin que haya cáncer. |
| `recaida` | entero | 0 o 1 | sí | Si el paciente recayó en algún momento del seguimiento. Igual en todas sus filas. |
| `mes_de_la_recaida` | decimal | meses | sí | Meses desde la cirugía hasta la confirmación. Vacío si no recayó. |
| `ultimo_seguimiento_meses` | decimal | meses | sí | Último mes con información. Es lo que permite distinguir «no recayó» de «no se sabe». |

**Reglas de coherencia** que el código comprueba antes de calcular nada:

- Ningún CEA negativo y ninguna visita con tiempo negativo.
- Un paciente con `recaida = 1` tiene `mes_de_la_recaida`. Uno con `recaida = 0`
  no lo tiene. Los dos son los errores de captura más frecuentes en una base de
  seguimiento, y los dos producen análisis de aspecto normal.
- Al menos dos visitas por paciente.

**Qué no se recoge.** Ni nombre, ni fecha de nacimiento, ni número de historia,
ni ninguna fecha de calendario: todos los tiempos son meses desde la cirugía.
Con eso basta para el análisis y no hace falta nada más.

## 4. Plan de análisis estadístico

Se declara antes de ver los datos. Cualquier cambio posterior va a
`docs/enmiendas_v2.md` como enmienda, con su fecha y su motivo.

### 4.1 Variables de cada modelo

Se comparan tres modelos sobre las **mismas** visitas:

1. **Umbral**: `supera_el_umbral`, es decir si el CEA de esta visita pasa de
   5 ng/ml. Es lo que se hace hoy.
2. **Valor**: `log2_cea`, el valor como número y sin umbral.
3. **Trayectoria**: `log2_cea`, más `cambio_anual` (duplicaciones por año frente
   a la visita anterior), `pendiente` (recta por mínimos cuadrados sobre las tres
   últimas visitas) y `cambio_sobre_el_basal` (frente a la primera visita de ese
   mismo paciente).

El logaritmo en base 2 está ahí porque el CEA sube de forma multiplicativa: una
subida de 2 a 4 ng/ml y otra de 20 a 40 son el mismo fenómeno y deben contar
igual. El `cambio_anual` se divide por el tiempo transcurrido para que dos
visitas separadas seis meses y dos separadas tres no den el mismo número por el
mero hecho de estar más juntas.

### 4.2 Variable de resultado y censura

Para cada visita, el resultado es si el paciente recae en los **12 meses**
siguientes. Doce meses es el horizonte en el que una detección cambia el manejo.

Una visita que no llega a tener doce meses de seguimiento y en la que no ha
pasado nada **no entra en el análisis**. No se sabe si el paciente recayó
después, y meterla como «no recayó» inventaría información. Eso es censura, y la
forma de tratarla es excluir la visita, no adivinarla. El número de visitas
excluidas por este motivo se publica.

### 4.3 Comparación

- El AUC de cada modelo se calcula **fuera de muestra**, con validación cruzada
  de cinco pliegues **partidos por paciente**. Partir por visita dejaría visitas
  del mismo paciente a los dos lados y el modelo vería en entrenamiento a la
  misma persona que luego evalúa: el AUC saldría más alto y la causa sería la
  fuga, no el modelo.
- La comparación es de **diferencias pareadas**, no de intervalos que se solapan.
  En cada réplica del bootstrap se remuestrean **pacientes**, se recalculan los
  tres AUC sobre los mismos pacientes y se guarda la diferencia. Dos intervalos
  que se solapan pueden corresponder a una diferencia que no incluye el cero, así
  que lo que decide es el intervalo de la diferencia.
- **La unidad de remuestreo es el paciente y no la visita.** Un paciente aporta
  varias visitas correlacionadas entre sí; remuestrear visitas daría intervalos
  más estrechos de lo que corresponde.
- 2 000 réplicas. La semilla es 20260914 + 909, derivada de la semilla base del
  preregistro.

### 4.4 Regla de decisión

La hipótesis se considera apoyada si el intervalo del 95 % de
`AUC(trayectoria) − AUC(umbral)` queda **enteramente por encima de cero**. Se
declara además la comparación de la trayectoria contra el valor solo, que es la
que dice si la trayectoria aporta algo más que el número.

### 4.5 Tamaño necesario

Con una proporción de visitas con recaída del 11 %, cinco visitas por paciente,
potencia del 80 % y alfa del 5 %, y suponiendo una correlación entre visitas del
mismo paciente de 0,5:

| Diferencia de AUC a detectar | Pacientes |
| --- | --- |
| 0,10 | 330 |
| 0,05 | 1340 |

Es una aproximación: la varianza del AUC es la de Hanley y McNeil corregida por
el efecto de diseño del conglomerado. El cálculo exacto exigiría simular. Lo que
hace falta aquí es el orden de magnitud para poder pedir el estudio, no la cifra
al paciente.

## 5. Lo que este estudio no podría concluir

- **No diría si conviene hacer más CEA.** Detectar antes una recaída no es lo
  mismo que vivir más. Eso es un ensayo con supervivencia como resultado, no una
  cohorte de precisión diagnóstica.
- **No daría un umbral para usar en consulta.** Un AUC es una medida global; el
  punto de corte depende de qué cuesta un falso positivo y qué cuesta un falso
  negativo, y eso no sale de estos datos.
- **No valdría para quien tenga pocas determinaciones.** El modelo de trayectoria
  necesita al menos dos, y con dos la pendiente es simplemente una diferencia.
- **No se traslada a otra población sin comprobarlo.** El CEA basal depende del
  tabaquismo y de la función hepática, y la proporción de fumadores varía mucho
  entre cohortes.

## 6. Ética y datos

- El estudio **requiere aprobación de un comité de ética** antes de recoger nada.
  Mientras no la tenga, el estado de este protocolo es `pendiente`.
- Los datos serían seudonimizados en origen. No se registra ninguna fecha de
  calendario: todos los tiempos son meses desde la cirugía.
- No hay intervención. El CEA seriado ya forma parte del seguimiento habitual.
- **Ningún dato clínico de ninguna persona entra en este repositorio.** La
  cohorte de `tests/SINTETICO_cohorte_cea.csv` es generada por
  `scripts/literatura/generar_sintetico_cohorte_cea.py` con una semilla fija y no
  corresponde a nadie. El prefijo `SINTETICO_` y la carpeta `tests/` están ahí
  para que no se confunda nunca con datos reales.

## 7. Código

| Qué | Dónde |
| --- | --- |
| Análisis | `src/oncogo_crc/v2/cohorte_cea.py` |
| Cohorte sintética | `tests/SINTETICO_cohorte_cea.csv` |
| Generador de la cohorte | `scripts/literatura/generar_sintetico_cohorte_cea.py` |
| Tests | `tests/test_cohorte_cea_d9.py` |

El análisis corre entero sobre la cohorte sintética. Sobre esos 300 pacientes
sintéticos, el modelo del umbral da un AUC de 0,60, el del valor 0,70 y el de la
trayectoria 0,70; la diferencia entre el valor y el umbral no incluye el cero y la
diferencia entre la trayectoria y el valor sí.

**Esos números no dicen nada sobre la realidad.** Son los del mecanismo que se
programó en el generador. Lo único que demuestran es que el código detecta una
diferencia cuando existe y no la fabrica cuando no existe, que es exactamente lo
que hay que saber antes de tener datos de verdad.

---

Ningún resultado de esta tesis sirve para decidir sobre un paciente.
