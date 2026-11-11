# Despliegue de la web documental

Sitio: https://oncogo-crc-tesis.web.app

Proyecto de Firebase Hosting: `oncogo-crc-tesis`.

Versión de Hosting: `161f72104d6d9530`.

Hora de la versión en vivo: 2026-11-11T16:19:47.041Z (11:19:47 en UTC−5). `createTime` 2026-11-11T16:19:39.993489Z. `finalizeTime` 2026-11-11T16:19:47.252250Z.

El comando fue `firebase deploy --only hosting --project oncogo-crc-tesis -m "web documental, resultados v3.7-final 595118a, sitio e77d0e1"`. El recurso de la versión no guarda ese mensaje: `labels` solo tiene `deployment-tool=cli-firebase`. La configuración publicada tiene 16 redirecciones. El directorio público sigue siendo `sitio`. No hay React.

El commit desplegado es `e77d0e16ea8a14b0b5a8e4963b6f90ea6bc448dc`. La etiqueta `v3.7-final` no se movió y sigue en `595118ae09ea6006c744c86b28cb0a10d42e0442`.

La versión anterior de esta web documental, `0fe6e84d75d1c6df`, quedó liberada el 2026-11-11T15:42:46.751Z con el sitio `e4b8c81`. Sigue en el historial de Hosting. No se copió a otra carpeta. La v3.7 anterior sigue archivada en `/archivo/v3_7`.

Antes de publicar, en el emulador `127.0.0.1:5002` y con los mismos ficheros, Lighthouse 13.6.0 dio rendimiento, accesibilidad y buenas prácticas: `/` escritorio 100/100/96 y móvil 100/100/96; `/metodos/proteomica` 100/100/96 en los dos; `/proteinas/P02749` 100/100/96 en los dos. La portada, ya con `v4.a1fc0bd2e5.css`, no tuvo violaciones serious ni critical de axe. Esas medidas son del emulador, no de la URL pública.

## Comprobación en vivo, con `?v=4`

El aviso es «Ningún resultado de esta tesis sirve para decidir sobre un paciente.» La hoja es `/css/v4.`.

| Ruta | Código | Aviso | Hoja v4 |
|---|---|---|---|
| `/` | 200 | 2 | sí |
| `/investigacion/contexto` | 200 | 2 | sí |
| `/investigacion/historia` | 200 | 2 | sí |
| `/datos` | 200 | 2 | sí |
| `/datos/bases-de-datos` | 200 | 2 | sí |
| `/metodos/canalizacion` | 200 | 2 | sí |
| `/resultados/suero` | 200 | 2 | sí |
| `/resultados/tejido` | 200 | 2 | sí |
| `/resultados/aprendizaje` | 200 | 2 | sí |
| `/resultados/preanalitica` | 200 | 2 | sí |
| `/resultados/hipotesis` | 200 | 2 | sí |
| `/resultados/priorizacion` | 200 | 2 | sí |
| `/discusion` | 200 | 2 | sí |
| `/conclusiones` | 200 | 2 | sí |
| `/entregables` | 200 | 2 | sí |
| `/archivo` | 200 | 2 | sí |
| `/ayuda` | 200 | 2 | sí |
| `/glosario` | 200 | 2 | sí |
| `/bibliografia` | 200 | 2 | sí |
| `/proteinas` | 200 | 2 | sí |
| `/proteinas/P02749` | 200 | 2 | sí |
| `/proteinas/candidatas` | 200 | 2 | sí |
| `/proteinas/comparar?p=P02749,P05155` | 200 | 2 | sí |
| `/buscar` | 200 | 2 | sí |
| `/no-existe` | 404 | 2 | sí |

`/bibliografia` contiene «Rai» y `rai2025`. El índice de `/buscar` trae `P02749 APOH` y `TOST`; no lista el apellido. `/proteinas/comparar?p=P02749,P05155` trae `vista-resumen`, `vista-go`, `vista-secuencia`, `vista-dominios` y `vista-estructura`, y un `data-modelos` local.

| Ruta vieja | Código | Destino |
|---|---|---|
| `/lo-prometido` | 301 | `/investigacion/historia` |
| `/historial` | 301 | `/investigacion/historia` |
| `/tejido` | 301 | `/resultados/tejido` |
| `/suero` | 301 | `/resultados/suero` |
| `/aprendizaje-automatico` | 301 | `/resultados/aprendizaje` |
| `/preanalitica` | 301 | `/resultados/preanalitica` |
| `/hipotesis` | 301 | `/resultados/hipotesis` |
| `/priorizacion` | 301 | `/resultados/priorizacion` |
| `/desviaciones` | 301 | `/entregables` |
| `/reproducir` | 301 | `/investigacion/reproducibilidad` |
| `/como-se-hizo` | 301 | `/investigacion/reproducibilidad` |
| `/aprende` | 301 | `/datos/bases-de-datos` |
| `/preguntas` | 301 | `/ayuda` |
| `/v3_3` | 301 | `/archivo/v3_3` |
| `/v3_4` | 301 | `/archivo/v3_4` |
| `/lectura-31-oct` | 301 | `/archivo/lectura-31-oct` |

La petición llevaba `?v=4`, y el `Location` lo conservó. `/archivo/v3_7` y `/archivo/v3_6` respondieron 200. Con barra final, 301 a la ruta sin barra.

Cinco fichas, muestra `random.Random(20261111)` sobre `suero_k4_crc_cn.tsv`. El `log2fc` de CP frente a CN coincidió con la celda:

| protein_id | log2fc |
|---|---|
| `P49747;P49747-2` | -0.0401515185583427 |
| `P05156` | -0.0381274240177411 |
| `P0DOX7` | -0.0405369755463027 |
| `P01743` | 0.4058961730010573 |
| `P59665;P59666` | -0.1436413247840597 |

`/proteinas/candidatas` nombra los 20 genes con `ultimo_paso` 5 en `priorizacion.tsv`. No falta ninguno.

La portada trae `no localiza un tumor` seis veces, `Las dos pistas no se unen` una vez, el embudo y las dos matrices. Los n del embudo están en `value`: 688, 135, 118, 20 y 0. Las cinco fichas de arriba incluyen «Resultado calculado».
