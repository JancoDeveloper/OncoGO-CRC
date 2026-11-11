# Trazabilidad del andamiaje

Cada cifra del PDF sale de `tesis-ucaldas/latex/cifras.tex`, que escribe `scripts/v3/exportar_cifras_v3_7.py` desde `results/v3_7/`. No se editó ese archivo a mano. No se modificaron `results/` ni los scripts que calculan.

El sello de la corrida es `71b310a3fcaebe9267fac971156df3aab491e009`, 2026-11-11T00:48:45Z, en `results/v3_7/certificado.json`. La etiqueta `v3.7-final` está en `595118ae09ea6006c744c86b28cb0a10d42e0442`.

## Cifras

| Qué aparece | Comando o valor | Archivo |
|---|---|---|
| Celda 284, 0, 166, 118 | `\SueroN` `\SueroDif` `\SueroEq` `\SueroInd` | `suero_embudo.json` `celda_preregistrada` |
| Embudo 897, 329, 312 | `\EmbudoGrupos` `\EmbudoQ` `\EmbudoHumanos` | `suero_embudo.json` `controles` |
| Signo 147 de 147 | `\TejidoNum` `\TejidoDen` | `k1_etiquetas.json` |
| Metástasis 116 de 143 | `\TejidoMetaNum` `\TejidoMetaDen` | `k1_etiquetas.json` |
| Cobertura 147 de 688 y H1 | `\CoberturaNum` `\CoberturaDen` `\HunoVeredicto` | `k1_etiquetas.json` |
| AUC 0,880, 0,915, 0,647, p 0,12 | `\MlTodas` `\MlCalidad` `\MlSinCalidad` `\MlSinCalidadP` | `ml_suero_permutacion.json` |
| AUC 0,947, 0,818, 0,887 | `\MlHumana` `\MlCoagulacion` `\MlPca` | el mismo JSON, M1, M4 y M6 |
| H5 a H10, estados | `\HcincoEstado` … `\HdiezEstado` | `h5` … `h10` JSON |
| H5, universo 29 | `\HcincoUniverso` | `h5_especificidad_organo.json` |
| H5, fracción 0,119 y 624 proteínas | valor leído del JSON, no hay comando | `fraccion_pan_carcinoma` 0,11859; `n_clasificadas` 624 |
| H7, AUC 0,483 y 0,464 | `\HsieteAlvez` `\HsieteHpa` | `h7_visibilidad_en_sangre.json` |
| H7, intervalos 0,381–0,584 (53 y 73) y 0,353–0,573 (53 y 63) | leídos del JSON `ic95_bajo`, `ic95_alto`, `n_pos`, `n_neg` | el mismo archivo, `por_recurso` |
| H7 secundario 0,474; 668 y 20 | `\HsieteSec` `\HsieteNEsm` `\HsieteOmitidas` | `modelo_secundario` |
| H8, 0 de 53 | `\HochoTumorales` `\HochoUniverso` | `h8_carga_tumoral.json` |
| H10, 5 proteínas | `\HdiezN` | `h10_especificidad.json` |
| Preanalítica, estado, 60 y 60 | `\PreEstado` `\PreCrudos` `\PreFecha` | `preanalitica.json` |
| AUC del orden 0,000 | `\PreAuc` | `preanalitica.json` `calculado_con_legibles.auc.valor` |
| p de Mann-Whitney \(1,5\times 10^{-11}\) | el JSON da 1,4508889103849684e-11; `\PreP` lo escribe con un exponente que la fuente no tiene | `preanalitica.json` |
| Embudo de priorización 688, 135, 135, 118, 20, 0, 0 | `\PrioriPasoUno` … `\PrioriPasoSiete` | `priorizacion_embudo.json` |
| 18 bajan y 2 suben | `\PrioriBajan` `\PrioriSuben` | `priorizacion.tsv` |
| HC: 3 y 1 | `\HcCpDif` `\HcCnDif` | `suero_k4_crc_sanos.tsv`, `suero_k4_cn_sanos.tsv` |
| GO 16 de 16; CPTAC 787, 17, 63, 64; Tanaka 101, 0, 0, 21 | comandos `\Go…` | `go_tejido.json` |
| Multiverso 288, 72, 13 | `\MultiSuero` `\MultiMl` `\MultiDif` | `multiverso_resumen.json` |
| Semilla 20260914 y paquetes | campo `sello` del certificado; también `semilla` en los JSON de ML y de H7 | `certificado.json` |
| Alfa 0,05, margen 0,585, mínimos 3 y 8, fracción 0,70 | no son una cifra de resultado | `config/preregistro_v3.yaml` `inferencia` |
| PDC000116, PXD034575, PXD055218, PXD064655 | accesiones | `config/datasets.yaml`, `config/preregistro_v3.yaml` |

## Figuras

El script `scripts/v3/exportar_sitio_v3_7.py` escribe el SVG en `sitio/figuras/`. Esas nueve figuras se copiaron a `tesis-ucaldas/figures/` y a `tesis-ucaldas/latex/figures/`. El PDF de cada una se obtuvo del SVG con svglib y reportlab, fuera del entorno del proyecto. No hay PNG a 300 dpi: [PENDIENTE: falta el backend rlPyCairo].

| Figura | Archivo | Qué muestra | Datos |
|---|---|---|---|
| f01 | `f01.svg`, `f01.pdf` | Mapa de tejido, suero, hipótesis y priorización | El trazo es el orden de lectura, no una proteína que sale del tumor |
| f02 | `f02.svg`, `f02.pdf` | Embudo del suero | `suero_embudo.json` |
| f03 | `f03.svg`, `f03.pdf` | Celda CP frente a CN y contrastes con HC | `suero_embudo.json`, `suero_k4_crc_sanos.tsv`, `suero_k4_cn_sanos.tsv` |
| f04 | `f04.svg`, `f04.pdf` | Rejilla de 288 celdas | `multiverso_resumen.json` |
| f05 | `f05.svg`, `f05.pdf` | log2FC de CPTAC y de Tanaka | Tabla de H7 y Tanaka |
| f06 | `f06.svg`, `f06.pdf` | Controles de tejido no medibles | [PENDIENTE: verificar el veredicto de CA1 en el fichero antes de escribirlo. AQP8, GUCA2A y SLC26A3 salen del título de la figura.] |
| f07 | `f07.svg`, `f07.pdf` | NES de CPTAC frente al de Tanaka | `go_tejido.json` |
| f08 | `f08.svg`, `f08.pdf` | AUC de M1 a M6 | `ml_suero_permutacion.json` |
| f09 | `f09.svg`, `f09.pdf` | SHAP de M2 y de M3 | El SVG lo produce el mismo script |

`Fig4-1`, `Fig4-2` y `Fig4-3` siguen en `latex/figures/` de la versión anterior. `Fig4-2_orden_adquisicion.png.sello.json` apunta a `results/v3_7/preanalitica_orden.tsv`. No se insertaron en este andamiaje: la figura de la corrida que entra en el PDF es `f08` para los modelos y los números de preanalítica van en texto. `F1` a `F14` tienen sello de la etapa `v3_3_figuras` (2026-11-09) y no se usan como cifras de la corrida firmada.

## Tablas

| Tabla | Dónde | Fuente |
|---|---|---|
| Herramientas de IA | `tab:ia` | Declaración del autor. No es un resultado. Cursor 3.28.7, leído de la aplicación instalada. SciSpace: consultado entre 2026-03 y 2026-11. |
| Cohortes | `tab:cohortes` | `config/datasets.yaml`, `config/preregistro_v3.yaml` |
| Celda del suero | `tab:celda` | `suero_embudo.json` |
| Hipótesis | `tab:hipotesis` | `h5` a `h10` JSON |
| Correspondencia | anexo, `\input{correspondencia}` | `config/correspondencia_propuesta.tsv`, vía `scripts` que escriben `correspondencia.tex` |
| Abreviaturas | preliminares | Definiciones usadas en el andamiaje. HC no se llama el contraste que decide la pregunta |

## [PENDIENTE]

- Veredicto de CA1 en el fichero de controles, antes de escribirlo.
- Si una frase sobre población latinoamericana va a entrar, hace falta el fichero que la sostenga. Ahora no está.
- Desplazamientos 720, 721, 840 y 841: están escritos en el preregistro. El JSON de la corrida guarda la semilla 20260914. Confirmar en el código que el desplazamiento se aplica antes de citar el entero.
- Claves de la bibliografía que aún no están en una sección: ver `GUIA-DE-REDACCION.md`.
