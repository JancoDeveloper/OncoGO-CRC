# Catalogo de datasets y estado operativo

Este documento resume, para cada cohorte del proyecto, **de donde se obtiene**, **que informacion aporta**, **para que sirve en la tesis** y **si ya esta funcional en descarga remota real**.

## Criterios de estado

- **Funcional (remoto OK):** `download_real` resuelve URLs HTTP/HTTPS utiles (archivo de datos), no assets genericos.
- **Parcial:** hay inventario/fuente, pero falta URL final de matriz o metadata.
- **No funcional:** fuente sin resolver automaticamente o URL invalida para uso analitico.
- **Contextual:** dataset de soporte narrativo (no bloquea pipeline nuclear).

## Resumen ejecutivo

| Dataset | Proveedor | Rol en proyecto | Estado remoto actual |
|---|---|---|---|
| `cptac_colorectal` | PDC | Nucleo biologico (tumor vs normal) | **No funcional remoto definitivo** (URL autodetectada no util de archivo) |
| `pride_crc_meta` | PRIDE | Validacion externa en tejido | **Parcial** (accesion definida; falta URL final de matriz/metadata) |
| `plasma_crc_2024` | PX/iProX | Capa clinica/traslacional | **No funcional remoto definitivo** (URL autodetectada generica) |
| `crcdb_context` | Manual/CRCDB | Contexto EOCRC y anotacion | **Contextual** (manual-curated, no descarga automatica de matriz proteomica) |
| `pharmacoproteomics` | PRIDE | Extension de respuesta a tratamiento | **Parcial/No funcional** (accesion de referencia requiere verificacion final) |
| `metastatic_crc_optional` | PX | Extension metastasis | **Parcial/No funcional** (accesion definida; falta URL final util) |

## Detalle por dataset

### 1) `cptac_colorectal`
- **Fuente:** PDC (CPTAC-2 / Colon Cancer, via GraphQL).
- **Que trae:** abundancia proteica en tejido y metadata clinica por caso/muestra.
- **Utilidad:** base del analisis principal (`cptac_tumor_vs_normal`), GO, red, ML, traduccion.
- **Estado:** no funcional remoto definitivo al momento (la autodeteccion puede elegir URLs HTTP no tabulares); requiere fijar URL de archivo final o mejorar selector PDC por tipo de archivo.

### 2) `pride_crc_meta`
- **Fuente:** PRIDE (`PXD015163` configurado en `datasets.yaml`).
- **Que trae:** archivos de proyecto (matrices/metadata si el estudio los publica como procesados).
- **Utilidad:** validacion externa de hallazgos del nucleo.
- **Estado:** parcial; la API de archivos esta integrada, pero falta resolver de forma deterministica el par `matrix_url`/`metadata_url` final.

### 3) `plasma_crc_2024`
- **Fuente:** ProteomeXchange / iProX (`PXD047207`, fallback iProX configurado).
- **Que trae:** potencial de matriz de plasma y variables clinicas.
- **Utilidad:** capa traslacional (CRC vs healthy, respuesta terapeutica).
- **Estado:** no funcional remoto definitivo; la resolucion actual puede caer en URL de pagina/endpoint generico y no archivo final.

### 4) `crcdb_context`
- **Fuente:** CRCDB (manual curated).
- **Que trae:** contexto, anotaciones y soporte de discusion (no necesariamente matriz proteomica estandar).
- **Utilidad:** fortalecer discusion y contexto EOCRC.
- **Estado:** contextual; no bloquea pipeline nuclear, pero no opera como descarga automatica de matriz.

### 5) `pharmacoproteomics`
- **Fuente:** PRIDE (`PXD019234` configurado).
- **Que trae:** potencial de perfiles asociados a tratamiento/respuesta.
- **Utilidad:** extension para biomarcadores funcionales y priorizacion traslacional.
- **Estado:** parcial/no funcional por ahora; requiere confirmar accession util y archivos procesados reales.

### 6) `metastatic_crc_optional`
- **Fuente:** PX (`PXD021238`).
- **Que trae:** posible cohorte metastasica para contraste biologico adicional.
- **Utilidad:** extension de progresion/metastasis y estabilidad de hallazgos.
- **Estado:** parcial/no funcional remoto definitivo; falta URL de matriz/metadata final util.

## Que informacion se descarga y por que es util

Cada dataset idealmente aporta:
1. **Matriz proteomica** (`protein_id` + columnas `sample_id`) para analisis cuantitativo.
2. **Metadata de muestra** (`sample_id`, `condition`, y covariables clinicas) para comparaciones y estratificacion.
3. **Evidencia de procedencia** (manifest/provider manifest) para reproducibilidad.

Esto alimenta:
- diferencial (proteinas alteradas),
- GO (procesos funcionales),
- redes/hubs,
- ML explicable,
- ranking traslacional.

## Archivos de evidencia

- `data/raw/download_real_manifest.tsv`
- `data/manifests/_summary.json`
- `results/validation/input_validation_report.json`
- `reports/web/datasets.md`
- `reports/web/clinical_coverage.md`

## Validacion de que la consulta es correcta (y por que)

La validacion no solo revisa formato, tambien comprueba si los datos descargados sirven para la **pregunta analitica** de cada dataset:

1. **Estructura**  
   - matriz con `protein_id` + `sample_id`  
   - metadata con `sample_id` y `condition`

2. **Coherencia de muestras**  
   - matching exacto entre columnas de matriz y `sample_id` de metadata

3. **Fuentes remotas utiles**  
   - URL `http(s)` real  
   - no `file://`  
   - no assets/web genricos (css/js/fonts/cgi sin archivo final)

4. **Soporte de comparaciones nucleares**  
   Para cada comparacion no-opcional en `config/comparisons.yaml`:
   - existe la columna de grupo (`condition_column`/`group_col`)
   - existen los valores `case/control` definidos
   - hay solapamiento real de muestras matriz-metadata

Esto queda reportado en:
- `results/validation/input_validation_report.json`  
  - bloque `datasets.<dataset>.intent_validation`  
  - bloque `remote_sources`

Si `fail_on_validation_error: true`, la corrida falla cuando:
- hay errores estructurales,
- no hay matching exacto,
- hay `remote_sources_invalid_url`,
- o alguna comparacion nuclear no esta soportada por los datos.

## Checklist tecnico rapido por dataset

Para marcar un dataset como **funcional**, deben cumplirse todos:

- [ ] Fuente remota resoluble y trazable (accession/study + manifest)
- [ ] `matrix_url` y `metadata_url` utiles (no genericos)
- [ ] Matriz parseable y estandar (`protein_id`, `sample_id`)
- [ ] Metadata parseable con `sample_id`, `condition`
- [ ] Matching exacto matriz-metadata
- [ ] Comparaciones nucleares del dataset soportadas (case/control presentes)
- [ ] Estado final en validacion: sin errores bloqueantes

## Siguiente paso recomendado para dejar todo funcional remoto

1. Endurecer seleccion por proveedor para aceptar solo archivos tabulares finales.
2. Fijar `matrix_url` y `metadata_url` deterministas por dataset activo.
3. Re-ejecutar:
   - `oncogo-crc run --stage download_real`
   - `oncogo-crc run --stage validate_inputs`
   - `oncogo-crc run --stage all`
4. Confirmar estado **Funcional (remoto OK)** en todos los datasets activos.

## Modo arranque (para iniciar ya)

Para no bloquear el inicio del trabajo mientras se cierran URLs finales:

1. Mantener activos solo datasets nucleo (`cptac_colorectal`, `pride_crc_meta`, `plasma_crc_2024`).
2. Usar `download_real.fail_on_validation_error: false` durante bootstrap.
3. Ejecutar:
   - `oncogo-crc run --stage download_real`
   - `oncogo-crc run --stage validate_inputs`
   - `oncogo-crc run --stage all`
4. Revisar manifests y convertir warnings en tareas de cierre.
5. Cuando nucleo este limpio, volver a modo estricto:
   - `download_real.fail_on_validation_error: true`
