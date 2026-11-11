# Inventario de datos científicos

Fecha: 2026-11-11. La web no copia estos ficheros a otro formato para «mejorar» un gráfico. `dato_html` solo acepta rutas bajo `results/v3_7/` y `config/`. Si la ruta no existe, el generador se detiene.

## Corrida

`results/v3_7/`. Certificado publicable. Sello `71b310a3fcaebe9267fac971156df3aab491e009`, fecha 2026-11-11T00:48:45Z. Etiqueta `v3.7-final` en `595118ae09ea6006c744c86b28cb0a10d42e0442`. Hay ficheros TSV y JSON de tejido (k1), suero (k4), hipótesis H5–H10, aprendizaje, preanalítica, priorización y el embudo del suero. Los `*.sello.json` acompañan a esos resultados. No se volvieron a calcular.

## Suero y tejido, por separado

Suero, celda preregistrada, `suero_embudo.json`: n 284, diferentes 0, equivalentes 166, indeterminadas 118, contraste CP frente a CN. La matriz principal tiene `n_grupos_en_matriz_principal` 312.

Tejido: `k1_etiquetas.json`, contraste primario frente a mucosa, numerador y denominador de la concordancia de signo. Las tablas `k1_*.tsv` son tejido. Las tablas `k4_*.tsv` y `suero_k4_*.tsv` son suero. No se fusionan.

## Priorización

`priorizacion_embudo.json`, siete pasos, con el nombre y el n guardados en el JSON:

| Paso | n |
|---|---|
| 1. diferencial certificada en CPTAC (TREAT) | 688 |
| 2. medida en Tanaka P frente a N y con el mismo signo | 135 |
| 3. clase H5 (se informa; no elimina) | 135 |
| 4. HPA: secretada a sangre o con concentración plasmática conocida | 118 |
| 5. medida en la matriz principal del suero PXD055218 | 20 |
| 6. K4 CP frente a CN: clase diferente con el mismo signo que el tejido | 0 |
| 7. H10: especifica_crc en el HPA Disease Blood Atlas | 0 |

`nadie_pasa_el_paso_6` es verdadero. La lista `hipotesis_para_cohortes_futuras` tiene 20 genes. La lectura del propio JSON dice que no son un hallazgo para un paciente.

## Hipótesis

Estados en JSON: H5 `sostiene`, H6 `no_sostiene`, H7 `refuta`, H8 `no_evaluable`, H9 `exploratorio`, H10 `descriptiva`. H7 guarda además el AUC del modelo. Esos estados no se promedian.

## Anotación, fuera de la corrida

UniProt `2026_04`. Gene Ontology `releases/2026-07-26`. AlphaFold: 348 respuestas con `latestVersion` 6 y cuatro respuestas vacías (`P04114`, `P22352`, `P49908`, `Q9Y6R7`). STRING: estado `ok`, 59 aristas, versión vacía. Es asociación previa, no un resultado de la tesis.

## Qué no está

No hay respuesta de API guardada para el instrumento, el tamaño ni los grupos de PXD055218, PXD034575, PDC000116, HPA ni Gene Ontology. LinkedOmicsKB y PXD064655 no están en `manifest_entradas.sha256`. No hay un campo de segmento del colon para las muestras. No se calculó la similitud de Wang ni una superposición de estructuras. El mapa de las 284 parejas no se dibuja.

Cuando una pieza de la interfaz necesita ese dato, el texto es: «Información no disponible en las fuentes del proyecto».
