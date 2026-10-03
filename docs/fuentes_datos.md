# Fuentes de datos de v2

Una ficha por fuente. El estado solo puede ser uno de tres:

- **usada**: descargada, parseada y en uso en algún resultado de `results/v2/`.
- **no_disponible**: la descarga falló o el recurso no responde. Se anota el
  motivo exacto y qué análisis queda sin evaluar por su causa.
- **requiere_acceso**: existe, pero hace falta una solicitud, un acuerdo de uso
  o credenciales que este trabajo no tiene.

Ninguna fuente marcada `no_disponible` o `requiere_acceso` se sustituye por
datos simulados. Los datos sintéticos solo existen dentro de `tests/`, en
ficheros con el prefijo `SINTETICO_`.

El estado legible por máquina vive junto a cada descarga, en
`external/<fuente>/estado.json`. Esta página es la versión para leer.

---

## Cohortes del certificado (heredadas de v1, sin volver a descargar)

### PDC000116 — CPTAC cáncer de colon, cohorte de descubrimiento

- **Estado**: usada
- **Qué es**: proteoma de tumor y mucosa normal del mismo paciente, marcaje
  isobárico TMT, cuantificación del Proteomic Data Commons.
- **Qué aporta a v2**: los 94 pares tumor–normal del contraste principal.
- **Columnas**: se usa exclusivamente la familia `Unshared Log Ratio`, que son
  los péptidos no compartidos. La matriz **ya viene en log2**, así que no se
  vuelve a transformar. Cada espécimen aparece dos veces en el fichero, una por
  familia de columnas, y las dos familias nunca se mezclan.
- **Trampa conocida**: las filas `Mean`, `Median` y `StdDev` no son proteínas y
  se descartan.

### PXD034575 — Tanaka 2024, réplica independiente

- **Estado**: usada
- **Qué es**: LFQ con MaxQuant de tumor, mucosa normal y pólipo del mismo
  paciente, 16 pacientes.
- **Qué aporta a v2**: la réplica de dirección (R1), el brazo pareado del nulo
  exacto de C1, el brazo de pólipo de C5 y las 16 unidades de validación de C6.
- **Nombres de columna**: `LFQ intensity CRC<n><M|N|P>`, donde **M es el tumor**
  (maligno), N el tejido normal y P el pólipo. No es `T`.
- **Trampa conocida**: el cero de MaxQuant significa ausencia de identificación,
  no abundancia cero. No se imputa.

### PXD055218 — suero DIA, cáncer frente a sano

- **Estado**: usada
- **Qué es**: 20 sueros de cáncer colorrectal y 20 de control sano, DIA sin
  enriquecimiento.
- **Qué aporta a v2**: el lado de circulación de la brecha tejido→plasma y de la
  ventana de detección de C8.
- **Limitación estructural**: no hay ningún paciente con tejido y sangre a la
  vez. Toda comparación tejido→plasma es entre cohortes, nunca intraindividual.

---

## Fuentes nuevas de v2

### Gene Ontology — Biological Process

- **Estado**: usada
- **Fichero**: `metadata/gene_ontology/go_gene_long.tsv`, el mismo de v1
- **Qué aporta**: una de las tres colecciones de la cláusula de vías (C4).
  11 248 términos BP.
- **Resultado en v2**: de 132 términos significativos en CPTAC, solo 6 llegan a
  ser evaluables en la réplica, por debajo del mínimo preregistrado de 10. La
  colección queda **no evaluable**.

### Reactome

- **Estado**: usada
- **Script**: `scripts/download/colecciones_vias.py`
- **Fichero**: `external/vias/reactome_pathways.gmt`
- **Origen**: `https://reactome.org/download/current/ReactomePathways.gmt.zip`
- **Licencia**: CC BY 4.0
- **Cita**: Milacic M et al. The Reactome Pathway Knowledgebase 2024. *Nucleic
  Acids Research*.
- **Contenido**: 2 868 términos, 142 657 pares término-gen.
- **Resultado en v2**: la única de las tres colecciones evaluable. 449 términos
  significativos en CPTAC, 51 evaluables en la réplica, concordancia 47/51.

### MSigDB Hallmark

- **Estado**: usada
- **Script**: `scripts/download/colecciones_vias.py`
- **Fichero**: `external/vias/msigdb_hallmark.gmt`
- **Licencia**: CC BY 4.0 para la colección Hallmark. El resto de MSigDB tiene
  condiciones propias y no se usa.
- **Cita**: Liberzon A et al. The Molecular Signatures Database Hallmark Gene Set
  Collection. *Cell Systems* 2015.
- **Contenido**: 50 términos, 7 322 pares término-gen.
- **Resultado en v2**: 7 términos significativos en CPTAC, 2 evaluables en la
  réplica. **No evaluable**.

### Human Protein Atlas — concentraciones en sangre

- **Estado**: usada
- **Script**: `scripts/download/human_protein_atlas.py`
- **Fichero**: `external/hpa/hpa_concentracion_sangre.tsv`
- **Origen**: `https://www.proteinatlas.org/download/proteinatlas.tsv.zip`
- **Licencia**: CC BY-SA 4.0
- **Cita**: Uhlén M et al. Tissue-based map of the human proteome. *Science*
  2015.
- **Contenido**: 20 162 genes; 4 285 con concentración por espectrometría de
  masas y 457 por inmunoensayo. Del TSV completo solo se guardan las columnas
  necesarias; el resto no entra en el workspace.
- **Qué aporta**: la corrección C8. Permite estimar la probabilidad de detección
  en función de la concentración y medir la concordancia solo dentro de la
  ventana de detección.
- **Limitación**: las concentraciones son de sujetos de referencia sanos, no de
  los pacientes de PXD055218. Sirven para el orden de magnitud, no como medida
  del individuo. Además solo una fracción del proteoma tiene concentración
  publicada, y está sesgada hacia proteínas abundantes y ya estudiadas.

---

## Recordatorio

Ningún resultado de esta tesis sirve para decidir sobre un paciente.
