# Auditoría para la lectura visual

Fecha: 2026-11-11. La ciencia sigue siendo la v3.7-final, etiqueta `595118ae09ea6006c744c86b28cb0a10d42e0442`. No se modifican `results/v3_7/` ni la tesis.

El sitio se genera con Jinja en `scripts/sitio/construir.py` y Firebase Hosting sirve `sitio/`. No hay `package.json`. React, React Flow, ECharts y Lucide no entran: un diagrama o una barra de esta tesis se resuelve con SVG y CSS, y una librería nueva no aporta un dato. BioRender no se usa: sus ilustraciones no son de esta investigación y no hay permiso de publicación revisado. Inter, Manrope e IBM Plex no se cargan: el sitio ya aloja Source Sans 3 y Source Serif 4, y no descarga fuentes de terceros.

## Qué se puede mostrar

| Elemento | Fuente | Tipo | ¿Se visualiza? |
|---|---|---|---|
| Celda 0 / 166 / 118 / 284 | `results/v3_7/suero_embudo.json#celda_preregistrada` | Dato real | Sí. 0 + 166 + 118 = 284 |
| Diferente | Glosario, TREAT | Definición | Sí, como clase, no como correspondencia |
| Equivalente | Glosario, TOST, margen de 1,5 veces | Definición | Sí |
| Indeterminada | Glosario, potencia | Definición | Sí. No es «no clasificada por falta de ficha» |
| Embudo, siete pasos | `priorizacion_embudo.json` | Dato real | Sí. El paso 6 y el 7 valen 0 |
| Concordancia de signo 147 / 147 | `k1_etiquetas.json` | Dato real | Sí, como tejido |
| H7 `refuta` | `h7_visibilidad_en_sangre.json#estado` | Estado declarado | Sí, como estado, no como hallazgo positivo |
| 284 fichas | `suero_k4_crc_cn.tsv` y UniProt guardado | Dato real | Sí. Nombre y gen solo si la descarga los trae |
| Paso 5 | `priorizacion.tsv` `ultimo_paso` 5 | Hipótesis | Sí, con la etiqueta ya fijada |
| Colon por segmentos | No hay campo de localización | Limitación | Sí, como esquema que dice que no localiza |
| Redes | Estado declarado, no calculadas | Limitación | No como red |
| Cohorte nueva | No medida en esta tesis | Limitación | No como resultado |
| Tránsito tumor → sangre | No está en los resultados | Limitación | No como animación |

CN sigue siendo pacientes sintomáticos con otra enfermedad intestinal. El aviso exacto permanece dos veces en cada página.
