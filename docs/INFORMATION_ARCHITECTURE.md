# Arquitectura de la información

Fecha: 2026-11-11.

## Menú principal

| Etiqueta | Ruta | Qué había ahí |
|---|---|---|
| Inicio | `/` | La portada |
| La investigación | `/investigacion/contexto` | Contexto, y desde ahí la pregunta, la historia y la reproducibilidad |
| ¿Cómo funciona? | `/metodos/proteomica` | Proteómica y el resto de métodos |
| Explorar resultados | `/resultados` | Tejido, suero, aprendizaje, preanalítica, hipótesis, priorización |
| Aprende | `/ayuda` | Cómo leer una cifra, una cita, y tres comprobaciones |

`/aprende` no es una página nueva. Esa ruta ya redirige a `/datos/bases-de-datos` desde la web v3.7. «Aprende» apunta a `/ayuda` para no pisar el 301.

## Pie

Datos, Métodos (si su primera página no es ya el menú principal), Resultados cuando no esté repetido, Proteínas, Discusión, Conclusiones, Entregables, Archivo, Glosario y Bibliografía. Ninguna ruta de `sitio_fuente/contenido/` se borró.

## Portada

1. Pregunta, tal como está escrita en `/investigacion/pregunta-e-hipotesis`.
2. Dos llamadas: la investigación y los resultados.
3. La celda: diferentes, equivalentes, indeterminadas y n.
4. Definición de CN.
5. El problema, con las citas de tejido y el acceso del suero.
6. H7 y un enlace al resto de hipótesis.
7. D1 y el embudo de `priorizacion_embudo.json`.
8. Concordancia de signo en tejido y el cero del suero, en párrafos distintos.
9. Esquema del colon que no localiza muestras.
10. Seis accesos: proteómica, resultados, proteínas, aprendizaje automático, conclusiones, bibliografía.
11. Alcance.

## Qué no se añadió

Una red nueva: la de STRING sigue en `/proteinas/candidatas` y se presenta como asociación previa. Un cuestionario clínico. Una ilustración de proteínas saliendo del tumor. Filtros de matriz que no estén ya en las fichas.
