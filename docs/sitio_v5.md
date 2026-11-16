# Sitio v5

Sitio estático generado por `scripts/sitio/construir.py`. Lo que se publica sigue siendo `sitio/`. No se migró el generador.

## Qué cambió

- Bibliografía de la tesis v2.7: 86 entradas APA en `sitio_fuente/referencias.bib` y `sitio_fuente/referencias_apa.yaml`. El generador ya no lee `tesis-ucaldas/latex/references.bib`.
- Citas en el texto: «(Apellido, año)», «(A y B, año)» o «(A et al., año)», con enlace a `/bibliografia#clave` y un recuadro con autores, año, título, revista y DOI.
- Hipótesis copiadas de `sitio_fuente/tesis_v2_7/hipotesis-objetivos.md`. El estado va en un chip. H4 muestra el umbral 0,90 y no un veredicto. H3 conserva la frase de la tesis y, aparte, el chip de la enmienda de cierre.
- `sitio_fuente/etiquetas.yaml` traduce `true`, `rechaza_h1`, `no_evaluable` y el resto de tokens. El build se detiene si la prosa visible (fuera de `<details>`, el pie y los `<script>`) contiene esos tokens o una ruta `results/`, `config/`, `.tsv`, `.json` o `.yaml`.
- Los enlaces al repositorio usan la etiqueta `v3.7-final`. El sello de cada figura sigue guardando el commit `595118ae09ea6006c744c86b28cb0a10d42e0442`.
- Cabecera con Aprende, Datos, Métodos, Resultados y Proteínas. Glosario y Bibliografía quedan a la derecha. En pantallas de menos de 900 px el menú principal va a un cajón.
- Aviso ético: bloque «Para qué sirve y para qué no» en la portada, una línea en el pie de todas las páginas, y un recuadro en Resultados y en las páginas de aprendizaje automático. No hay banner rosado (`class="aviso"`).
- Recorrido `/aprende/` de nueve páginas, con anterior y siguiente. `/aprende` sigue redirigiendo a `/datos/bases-de-datos`.
- Glosario: 145 términos. `{{termino:X}}` es un botón que abre el texto llano, un enlace al glosario y los relacionados. Se comprobó con el término CN en la portada.
- Diagramas D1–D11 dibujados en el build como SVG con ramas, rótulos y nodos enlazados. No se añadió Mermaid.
- Gráficos: se mantienen los SVG de la corrida, con «Cómo leer este gráfico», tabla y descarga. El tablero n11 ya no dice que el preregistro no tiene H3.
- Buscador: `sitio/buscar/indice.json`. Encuentra KNG1 (P01042 y P01042-2), p02749, apoa4 (P06727) y S100A8;S100A9 (P05109 y P06702).
- Fichas: resumen en lenguaje llano al inicio. Mol* sigue cargándose solo al pulsar el botón.
- `sitemap.xml`, `robots.txt`, descripción, Open Graph de título y JSON-LD `ScholarlyArticle` en la portada.
- Páginas `/como-citar` y `/sobre`. `CITATION.cff` declara la licencia MIT del repositorio.

## Pruebas

`pytest -q tests/test_sitio_v4.py tests/test_sitio_v3_7.py`: 59 passed.

El emulador de hosting quedó en el puerto 5002. En la portada, a 1440 px de ancho, hay una sola cabecera, un solo título y el documento no desborda en horizontal. El botón CN abre el recuadro del glosario. `/aprende/cancer-colorrectal` enlaza con «Siguiente: Del tejido a la sangre» y conserva la frase de que la tesis no asigna muestras a un segmento del colon.

## Qué no se midió

Lighthouse y axe no se ejecutaron en esta pasada: `lighthouse` no está instalado en la máquina. No hay cifras nuevas de rendimiento ni de accesibilidad automática. Las de `docs/sitio_v4.md` no se reutilizan como si fueran de esta versión. Tampoco dejé capturas en `docs/sitio_v5_capturas/`: la superficie de captura repetía la cabecera fija y no sirve como registro visual.

No hay imagen Open Graph propia. El JavaScript de la portada, ya incluido en la página, pesa unos 17 KB. No se cargó Vega.

## Qué no se hizo, y por qué

- Los gráficos no pasaron a Vega-Lite. `vega-embed` no cabe en el presupuesto de la portada y las páginas no pueden llevar `<script src>` ni una CDN. Sin JavaScript cada figura sigue siendo el SVG y su tabla.
- Pagefind no se indexó. El comando de despliegue lo incluye. El buscador que ya funciona es el índice JSON.
- No hay un grafo distinto por dataset en `/datos/<id>`. Las páginas de datos siguen sin inventar instrumento, n o grupos. Esos datos no están en `config/` ni en `sitio_fuente/fuentes/`; el detalle está en `docs/sitio_v4_omitido.md`. PXD064655 sigue nombrado y no evaluado. H8 sigue en «No se pudo evaluar».
- No se reescribió, página por página, la prosa antigua de las 39 páginas previas. Las cifras de esas páginas salen de `{{dato:}}`.
- No se desplegó.

## Despliegue

```bash
python scripts/sitio/construir.py && npx pagefind --site sitio && firebase deploy --only hosting
```
