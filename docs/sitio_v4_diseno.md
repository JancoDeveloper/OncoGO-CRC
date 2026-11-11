# Diseño de la web documental

La distribución se miró en UniProt el 11 de noviembre de 2026. No se copia su logo, su marca, su código ni sus hojas de estilo. Las capturas propias están en `docs/sitio_v4_capturas/`.

## Patrones que se adoptan

La portada (`uniprot-portada.png`) tiene una barra fija con el nombre y el menú, un buscador grande y, debajo, tarjetas de acceso a los recursos. La web de OncoGO usa esa misma distribución: barra con buscador, la pregunta de la tesis en una frase, tarjetas y un bloque «Empieza aquí».

La entrada de P02749 (`uniprot-entrada.png`) confirma el resto: índice lateral con Function, Names & Taxonomy, Subcellular Location, Disease & Variants, PTM, Expression, Interaction, Structure, Family & Domains, Sequence, Similar Proteins y Homologs; título con el acceso y el nombre; ficha de campos cortos; pestañas; un párrafo de función; y la tabla de Gene Ontology con aspecto, término y evidencia. La ayuda (`uniprot-ayuda.png`) es una sola columna de preguntas como encabezados y respuestas en prosa. La búsqueda (`uniprot-busqueda.png`) es una tabla de accesiones con filtros a la izquierda; un aviso pide elegir tarjetas o tabla, y la tabla se ve detrás. Aquí la búsqueda lista todas las páginas sin JavaScript.

## Lo que no se adopta

No se adopta el azul institucional, el logotipo de puntos, la cesta, BLAST ni las noticias del consorcio. El color, la tipografía y los componentes están en `sitio_fuente/estilos/tokens.css`. La tipografía es Source Sans 3 y Source Serif 4, con licencia OFL, alojadas en el sitio.

## Esta fase

La fase 1 deja la portada, la página de proteómica, el glosario de los términos que esas páginas usan, la bibliografía de esas citas y el diagrama D1. El sitio anterior queda en `sitio/archivo/v3_7/`.

## Umbral del panorama de las candidatas

Quedó fijado el 2026-11-11, antes de calcular el mapa.

La vista directa es el conjunto de términos GO con al menos una evidencia distinta de IEA, juntando los tres aspectos. La vista propagada añade los ancestros `is_a` de ese conjunto. Si los dos conjuntos están vacíos, el índice de Jaccard es 0. Una arista de la red existe cuando el Jaccard de la vista directa es mayor o igual que 0,15. Un BinaryCIF de muestra (P02749, versión que devolvió la API) pesa 160133 bytes. Con ese tamaño, el corte de 40 MB que se había apuntado dejaría fuera el conjunto. Antes de descargar el resto, el corte pasa a 80 MB: si la suma de los BinaryCIF queda por debajo, se alojan; si no, no entran al sitio y queda el motivo. El mapa de calor ordena las hojas con enlace promedio (UPGMA) y distancia 1 menos ese Jaccard. Si dos hojas empatan, manda el acceso en orden alfabético.
