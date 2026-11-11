# Sistema de diseño

Fecha: 2026-11-11.

## Tecnología

Siguen `sitio_fuente/estilos/tokens.css` y `componentes.css`. El generador publica una sola hoja `v4.` más diez hexadecimales del contenido. No hay fuentes de Google ni de un CDN: la prueba `test_el_css_lleva_su_hash_y_no_hay_terceros` lo impide.

Manrope, Plus Jakarta Sans e Inter no se cargan. Exigirían una descarga externa o ficheros nuevos. El sitio ya aloja Source Sans 3 (texto) y Source Serif 4 (títulos), con licencia OFL. Los identificadores usan `ui-monospace` del sistema, sin fichero añadido.

## Color

El texto y el fondo tienen que superar 4,5:1 en claro y en oscuro. Los turquesas y el coral claros de una paleta decorativa no llegan a ese contraste como texto, así que no se usan para leer.

Claro:

| Token | Hex | Uso |
|---|---|---|
| `--papel` | `#f5f7fa` | Fondo |
| `--tinta` | `#142d4e` | Texto |
| `--tinta-suave` | `#34506e` | Texto secundario |
| `--marca` | `#0e4d6e` | Enlaces, barras del embudo, botones |
| `--sobre-marca` | `#f5f7fa` | Texto sobre el botón |
| `--aviso-fondo` | `#f8ebe3` | Fondo del aviso |
| `--aviso-tinta` | `#6a2c18` | Texto del aviso |
| `--tarjeta` | `#ffffff` | Superficie |

El tema oscuro conserva los tokens que ya pasaban AA. El conmutador sigue en el botón «Tema».

Azul (`--tinta`, `--marca`) es la información y el proceso. El coral oscuro del aviso es la advertencia. No se pinta un resultado de más o de menos con un color que no esté en la leyenda de la figura ya existente.

## Tipo y espacio

El título principal usa `clamp(1.8rem, 4vw, 2.7rem)`. El cuerpo sigue en 1,05 rem. Las cifras de la celda van en una lista `.cifras`. Las barras del embudo no tienen un mínimo de ancho: un n de 0 mide 0.

## Movimiento

No hay animación decorativa. `prefers-reduced-motion: reduce` anula la transición que pudiera heredar una barra.

## Componentes

`.boton` y `.boton-secundario` para las dos llamadas de la portada. `.embudo` para los siete pasos. `.colon` para el esquema. `.secundario` para el pie. `.tarjeta` se reutiliza en «Explorar».
