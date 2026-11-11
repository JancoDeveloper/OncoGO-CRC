# Web v3.7, pendiente de publicación

El sitio se genera con `scripts/v3/exportar_sitio_v3_7.py` y los tests de `tests/test_sitio_v3_7.py` pasan. Las páginas de la v3.6 están en `sitio/archivo/v3_6/`. No se ha desplegado. La página pública sigue en la v3.6.

`scripts/v3/verificar_enlaces.py` no termina en cero. No se publica un enlace que no responda 200 con el título esperado, y no se despliega mientras esa comprobación falle.

Siguen sin poder publicarse:

- `https://europepmc.org/article/MED/40652393` responde 403.
- `https://doi.org/10.1021/acs.jproteome.5c00285` termina en 403.
- `https://doi.org/10.1021/acs.jproteome.2c00145` termina en 403.
- `https://doi.org/10.1016/j.celrep.2024.113810` y `https://doi.org/10.1016/j.cell.2019.03.030` responden 200 con el título «Redirecting», no con el título del artículo. El sitio del editor responde 403 a este cliente.

El resto de `config/enlaces_externos.tsv` respondió 200 y sí entra en la web local.
