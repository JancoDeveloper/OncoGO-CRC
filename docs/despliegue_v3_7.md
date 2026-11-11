# Despliegue de la web v3.7

Sitio: https://oncogo-crc-tesis.web.app

Proyecto de Firebase Hosting: `oncogo-crc-tesis`.

Versión de Hosting: `d5ce6e8f744c7937`.

Hora de la versión en vivo: 2026-11-11T11:29:29.207Z (06:29:29 en UTC−5).

Mensaje: `v3.7-final 595118a, sitio 6e1070d`.

El commit del sitio desplegado es `6e1070d`. La etiqueta `v3.7-final` no se movió y sigue en `595118ae09ea6006c744c86b28cb0a10d42e0442`.

## Comprobación en vivo, con `?v=37`

| Ruta | Código | Aviso | Franja publicable y etiqueta | `results/v3_6` fuera de `/archivo/` | CN como pacientes sintomáticos |
|---|---|---|---|---|---|
| `/` | 200 | 2 | sí | no | sí, en la primera frase |
| `/lo-prometido` | 200 | 2 | sí | no | dice que CN no son personas sanas; la definición está en el inicio |
| `/datos` | 200 | 2 | sí | no | sí |
| `/tejido` | 200 | 2 | sí | no | no nombra CN |
| `/suero` | 200 | 2 | sí | no | sí |
| `/aprendizaje-automatico` | 200 | 2 | sí | no | no nombra CN |
| `/preanalitica` | 200 | 2 | sí | no | no nombra CN |
| `/hipotesis` | 200 | 2 | sí | no | no nombra CN |
| `/priorizacion` | 200 | 2 | sí | no | nombra el contraste CP frente a CN, sin llamarlos sanos |
| `/desviaciones` | 200 | 2 | sí | no | no nombra CN |
| `/reproducir` | 200 | 2 | sí | no | no nombra CN |
| `/como-se-hizo` | 200 | 2 | sí | no | no nombra CN |
| `/aprende` | 200 | 2 | sí | no | no nombra CN |
| `/glosario` | 200 | 2 | sí | no | sí |
| `/preguntas` | 200 | 2 | sí | no | sí |
| `/historial` | 200 | 2 | sí | no | no nombra CN |
| `/archivo/` | 200 | 2 | sí | no | no nombra CN |
| `/archivo/v3_6/` | 200 | 2 | es el sitio anterior | la ruta está dentro de `/archivo/` | no aplica |
| `/no-existe` | 404 | 2 | sí | no | no nombra CN |

La portada en vivo muestra el certificado publicable, la etiqueta `v3.7-final` y, lo primero, 0/166/118 frente a CN.
