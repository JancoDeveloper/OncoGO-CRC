# Web v3.7 pendiente

La corrida B quedó en `results/v3_7/` con certificado publicable. La web no se reconstruyó en esta sesión y no se desplegó.

Queda por hacer, en este orden:

- `scripts/v3/exportar_sitio_v3_7.py`, sin modificar el generador de v3.6.
- Mover las páginas actuales a `sitio/archivo/v3_6/` sin borrarlas.
- Franja del certificado leída de `results/v3_7/certificado.json`, visible en todas las páginas.
- CN descrito como pacientes sintomáticos con otra enfermedad intestinal, no como sanos.
- Página `lo-prometido` desde `config/correspondencia_propuesta.tsv`.
- Gráficos 2, 3, 5, 8 y 10 como mínimo: embudo del suero, celda 0/166/118, R1, bosque de AUC y orden de las 60 muestras.

No se crea la etiqueta `v3.7-final`. Las dos corridas coinciden en los números principales, pero no comparten sello: la primera no declaraba `results/certified/inference_audit.json` y su certificado quedó `no_publicable`. Publicar en producción es decisión del autor.
