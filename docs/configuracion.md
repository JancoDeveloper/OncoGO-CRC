# Configuración

La corrida lee estos ficheros. `config/preregistro_v3.yaml` conserva el
nombre: su SHA-256 está en el sello (`34603ba0…`) y cambiar el nombre
rompería la comprobación.

| Fichero | Qué fija |
| --- | --- |
| `config/preregistro_v3.yaml` | Alfa, margen TREAT, semillas, paneles, multiverso y enmiendas fechadas |
| `config/config.yaml` | Semilla global 42 y rutas. `paths.data_workspace` puede apuntar a un volumen; la variable `WORKSPACE` o `--data-workspace` ganan |
| `config/datasets.yaml` | Las cuatro cohortes: PDC000116, PXD034575, PXD055218 y PXD064655 |
| `config/entradas_v3_7.yaml` | Excepción del vigilante para `results/certification/certificado_resumen.json` |
| `config/esperado_v3_7.tsv` | 35 métricas y su tolerancia (0 en conteos, 1e-9 en AUC y p) |
| `config/comparisons.yaml` | Lo carga `load_settings`. La etapa v3.7 no lo reescribe |

`results/certified/` ya está permitido por el vigilante. No hace falta
listarlo en `entradas_v3_7.yaml`.

Para mover los datos de sitio:

```bash
export WORKSPACE=/ruta/con/espacio
# o
oncogo-crc run --data-workspace "$WORKSPACE" --stage v3_7
```

No hace falta editar el YAML. `config/config.yaml` sigue llevando la ruta
que usó la corrida sellada para que su huella no cambie.
