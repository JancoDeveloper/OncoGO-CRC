# Diagnóstico del corte de `make v3-7`

Fecha de la prueba: 2026-11-10, 14:44–14:59 (UTC−5). Árbol en `a3c5636`, que solo añade `logs/` a `.gitignore` respecto de `9af0d01`.

## Qué se ejecutó

La orden fue la de la etapa, con el mismo entorno que el Makefile (`PYTHONHASHSEED=0` y un hilo para BLAS), más `caffeinate -dims`, `/usr/bin/time -l` y `python -X faulthandler`. La salida está en `logs/v3_7_corrida.log`. El RSS del proceso `Python -X faulthandler` se muestreó cada 30 s en `logs/v3_7_memoria.tsv`.

Hubo un intento anterior, a las 14:44, que salió en 2,6 s con `codigo_python=1`. El árbol tenía `?? logs/` y el sello se negó:

```
oncogo_crc.utils.sello.ArbolSucio: el árbol de git no está limpio; no se firma la corrida
```

Ese intento no llegó a calcular ninguna pieza. El registro quedó en `logs/v3_7_intento_arbol_sucio.log`.

## Causa demostrada de la parada

El segundo intento, de un solo proceso, pasó H9 y H10 y se detuvo en la preanalítica. El código del proceso es 1. La traza es:

```
[v3.7] h9 en 9.4 s
[v3.7] h10 en 1.2 s
[v3.7] report_y_union en 2.8 s
[v3.7] preanalitica
...
File "src/oncogo_crc/stages/v3_5_preanalitica.py", line 158, in ejecutar
    viejo = pd.read_csv(previo, sep="\t")
...
pandas.errors.EmptyDataError: No columns to parse from file
codigo_python=1
```

`/usr/bin/time -l` cierra así: `858.27 real`, `maximum resident set size` 2 076 917 760 bytes, `peak memory footprint` 23 195 191 888. El disco `/Volumes/JANCODATA` se mantuvo en 400 866 560 KiB libres en todas las muestras.

El fichero que no se puede leer es `results/v3_7/preanalitica_intentos.tsv` en el volumen. Pesa 1 byte y su contenido es un salto de línea (`0a`). Lo dejó la corrida parcial anterior, cuando no hubo intentos nuevos y `volcar_tsv` escribió un marco vacío. `ejecutar` lo abre solo porque existe.

## Qué queda sin demostrar

Los cortes anteriores, con código 2 de `make` y sin traza, al anunciar H9 o H10, no se repitieron. Esta vez H9 terminó en 9,4 s y H10 en 1,2 s. Para esos cortes la causa sigue siendo no demostrada: el proceso de entonces no dejó un código de Python.

En esta prueba se descarta, para el proceso 28570:

- memoria agotada del cálculo: el residente máximo fue 2,1 GB y el proceso salió con la traza de pandas, no con 137 ni 9;
- disco lleno: quedaban unos 382 GB;
- desmontaje del volumen: las líneas `Bad volume: /Volumes/JANCODATA` son de Finder y `cache_delete` a las 14:59:43, después de la salida del proceso, y no coinciden con un fallo de escritura;
- un jetsam de este Python: las líneas `memorystatus` del mismo minuto nombran `SAMaintenance`, `AssetMetrics` y otras extensiones, no el pid 28570.

`ulimit -a` deja sin tope el tiempo de CPU, el tamaño de fichero y el espacio de direcciones. La pila es de 8176 KB. El núcleo de volcado es 0.
