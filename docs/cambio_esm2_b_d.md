# ESM-2, de B a D

B está en `results/archivo/v3_7_final_b/` del volumen. D está en `results/archivo/v3_7_final_d/`. Los números principales coinciden. Cambian el secundario de H7 y los sellos, estos por el commit y por el manifiesto.

| | AUC secundario | Incremento sobre el principal | Proteínas con ESM-2 |
| --- | --- | --- | --- |
| B | 0,4773843370379944 | −0,005427759110881281 | 181 de 668 con secuencia |
| D | 0,47376583096407343 | −0,009046265184802227 | 668 |

El AUC secundario de D menos el de B es −0,00361850607392095. El principal sigue en 0,48281209614887566 y el veredicto sigue en refuta. Fuentes: `h7_visibilidad_en_sangre.json` de cada corrida y `h7_tabla_proteinas.tsv`.

En B, `logs/v3_7_final_b.log` mide un pico de huella de 23 209 986 616 bytes. En D, `logs/v3_7_final_d.log` mide 1 431 275 728. La corrida ya no carga el modelo.

D y E son idénticas en los números. El detalle está en `docs/determinismo_d_e.md`.
