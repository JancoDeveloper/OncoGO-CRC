# Determinismo D frente a E, regla estricta

D: `/Volumes/JANCODATA/OncoGO-CRC/results/archivo/v3_7_final_d`. E: `/Volumes/JANCODATA/OncoGO-CRC/results/v3_7`.

Se quita el bloque `sello` y las claves `fecha_utc`, `pid`, `inicio` y `fin`. El resto se compara tal cual, sha256 incluidos. Si un sha256 difiere, el fichero al que apunta solo puede diferir en esas claves o en otro sha256 explicado de la misma forma. `gsea_tejido.tsv` se compara ordenado por cohorte y término. En `tests.txt` solo se normalizan los tiempos y los pid.

Diferencias explicadas por fecha o proceso: 165.

Diferencias sin explicar: ninguna.

D y E coinciden con esta regla.
