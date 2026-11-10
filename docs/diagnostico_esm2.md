# Diagnóstico de ESM-2

Lectura de `results/v3_7/h7_tabla_proteinas.tsv` (corrida B), de `results/archivo/v3_7_final_c/h7_tabla_proteinas.tsv` en el volumen (corrida C) y del FASTA que usa el código, `external/uniprot/human_reviewed_isoforms.fasta.gz`. Una proteína «tiene ESM-2» si `esm2_000` no está vacío. La longitud sale de `secuencias_por_gen`.

`h7_visibilidad_en_sangre.json` de B dice 668 proteínas con secuencia. La tabla tiene 688 filas. Hay secuencia para 668 y no la hay para 20, en las dos corridas.

| Corrida | Con vector | Sin vector, con secuencia | Sin secuencia |
| --- | --- | --- | --- |
| B | 181 | 487 | 20 |
| C | 196 | 472 | 20 |

En B se perdieron 487 de las 668 con secuencia. No es un puñado de proteínas largas.

En C y no en B, 15: ACTN1, ADH5, CKB, COL6A2, FABP2, FABP4, FOLR2, HMGA1, MTHFD2, NCAM1, P4HA1, PLA2G2A, SCG2, SLC25A4, SLC2A4. En B y no en C, ninguna.

Longitud de las que tienen secuencia:

| Grupo | n | Mediana | Mínimo | Máximo |
| --- | --- | --- | --- | --- |
| B, con vector | 181 | 444 | 61 | 3695 |
| B, sin vector | 487 | 452 | 61 | 5890 |
| C, con vector | 196 | 435,5 | 61 | 3695 |
| C, sin vector | 472 | 455 | 61 | 5890 |

La causa está en `embeddings_esm2`: el modelo va a MPS y un `RuntimeError` se salta con `continue`, igual que un desajuste de `oculto.shape[0]`. No queda registro.

Una sola ausente, sola, sí corre en MPS. Se probó TAGLN3 (199), ILK (452) y AHNAK (5890). El error aparece al encadenarlas. Tras AHNAK (5890) y AHNAK2 (5795), FCGBP (5405) falla:

```
RuntimeError: MPS backend out of memory (MPS allocated: 3.77 GiB, other allocations: 14.55 GiB, max allowed: 20.13 GiB). Tried to allocate 2.19 GiB on shared pool.
```

A partir de ahí fallan también las cortas (ILK, 452; HBB, 147), con el mismo mensaje y una petición de unos pocos KiB. El `continue` las deja fuera. Por eso B y C no traen el mismo subconjunto: depende de cuándo se llena la memoria.

## Cálculo en CPU

El script `scripts/v3/esm2_embeddings.py` corre en CPU con 4 hilos. Parte cada secuencia en ventanas de hasta 1022 residuos. El TSV queda en el volumen, `external/modelos/h7_esm2_embeddings.tsv`. El repositorio guarda el sello y el sha256.

Los dos procesos dan el mismo sha256 del TSV: `2b1e48e3d90476bd1cd9f2c8a2091b18be9b0dd119c61c87ac912729bd64e6b7`. No hizo falta repetir con 1 hilo.

| | Calculadas | Omitidas | Ventanas | Tiempo real | Memoria residente máxima | Pico de huella |
| --- | --- | --- | --- | --- | --- | --- |
| 1, `logs/esm2_1.log` | 668 | 20 | 832 | 106,19 s | 1 389 969 408 | 1 074 432 064 |
| 2, `logs/esm2_2.log` | 668 | 20 | 832 | 103,53 s | 1 488 437 248 | 1 235 912 840 |

Las 20 omisiones tienen el mismo motivo, «sin secuencia en el FASTA revisado»: C10ORF99, C14ORF159, C15ORF48, C15ORF52, C1ORF123, C2ORF88, CENPJ, FAM129A, H1F0, H3F3A, HIST1H2BC, HIST1H2BO, IGKV4-1, MARC2, METTL7A, MRVI1, SAA2-SAA4, SEPT10, SEPT6, WISP2. Ninguna se quedó fuera por un error del modelo.

`hw.memsize` es 17 179 869 184. La residente se quedó cerca de 1,5 GB, por debajo de la mitad. Torch 2.15.0 y transformers 5.22.2. El sello del segundo cálculo, que es el que lee la corrida, está en `external/modelos/h7_esm2_embeddings.tsv.sello.json`.
