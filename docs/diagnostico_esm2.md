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
