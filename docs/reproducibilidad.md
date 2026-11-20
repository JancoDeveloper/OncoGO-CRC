# Reproducibilidad

La corrida de referencia quedó sellada el 11 de noviembre de 2026.

| Campo | Valor |
| --- | --- |
| Commit del sello | `71b310a3fcaebe9267fac971156df3aab491e009` |
| Etiqueta | `v3.7-final` → `595118ae09ea6006c744c86b28cb0a10d42e0442` |
| Python | 3.14.3 |
| numpy, pandas, scipy, scikit-learn | 2.4.4, 3.0.2, 1.17.1, 1.8.0 |
| Inicio y fin | 2026-11-11T00:48:45Z a 2026-11-11T01:04:00Z |
| Certificado | `publicable`, 22 casillas |
| Semilla base | 20260914 |

`PYTHONHASHSEED=0` y un hilo en OpenMP, OpenBLAS y MKL los fija el Makefile.
La etapa vuelve a evaluar el modelo M2 y exige que el AUC coincida a 1e-12.
Eso queda en `results/v3_7/determinismo.json`.

Dos ejecuciones completas en el equipo del autor se comparan en
`docs/reproducibilidad/determinismo_d_e.md` y en el JSON del mismo directorio.
`make verificar` repite esa comparación contra `results/v3_7`:

- las 35 filas de `config/esperado_v3_7.tsv`, con tolerancia 0 en conteos y
  1e-9 en AUC y en p;
- SHA-256 archivo por archivo;
- se admiten los sellos (`*.sello.json`), `corrida.json`, `tests.txt`,
  `procedencia.json`, `sello_corrida.json`, `lecturas_v3_7.json` y el orden
  de filas de `gsea_tejido.tsv`.

`results/v3_7/tests.txt` conserva rutas del equipo en el que se firmó la
corrida. No se edita: forma parte de la referencia. Una corrida nueva escribe
su propio `tests.txt` en el workspace.

`make descargar` trae, desde el depósito
[10.5281/zenodo.24398807](https://doi.org/10.5281/zenodo.24398807)
(CC-BY-4.0), y comprueba el SHA-256:

1. La matriz TMT `data/raw/cptac/cptac_crc_matrix.tsv`, derivada del estudio
   de acceso abierto PDC000116 (Vasaikar et al. 2019) y redistribuida con
   atribución. Huella `a1272baf…`.
2. Los tres ficheros certificados de v1. `make certified` también puede
   regenerarlos cuando esa matriz ya está.
3. El TSV de embeddings ESM-2. El modelo principal no lo usa.

Lo que todavía no puede repetirse solo con este árbol:

1. Las fechas de adquisición de los 60 crudos, que piden los `.raw` y
   ThermoRawFileParser.

El registro `corrida.json` de la referencia no guarda CPU, RAM ni sistema
operativo. Esta versión no altera ese archivo sellado.
