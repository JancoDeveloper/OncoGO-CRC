# Entrega del 11 de noviembre de 2026

El archivo para entregar es `tesis-ucaldas/word/tesis-maestria-ucaldas.docx`. La hoja «Qué completar antes de exportar el PDF» está después de la portada. Los pasos de formato están en `tesis-ucaldas/word/INSTRUCCIONES-WORD.md`.

## Lo que solo puede completar el autor

En la portada siguen en gris:

- Director.
- Codirector.
- Grupo de investigación.
- Línea de investigación.

El repositorio no tiene esos cuatro nombres. La dedicatoria y los agradecimientos siguen fuera: la guía de 2022 los marca como opcionales.

No exporte `JuanDiego_PolancoFranco_2026.pdf` mientras un corchete siga visible. Después de escribirlos, actualice la tabla de contenido, la lista de tablas y la lista de figuras.

## Cifras de la corrida firmada

Salen de `results/v3_7/`. El certificado es publicable. El sello es `71b310a3fcaebe9267fac971156df3aab491e009`, 2026-11-11T00:48:45Z. La etiqueta `v3.7-final` está en `595118ae09ea6006c744c86b28cb0a10d42e0442` y no se movió.

| Qué | Cifra |
|---|---|
| Celda del suero, CP frente a CN | 284 proteínas: 0 diferentes, 166 equivalentes, 118 indeterminadas |
| Signo del tumor primario, Tanaka y CPTAC | 147 de 147 |
| Modelo, todas las proteínas | AUC 0,880 |
| Modelo, solo calidad de la muestra | AUC 0,915 |
| Modelo, sin la calidad | AUC 0,647; p de permutación 0,12 |
| H5 | sostiene; 624 proteínas; universo 29; fracción pan-carcinoma 0,119 |
| H6 | no sostiene |
| H7 | refuta. Álvez 2023: AUC 0,483 (0,381 a 0,584; 53 y 73). HPA: AUC 0,464 (0,353 a 0,573; 53 y 63). Papier 2024: no evaluable |
| H7, modelo secundario | AUC 0,474; incremento −0,009. No entra en el veredicto. Vector en 668; se omiten 20 |
| H8 | no evaluable: 0 tumorales epiteliales en un universo de 53 |
| H9 | exploratorio |
| H10 | descriptiva, 5 proteínas |
| Preanalítica | calculado; 60 crudos con fecha; CP, CN y HC en bloques |
| Priorización | 688, 135, 135, 118, 20, 0 y 0. Nadie pasa el paso 6 |

Esas cifras están en el cuerpo y en el apartado 6.7 del Word. `results/v3/` a `results/v3_6/` quedan como historial.

## Lo que no es un resultado

No se ejecutaron las vías del suero, la cláusula de vías de H1, la unión de HC por XIC, la deriva del PRTC, las redes de coexpresión, el catálogo de candidatas ni los neoantígenos. H8 no evaluable no dice que el estadio IV no lleve proteína del tumor a la sangre. Cero proteínas diferentes frente a CN no dice que el tumor no cambie en el tejido. Ningún resultado de esta tesis sirve para decidir sobre un paciente.
