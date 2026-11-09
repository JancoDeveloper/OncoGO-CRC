# Para sustentar OncoGO-CRC

La página con gráficas y fichas está en [divulgacion.html](divulgacion.html). Esta nota es lo que hay que poder decir.

## El resultado, en una frase

Entre las 688 proteínas certificadas, el primario coincide en signo con CPTAC en 89 de 90. En el suero oficial, de 205 proteínas, el corte de 1,5 veces no separa el cáncer de la colonoscopia negativa: 105 tienen ese cambio descartado y 100 quedan indeterminadas.

Sirve como tesis, porque la pregunta queda respondida con un control de otra enfermedad intestinal. No sirve como prueba de un marcador para la consulta.

## Tres cohortes

| Cohorte | Qué es | Qué responde |
|---|---|---|
| PDC000116, Vasaikar, *Cell* 2019 | Tejido tumoral y mucosa, TMT, ya en log2 | El ancla. El tumor tiene un proteoma distinto. |
| PXD034575, Tanaka, *Cell Reports* 2024 | 16 pacientes. P primario, M metástasis hepática, N mucosa | Si ese proteoma se repite en otro tejido. En las 688 certificadas, corrida v3.2: 89 de 90 con el primario, 69 de 82 con la metástasis. El cierre del 31 de octubre, en otro universo, es 687 de 874 y correlación 0,83. |
| PXD055218, Rai, *J Proteome Res* 2025 | Suero DIA, 20 cáncer, 20 colonoscopia negativa, 20 sanos. Contraste oficial: 205 proteínas | Frente a la colonoscopia negativa el corte de 1,5 veces deja 0 proteínas. De 205, 105 son equivalentes y 100 indeterminadas. |

## Orden de una tarde

1. Decir el resultado en voz alta, sin la hoja.
2. Dibujar las tres cohortes y la pregunta de cada una.
3. Explicar LBP (pasa el FDR contra sanos, no contra la enfermedad) y el AUC vigente: SERPINA3, 0,865, p = 0,0464, clase indeterminada. BAIAP3 (AUC 0,82, FDR 0,20) es el catálogo del 31 de octubre y no se recalculó.
4. Explicar P, M y N, y por qué el cero de MaxQuant es una ausencia.
5. Explicar que 89/90, 687/874 y 0,104 son tres cuentas distintas. Ninguna mide tejido contra sangre.
6. Cerrar con el límite: 20 por grupo, sin el inmunoensayo de CEA, sin secuencia de péptido mutado.

## Palabras

- **log2FC.** Cambio en logaritmo base 2. 1 es el doble. El umbral 0,585 es log2(1,5).
- **FDR.** Fracción esperada de falsos entre lo que se declara distinto. Corte 0,05, Benjamini–Hochberg.
- **treat.** La t moderada pregunta si el cambio supera 0,585, no si es distinto de cero.
- **AUC.** Si la proteína ordena al cáncer por encima de la otra condición. 0,50 es azar. No dice si sube o baja y no sustituye al FDR.
- **Pareado.** En Tanaka, tumor y mucosa son del mismo paciente. En el suero, los grupos son personas distintas.
- **Grupo proteico.** La fila del DIA no guarda la secuencia mutada de un péptido.
- **Fase aguda.** Respuesta del organismo. C9, LRG1, SERPINA1, CRP y HP ya estaban marcadas así.

## Si lo preguntan

**¿Cero proteínas es un fracaso?** Es el resultado. El umbral se fijó antes de mirar el suero.

**¿Por qué no basta el sano?** Contra sanos suben SAA1 y LBP. Eso mide respuesta del organismo. La colonoscopia negativa es el control que decide la especificidad.

**¿SAA1 es un marcador?** Separa el cáncer de los sanos. No separa el cáncer de la otra enfermedad intestinal.

**¿BAIAP3?** En el catálogo del 31 de octubre, candidata de Rai: AUC 0,82 y FDR 0,20. Esa tabla no se recalculó en v3.2. No es un hallazgo de esta tesis.

**¿CHMP4A?** No entra en las 312 de v3.2: su q global de grupo es 0,037. El log2FC 5,42 y las 16 de 40 observaciones son la matriz sin filtrar del 31 de octubre. No se presenta como marcador.

**¿Neoantígenos?** Esta matriz no puede confirmarlos. Rai, Xiang y Hembrough describen péptidos. Aquí cada fila es un grupo de proteínas.

**¿Por qué no se fusionaron las cohortes?** Tejido TMT, tejido LFQ y suero DIA no son la misma medición ni las mismas personas.

## No se dice

- Que hay un marcador de sangre listo para la consulta.
- Que la columna P es un pólipo.
- Que 0,83 sustituye a 0,104.
- Que las vías del tumor se ven en este suero.
