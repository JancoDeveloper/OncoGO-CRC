# v3.6 frente a v3.7

Cada valor lo formatea `src/oncogo_crc/formato.py` desde el fichero de esa corrida. Una fila `cambia` es una diferencia real. El resto se reprodujo.

| Macro | v3.6 | v3.7 | |
| --- | --- | --- | --- |
| \SueroN | 284 | 284 | igual |
| \SueroDif | 0 | 0 | igual |
| \SueroEq | 166 | 166 | igual |
| \SueroInd | 118 | 118 | igual |
| \EmbudoGrupos | 897 | 897 | igual |
| \EmbudoQ | 329 | 329 | igual |
| \EmbudoHumanos | 312 | 312 | igual |
| \MlTodas | 0,880 | 0,880 | igual |
| \MlCalidad | 0,915 | 0,915 | igual |
| \MlSinCalidad | 0,647 | 0,647 | igual |
| \MlSinCalidadP | 0,12 | 0,12 | igual |
| \MlHumana | 0,947 | 0,947 | igual |
| \MlCoagulacion | 0,818 | 0,818 | igual |
| \MlPca | 0,887 | 0,887 | igual |
| \MlContaminantes | 0,942 | 0,942 | igual |
| \TejidoNum | 147 | 147 | igual |
| \TejidoDen | 147 | 147 | igual |
| \TejidoMetaNum | 116 | 116 | igual |
| \TejidoMetaDen | 143 | 143 | igual |
| \HcincoEstado | Sostiene | Sostiene | igual |
| \HcincoUniverso | 29 | 29 | igual |
| \HcincoDif | 0,069 | 0,069 | igual |
| \HseisEstado | No sostiene | No sostiene | igual |
| \HsieteEstado | Refuta | Refuta | igual |
| \HsieteAlvez | 0,483 | 0,483 | igual |
| \HsieteHpa | 0,464 | 0,464 | igual |
| \HochoEstado | No evaluable | No evaluable | igual |
| \HochoUniverso | 33 | 53 | cambia |
| \HochoTumorales | 0 | 0 | igual |
| \HochoPlaquetasP | 0,00095 | 0,00095 | igual |
| \HnueveEstado | Exploratorio | Exploratorio | igual |
| \HdiezEstado | Descriptiva | Descriptiva | igual |
| \HdiezN | 5 | 5 | igual |
| \PreEstado | No evaluable | Calculado | cambia |
| \PreAuc | 0,000 | 0,000 | igual |
| \PreP | 2,3 × 10⁻¹⁰ | 1,5 × 10⁻¹¹ | cambia |
| \PreRhoPiel | −0,78 | −0,76 | cambia |
| \PreFecha | 54 | 60 | cambia |
| \CoberturaNum | 147 | 147 | igual |
| \CoberturaDen | 688 | 688 | igual |
| \HunoVeredicto | Rechaza H1 | Rechaza H1 | igual |
| \BhPrecursoresCn | no medido | 0,027 | cambia |
| \BhPielCn | no medido | 0,056 | cambia |
| \GoConNum | no medido | 16 | cambia |
| \GoConDen | no medido | 16 | cambia |
| \Priori5 | no medido | 20 | cambia |
| \Priori6 | no medido | 0 | cambia |
| \HsieteSec | 0,48591 | 0,47738 | cambia |
| \PreCrudos | 60 | 60 | igual |
| \CalidadCoagP | 0,00099 | 0,00099 | igual |
| \CalidadPielP | 3,4 × 10⁻⁶ | 3,4 × 10⁻⁶ | igual |
| \CertificadoEstado | Publicable | Publicable | igual |
| \MultiSuero | 288 | 288 | igual |
| \MultiMl | 72 | 72 | igual |
| \MultiDif | 13 | 13 | igual |
| n de Spearman, piel | 36 | 40 | cambia |
| rho de piel | −0,78 | −0,76 | cambia |
| n de Spearman, coagulacion | 36 | 40 | cambia |
| rho de coagulacion | −0,65 | −0,53 | cambia |

La causa de los cambios está en `docs/estado_v3_7.md`.
