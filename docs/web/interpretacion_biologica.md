> Documento de la v1. No forma parte de los resultados certificados de la tesis; ver results/v3_7.

# Interpretación biológica de los módulos detectados

Esta página **traduce los términos GO/GSEA significativos a contexto biológico de cáncer colorrectal**, agrupando los _leading edges_ por familia funcional y emparejándolos con un catálogo curado a partir de literatura proteogenómica reciente (Cancer Cell 2023, Nature 2014, PLOS Comp Biol 2024, García-Foncillas et al. 2026).

_Los textos siguientes son **hipótesis generadas automáticamente**, no diagnóstico. La discusión final de la tesis debe contrastar cada módulo con la literatura específica._

Para definiciones de NES, leading edge, ORA vs GSEA, etc. ver el [Glosario](glosario.md).

---

### Tumor vs tejido adyacente normal (CPTAC)

**Módulos biológicos detectados** (leading edges con FDR < 0.25 agregadas y emparejadas con catálogo curado de CRC):

- **Traducción mitocondrial (ribosoma 39S/28S, MRPL/MRPS)** — `64` genes en leading edges. Ejemplos: _DAP3, GADD45GIP1, GATC, MRPL10, MRPL11, MRPL13, MRPL14, MRPL15, MRPL16, MRPL17, MRPL19, MRPL2_.
  > Aumento de la maquinaria de **traducción mitocondrial** (subunidades MRPL/MRPS) en tumor se ha asociado a reprogramación bioenergética en CCR: las células tumorales preservan o amplifican la biogénesis mitocondrial incluso cuando bajan los complejos OXPHOS terminales (NDUF/UQCR). Es coherente con el modelo de **plasticidad metabólica** (no Warburg estricto) descrito en proteogenómica reciente del CPTAC.

- **OXPHOS / cadena respiratoria (complejos I–V, NDUF/UQCR/COX/ATP5)** — `29` genes en leading edges. Ejemplos: _ATP5F1C, ATP5F1E, ATP5ME, ATP5MG, ATP5PB, ATP5PD, ATP5PF, ATP5PO, NDUFA1, NDUFA11, NDUFA12, NDUFA13_.
  > Disminución del **complejo I (NDUFs)**, complejo III (UQCR), citocromo c oxidasa (COX) y ATP sintasa mitocondrial en tumor es un patrón clásico en CCR: indica desacople parcial de fosforilación oxidativa y favorece glucólisis aeróbica (efecto Warburg parcial). Particularmente prominente en CMS3 metabólico (Guinney et al., 2015).

- **Proteasoma 26S (subunidades PSMA/PSMB/PSMC/PSMD, reguladores)** — `25` genes en leading edges. Ejemplos: _PSMA1, PSMA3, PSMA4, PSMA5, PSMA6, PSMA7, PSMB10, PSMB7, PSMC1, PSMC2, PSMC3, PSMC4_.
  > Sobre-expresión de **subunidades del proteasoma 26S** (PSMA/PSMB catalíticas, PSMC/PSMD regulatorias) en tumor refleja mayor recambio proteico para sostener proliferación. Es una vulnerabilidad farmacológica documentada (bortezomib y análogos) pero el éxito en CCR es limitado en monoterapia — útil como hipótesis de combinación, **nunca como recomendación** clínica desde un estudio computacional.

- **Matriz extracelular / mesénquima / músculo liso (CMS4 mesenquimal)** — `23` genes en leading edges. Ejemplos: _AHNAK, AHNAK2, CAV1, CAV2, COL4A5, COL6A2, DAG1, DES, DMD, DTNA, FBN1, FLNA_.
  > Disminución de proteínas de **matriz extracelular**, **citoesqueleto muscular liso** y **complejo distrofina-glicoproteína** en tumor refleja pérdida del compartimento estromal y epitelio de criptas en favor del epitelio tumoral expandido. Cuando estas vías están **aumentadas** en tumor son típicas de **CMS4 mesenquimal** (TGF-β, angiogénesis), con peor pronóstico (Guinney et al., 2015; García-Foncillas et al., 2026). Revisar siempre el signo del NES antes de interpretar.

- **Procesamiento de mRNA (poliadenilación 3', factores CPSF/CSTF)** — `21` genes en leading edges. Ejemplos: _CDC73, CLP1, CPSF1, CPSF2, CPSF3, CPSF4, CPSF6, CPSF7, CSTF1, CSTF2, MTPAP, NCBP1_.
  > Activación coordinada del **complejo de poliadenilación 3'** (CPSF/CSTF) en tumor sugiere regulación alternativa del extremo 3' UTR — un mecanismo conocido de **APA (alternative polyadenylation)** que en cáncer suele acortar 3' UTRs y eliminar sitios miRNA, aumentando estabilidad de mRNAs oncogénicos. Reportado en proteogenómica CPTAC pan-cancer.

- **Citoesqueleto de actina y motilidad (filaminas, septinas, ARPC)** — `8` genes en leading edges. Ejemplos: _FLNA, NEXN, PALLD, SYNPO, SYNPO2, TPM1, TPM2, TPM4_.
  > Reorganización del **citoesqueleto de actina** (filaminas, septinas, complejo ARP2/3) es central en invasión y migración celular. Cambios direccionales aquí son hipótesis razonables de **fenotipo invasivo / EMT parcial**, especialmente si coocurren con ECM/CMS4.

- **Chaperonas y proteostasis (HSP70/HSP90, CCT/TRiC, PDI/ERP)** — `4` genes en leading edges. Ejemplos: _BAG3, BAG6, DNAJA3, DNAJB4_.
  > Activación de **chaperonas citosólicas (HSP70/HSP90, complejo CCT/TRiC)** y de la **PDI/ERp** del retículo endoplásmico en tumor indica estrés proteostático sostenido — coherente con alta tasa de síntesis de proteínas oncogénicas y respuesta a unfolded protein response (UPR).

_Estos módulos son **hipótesis cuantitativas generadas por el pipeline**: requieren validación experimental y deben citarse junto a la literatura específica (Cancer Cell 2023, PLOS Comp Biol 2024, Nature 2014, García-Foncillas et al. 2026) antes de cualquier afirmación clínica._

### Plasma CRC vs control sano

**Módulos biológicos detectados** (leading edges con FDR < 0.25 agregadas y emparejadas con catálogo curado de CRC):

- **Proteínas de fase aguda / cascada del complemento / lipoproteínas (plasma)** — `28` genes en leading edges. Ejemplos: _APOB, APOE, C1RL, C5, C6, C9, CFB, CFH, CRP, F11, F9, FGA_.
  > Elevación de **proteínas de fase aguda** (LBP, SAA1/2, CRP, SERPINAs, ORM, HP) y activación de **complemento** (C3, C9, CFB, CFH) en plasma CRC es una **respuesta inflamatoria sistémica** bien documentada (Nature Commun 2024; Cell Rep Med 2025). No es específica de CRC pero es una **firma robusta de cáncer activo** útil para discusión exploratoria de biopsia líquida. Las apolipoproteínas (APOA1/4) frecuentemente **bajan** por reprogramación lipídica hepática.

- **Matriz extracelular / mesénquima / músculo liso (CMS4 mesenquimal)** — `4` genes en leading edges. Ejemplos: _AHNAK, COL4A2, FN1, VCAM1_.
  > Disminución de proteínas de **matriz extracelular**, **citoesqueleto muscular liso** y **complejo distrofina-glicoproteína** en tumor refleja pérdida del compartimento estromal y epitelio de criptas en favor del epitelio tumoral expandido. Cuando estas vías están **aumentadas** en tumor son típicas de **CMS4 mesenquimal** (TGF-β, angiogénesis), con peor pronóstico (Guinney et al., 2015; García-Foncillas et al., 2026). Revisar siempre el signo del NES antes de interpretar.

_Estos módulos son **hipótesis cuantitativas generadas por el pipeline**: requieren validación experimental y deben citarse junto a la literatura específica (Cancer Cell 2023, PLOS Comp Biol 2024, Nature 2014, García-Foncillas et al. 2026) antes de cualquier afirmación clínica._

### PRIDE validación externa: tumor vs normal

**Módulos biológicos detectados** (leading edges con FDR < 0.25 agregadas y emparejadas con catálogo curado de CRC):

- **OXPHOS / cadena respiratoria (complejos I–V, NDUF/UQCR/COX/ATP5)** — `21` genes en leading edges. Ejemplos: _ATP5F1D, ATP5F1EP2, ATP5ME, ATP5MG, ATP5PD, ATP5PO, COX4I1, COX5A, COX5B, NDUFA13, NDUFA2, NDUFA8_.
  > Disminución del **complejo I (NDUFs)**, complejo III (UQCR), citocromo c oxidasa (COX) y ATP sintasa mitocondrial en tumor es un patrón clásico en CCR: indica desacople parcial de fosforilación oxidativa y favorece glucólisis aeróbica (efecto Warburg parcial). Particularmente prominente en CMS3 metabólico (Guinney et al., 2015).

- **Chaperonas y proteostasis (HSP70/HSP90, CCT/TRiC, PDI/ERP)** — `18` genes en leading edges. Ejemplos: _AHSA1, BAG3, CANX, CCT6B, CCT8, CLU, DNAJA1, DNAJB5, ERP44, GRPEL1, HSP90AA1, HSP90AB2P_.
  > Activación de **chaperonas citosólicas (HSP70/HSP90, complejo CCT/TRiC)** y de la **PDI/ERp** del retículo endoplásmico en tumor indica estrés proteostático sostenido — coherente con alta tasa de síntesis de proteínas oncogénicas y respuesta a unfolded protein response (UPR).

- **Citoesqueleto de actina y motilidad (filaminas, septinas, ARPC)** — `17` genes en leading edges. Ejemplos: _ACTG2, ARPC5, FLNA, GSN, MARCKS, MYO1A, MYO1C, NEXN, PALLD, PFN1, SCIN, SEPTIN2_.
  > Reorganización del **citoesqueleto de actina** (filaminas, septinas, complejo ARP2/3) es central en invasión y migración celular. Cambios direccionales aquí son hipótesis razonables de **fenotipo invasivo / EMT parcial**, especialmente si coocurren con ECM/CMS4.

- **Proteasoma 26S (subunidades PSMA/PSMB/PSMC/PSMD, reguladores)** — `11` genes en leading edges. Ejemplos: _ADRM1, PSMA1, PSMA4, PSMB10, PSMB3, PSMB5, PSMD1, PSMD11, PSMD2, PSMD6, PSMD7_.
  > Sobre-expresión de **subunidades del proteasoma 26S** (PSMA/PSMB catalíticas, PSMC/PSMD regulatorias) en tumor refleja mayor recambio proteico para sostener proliferación. Es una vulnerabilidad farmacológica documentada (bortezomib y análogos) pero el éxito en CCR es limitado en monoterapia — útil como hipótesis de combinación, **nunca como recomendación** clínica desde un estudio computacional.

- **Matriz extracelular / mesénquima / músculo liso (CMS4 mesenquimal)** — `11` genes en leading edges. Ejemplos: _ACTG2, AHNAK, COL14A1, COL6A2, FLNA, FLNC, SGCD, SYNM, TAGLN3, TLN1, TLN2_.
  > Disminución de proteínas de **matriz extracelular**, **citoesqueleto muscular liso** y **complejo distrofina-glicoproteína** en tumor refleja pérdida del compartimento estromal y epitelio de criptas en favor del epitelio tumoral expandido. Cuando estas vías están **aumentadas** en tumor son típicas de **CMS4 mesenquimal** (TGF-β, angiogénesis), con peor pronóstico (Guinney et al., 2015; García-Foncillas et al., 2026). Revisar siempre el signo del NES antes de interpretar.

_Estos módulos son **hipótesis cuantitativas generadas por el pipeline**: requieren validación experimental y deben citarse junto a la literatura específica (Cancer Cell 2023, PLOS Comp Biol 2024, Nature 2014, García-Foncillas et al. 2026) antes de cualquier afirmación clínica._

### Estadio avanzado vs temprano (si stage está disponible)

**Módulos biológicos detectados** (leading edges con FDR < 0.25 agregadas y emparejadas con catálogo curado de CRC):

- **Traducción mitocondrial (ribosoma 39S/28S, MRPL/MRPS)** — `27` genes en leading edges. Ejemplos: _GADD45GIP1, MRPL11, MRPL12, MRPL15, MRPL17, MRPL2, MRPL20, MRPL22, MRPL27, MRPL3, MRPL30, MRPL32_.
  > Aumento de la maquinaria de **traducción mitocondrial** (subunidades MRPL/MRPS) en tumor se ha asociado a reprogramación bioenergética en CCR: las células tumorales preservan o amplifican la biogénesis mitocondrial incluso cuando bajan los complejos OXPHOS terminales (NDUF/UQCR). Es coherente con el modelo de **plasticidad metabólica** (no Warburg estricto) descrito en proteogenómica reciente del CPTAC.

- **Matriz extracelular / mesénquima / músculo liso (CMS4 mesenquimal)** — `15` genes en leading edges. Ejemplos: _COL10A1, COL14A1, COL15A1, COL16A1, COL17A1, COL18A1, COL22A1, COL3A1, COL4A5, COL8A1, COLEC12, DAG1_.
  > Disminución de proteínas de **matriz extracelular**, **citoesqueleto muscular liso** y **complejo distrofina-glicoproteína** en tumor refleja pérdida del compartimento estromal y epitelio de criptas en favor del epitelio tumoral expandido. Cuando estas vías están **aumentadas** en tumor son típicas de **CMS4 mesenquimal** (TGF-β, angiogénesis), con peor pronóstico (Guinney et al., 2015; García-Foncillas et al., 2026). Revisar siempre el signo del NES antes de interpretar.

_Estos módulos son **hipótesis cuantitativas generadas por el pipeline**: requieren validación experimental y deben citarse junto a la literatura específica (Cancer Cell 2023, PLOS Comp Biol 2024, Nature 2014, García-Foncillas et al. 2026) antes de cualquier afirmación clínica._

### Hispano/latino vs no hispano/latino (etnicidad NIH; contexto US — ver clinical_demographics.md)

**Módulos biológicos detectados** (leading edges con FDR < 0.25 agregadas y emparejadas con catálogo curado de CRC):

- **Proteínas de fase aguda / cascada del complemento / lipoproteínas (plasma)** — `23` genes en leading edges. Ejemplos: _APOB, APOE, C1R, C2, C3, C5, C6, C7, C9, CFB, CFH, CFHR3_.
  > Elevación de **proteínas de fase aguda** (LBP, SAA1/2, CRP, SERPINAs, ORM, HP) y activación de **complemento** (C3, C9, CFB, CFH) en plasma CRC es una **respuesta inflamatoria sistémica** bien documentada (Nature Commun 2024; Cell Rep Med 2025). No es específica de CRC pero es una **firma robusta de cáncer activo** útil para discusión exploratoria de biopsia líquida. Las apolipoproteínas (APOA1/4) frecuentemente **bajan** por reprogramación lipídica hepática.

- **Traducción mitocondrial (ribosoma 39S/28S, MRPL/MRPS)** — `21` genes en leading edges. Ejemplos: _DAP3, MRPL1, MRPL14, MRPL16, MRPL20, MRPL24, MRPL28, MRPL3, MRPL32, MRPL51, MRPL55, MRPS11_.
  > Aumento de la maquinaria de **traducción mitocondrial** (subunidades MRPL/MRPS) en tumor se ha asociado a reprogramación bioenergética en CCR: las células tumorales preservan o amplifican la biogénesis mitocondrial incluso cuando bajan los complejos OXPHOS terminales (NDUF/UQCR). Es coherente con el modelo de **plasticidad metabólica** (no Warburg estricto) descrito en proteogenómica reciente del CPTAC.

- **Chaperonas y proteostasis (HSP70/HSP90, CCT/TRiC, PDI/ERP)** — `6` genes en leading edges. Ejemplos: _DNAJA1, HSP90B1, HSPA14, HSPA1A, HSPA1L, PPIA_.
  > Activación de **chaperonas citosólicas (HSP70/HSP90, complejo CCT/TRiC)** y de la **PDI/ERp** del retículo endoplásmico en tumor indica estrés proteostático sostenido — coherente con alta tasa de síntesis de proteínas oncogénicas y respuesta a unfolded protein response (UPR).

- **Matriz extracelular / mesénquima / músculo liso (CMS4 mesenquimal)** — `5` genes en leading edges. Ejemplos: _ACTG2, FLNA, FN1, MMP14, TLN1_.
  > Disminución de proteínas de **matriz extracelular**, **citoesqueleto muscular liso** y **complejo distrofina-glicoproteína** en tumor refleja pérdida del compartimento estromal y epitelio de criptas en favor del epitelio tumoral expandido. Cuando estas vías están **aumentadas** en tumor son típicas de **CMS4 mesenquimal** (TGF-β, angiogénesis), con peor pronóstico (Guinney et al., 2015; García-Foncillas et al., 2026). Revisar siempre el signo del NES antes de interpretar.

- **Citoesqueleto de actina y motilidad (filaminas, septinas, ARPC)** — `5` genes en leading edges. Ejemplos: _ACTA1, ACTB, ACTG2, FLNA, TLN1_.
  > Reorganización del **citoesqueleto de actina** (filaminas, septinas, complejo ARP2/3) es central en invasión y migración celular. Cambios direccionales aquí son hipótesis razonables de **fenotipo invasivo / EMT parcial**, especialmente si coocurren con ECM/CMS4.

_Estos módulos son **hipótesis cuantitativas generadas por el pipeline**: requieren validación experimental y deben citarse junto a la literatura específica (Cancer Cell 2023, PLOS Comp Biol 2024, Nature 2014, García-Foncillas et al. 2026) antes de cualquier afirmación clínica._

