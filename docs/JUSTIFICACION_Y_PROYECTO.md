# OncoGO-CRC — justificación, qué hace y qué busca

> **Encuadre (invariante).** OncoGO-CRC es una **metodología/herramienta de análisis computacional
> reproducible**. Su producto son **hipótesis priorizadas con evidencia convergente**, no dianas
> terapéuticas validadas ni afirmaciones clínicas. Todo hallazgo es **exploratorio** y requiere validación
> experimental independiente (fuera del alcance de una tesis de maestría computacional).

Este documento explica **por qué existe el proyecto**, **qué hace** y **qué busca**, con soporte
bibliográfico. Es la base narrativa para el marco teórico y la introducción de la tesis.

---

## 1. Motivación

### 1.1 Relevancia del problema
El cáncer colorrectal (CCR) es una de las principales causas de mortalidad oncológica a nivel mundial y muestra
una **incidencia creciente en población joven** (early-onset, EOCRC), con disparidades relevantes en poblaciones
latinoamericanas. Esta motivación tiene además un componente personal: el impacto directo de la enfermedad en la
vida del autor, que orienta el trabajo hacia una **utilidad real** para la investigación de la enfermedad, no
hacia un ejercicio meramente académico.

### 1.2 Vacío científico que aborda
La literatura reciente identifica cuatro problemas que este proyecto ataca de forma computacional:

1. **Heterogeneidad molecular del CCR** (subtipos CMS) no capturada por biomarcadores genómicos rutinarios
   (Guinney et al. 2015, *Nat Med*).
2. **Desacople proteoma–transcriptoma**: la proteómica aporta una capa funcional que el ARN no siempre refleja
   (Zhang et al. 2014, *Nature*; consorcio CPTAC).
3. **Fragmentación de los datos públicos**: cohortes proteómicas dispersas (CPTAC, PRIDE) difíciles de integrar
   por heterogeneidad de plataforma (TMT/LFQ/DIA) e identificadores.
4. **Necesidad de reuso riguroso**: la integración con armonización de identificadores y control de
   reproducibilidad permite recuperar firmas consistentes entre cohortes (Robles et al. 2024, *PLOS Comput Biol*).

### 1.3 Contexto traslacional (Latinoamérica / Colombia)
La revisión de García-Foncillas et al. (2026) subraya la heterogeneidad clínica del CCR, el aumento en población
joven, la relevancia de MSI y subtipos moleculares, y el potencial —pero también las limitaciones— de la biopsia
líquida, desde una perspectiva prudente y contextualizada para la región.

---

## 2. Qué hace (capacidades)

OncoGO-CRC integra proteómica pública de CCR y ejecuta, de forma reproducible y trazable:

- **Ingesta y control de calidad** de matrices proteómicas (CPTAC tejido TMT10, PRIDE tejido LFQ, plasma DIA-NN exploratorio).
- **Armonización de identificadores** proteína→gen con UniProt (mapa maestro + bulk humano revisado con alias).
- **Análisis diferencial** (t-test de Welch/pareado + control de FDR por Benjamini–Hochberg).
- **Reproducibilidad multi-cohorte**: firma núcleo por puente UniProt con concordancia de signo entre cohortes.
- **Enriquecimiento funcional**: GO (ORA + GSEA pre-ranked con anotación GOA humana real) y KEGG.
- **Redes de coexpresión** y detección de hubs.
- **ML interpretable**: L1/ElasticNet/RandomForest con validación cruzada anidada, SHAP y bootstrap de estabilidad.
- **Score de priorización traslacional (TPS)** y evidencia externa (concordancia TCGA RNA-seq / Cox).
- **Trazabilidad completa**: `run_manifest.json` con hashes SHA-256 y (2.0) versión de cada base de datos.
- **Publicación web** reproducible (MkDocs) y (2.0) dashboard interactivo.

La evolución **2.0** añade tres capas metodológicas: homogeneización/integración multi-ómica, análisis de red
(módulos), y priorización por evidencia convergente. Ver `ARQUITECTURA_ONCOGO_2.0.md` y
`DISENO_TECNICO_ONCOGO_2.0.md`.

---

## 3. Qué busca (objetivos)

### 3.1 Objetivo general
Desarrollar una **herramienta metodológica reproducible y reutilizable** que integre, homogeneíce y analice
proteómica pública de CCR y **priorice candidatos como hipótesis con evidencia convergente multi-fuente**, con
interpretabilidad y trazabilidad FAIR.

### 3.2 Objetivos específicos
1. Armonizar identificadores y **homogeneizar** cohortes heterogéneas (TMT/LFQ/DIA) evaluando la corrección con
   métricas de mezcla de lotes y **preservación de señal biológica**.
2. Identificar una **firma proteica reproducible** entre cohortes independientes (concordancia de signo/magnitud).
3. Interpretar la firma en **contexto de red/módulos** (biología de sistemas), no como proteínas aisladas.
4. **Priorizar candidatos** integrando evidencia pública objetiva (esencialidad, tractabilidad, especificidad
   tumoral, contexto genómico) como **marco de convergencia**, no como recomendación clínica.
5. Garantizar **reproducibilidad e interpretabilidad** (manifest, versiones, SHAP/explicabilidad, estándares
   FAIR/DOME/TRIPOD+AI) y **publicar los procesos y resultados en la web** para consulta abierta.

---

## 4. Alcance y NO-objetivos (prudencia traslacional)

**Sí busca:** una metodología, una firma reproducible, un marco de priorización de hipótesis, y una plataforma
consultable y reutilizable.

**No busca (explícitamente):**
- Proponer biomarcadores o dianas listos para uso clínico.
- Sustituir guías terapéuticas ni realizar predicciones clínicas individuales.
- Afirmar causalidad biológica a partir de asociaciones computacionales.

Los resultados se interpretan como **hipótesis generadoras**, en línea con la prudencia traslacional recomendada
en ómicas (Nature 2013; estándares TRIPOD+AI, REMARK, DOME).

---

## 5. Contribución esperada

1. **Metodológica:** un pipeline FAIR con trazabilidad (manifest + versiones de datos), armonización canónica y
   evaluación honesta de reproducibilidad.
2. **De evidencia:** una firma reproducible con concordancia externa cuantificada (no una sola cohorte).
3. **Interpretativa:** mapeo de la firma a biología de subtipo (p. ej. CMS4/estroma) con literatura 2024–2026.
4. **De utilidad:** una **tabla de hipótesis priorizadas** con evidencia convergente y una **web consultable**,
   como base para futuras líneas de investigación (validación experimental, cohortes regionales).

---

## 6. Referencias base

- Guinney J. et al. (2015). The consensus molecular subtypes of colorectal cancer. *Nat Med* 21:1350.
- Zhang B. et al. (2014). Proteogenomic characterization of human colon and rectal cancer. *Nature* 513:382.
- Robles J., Prakash A., Vizcaíno J. A. y Casal J. I. (2024). Integrated meta-analysis of colorectal cancer public proteomic datasets. *PLOS Comput Biol*. doi:10.1371/journal.pcbi.1011828
- Vasaikar S. et al. (2019). Proteogenomic analysis of human colon cancer. *Cell* 177:1035. doi:10.1016/j.cell.2019.03.030. Es la ficha de PDC000116. El pan-cáncer de CPTAC en *Cancer Cell* 2023 (doi:10.1016/j.ccell.2023.08.014) es otro artículo.
- García-Foncillas J. (2026). Biología molecular del cáncer colorrectal. *Rev. colomb. hematol. oncol.* doi:10.51643/22562915.853. Revisión narrativa de vías y biopsia líquida de ADN. No mide el proteoma ni define un subtipo CMS en estas muestras.
- Colección metodológica y de bases de datos: ver `HERRAMIENTAS_Y_TECNICAS.md` y `DISENO_TECNICO_ONCOGO_2.0.md`.
