#!/usr/bin/env python3
"""Genera el mazo de Anki y los TSV de respaldo de las tarjetas de estudio.

Salida, en esta misma carpeta:
  OncoGO-CRC_glosario_y_banco.apkg
  glosario.tsv
  banco_preguntas.tsv

El texto sale del banco de la defensa (secciones 11.3, 12 y 14.13 del plan
de estudio). Para regenerar: pip install genanki && python3 generar_anki.py
"""

from __future__ import annotations

import html
import sqlite3
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
APKG = ROOT / "OncoGO-CRC_glosario_y_banco.apkg"
GLOSARIO_TSV = ROOT / "glosario.tsv"
BANCO_TSV = ROOT / "banco_preguntas.tsv"

# Fijo a propósito: el mismo texto produce el mismo paquete.
TIMESTAMP = 1_762_387_200

MODEL_GLOSARIO_ID = 1850264101
MODEL_DEFENSA_ID = 1850264102
DECK_GLOSARIO_ID = 1850264201
DECK_BANCO_ID = 1850264202

DECK_GLOSARIO = "OncoGO-CRC::01 Glosario"
DECK_BANCO = "OncoGO-CRC::02 Banco de preguntas de la defensa"

BLOQUES = {
    "A": "Problema y diseño",
    "B": "Proteómica",
    "C": "Estadística",
    "D": "Biología y sangre",
    "E": "Bioinformática, ML y DL",
    "F": "Reproducibilidad y futuro",
}

CSS = """
.card {
  font-family: "Liberation Sans", "Segoe UI", sans-serif;
  font-size: 20px;
  line-height: 1.45;
  text-align: left;
  color: #1c1917;
  background-color: #fafaf9;
  padding: 6px 2px;
}
.card.nightMode {
  color: #f5f5f4;
  background-color: #1c1917;
}
.meta, .etiqueta {
  display: block;
  font-size: 12px;
  font-weight: 700;
  letter-spacing: 0.04em;
  text-transform: uppercase;
  color: #57534e;
}
.card.nightMode .meta, .card.nightMode .etiqueta { color: #a8a29e; }
.meta { margin-bottom: 0.7em; }
.etiqueta { margin-top: 0.95em; margin-bottom: 0.25em; }
.pregunta, .termino { font-size: 22px; font-weight: 600; }
.definicion-frente { font-size: 22px; font-weight: 600; }
.pista { margin-top: 1.3em; font-size: 15px; color: #78716c; }
.card.nightMode .pista { color: #a8a29e; }
.corta, .desarrollada, .apoyo {
  margin: 0.15em 0 0.4em 0;
  padding-left: 0.75em;
}
.corta { border-left: 4px solid #2f7d57; }
.desarrollada { border-left: 4px solid #3d5a80; }
.apoyo { border-left: 4px solid #a16207; }
"""

GLOSARIO_QFMT = """
<div class="meta">Glosario · {{Letra}}</div>
<div class="termino">{{Término}}</div>
"""

GLOSARIO_AFMT = """
{{FrontSide}}
<hr id="answer">
<div class="definicion">{{Definición}}</div>
"""

GLOSARIO_REV_QFMT = """
<div class="meta">Glosario · {{Letra}} · di el término</div>
<div class="definicion-frente">{{Definición}}</div>
<div class="pista">Di el término a partir de la definición.</div>
"""

GLOSARIO_REV_AFMT = """
{{FrontSide}}
<hr id="answer">
<div class="etiqueta">Término</div>
<div class="termino">{{Término}}</div>
"""

DEFENSA_QFMT = """
<div class="meta">{{Bloque}} · {{Numero}} de 55</div>
<div class="pregunta">{{Pregunta}}</div>
<div class="pista">En voz alta: frase directa, luego el número o la fuente, y al final el límite.</div>
"""

DEFENSA_AFMT = """
{{FrontSide}}
<hr id="answer">
<div class="etiqueta">Respuesta corta</div>
<div class="corta">{{Respuesta corta}}</div>
<div class="etiqueta">Respuesta desarrollada</div>
<div class="desarrollada">{{Respuesta desarrollada}}</div>
{{#Cierre}}
<div class="apoyo">{{Cierre}}</div>
{{/Cierre}}
"""

# (letra, término, definición) — 108, en el orden del glosario.
GLOSARIO: list[tuple[str, str, str]] = [
    ("A", "ABSOLUTE", "Método que estima la pureza tumoral desde el ADN (Carter 2012)"),
    ("A", "Accesión", "Identificador único: P02768 en UniProt, PXD055218 en ProteomeXchange"),
    ("A", "Adenoma", "Pólipo benigno que puede progresar a cáncer"),
    ("A", "AlphaDIA", "Motor abierto para analizar DIA con aprendizaje profundo (Wallmann 2026)"),
    ("A", "Aleatorización mendeliana (MR)", "Inferencia causal con variantes genéticas como instrumentos"),
    ("A", "AUC", "Área bajo la curva ROC: probabilidad de que un caso puntúe más que un control"),
    ("B", "Bayes empírico", 'Estimar parámetros "prestando" información de todas las proteínas (limma)'),
    ("B", "BH", "Benjamini-Hochberg: corrección que controla el FDR"),
    ("B", "Bootstrap", "Remuestrear con reemplazo para estimar la incertidumbre"),
    ("C", "camera", "Prueba de conjuntos competitiva que corrige la correlación (Wu y Smyth 2012)"),
    ("C", "CEA", "Antígeno carcinoembrionario (CEACAM5), marcador clínico de seguimiento"),
    ("C", "CI", "Integración continua: los tests corren solos en cada cambio"),
    ("C", "Clopper-Pearson", "Intervalo exacto para una proporción"),
    ("C", "CMS", "Subtipos moleculares de consenso del CCR (1 a 4)"),
    ("C", "Colocalización", "Probar si la misma variante causa dos señales (PP.H4)"),
    ("C", "Confusión", "Una tercera variable que explica la diferencia observada"),
    ("C", "CON__ / REV__", "Prefijos de MaxQuant para contaminantes y decoys"),
    ("C", "CPTAC", "Consorcio de proteogenómica del cáncer (NCI)"),
    ("C", "CV anidada", "Validación cruzada con un bucle interno para elegir hiperparámetros"),
    ("D", "Dalton (Da)", "Unidad de masa molecular; ~ la masa de un protón"),
    ("D", "DDA / DIA", "Adquisición dependiente o independiente de datos"),
    ("D", "Decoy", "Secuencia falsa (invertida) para estimar errores"),
    ("D", "Deming", "Regresión con error en ambas variables"),
    ("D", "DeLong", "Método para el IC de un AUC y para comparar dos AUC en las mismas personas"),
    ("D", "Dermcidina (DCD)", "Proteína del sudor; contaminante típico"),
    ("D", "DIA-NN", "Software de análisis DIA usado por los autores de PXD055218"),
    ("E", "Efecto de lote", "Diferencia técnica entre tandas de medición"),
    ("E", "EII", "Enfermedad inflamatoria intestinal (Crohn, colitis ulcerosa)"),
    ("E", "Embedding", "Vector que resume un objeto (una proteína)"),
    ("E", "Enmienda", "Cambio fechado y justificado a un preregistro"),
    ("E", "Entrapment", "Secuencias que no pueden estar en la muestra, para medir el FDR real"),
    ("E", "Equivalencia", "Demostrar que un efecto es menor que un margen (TOST)"),
    ("E", "ESM-2", "Modelo de lenguaje de proteínas de Meta (Lin 2023)"),
    ("E", "Estroma", "Tejido de soporte: fibroblastos, vasos, matriz"),
    ("F", "Fase aguda", "Respuesta inflamatoria sistémica (SAA1, CRP, LBP)"),
    ("F", "FDR", "Tasa de falsos descubrimientos"),
    ("F", "FIT", "Prueba inmunoquímica fecal de sangre oculta"),
    ("F", "Fisher (prueba exacta)", "Prueba de asociación en tablas 2 × 2 (presencia y ausencia)"),
    ("F", "Fuga de datos", "Información de la prueba que se cuela en el entrenamiento"),
    ("F", "FWER", "Probabilidad de al menos un falso positivo"),
    ("G", "GO", "Gene Ontology: proceso, función y componente"),
    ("G", "Gradient boosting", "Conjunto de árboles que corrigen errores sucesivos"),
    ("G", "Grupo proteico", "Conjunto de proteínas que los péptidos no permiten distinguir"),
    ("G", "GSEA", "Enriquecimiento sobre una lista ordenada completa"),
    ("G", "GWAS", "Estudio de asociación del genoma completo"),
    ("H", "Hemólisis", "Rotura de eritrocitos; libera hemoglobina y otras proteínas"),
    ("H", "HLA", "Moléculas que presentan péptidos al sistema inmune"),
    ("H", "Hodges-Lehmann", "Mediana de las diferencias entre grupos; tamaño del efecto de Mann-Whitney"),
    ("H", "HPA", "Human Protein Atlas"),
    ("I", "I²", "Fracción de la variación entre estudios debida a heterogeneidad real"),
    ("I", "Imputación", "Rellenar valores ausentes; la tesis no lo hace"),
    ("I", "Isoforma", "Variante de una proteína del mismo gen"),
    ("I", "IVW", "Media ponderada por la inversa de la varianza (MR)"),
    ("K", "Kruskal-Wallis", "Prueba no paramétrica para comparar tres o más grupos"),
    ("L", "LFQ", "Cuantificación sin marcaje"),
    ("L", "limma", "Paquete de modelos lineales con varianza moderada"),
    ("L", "Lin (CCC)", "Coeficiente de concordancia sobre y = x"),
    ("L", "log2FC", "Logaritmo en base 2 del cociente; 1 = doble, 0,585 = 1,5 veces"),
    ("M", "Mann-Whitney", "Prueba no paramétrica para dos grupos"),
    ("M", "Manifiesto", "Lista de ficheros con su huella SHA-256"),
    ("M", "MaxLFQ", "Algoritmo de LFQ basado en razones de péptidos (Cox 2014)"),
    ("M", "MAR / MCAR / MNAR", "Tipos de ausencia: al azar condicional, completamente al azar, no al azar"),
    ("M", "m/z", "Masa dividida por carga: lo que mide el espectrómetro"),
    ("M", "MDE", "Efecto mínimo detectable con la potencia fijada"),
    ("M", "MHCflurry", "Predictor abierto de presentación en HLA de clase I"),
    ("M", "MOFA+", "Factorización para integrar varias ómicas"),
    ("M", "MS1 / MS2", "Espectro de iones intactos / de fragmentos"),
    ("M", "MSI", "Inestabilidad de microsatélites; tumores con muchas mutaciones"),
    ("M", "Multiverso", "Repetir el análisis con todas las decisiones razonables"),
    ("N", "Neoantígeno", "Péptido mutado del tumor que el sistema inmune puede reconocer"),
    ("N", "Nulo del máximo", "Distribución del mejor resultado cuando no hay señal"),
    ("O", "Olink (PEA)", "Plataforma de afinidad por anticuerpos pareados, no MS"),
    ("P", "p-valor", "Probabilidad de un resultado tan extremo si H0 fuera cierta"),
    ("P", "Péptido proteotípico", "Péptido que pertenece a una sola proteína"),
    ("P", "Péptido señal", "Secuencia que dirige una proteína a secreción"),
    ("P", "Permutación", "Barajar etiquetas para construir un nulo"),
    ("P", "Pivotar", "Pasar una tabla de formato largo a ancho"),
    ("P", "Plasma / suero", "Sangre sin células con anticoagulante / tras coagular"),
    ("P", "Potencia", "Probabilidad de detectar un efecto real"),
    ("P", "pQTL", "Variante genética asociada al nivel de una proteína"),
    ("P", "Precursor", "Ion de un péptido con una carga dada"),
    ("P", "Preanalítica", "Todo lo que ocurre antes de medir la muestra"),
    ("P", "Preregistro", "Plan fijado antes de ver los datos"),
    ("P", "ProtT5", "Modelo de lenguaje de proteínas (Elnaggar 2022)"),
    ("P", "Pureza tumoral", "Fracción de células tumorales en la muestra"),
    ("Q", "Q de Cochran", "Suma ponderada de las distancias de cada estudio al efecto común; base de I²"),
    ("Q", "q-valor", "El menor FDR al que una identificación o prueba sería aceptada"),
    ("R", "R1 a R5", "Índices de la tesis: concordancia, detectabilidad, homogeneidad, composición y nulo"),
    ("R", "Regresión logística", "Modelo de probabilidad de clase"),
    ("R", "Regularización", "Penalización para evitar el sobreajuste (L1, L2, elastic-net)"),
    ("R", "Reitsma", "Meta-análisis bivariado de sensibilidad y especificidad"),
    ("R", "Rescoring", "Reordenar identificaciones con predicciones de espectro y tiempo"),
    ("S", "Sello de corrida", "JSON con commit, hash del preregistro y fecha"),
    ("S", "Semilla", "Número que fija el azar para reproducir resultados"),
    ("S", "SHA-256", "Huella digital de un fichero"),
    ("S", "SHAP", "Reparto de una predicción entre variables; explica el modelo, no la biología"),
    ("S", "Sobreajuste", "Memorizar el ruido del entrenamiento"),
    ("T", "Tamizaje", "Buscar enfermedad en personas sin síntomas"),
    ("T", "Target-decoy", "Estrategia de búsqueda con secuencias falsas"),
    ("T", "τ² (tau cuadrado)", "Varianza del efecto verdadero entre estudios en un meta-análisis"),
    ("T", "TMbed", "Predictor abierto de transmembrana y péptido señal"),
    ("T", "TMT", "Marcaje isobárico para multiplexar muestras"),
    ("T", "TOST", "Dos pruebas unilaterales para equivalencia"),
    ("T", "treat", "Prueba de cambio por encima de un umbral (McCarthy y Smyth 2009)"),
    ("T", "TRIPOD+AI", "Guía para reportar modelos de predicción"),
    ("U", "Universo", "Conjunto de fondo en un enriquecimiento"),
    ("V", "Validación externa", "Probar el modelo en datos de otro recurso"),
    ("V", "VPP", "Valor predictivo positivo: probabilidad de enfermedad si la prueba sale positiva; depende de la prevalencia"),
]

# (número, bloque, pregunta, respuesta corta, respuesta desarrollada, cierre)
PREGUNTAS: list[tuple[int, str, str, str, str, str]] = [
    (
        1,
        "A",
        "¿Por qué no propones un biomarcador?",
        "Porque con 20 sueros por grupo y sin réplica externa no se puede. Lo que aporto es medir qué candidatas tienen opciones. (Parte 1.3.)",
        "Porque los datos no lo permiten. Con 20 sueros por grupo, el efecto mínimo detectable es de unas 4 veces, solo hay un estudio con un grupo de enfermedad intestinal, y no conozco la preanalítica. Lo que aporto es un método y un estándar para decidir qué candidatas merecen una validación cara.",
        "Apoyo: K4 (0 diferentes frente a colonoscopia negativa), Pepe 2001 fase 1.",
    ),
    (
        2,
        "A",
        "¿Por qué no fusionaste cohortes?",
        "Distintas plataformas y personas; los efectos de lote no se separan de la biología (Leek 2010).",
        "TMT, LFQ y DIA miden en escalas distintas y en personas distintas. Si junto las intensidades, la plataforma explica más varianza que la enfermedad, y no puedo separar una de otra. Por eso comparo efectos calculados dentro de cada cohorte: signo, log2FC, AUC.",
        "Apoyo: Leek 2010.",
    ),
    (
        3,
        "A",
        "¿Por qué es tan importante el grupo CN?",
        "Es el único control de otra enfermedad intestinal; sin él, la fase aguda parece señal de cáncer.",
        "Porque es lo más parecido al uso real: personas que llegan a colonoscopia y no tienen cáncer. Frente a sanos, cualquier marcador de inflamación parece de cáncer: SAA1 y LBP salen frente a sanos y desaparecen frente a CN.",
        "Apoyo: K4.",
    ),
    (
        4,
        "A",
        "¿Qué es exactamente CN?",
        "Responde con lo que verifiques en Rai 2025 (Parte 10.4, punto 3).",
        'Responde con lo que diga Rai 2025 (métodos y tabla de participantes). Si es colonoscopia negativa, di "personas con indicación de colonoscopia y resultado negativo". Si es EII, ajusta todas las frases del plan y explica que entonces el contraste es "cáncer frente a enfermedad inflamatoria intestinal". Hazlo antes de la defensa.',
        "",
    ),
    (
        5,
        "A",
        "¿Por qué no imputaste?",
        "En proteómica las ausencias son sobre todo MNAR; imputar crea diferencias falsas. La presencia se prueba aparte (K6).",
        "Porque en proteómica falta sobre todo lo que está por debajo del límite de detección (MNAR). Imputar crea o borra diferencias que nadie midió. En vez de eso, probé la presencia por separado con Fisher: 358 proteínas por debajo del 70 % de observaciones; 11 significativas en presencia.",
        "Apoyo: K6, Lazar 2016.",
    ),
    (
        6,
        "A",
        "¿Qué es un preregistro y por qué te sirve?",
        "Fijar hipótesis, umbrales y semillas antes de ver; evita el p-hacking (Nosek 2018).",
        "Un documento con hipótesis, umbrales, semillas y reglas de decisión, fechado en Git antes de calcular. Impide ajustar el análisis al resultado.",
        "Apoyo: config/preregistro_v3.yaml, Nosek 2018.",
    ),
    (
        7,
        "A",
        "¿Viste los datos antes de la enmienda?",
        "Sí, la corrida v3 estaba en disco. Lo declaro; v3.2 añade sello de corrida y declaración de transparencia.",
        "Sí, y lo declaro. El preregistro v3 se commiteó 14 minutos antes que las tablas, con la corrida en disco. Por eso la v3.2 añade un sello en cada resultado y una declaración de transparencia, y mantiene el análisis original como referencia. No lo minimices.",
        "",
    ),
    (
        8,
        "A",
        "¿Qué es H3?",
        "No existe en ningún preregistro. Lo digo así.",
        "No existe. Ni el preregistro v2 ni el v3 tienen una H3. Lo dejo escrito así para que nadie piense que se escondió.",
        "Apoyo: hipotesis_v3.tsv.",
    ),
    (
        9,
        "A",
        "¿Por qué H7 es la confirmatoria?",
        "Es la pregunta propia de la tesis y la única con validación fuera de muestra en varios recursos.",
        "Porque es la pregunta propia de la tesis: si se puede predecir de antemano qué cambios del tumor se ven en sangre. Es la única con validación fuera de muestra en varios recursos independientes, y con una regla de decisión numérica.",
        "",
    ),
    (
        10,
        "A",
        "¿Por qué Holm entre H5, H8 y H10?",
        "Son tres pruebas secundarias; Holm controla el error de familia y es más potente que Bonferroni.",
        "Son tres pruebas secundarias; si cada una usa 0,05, la probabilidad de al menos un falso positivo sube a casi 0,14. Holm la mantiene en 0,05 y es más potente que Bonferroni.",
        "Ejemplo: con p = 0,01, 0,03 y 0,04, se comparan con 0,0167, 0,025 y 0,05; pasa la primera, la segunda no (0,03 > 0,025) y ahí se para.",
    ),
    (
        11,
        "B",
        "¿Qué es el FDR target-decoy?",
        "Se buscan también secuencias invertidas; los aciertos falsos con decoy estiman los falsos con target (Elias y Gygi 2007).",
        "Se buscan los espectros contra proteínas reales y contra secuencias invertidas que no pueden existir. Cuántas invertidas pasan un umbral estima cuántas reales falsas lo pasan. Se acepta hasta que la fracción estimada de falsos es 1 %.",
        "Apoyo: Elias y Gygi 2007.",
    ),
    (
        12,
        "B",
        "¿Por qué el FDR de precursor no basta?",
        "Los errores no se heredan; el FDR de proteína puede ser mucho mayor (Savitski 2015).",
        "Porque las proteínas falsas suelen tener un solo precursor y se acumulan al agregar. Con 1 % de precursores falsos, la fracción de proteínas falsas puede ser varias veces mayor. Hay que filtrar por q de proteína.",
        "Apoyo: Savitski 2015; en el suero, 432 grupos nunca pasan PG.Q.Value ≤ 0,01.",
    ),
    (
        13,
        "B",
        "¿Qué diferencia hay entre PG.Q.Value y Global.PG.Q.Value?",
        "El primero es por corrida; el segundo, sobre todo el experimento.",
        "El primero se calcula en cada corrida; el segundo con todas las corridas juntas. Uso los dos: el global asegura que el grupo existe en el experimento; el de corrida, que el valor de esa muestra es fiable.",
        "",
    ),
    (
        14,
        "B",
        "¿Qué es MaxLFQ?",
        "Cuantificación sin marcaje que usa razones de péptidos compartidos entre muestras (Cox 2014).",
        "Un algoritmo que estima la abundancia relativa de una proteína en cada muestra usando los cocientes entre muestras de los péptidos que comparten. Evita que un péptido muy intenso o que falta en unas muestras distorsione el valor.",
        "Apoyo: Cox 2014.",
    ),
    (
        15,
        "B",
        "¿Qué es la compresión de razones del TMT?",
        "Los iones que se cofragmentan aplastan las razones hacia 1; por eso Lin se ve afectado y T4 usa Deming.",
        "En TMT, otros péptidos que se aíslan a la vez mezclan sus reporteros y aplastan los cocientes hacia 1. Por eso los log2FC de CPTAC son más pequeños que los de Tanaka y uso Deming y Spearman además de Lin.",
        "Apoyo: Savitski 2013.",
    ),
    (
        16,
        "B",
        "¿Qué diferencia hay entre DDA y DIA?",
        "DDA elige los iones más intensos; DIA fragmenta ventanas completas y es más reproducible.",
        "DDA fragmenta solo los iones más intensos de cada momento: muestrea al azar lo poco abundante y deja ausencias. DIA fragmenta todo por ventanas: más completo y reproducible, pero necesita software que separe espectros mezclados.",
        "",
    ),
    (
        17,
        "B",
        "¿Cómo sabes que una proteína es humana y no un artefacto de la base?",
        "S2 (especie y ortólogo compartido) y S3 (entrapment).",
        "Miro todos los miembros del grupo: si ninguno es _HUMAN, el grupo sale. Además compruebo si sus péptidos podrían ser también humanos (ortólogo compartido) y uso los no humanos para estimar el FDR real.",
        "Apoyo: S2 y S3.",
    ),
    (
        18,
        "B",
        "¿Las queratinas no son biología?",
        "Son el contaminante más común (Hodge 2013; Frankenfield 2022): piel y polvo, calculadas por muestra.",
        "Las del epitelio sí podrían serlo, pero KRT1, KRT9 y KRT10 son de piel, KRT9 de palmas y plantas, y van con la dermcidina del sudor. Son el contaminante número uno en proteómica. Las calculo por muestra como puntuación de contaminación.",
        "Apoyo: Hodge 2013; Frankenfield 2022.",
    ),
    (
        19,
        "B",
        "¿Qué es el entrapment y qué FDR real estimas?",
        "Proteínas que no pueden estar en la muestra; su tasa estima el FDR real (Wen 2025). El número sale de S3.",
        "Añadir a la búsqueda secuencias de especies que no pueden estar en la muestra; cuántas se identifican estima el FDR real. Aquí los 215 grupos no humanos funcionan así sin haberlo planeado. El número exacto lo da S3 como cota, porque no conozco la base que usaron los autores.",
        "Fórmula de Wen 2025: FDP ≈ N_e·(1 + 1/r)/(N_t + N_e), donde r es el tamaño de la base de entrapment dividido por el de la humana.",
    ),
    (
        20,
        "B",
        "¿Por qué no corriste DIA-NN?",
        "Su licencia no es libre (regla 15). Se cita porque los autores lo usaron. AlphaDIA es el trabajo futuro abierto.",
        "La regla del proyecto es software libre que corra en local en macOS y Linux. DIA-NN no cumple. Lo cito porque los autores lo usaron. Para repetir la búsqueda propongo AlphaDIA, que es abierto.",
        "",
    ),
    (
        21,
        "B",
        "¿Repetirías la búsqueda?",
        "Sí, con AlphaDIA contra FASTA humano más contaminantes (10.4 de v3.2).",
        "Sí, como trabajo futuro o extensión: AlphaDIA con FASTA humano más la lista de contaminantes de Frankenfield, FDR de proteína al 1 %, y comparar con el report de los autores.",
        "",
    ),
    (
        22,
        "C",
        "¿Qué es un p-valor?",
        "La probabilidad de ver algo tan extremo si H0 fuera cierta. No es la probabilidad de H0.",
        "La probabilidad de obtener un resultado al menos tan extremo como el observado si la hipótesis nula fuera cierta. No es la probabilidad de que la nula sea cierta, ni la probabilidad de que el resultado sea un error, ni mide el tamaño del efecto.",
        'Añade: "Por eso no digo \'hay un 2,45 % de probabilidad de que sea azar\' para K5, sino \'si no hubiera relación, un máximo así saldría en el 2,45 % de las permutaciones\'."',
    ),
    (
        23,
        "C",
        "¿Qué hace BH?",
        "Controla la fracción esperada de falsos entre los rechazos.",
        "Ordena los p y acepta los más pequeños hasta que la fracción esperada de falsos entre los aceptados es como mucho 5 %. No promete cero falsos; promete que en promedio no pasan del 5 %.",
        "",
    ),
    (
        24,
        "C",
        "¿Por qué el alfa efectivo era 7,14 × 10⁻⁵?",
        "876 pruebas, casi ninguna señal; BH se vuelve casi Bonferroni.",
        "Si casi ninguna proteína tiene señal, BH se comporta casi como Bonferroni: 0,05/876 = 5,7 × 10⁻⁵. La simulación con las varianzas observadas dio 7,14 × 10⁻⁵. Por eso es tan importante no inflar el denominador con proteínas que no pasan el FDR.",
        "",
    ),
    (
        25,
        "C",
        "¿Qué diferencia hay entre treat y una prueba t?",
        "treat prueba si el cambio supera un umbral (0,585 log2, 1,5 veces); la t prueba si es distinto de 0.",
        "La t pregunta si el cambio es distinto de 0; con suficiente gente, cualquier cambio diminuto lo es. treat pregunta si supera 1,5 veces. Un cambio de 0,9 log2 con error 0,2 es muy significativo con la t (p ≈ 6 × 10⁻⁵) pero no con treat (p ≈ 0,06).",
        "",
    ),
    (
        26,
        "C",
        "¿Qué aporta TOST?",
        'Permite afirmar equivalencia: que el cambio es menor que ±0,585. Sin TOST, "no significativo" no es "igual".',
        'Permite decir "no cambia más de 1,5 veces" con evidencia. Sin TOST, "no significativo" solo significa "no lo sé".',
        "",
    ),
    (
        27,
        "C",
        "¿Cero diferentes es ausencia de efecto?",
        "No. Solo 145 son equivalentes; 731 son indeterminadas.",
        "No. Solo 145 proteínas son equivalentes: en ellas sí puedo decir que el cambio es pequeño. En 731 no tengo información suficiente.",
        "",
    ),
    (
        28,
        "C",
        "¿Por qué el máximo AUC es alto sin señal?",
        "Elegir el mejor de cientos infla el máximo; el nulo permutado del máximo lo corrige (media 0,843).",
        'Cada proteína tiene su propio AUC por azar alrededor de 0,5; al "plegarlo" (máximo entre AUC y 1 − AUC) siempre queda ≥ 0,5. Si eliges el mayor de cientos, eliges el extremo de la distribución del azar. En los datos reales el nulo permutado dio media 0,843; el observado (0,9375) apenas lo supera, y además era una queratina. Simulación sin ninguna señal, 20 contra 20: una sola proteína da un AUC plegado medio de 0,57; el máximo de 600 proteínas independientes, 0,79 de media.',
        "",
    ),
    (
        29,
        "C",
        "¿Qué mide Lin que no mide Spearman?",
        "Si los valores caen sobre y = x: acuerdo en magnitud, no solo en orden.",
        "El acuerdo sobre la diagonal: que los cambios tengan el mismo signo y el mismo tamaño. Lin 0,734 con el primario y 0,440 con la metástasis: el signo se conserva, pero el hígado distorsiona el tamaño.",
        "",
    ),
    (
        30,
        "C",
        "¿Por qué Clopper-Pearson?",
        "Es exacto para proporciones binomiales, también cerca de 1.",
        "Es el intervalo exacto de una proporción binomial. Cerca de 1 (116/118) los intervalos aproximados se salen del rango o son demasiado estrechos. Clopper-Pearson es conservador y nunca pasa de 1. 116/118 da 0,940 a 0,998.",
        "",
    ),
    (
        31,
        "C",
        "¿Por qué R3 necesita ≥ 3 cohortes independientes?",
        "Con 2, I² es muy inestable; con cohortes solapadas, la heterogeneidad es falsa.",
        "I² se calcula con Q de Cochran y sus grados de libertad (cohortes − 1). Con 2 cohortes hay 1 grado de libertad y la estimación de la heterogeneidad salta entre 0 y valores altos con cualquier cambio pequeño. Además, las cohortes tienen que ser independientes (K2).",
        "",
    ),
    (
        32,
        "C",
        "¿Qué son los errores tipo S y M?",
        "Con poca potencia, lo significativo puede tener el signo equivocado (S) o un tamaño exagerado (M) (Gelman y Carlin 2014).",
        "Con poca potencia, lo que sale significativo tiene que tener un efecto observado grande, aunque el real sea pequeño: se exagera (M) o incluso se invierte el signo (S). Por eso desconfío de efectos grandes en n pequeños, como CHMP4A con log2FC 5,42 sobre pocos observados.",
        "Apoyo: Gelman y Carlin 2014.",
    ),
    (
        33,
        "C",
        "¿Cuántas versiones del análisis probaste?",
        "Todas están en el multiverso, con la fracción que mantiene cada conclusión (Steegen 2016).",
        "Todas las razonables están en el multiverso: umbral treat, FDR, presencia mínima, filtros del suero, agregación y contraste del tejido. Informo qué fracción de combinaciones mantiene cada conclusión. El análisis principal es el preregistrado.",
        "Apoyo: Steegen 2016.",
    ),
    (
        34,
        "C",
        "¿Por qué cambió el resultado del suero entre v3 y v3.2?",
        "El filtro de FDR de proteína y la retirada de contaminantes cambian los denominadores (suero_resumen.json).",
        "Porque la matriz v3 incluía 876 grupos sin filtro de FDR de proteína, 215 no humanos y contaminantes de piel. Al filtrar quedan 307 grupos humanos fiables. Cambian el denominador de BH, la potencia y el nulo del máximo. La conclusión cualitativa se comprueba en suero_resumen.json.",
        "",
    ),
    (
        35,
        "D",
        "¿Por qué la fase aguda confunde?",
        "SAA1, CRP, LBP suben con cualquier inflamación; no distinguen cáncer de otra enfermedad.",
        "El hígado sube SAA1, CRP, LBP o SERPINA3 ante cualquier inflamación: infección, EII, cirugía o un tumor. Son reales, pero no dicen qué enfermedad hay. Por eso H10 las contrasta con 59 enfermedades.",
        "",
    ),
    (
        36,
        "D",
        "¿Por qué FGB alto en suero es sospechoso?",
        "El fibrinógeno se consume al coagular; si queda mucho, la coagulación fue distinta.",
        "El suero se obtiene dejando coagular; el fibrinógeno se gasta en el coágulo. Si queda mucho, el coágulo no se formó bien o el tiempo y la temperatura fueron distintos. Es una señal del proceso, no de la enfermedad (FGB +1,40 en la revisión; coagulación p = 0,0063 en K7).",
        "",
    ),
    (
        37,
        "D",
        "¿Y si cáncer y controles se procesaron distinto?",
        "S4 y S8: paneles de Geyer 2019, ajuste por coagulación, profundidad por corrida (igual entre grupos, revisión).",
        "Es la explicación más probable de parte de la señal. Lo mido: profundidad igual (precursores, p = 0,86), pero fracción de alta abundancia (p = 1,23 × 10⁻⁵), coagulación (p = 0,0063) y piel distintas. Por eso declaro confusión preanalítica y ajusto por coagulación como sensibilidad.",
        "Apoyo: K7, Geyer 2019.",
    ),
    (
        38,
        "D",
        "¿Qué es el CEA y por qué no basta?",
        "Antígeno carcinoembrionario (CEACAM5); poco sensible en estadio temprano. Sirve de seguimiento.",
        "El antígeno carcinoembrionario, producto de CEACAM5. Se usa para seguir pacientes ya tratados. No sirve para tamizaje porque en estadios tempranos sale normal con frecuencia, y también sube en fumadores e inflamación.",
        "",
    ),
    (
        39,
        "D",
        "¿Por qué la metástasis hepática no sustituye al primario?",
        "El hígado aporta su proteoma; Lin 0,440 frente a 0,734.",
        'La muestra M contiene tumor y tejido de hígado. Las proteínas de hepatocito (albúmina, CPS1, ARG1, CYP2E1...) aparecerían "subidas en el tumor" cuando en realidad vienen del hígado. Por eso Lin bajó a 0,440 y el contraste correcto es P frente a N.',
        "",
    ),
    (
        40,
        "D",
        "¿Qué es CMS4 y por qué importa?",
        "Subtipo mesenquimal, rico en estroma; cambia la composición del tejido (Guinney 2015).",
        "El subtipo mesenquimal: mucho estroma (fibroblastos), peor pronóstico. Importa porque un tumor con mucho estroma cambia el proteoma del tejido aunque el epitelio no cambie, y eso afecta a R4 y H6.",
        "Apoyo: Guinney 2015.",
    ),
    (
        41,
        "E",
        "¿Por qué el universo de GO debe ser lo detectado?",
        "Si usas el genoma, lo abundante parece enriquecido.",
        'El espectrómetro solo ve lo abundante (proteínas extracelulares, metabólicas, estructurales). Si comparas con el genoma completo, esos términos salen "enriquecidos" solo porque son los que se pueden medir. El fondo correcto es lo que la cohorte detectó, y mejor aún, emparejado por abundancia.',
        "",
    ),
    (
        42,
        "E",
        "¿Por qué no interpretas redes (H2)?",
        "Falta un nulo que conserve el grado; los nodos muy conectados siempre salen.",
        "Las proteínas muy conectadas en STRING salen en cualquier lista. Sin comparar con al menos 1000 redes aleatorias que conserven el número de conexiones de cada proteína, no puedo distinguir una red real de un artefacto de conectividad. Por eso H2 queda no evaluada.",
        "",
    ),
    (
        43,
        "E",
        "¿Por qué no entrenaste una red neuronal con los sueros?",
        "60 muestras y cientos de variables: sobreajuste seguro (Vabalas 2019).",
        "60 muestras con cientos de proteínas. Una red tiene miles de parámetros: memorizaría las muestras. Además, las diferencias que hay son de contaminación y proceso; una red las aprendería encantada.",
        "",
    ),
    (
        44,
        "E",
        "¿Cómo evitaste fugas?",
        "Todo (normalizar, seleccionar) dentro de cada pliegue; CV anidada (Kapoor y Narayanan 2023).",
        "Todo lo que aprende de los datos (normalizar, seleccionar, ajustar) ocurre dentro de cada pliegue de entrenamiento; hiperparámetros en un bucle interno; y en H7 dejo fuera un recurso de sangre completo.",
        "Apoyo: Kapoor y Narayanan 2023.",
    ),
    (
        45,
        "E",
        "¿Qué es la prueba de permutación del clasificador?",
        "Repetir toda la validación con etiquetas barajadas; el AUC real debe superar al barajado (Ojala y Garriga 2010).",
        'Responde "¿el rendimiento es mejor que el que obtendría el mismo procedimiento con etiquetas sin sentido?". Con 20 contra 20, un AUC de 0,70 puede salir por azar; la permutación da la distribución de referencia.',
        "",
    ),
    (
        46,
        "E",
        "¿Qué aporta ESM-2 que no aporte la abundancia?",
        "Es justo lo que mide H7: el aporte incremental sobre el modelo base.",
        "Esa es la pregunta de H7: comparo un modelo base con solo abundancia contra uno completo que añade localización, secreción, tipo celular y variables de secuencia. Si la diferencia de AUC tiene IC por encima de 0, la secuencia aporta; si no, la conclusión es que manda la abundancia.",
        "",
    ),
    (
        47,
        "E",
        "¿Por qué no SignalP, DeepLoc o NetMHCpan?",
        "Exigen licencia académica; TMbed, Light Attention y MHCflurry son abiertos y se compararon con UniProt.",
        "La regla 15 exige software libre que corra en local en macOS y Linux. SignalP y DeepLoc requieren licencia académica. TMbed (Apache-2.0) y Light Attention (MIT) son abiertos y se comparan con UniProt para mostrar que dan resultados razonables.",
        'Añade: "Comparé sus sustitutos con la anotación de UniProt (h7_concordancia_predictores.tsv)."',
    ),
    (
        48,
        "E",
        "¿Los neoantígenos de Rai se presentan en HLA?",
        "Exploratorio con MHCflurry; una matriz de grupos proteicos no puede confirmarlos.",
        "Lo exploro con MHCflurry en alelos frecuentes, pero predecir unión no es demostrar presentación ni detección. Y con grupos proteicos no puedo confirmar que el péptido mutado esté en el suero.",
        "",
    ),
    (
        49,
        "E",
        "¿Qué es la aleatorización mendeliana?",
        "Un ensayo natural con alelos repartidos al azar; exploratoria (H9).",
        'Al nacer, cada persona recibe al azar una de las dos copias de cada gen de sus padres. Si una variante sube el nivel de la proteína X, las personas con esa variante son como un grupo "tratado" con más X desde el nacimiento, sin que influyan la dieta o el tabaco. Si ese grupo tiene más cáncer colorrectal, y la variante solo actúa a través de X, X causa riesgo.',
        'Añade: "Es exploratoria: no decide la tesis."',
    ),
    (
        50,
        "E",
        "¿Por qué no mezclar el GWAS multiancestral con pQTL europeos?",
        "Las frecuencias alélicas y el desequilibrio de ligamiento difieren; se sesga el efecto.",
        "Las frecuencias de las variantes y su correlación con las vecinas cambian entre poblaciones. Si el efecto sobre la proteína se estimó en europeos y el efecto sobre el cáncer en una mezcla, el cociente queda sesgado. Uso el GWAS sin UK Biobank (GCST90255676) con los pQTL de UKB-PPP para evitar además el solapamiento de personas.",
        "",
    ),
    (
        51,
        "F",
        "¿Cualquiera puede repetir tu análisis?",
        "Sí: software libre, macOS y Linux, CI en ambos, un comando, sello y manifiestos SHA-256.",
        "Sí. Código abierto en GitHub, solo software libre, un comando por etapa, tests en macOS y Linux, semillas fijas, manifiestos SHA-256 de entradas y salidas, y un sello en cada resultado.",
        "",
    ),
    (
        52,
        "F",
        "¿Cómo pruebas que no tocaste los datos de entrada?",
        "Los manifiestos SHA-256 cambian si cambia un bit.",
        "Demuestra que los ficheros son los mismos que se registraron (no se tocaron los datos ni los resultados congelados). No demuestra que el análisis sea correcto ni que se hiciera antes de ver los datos. Añade la huella del report: 213 313 241 bytes, SHA-256 A8A5AADD...59E8.",
        "",
    ),
    (
        53,
        "F",
        "¿Qué recomendarías a un laboratorio en Colombia?",
        "Protocolo preanalítico único, grupo de enfermedad intestinal, tamaño muestral calculado con el MDE de esta tesis y CEA medido.",
        "Cuatro grupos (cáncer, colonoscopia negativa, EII y sanos), mismo protocolo de toma para todos (tubo, tiempo a centrifugar, temperatura, orden de corrida al azar), paneles de calidad por muestra, CEA medido, y tamaño calculado con el MDE de esta tesis: con desviación 1,5 log2, para detectar 1,5 veces hacen falta muchas más de 20 personas por grupo.",
        "",
    ),
    (
        54,
        "F",
        "¿Qué es lo más débil de tu tesis?",
        "El suero: 20 por grupo, un único estudio, preanalítica desconocida. Por eso no propongo biomarcador.",
        "El suero: un solo estudio con grupo de enfermedad, 20 por grupo y preanalítica desconocida. Lo digo yo primero, y por eso la conclusión del suero es sobre lo que no se puede afirmar.",
        "",
    ),
    (
        55,
        "F",
        "¿Qué harías con más presupuesto?",
        "Tejido y sangre de la misma persona, y proteómica dirigida de las candidatas que resistan H10.",
        "Tejido y sangre de la misma persona el mismo día, para seguir una proteína del tumor a la sangre; y proteómica dirigida de las pocas candidatas que sobrevivan H10.",
        "",
    ),
]


def comprobar_fuente() -> None:
    if len(GLOSARIO) != 108:
        raise SystemExit(f"El glosario tiene {len(GLOSARIO)} términos; esperaba 108.")
    if len(PREGUNTAS) != 55:
        raise SystemExit(f"El banco tiene {len(PREGUNTAS)} preguntas; esperaba 55.")
    numeros = [n for n, *_ in PREGUNTAS]
    if numeros != list(range(1, 56)):
        raise SystemExit(f"La numeración del banco no es 1..55: {numeros}")
    terminos = [t for _, t, _ in GLOSARIO]
    if len(terminos) != len(set(terminos)):
        raise SystemExit("Hay términos repetidos en el glosario.")
    definiciones = [d for *_, d in GLOSARIO]
    if len(definiciones) != len(set(definiciones)):
        raise SystemExit("Hay definiciones repetidas; la tarjeta inversa sería ambigua.")
    preguntas = [p for _, _, p, *_ in PREGUNTAS]
    if len(preguntas) != len(set(preguntas)):
        raise SystemExit("Hay preguntas repetidas.")
    for n, bloque, pregunta, corta, desarrollada, _cierre in PREGUNTAS:
        if bloque not in BLOQUES:
            raise SystemExit(f"Bloque desconocido en la pregunta {n}: {bloque}")
        if not pregunta.strip() or not corta.strip() or not desarrollada.strip():
            raise SystemExit(f"La pregunta {n} tiene un campo vacío.")
    agujas = {
        1: "20 sueros",
        4: "Rai 2025",
        8: "No existe",
        10: "0,025",
        19: "Wen 2025",
        24: "7,14 × 10⁻⁵",
        27: "731",
        28: "0,843",
        34: "307 grupos",
        37: "1,23 × 10⁻⁵",
        47: "h7_concordancia_predictores.tsv",
        52: "A8A5AADD",
        55: "H10",
    }
    por_numero = {n: " ".join(resto) for n, *resto in PREGUNTAS}
    for n, aguja in agujas.items():
        if aguja not in por_numero[n]:
            raise SystemExit(f"La pregunta {n} no contiene {aguja!r}.")
    if not any(t == "VPP" for _, t, _ in GLOSARIO):
        raise SystemExit("Falta VPP en el glosario.")
    if not any(t == "τ² (tau cuadrado)" for _, t, _ in GLOSARIO):
        raise SystemExit("Falta τ² en el glosario.")


def texto_html(texto: str) -> str:
    return html.escape(texto, quote=False)


def escribir_tsv(path: Path, directivas: list[str], filas: list[tuple[str, ...]]) -> None:
    lineas = list(directivas)
    for fila in filas:
        for campo in fila:
            if "\t" in campo or "\n" in campo or "\r" in campo:
                raise SystemExit(f"Campo con tabulador o salto de línea en {path.name}: {campo[:80]!r}")
        lineas.append("\t".join(fila))
    path.write_text("\n".join(lineas) + "\n", encoding="utf-8")


def escribir_glosario_tsv() -> None:
    filas = []
    for letra, termino, definicion in GLOSARIO:
        filas.append((termino, definicion, f"letra::{letra} glosario"))
    escribir_tsv(
        GLOSARIO_TSV,
        [
            "#separator:tab",
            "#html:false",
            f"#deck:{DECK_GLOSARIO}",
            "#tags column:3",
        ],
        filas,
    )


def escribir_banco_tsv() -> None:
    filas = []
    for n, bloque, pregunta, corta, desarrollada, cierre in PREGUNTAS:
        frente = f"Bloque {bloque} · {BLOQUES[bloque]} · {n}. {pregunta}"
        partes = [
            "<b>Respuesta corta</b><br>",
            texto_html(corta),
            "<br><br><b>Respuesta desarrollada</b><br>",
            texto_html(desarrollada),
        ]
        if cierre:
            partes.extend(["<br><br>", texto_html(cierre)])
        reverso = "".join(partes)
        filas.append((frente, reverso, f"bloque::{bloque} defensa"))
    escribir_tsv(
        BANCO_TSV,
        [
            "#separator:tab",
            "#html:true",
            f"#deck:{DECK_BANCO}",
            "#tags column:3",
        ],
        filas,
    )


def escribir_apkg() -> None:
    try:
        import genanki
    except ImportError:
        raise SystemExit("Falta genanki. Instálalo con: pip install genanki") from None

    modelo_glosario = genanki.Model(
        MODEL_GLOSARIO_ID,
        "OncoGO Glosario",
        fields=[{"name": "Término"}, {"name": "Definición"}, {"name": "Letra"}],
        templates=[
            {"name": "Término → definición", "qfmt": GLOSARIO_QFMT, "afmt": GLOSARIO_AFMT},
            {"name": "Definición → término", "qfmt": GLOSARIO_REV_QFMT, "afmt": GLOSARIO_REV_AFMT},
        ],
        css=CSS,
    )
    modelo_defensa = genanki.Model(
        MODEL_DEFENSA_ID,
        "OncoGO Defensa",
        fields=[
            {"name": "Pregunta"},
            {"name": "Bloque"},
            {"name": "Numero"},
            {"name": "Orden"},
            {"name": "Respuesta corta"},
            {"name": "Respuesta desarrollada"},
            {"name": "Cierre"},
        ],
        templates=[{"name": "Pregunta", "qfmt": DEFENSA_QFMT, "afmt": DEFENSA_AFMT}],
        css=CSS,
        sort_field_index=3,
    )

    mazo_glosario = genanki.Deck(
        DECK_GLOSARIO_ID,
        DECK_GLOSARIO,
        description=(
            "108 términos. Cada nota tiene dos tarjetas: término → definición, "
            "y definición → término. Con «enterrar relacionadas» activado, la inversa "
            "espera al día siguiente. Para una sola dirección, en el buscador: "
            'card:"Término → definición" o card:"Definición → término".'
        ),
    )
    mazo_banco = genanki.Deck(
        DECK_BANCO_ID,
        DECK_BANCO,
        description=(
            "55 preguntas de la defensa, en seis bloques. Di la frase directa, "
            "luego el número o la fuente, y al final el límite. Después elige "
            "Otra vez, Difícil, Bien o Fácil. El primer día sube el límite de "
            "tarjetas nuevas de este mazo para poder ver las 55."
        ),
    )

    for i, (letra, termino, definicion) in enumerate(GLOSARIO, start=1):
        mazo_glosario.add_note(
            genanki.Note(
                model=modelo_glosario,
                fields=[texto_html(termino), texto_html(definicion), texto_html(letra)],
                tags=[f"letra::{letra}", "glosario"],
                guid=genanki.guid_for("oncogo-glosario", termino),
                due=i,
            )
        )

    for n, bloque, pregunta, corta, desarrollada, cierre in PREGUNTAS:
        mazo_banco.add_note(
            genanki.Note(
                model=modelo_defensa,
                fields=[
                    texto_html(pregunta),
                    texto_html(f"Bloque {bloque} · {BLOQUES[bloque]}"),
                    str(n),
                    f"{n:02d}",
                    texto_html(corta),
                    texto_html(desarrollada),
                    texto_html(cierre),
                ],
                tags=[f"bloque::{bloque}", "defensa"],
                guid=genanki.guid_for("oncogo-defensa", f"{n:02d}"),
                due=n,
            )
        )

    paquete = genanki.Package([mazo_glosario, mazo_banco])
    paquete.write_to_file(APKG, timestamp=TIMESTAMP)


def verificar_apkg() -> None:
    with zipfile.ZipFile(APKG) as paquete:
        nombres = set(paquete.namelist())
        if "collection.anki2" not in nombres or "media" not in nombres:
            raise SystemExit(f"El paquete no tiene la estructura de Anki: {sorted(nombres)}")
        db_bytes = paquete.read("collection.anki2")
    with tempfile.NamedTemporaryFile(suffix=".anki2") as tmp:
        tmp.write(db_bytes)
        tmp.flush()
        con = sqlite3.connect(tmp.name)
        notas = con.execute("select count(*) from notes").fetchone()[0]
        tarjetas = con.execute("select count(*) from cards").fetchone()[0]
        mazos = con.execute("select decks from col").fetchone()[0]
        con.close()
    if notas != 163 or tarjetas != 271:
        raise SystemExit(f"Conteo inesperado: {notas} notas y {tarjetas} tarjetas.")
    if DECK_GLOSARIO not in mazos or DECK_BANCO not in mazos:
        raise SystemExit("Falta uno de los dos mazos dentro del paquete.")


def main() -> None:
    comprobar_fuente()
    escribir_glosario_tsv()
    escribir_banco_tsv()
    escribir_apkg()
    verificar_apkg()
    print(f"Glosario: {len(GLOSARIO)} términos → {GLOSARIO_TSV.name}")
    print(f"Banco: {len(PREGUNTAS)} preguntas → {BANCO_TSV.name}")
    print(f"Paquete: 163 notas, 271 tarjetas → {APKG.name}")


if __name__ == "__main__":
    main()
