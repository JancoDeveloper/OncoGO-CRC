"""Interpretación biológica automática y glosario para los Markdown del sitio MkDocs.

Las funciones leen los TSV de GO/GSEA/diferencial y producen secciones Markdown que **traducen
las observaciones cuantitativas** a contexto biológico de CRC, citando referencias clave. Es
intencionalmente conservadora: el lenguaje es de hipótesis, no de afirmación clínica.

Estrategia para los leading edges GSEA:

1. Normalizar nombres de genes (uppercase).
2. Detectar familias funcionales recurrentes vía prefijos / patrones (p. ej. ``^MRPL|^MRPS$``).
3. Si una familia aparece como leading edge en ≥ N términos significativos, generar un párrafo
   resumido con sentido biológico para CRC (Warburg, ECM/EMT, proteasoma, inflamación, etc.).

Esta capa es **complementaria** a la narrativa del autor: nunca reemplaza la discusión de tesis.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import pandas as pd


# ---------------------------------------------------------------------------
# Catálogo curado de "temas biológicos" a partir de leading edges
# ---------------------------------------------------------------------------
#
# Cada entrada tiene:
#   - ``patterns``: regex sobre símbolo de gen (HGNC, uppercase). Se evalúan con re.match.
#   - ``label``: nombre humano del módulo.
#   - ``crc_note``: por qué importa en CRC, basado en literatura proteogenómica reciente.
#
# Referencias generales:
#   - CPTAC proteogenómica pan-cancer (Cancer Cell 2023, 10.1016/j.ccell.2023.08.014)
#   - Zhang et al. 2014 (proteogenómica colon/rectal, Nature)
#   - García-Foncillas et al. 2026 (contexto traslacional Colombia, CMS/MSI)
#   - Robles et al. 2024 (meta-análisis reuso proteómica pública CRC, PLOS Comp Biol)

_BIO_THEMES: list[dict[str, object]] = [
    {
        "key": "mito_ribosome",
        "patterns": [r"^MRPL\d", r"^MRPS\d", r"^MTG\d", r"^DAP3$", r"^GADD45GIP1$", r"^GATC$"],
        "label": "Traducción mitocondrial (ribosoma 39S/28S, MRPL/MRPS)",
        "crc_note": (
            "Aumento de la maquinaria de **traducción mitocondrial** (subunidades MRPL/MRPS) en tumor "
            "se ha asociado a reprogramación bioenergética en CCR: las células tumorales preservan o "
            "amplifican la biogénesis mitocondrial incluso cuando bajan los complejos OXPHOS terminales "
            "(NDUF/UQCR). Es coherente con el modelo de **plasticidad metabólica** (no Warburg estricto) "
            "descrito en proteogenómica reciente del CPTAC."
        ),
    },
    {
        "key": "oxphos_complex_i",
        "patterns": [r"^NDUF[ABCS]\d", r"^UQCR[ABCHQ]?\d?$", r"^COX\d", r"^ATP5"],
        "label": "OXPHOS / cadena respiratoria (complejos I–V, NDUF/UQCR/COX/ATP5)",
        "crc_note": (
            "Disminución del **complejo I (NDUFs)**, complejo III (UQCR), citocromo c oxidasa (COX) y "
            "ATP sintasa mitocondrial en tumor es un patrón clásico en CCR: indica desacople parcial de "
            "fosforilación oxidativa y favorece glucólisis aeróbica (efecto Warburg parcial). "
            "Particularmente prominente en CMS3 metabólico (Guinney et al., 2015)."
        ),
    },
    {
        "key": "proteasome",
        "patterns": [r"^PSM[ABCDE]\d", r"^PSME\d", r"^ADRM1$", r"^UCHL5$"],
        "label": "Proteasoma 26S (subunidades PSMA/PSMB/PSMC/PSMD, reguladores)",
        "crc_note": (
            "Sobre-expresión de **subunidades del proteasoma 26S** (PSMA/PSMB catalíticas, PSMC/PSMD "
            "regulatorias) en tumor refleja mayor recambio proteico para sostener proliferación. Es una "
            "vulnerabilidad farmacológica documentada (bortezomib y análogos) pero el éxito en CCR es "
            "limitado en monoterapia — útil como hipótesis de combinación, **nunca como recomendación** "
            "clínica desde un estudio computacional."
        ),
    },
    {
        "key": "mrna_processing",
        "patterns": [r"^CPSF\d", r"^CSTF\d?$", r"^SYMPK$", r"^WDR33$", r"^NUDT21$", r"^CLP1$", r"^PAPOLA$",
                     r"^NCBP\d?$", r"^CDC73$", r"^PAF1$", r"^RPRD\d", r"^MTPAP$", r"^ZCCHC8$"],
        "label": "Procesamiento de mRNA (poliadenilación 3', factores CPSF/CSTF)",
        "crc_note": (
            "Activación coordinada del **complejo de poliadenilación 3'** (CPSF/CSTF) en tumor sugiere "
            "regulación alternativa del extremo 3' UTR — un mecanismo conocido de **APA (alternative "
            "polyadenylation)** que en cáncer suele acortar 3' UTRs y eliminar sitios miRNA, aumentando "
            "estabilidad de mRNAs oncogénicos. Reportado en proteogenómica CPTAC pan-cancer."
        ),
    },
    {
        "key": "ecm_emt_muscle",
        "patterns": [r"^COL\d?\w?\d?", r"^LAMA?\d?$", r"^LAMC?\d?$", r"^FBN\d?$", r"^VCAM\d?$", r"^FN\d?$",
                     r"^MMP\d", r"^FLNC$", r"^FLNA$", r"^DMD$", r"^ACTG\d?$", r"^DES$", r"^MYL?H?\d?$",
                     r"^TAGLN", r"^CNN\d?$", r"^TLN\d?$", r"^VCL$", r"^AHNAK\d?$", r"^SYNM$",
                     r"^DAG\d?$", r"^DTNA$", r"^CAV\d?$", r"^SGC[ABDE]$", r"^MATN\d?$"],
        "label": "Matriz extracelular / mesénquima / músculo liso (CMS4 mesenquimal)",
        "crc_note": (
            "Disminución de proteínas de **matriz extracelular**, **citoesqueleto muscular liso** y "
            "**complejo distrofina-glicoproteína** en tumor refleja pérdida del compartimento estromal "
            "y epitelio de criptas en favor del epitelio tumoral expandido. Cuando estas vías están "
            "**aumentadas** en tumor son típicas de **CMS4 mesenquimal** (TGF-β, angiogénesis), con peor "
            "pronóstico (Guinney et al., 2015; García-Foncillas et al., 2026). Revisar siempre el signo "
            "del NES antes de interpretar."
        ),
    },
    {
        "key": "chaperone_proteostasis",
        "patterns": [r"^HSP\d", r"^HSPA\d", r"^HSP90A[AB]\d?$", r"^DNAJ[ABCD]\d?$", r"^PPIA$",
                     r"^CCT\d", r"^TCP1$", r"^CLU$", r"^CANX$", r"^PDIA\d$", r"^ERP\d", r"^GRPEL\d$",
                     r"^PFDN\d?$", r"^BAG\d?$", r"^AHSA\d?$"],
        "label": "Chaperonas y proteostasis (HSP70/HSP90, CCT/TRiC, PDI/ERP)",
        "crc_note": (
            "Activación de **chaperonas citosólicas (HSP70/HSP90, complejo CCT/TRiC)** y de la **PDI/ERp** "
            "del retículo endoplásmico en tumor indica estrés proteostático sostenido — coherente con "
            "alta tasa de síntesis de proteínas oncogénicas y respuesta a unfolded protein response (UPR)."
        ),
    },
    {
        "key": "acute_phase_plasma",
        "patterns": [r"^LBP$", r"^SAA\d", r"^LRG\d?$", r"^CRP$", r"^SERPINA\d", r"^ORM\d?$",
                     r"^HP$", r"^C[1-9]R?L?$", r"^C9$", r"^CFB$", r"^CFH", r"^MBL\d?$",
                     r"^ITIH\d?$", r"^A2M$", r"^PROC$", r"^F\d?[0-9A-Z]?$", r"^FGA$", r"^FGB$", r"^FGG$",
                     r"^APOA\d?$", r"^APOB$", r"^APOE$", r"^TF$", r"^TFRC$"],
        "label": "Proteínas de fase aguda / cascada del complemento / lipoproteínas (plasma)",
        "crc_note": (
            "Elevación de **proteínas de fase aguda** (LBP, SAA1/2, CRP, SERPINAs, ORM, HP) y activación "
            "de **complemento** (C3, C9, CFB, CFH) en plasma CRC es una **respuesta inflamatoria sistémica** "
            "bien documentada (Nature Commun 2024; Cell Rep Med 2025). No es específica de CRC pero es "
            "una **firma robusta de cáncer activo** útil para discusión exploratoria de biopsia líquida. "
            "Las apolipoproteínas (APOA1/4) frecuentemente **bajan** por reprogramación lipídica hepática."
        ),
    },
    {
        "key": "cytoskeleton_actin",
        "patterns": [r"^ACT[ABG]\d?$", r"^MYO\d", r"^TPM\d?$", r"^TLN\d?$", r"^SEPTIN\d?$", r"^FLNA$",
                     r"^GSN$", r"^SCIN$", r"^PFN\d?$", r"^ARPC\d?$", r"^MARCKS$", r"^CTTN$",
                     r"^PALLD$", r"^NEXN$", r"^SVIL$", r"^SYNPO\d?$"],
        "label": "Citoesqueleto de actina y motilidad (filaminas, septinas, ARPC)",
        "crc_note": (
            "Reorganización del **citoesqueleto de actina** (filaminas, septinas, complejo ARP2/3) es "
            "central en invasión y migración celular. Cambios direccionales aquí son hipótesis razonables "
            "de **fenotipo invasivo / EMT parcial**, especialmente si coocurren con ECM/CMS4."
        ),
    },
    {
        "key": "immune_msi_like",
        "patterns": [r"^IGH[VJ]?\d?[A-Z]?\d?[A-Z]?$", r"^IGK", r"^IGL", r"^HLA[ABCDE]?", r"^CD\d{1,3}$",
                     r"^MSH\d?$", r"^MLH\d?$", r"^PMS\d?$", r"^ITGAL$", r"^PTPRC$"],
        "label": "Inmunidad / infiltrado / mismatch-repair (CMS1 MSI inmune)",
        "crc_note": (
            "Enriquecimiento de **cadenas de inmunoglobulinas** (IGHV/IGK/IGL), HLA y marcadores leucocitarios "
            "(CD\\*) sugiere **infiltrado inmune**. Combinado con deficiencia MSH/MLH/PMS apunta a perfil "
            "**MSI/CMS1**. En CRC colombiano y latinoamericano es relevante por su asociación con respuesta "
            "a immune checkpoint inhibitors (García-Foncillas et al., 2026)."
        ),
    },
    {
        "key": "membrane_trafficking",
        "patterns": [r"^CLTA$", r"^CLTC$", r"^EPN\d$", r"^AP[1-5][SAB]\d?$", r"^SORT\d?$", r"^DAB\d?$",
                     r"^DNM\d?L?$", r"^VPS\d", r"^GGA\d$", r"^SNX\d", r"^PICALM$", r"^EPS\d{1,2}$"],
        "label": "Endocitosis / tráfico vesicular (clatrina, AP-1/2, SNX, VPS)",
        "crc_note": (
            "Activación coordinada de **endocitosis mediada por clatrina** y tráfico endosómico es típica "
            "del reciclaje aumentado de receptores en cáncer (EGFR, TfR, etc.) — hipótesis útil para "
            "discutir el eje **MAPK/EGFR** sin afirmar mecanismo causal."
        ),
    },
]


def _matches_any(gene: str, patterns: Iterable[str]) -> bool:
    g = str(gene).strip().upper()
    return any(re.match(p, g) for p in patterns)


def _collect_leading_edge_genes(gsea_df: pd.DataFrame, *, fdr_threshold: float, top_n_terms: int) -> list[str]:
    if gsea_df is None or gsea_df.empty:
        return []
    fdr_col = next((c for c in ("FDR q-val", "FDR", "fdr") if c in gsea_df.columns), None)
    lead_col = next((c for c in ("Lead_genes", "Leading_genes", "Genes") if c in gsea_df.columns), None)
    nes_col = next((c for c in ("NES", "nes") if c in gsea_df.columns), None)
    if not (fdr_col and lead_col):
        return []
    sub = gsea_df.copy()
    sub["_fdr"] = pd.to_numeric(sub[fdr_col], errors="coerce")
    sub = sub.dropna(subset=["_fdr"])
    sub = sub[sub["_fdr"] < fdr_threshold]
    if sub.empty:
        return []
    if nes_col:
        sub["_abs_nes"] = pd.to_numeric(sub[nes_col], errors="coerce").abs()
        sub = sub.sort_values("_abs_nes", ascending=False)
    sub = sub.head(top_n_terms)
    out: list[str] = []
    for raw in sub[lead_col].astype(str):
        for tok in raw.replace(",", ";").split(";"):
            t = tok.strip().upper()
            if t and t != "NAN":
                out.append(t)
    return out


def interpret_gsea_modules(
    gsea_df: pd.DataFrame,
    *,
    fdr_threshold: float = 0.25,
    top_n_terms: int = 25,
    min_hits_per_theme: int = 4,
) -> list[dict[str, object]]:
    """Devuelve lista de temas detectados (label, hits, ejemplos, crc_note)."""
    genes = _collect_leading_edge_genes(gsea_df, fdr_threshold=fdr_threshold, top_n_terms=top_n_terms)
    if not genes:
        return []
    results: list[dict[str, object]] = []
    for theme in _BIO_THEMES:
        patterns = theme["patterns"]  # type: ignore[index]
        hits = [g for g in set(genes) if _matches_any(g, patterns)]
        if len(hits) >= min_hits_per_theme:
            hits_sorted = sorted(hits)
            results.append(
                {
                    "key": theme["key"],
                    "label": theme["label"],
                    "crc_note": theme["crc_note"],
                    "n_hits": len(hits_sorted),
                    "example_genes": hits_sorted[:12],
                }
            )
    results.sort(key=lambda r: int(r["n_hits"]), reverse=True)
    return results


def render_interpretation_markdown(
    comparison_label: str,
    comp_key: str,
    gsea_df: pd.DataFrame | None,
    *,
    fdr_threshold: float = 0.25,
) -> str:
    """Markdown con interpretación biológica de las leading edges para una comparación."""
    if gsea_df is None or gsea_df.empty:
        return (
            f"### {comparison_label}\n\n"
            f"_Sin GSEA con FDR < {fdr_threshold:g} en `{comp_key}` — no se puede inferir módulos "
            "funcionales en esta comparación._\n\n"
        )
    themes = interpret_gsea_modules(gsea_df, fdr_threshold=fdr_threshold)
    if not themes:
        return (
            f"### {comparison_label}\n\n"
            "_Se detectaron términos GSEA significativos, pero ninguno de los módulos curados "
            "(OXPHOS, traducción mitocondrial, proteasoma, ECM/CMS4, fase aguda, etc.) supera el "
            "umbral mínimo de hits en leading edges. Revisar `*_gsea_go.tsv` manualmente._\n\n"
        )
    parts = [f"### {comparison_label}\n\n"]
    parts.append(
        f"**Módulos biológicos detectados** (leading edges con FDR < {fdr_threshold:g} agregadas "
        "y emparejadas con catálogo curado de CRC):\n\n"
    )
    for t in themes:
        examples = ", ".join(t["example_genes"])  # type: ignore[arg-type]
        parts.append(
            f"- **{t['label']}** — `{int(t['n_hits'])}` genes en leading edges. "
            f"Ejemplos: _{examples}_.\n  > {t['crc_note']}\n\n"
        )
    parts.append(
        "_Estos módulos son **hipótesis cuantitativas generadas por el pipeline**: requieren "
        "validación experimental y deben citarse junto a la literatura específica (Cancer Cell 2023, "
        "PLOS Comp Biol 2024, Nature 2014, García-Foncillas et al. 2026) antes de cualquier afirmación "
        "clínica._\n\n"
    )
    return "".join(parts)


# ---------------------------------------------------------------------------
# Glosario
# ---------------------------------------------------------------------------

_GLOSSARY_ENTRIES: list[tuple[str, str]] = [
    (
        "log2FC (log<sub>2</sub> fold change)",
        "Logaritmo en base 2 del cociente de medias entre grupos. log2FC = +1 → grupo A duplica al B; "
        "log2FC = −1 → grupo A es la mitad. En matrices `log2norm`, las abundancias ya están en escala "
        "log2, así que el log2FC es directamente la diferencia de medias log2.",
    ),
    (
        "p-value vs. padj (FDR Benjamini–Hochberg)",
        "El **p-value** mide la probabilidad de observar la diferencia bajo la hipótesis nula. "
        "El **padj** (q-value) corrige por las miles de pruebas simultáneas (una por proteína) usando "
        "la tasa de descubrimiento falso (FDR-BH). En proteómica de alta dimensión el padj es el "
        "umbral correcto, no el p crudo. Convención habitual: padj < 0.05 = 'significativo'.",
    ),
    (
        "ORA (Over-Representation Analysis, Fisher)",
        "Pregunta: ¿hay **más genes de la lista significativa de los esperados por azar** en cada "
        "término GO? Usa lista binaria (sig vs no-sig). Sensible al tamaño de la lista: si casi todo el "
        "proteoma es significativo (p. ej. tumor vs normal en CPTAC), la ORA pierde poder porque el "
        "'foreground' se parece al universo. Por eso el pipeline siempre reporta también GSEA.",
    ),
    (
        "GSEA (Gene Set Enrichment Analysis, pre-ranked)",
        "Pregunta: ¿los genes de un término GO se **concentran en los extremos** del ranking (más "
        "arriba/abajo que por azar)? Usa el ranking continuo (t-stat o `signed_neglog10_p`), no listas "
        "binarias. Mucho más robusto que ORA cuando la señal es distribuida (caso típico de tumor vs "
        "normal proteómico). Subramanian et al., PNAS 2005.",
    ),
    (
        "NES (Normalized Enrichment Score)",
        "Versión normalizada del enrichment score de GSEA. **NES > 0** → término enriquecido en la cabeza "
        "del ranking (genes UP en el primer grupo). **NES < 0** → enriquecido en la cola (DOWN). "
        "|NES| > 1.5 con FDR q-val < 0.25 es la convención exploratoria; se reportan valores más "
        "permisivos que el padj diferencial porque el universo es de términos GO, no de genes.",
    ),
    (
        "Leading edge",
        "Subconjunto de genes de un término GO que aporta el grueso del enrichment score (genes "
        "más cerca del extremo del ranking dentro de ese término). Es el conjunto **más interpretable**: "
        "si ves MRPL/MRPS en el leading edge, ese término refleja realmente traducción mitocondrial.",
    ),
    (
        "Reducción semántica GO (`*_go_reduced.tsv`)",
        "Muchos términos GO se solapan (p. ej. 'mitochondrial translation' incluye 'mitoribosome'). "
        "El pipeline agrupa términos con Jaccard ≥ 0.45 sobre sus genes y conserva el representante con "
        "menor p — evita inflar listas con redundancia.",
    ),
    (
        "ROC AUC",
        "Área bajo la curva ROC. 0.5 = clasificador al azar; 1.0 = perfecto. En proteómica con N "
        "moderado y muchas features, AUC > 0.9 en hold-out puede ser optimista por azar — por eso el "
        "pipeline reporta nested CV y transferencia externa.",
    ),
    (
        "Nested cross-validation",
        "Doble bucle de CV: el **externo** estima el AUC de generalización; el **interno** elige "
        "hiperparámetros con GridSearchCV. Evita la fuga de datos típica de elegir hiperparámetros y "
        "evaluar en el mismo split (que infla el AUC). Es la métrica honesta para el capítulo de "
        "resultados ML.",
    ),
    (
        "Transferencia externa (CPTAC → PRIDE)",
        "Entrena en una cohorte (CPTAC), evalúa en otra independiente (PRIDE) **sin reajustar nada**. "
        "Es el experimento que pone a prueba la **transportabilidad** del modelo a dominio distinto "
        "(tecnología, lab, batch). Un AUC modesto aquí (p. ej. 0.65-0.80 frente a 0.95 dentro de CPTAC) "
        "es esperable y honesto: refleja el batch/tech shift, no fracaso biológico.",
    ),
    (
        "Calibración (sigmoid / isotonic)",
        "Las probabilidades de scikit-learn no siempre son fieles: un RF puede dar 0.9 cuando la "
        "frecuencia real es 0.7. La calibración (Platt sigmoid o isotonic) ajusta esa relación usando "
        "CV en train. En transfer, calibrar a veces baja el AUC pero mejora la utilidad clínica de "
        "umbrales — útil cuando se quiere decir 'probabilidad estimada de ser tumor', no solo ranking.",
    ),
    (
        "Bootstrap estabilidad de features (L1)",
        "Re-muestrea N veces el train, entrena un modelo L1 (sparse linear) y cuenta cuántas veces "
        "cada gen entra en el top-k por |coef|. Una **frecuencia alta** (p. ej. > 0.7) indica un gen "
        "robusto al re-muestreo; baja indica que su selección depende del split concreto.",
    ),
    (
        "TPS (Translational Priority Score)",
        "Score compuesto que combina: magnitud diferencial, FDR, centralidad en red, importancia ML "
        "(RF, L1, EN, SHAP), pertenencia a capas YAML (CMS, MAPK, módulos funcionales, tiers de "
        "evidencia literaria). **No es una decisión clínica**: es un **ranking explícito** de hipótesis "
        "investigativas, defendible ante un comité de tesis porque cada componente es trazable.",
    ),
    (
        "Firma núcleo CPTAC ∩ PRIDE",
        "Genes/proteínas con `padj < 0.05` y mismo signo de log2FC en **ambas cohortes** tras alinear "
        "por UniProt canónico. La fracción con mismo signo (`same_sign_fraction`) cuantifica "
        "reproducibilidad direccional; correlación Pearson/Spearman de log2FC cuantifica concordancia "
        "de magnitud (típicamente modesta entre TMT y LFQ, ~0.1–0.3).",
    ),
    (
        "TMT vs LFQ (¿por qué bajan las correlaciones entre cohortes?)",
        "**TMT** (Tandem Mass Tag, CPTAC) cuantifica multiplexado en un mismo MS; tiene buena "
        "precisión relativa pero compresión de ratios. **LFQ** (label-free, PRIDE/MaxQuant) cuantifica "
        "por intensidad de péptidos independientes; mayor varianza pero menos sesgo. Comparar log2FC "
        "entre ambos siempre baja la correlación absoluta — lo importante es la **concordancia de "
        "dirección** (signo) y la reproducibilidad de los **términos GO/GSEA**, no la magnitud bruta.",
    ),
    (
        "padj < 0.25 en GSEA (¿no es muy laxo?)",
        "Sí, intencionalmente. La convención original de Subramanian et al. (PNAS 2005) usa FDR q-val < "
        "0.25 como umbral **exploratorio** porque la corrección se hace sobre miles de términos GO y "
        "GSEA es conservador por diseño (permutaciones). Para módulos individuales se discute |NES| "
        "alto + leading edge biológicamente coherente, no el q-val aislado.",
    ),
    (
        "Catastrophic cancellation (warning SciPy)",
        "Aparece cuando dos números casi iguales se restan y se pierde precisión. En el pipeline ocurre "
        "para proteínas con varianza ~0 dentro de un grupo (todas las muestras con el mismo valor tras "
        "imputación). No afecta al contraste global; sí marca esas filas como poco informativas — el "
        "ranking las pone al final por padj=NaN o p~1.",
    ),
    (
        "Modo demo vs corrida real",
        "`--demo` genera matrices sintéticas con genes `GENE001…` y mapa GO mínimo para validar la "
        "**lógica** del pipeline sin descargar nada. Una corrida 'real' (sin `--demo`) lee `datasets.yaml` "
        "y `download_real` baja de PDC/PRIDE. El audit JSON (`*_gsea_input_*.json`) reporta el flag "
        "`contains_demo_gene_ids` para evitar publicar resultados demo como reales.",
    ),
    (
        "GOA (Gene Ontology Annotations)",
        "Anotación oficial gen→término GO. El pipeline usa `goa_human.gaf.gz` del Gene Ontology "
        "Consortium (~14 MiB, ~38 800 genes humanos, ~17 900 términos). Regenerar con "
        "`python scripts/fetch_go_human_annotation.py`. Sin GOA real (solo demo de 10 genes), ORA y "
        "GSEA producen 0 términos: es bug de cobertura, no biológico.",
    ),
    (
        "Mapping efectivo (`protein_mapping_master.tsv` + bulk UniProt)",
        "Tabla canónica gen↔UniProt usada para colapsar isoformas y para puentear CPTAC (genes) con "
        "PRIDE (UniProt). El bulk humano reviewed (~20 000 entradas) se descarga con "
        "`scripts/fetch_uniprot_human_reviewed_gene_map.py`. Sin él, la intersección CPTAC∩PRIDE se "
        "subestima drásticamente.",
    ),
    (
        "Workspace de datos (`paths.data_workspace`)",
        "Carpeta externa (típicamente disco USB/SSD) donde residen `data/`, `results/`, `reports/`, "
        "`metadata/` y `docs/web/`. Permite mantener el clon git ligero sin perder reproducibilidad. "
        "El pipeline sincroniza metadata desde el repo al workspace en cada `run` (mtime-based). "
        "Si el workspace no es accesible, la CLI cae a rutas del repo con un aviso.",
    ),
    (
        "run_manifest.json",
        "Manifest que se reinicializa al inicio de cada `oncogo-crc run`. Contiene `run_id`, timestamp, "
        "lista de etapas ejecutadas/omitidas, hashes SHA-256 de outputs clave y razones de skip. "
        "El publicador web (`publish_results`) **lee este manifest** y se niega a publicar resultados "
        "de etapas no ejecutadas en la corrida actual — evita publicar artefactos viejos como nuevos.",
    ),
    (
        "stage_bin (temprano vs avanzado)",
        "Agrupación binaria del estadio patológico AJCC: **early** (I–II) vs **advanced** (III–IV), "
        "derivada del texto de `stage` con tokens configurables en `clinical_filters.stage_mapping`. "
        "Alimenta la comparación `cptac_stage_advanced_vs_early`.",
    ),
    (
        "stage_group (I–IV)",
        "Grupo AJCC romano extraído automáticamente del campo `stage` (p. ej. «Stage IIIA» → `III`). "
        "Permite tabular estadios sin perder granularidad del Excel clínico PDC.",
    ),
    (
        "hispanic_latino_bin",
        "Estrato derivado de `ethnicity` según categorías NIH: `hispanic_latino`, `non_hispanic_latino` "
        "o `unknown`. **No** identifica nacionalidad colombiana ni residencia en Latinoamérica; en CPTAC "
        "refleja autodeclaración en cohortes reclutadas principalmente en EE.UU. Ver "
        "[Demografía clínica](clinical_demographics.md).",
    ),
    (
        "geographic_region_bin",
        "Región inferida desde `country` y, si falta, desde `enrollment_site`: `latin_america`, "
        "`united_states`, `other` o `unknown`. Útil para transparencia; en muchos estudios CPTAC el "
        "país no está informado y la variable queda mayoritariamente `unknown`.",
    ),
]


def render_glossary_markdown() -> str:
    parts = [
        "# Glosario — OncoGO-CRC\n\n",
        "Glosario de los términos estadísticos, bioinformáticos y traslacionales más usados en este sitio. "
        "Está pensado para acompañar la lectura de las páginas **Resultados base**, **Síntesis para tesis** "
        "e **Interpretación biológica** sin necesidad de saltar a textbooks externos.\n\n",
        "_Las definiciones priorizan precisión operativa sobre exhaustividad; las referencias completas "
        "están en `docs/RESULTADOS_Y_LIMITACIONES_TESIS.md`._\n\n",
    ]
    for term, body in _GLOSSARY_ENTRIES:
        parts.append(f"## {term}\n\n{body}\n\n")
    return "".join(parts)
