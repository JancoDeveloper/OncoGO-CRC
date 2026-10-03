"""Figuras y PDF explicativo de la corrida certificada. No recalcula el contraste primario."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image,
    ListFlowable,
    ListItem,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

REPO = Path(__file__).resolve().parents[1]
WORKSPACE_DEFAULT = Path("/Volumes/JANCODATA/OncoGO-CRC")
# Copia local de los insumos. Son tablas pequeñas: permiten rehacer el PDF sin el disco externo.
CACHE_DEFAULT = REPO / "docs" / "tesis_datos"
PDF_DEFAULT = REPO / "docs" / "OncoGO-CRC_tesis_explicada.pdf"
FONT = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT_B = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"

# Insumos, relativos a la raíz de datos. Nada de esto se recalcula aquí.
INPUTS = (
    "results/certified/inference_audit.json",
    "results/certified/cptac_unshared_treat.tsv",
    "results/certified/cptac_unshared_treat_stroma_adjusted.tsv",
    "results/certified/cptac_pride_overlap.tsv",
    "results/certification/certificado_resumen.json",
    "results/secondary/secondary_audit.json",
    "results/secondary/variables_inventario.tsv",
    "results/axes/axes_audit.json",
    "results/axes/panel_externo_plasma.tsv",
    "results/axes/pride_paired_tumor_vs_normal_treat.tsv",
    "results/axes/pride_paired_polyp_vs_normal_treat.tsv",
    "results/axes/tejido_a_plasma.tsv",
)

ROOT = CACHE_DEFAULT
CERT = ROOT / "results" / "certified"
SEC = ROOT / "results" / "secondary"
AXES = ROOT / "results" / "axes"
FIG = REPO / "docs" / "tesis_figuras"
PDF_PATH = PDF_DEFAULT


def set_paths(root: Path, pdf_path: Path, fig_dir: Path) -> None:
    global ROOT, CERT, SEC, AXES, FIG, PDF_PATH
    ROOT = root
    CERT = root / "results" / "certified"
    SEC = root / "results" / "secondary"
    AXES = root / "results" / "axes"
    FIG = fig_dir
    PDF_PATH = pdf_path


def missing_inputs(root: Path) -> list[str]:
    return [rel for rel in INPUTS if not (root / rel).exists()]


def sync_cache(source: Path, cache: Path) -> int:
    """Copia los insumos del espacio de datos a la caché del repo. No los modifica."""
    copied = 0
    for rel in INPUTS:
        src = source / rel
        if not src.exists():
            continue
        dst = cache / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        copied += 1
    return copied


def esc(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def comma(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}".replace(".", ",").replace("-", "\u2212")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def setup_font() -> None:
    font_manager.fontManager.addfont(FONT)
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["axes.unicode_minus"] = False
    pdfmetrics.registerFont(TTFont("Arial", FONT))
    bold = FONT_B if Path(FONT_B).exists() else FONT
    pdfmetrics.registerFont(TTFont("Arial-Bold", bold))


def fig_architecture(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.2, 5.6))
    ax.set_xlim(0, 10.2)
    ax.set_ylim(0, 5.6)
    ax.axis("off")
    boxes = [
        (0.3, 3.7, 2.2, 1.15, "Matrices públicas\nCPTAC, PRIDE, plasma\n(no se editan)", "#e7f0f8"),
        (2.9, 3.7, 2.2, 1.15, "Una medición\npor espécimen\n(péptido no compartido)", "#e7f0f8"),
        (5.5, 3.7, 2.2, 1.15, "Ausencias = ausencias\nsin imputar\nya en log2 si lo era", "#e7f0f8"),
        (8.1, 3.7, 1.8, 1.15, "treat\ndentro de\ncada cohorte", "#d9efe4"),
        (2.4, 1.45, 2.8, 1.15, "Cruce por UniProt\nsigno, vía, péptido", "#fff4d6"),
        (6.2, 1.45, 2.8, 1.15, "Certificado\nR1–R5 por separado\nH1 no se pliega", "#fde8d0"),
    ]
    for x, y, w, h, text, color in boxes:
        ax.add_patch(
            FancyBboxPatch(
                (x, y), w, h,
                boxstyle="round,pad=0.04,rounding_size=0.12",
                facecolor=color, edgecolor="#1f3a5f", linewidth=1.1,
            )
        )
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.5, color="#1a1a1a")
    for x1, x2 in ((2.5, 2.9), (5.1, 5.5), (7.7, 8.1)):
        ax.add_patch(FancyArrowPatch((x1, 4.27), (x2, 4.27), arrowstyle="-|>", mutation_scale=12, color="#1f3a5f"))
    ax.add_patch(FancyArrowPatch((9.0, 3.7), (3.8, 2.6), arrowstyle="-|>", mutation_scale=12, color="#1f3a5f"))
    ax.add_patch(FancyArrowPatch((5.2, 2.02), (6.2, 2.02), arrowstyle="-|>", mutation_scale=12, color="#1f3a5f"))
    ax.text(0.35, 0.35, "No entra al certificado: fusionar cohortes, imputar, segundo log2, ni un modelo con AUC interno de 1.", fontsize=8, color="#444444")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_volcano(path: Path) -> None:
    df = pd.read_csv(CERT / "cptac_unshared_treat.tsv", sep="\t")
    x = df["log2fc"].astype(float)
    q = df["padj"].astype(float).clip(lower=1e-300)
    y = -np.log10(q)
    sig = df["significant"].astype(bool) & x.notna()
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.scatter(x[~sig], y[~sig], s=8, c="#b0b7c3", linewidths=0, label="no pasa treat")
    ax.scatter(x[sig], y[sig], s=12, c="#0f6e6e", linewidths=0, label="q < 0,05 y |log2FC| > 0,585")
    ax.axvline(0.585, color="#888888", lw=0.7, ls="--")
    ax.axvline(-0.585, color="#888888", lw=0.7, ls="--")
    labels = df.loc[sig].nsmallest(6, "padj")
    for _, row in labels.iterrows():
        ax.annotate(str(row["protein_id"]), (row["log2fc"], -np.log10(max(row["padj"], 1e-300))), fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("log2FC tumor − normal (péptido no compartido)")
    ax.set_ylabel("−log10 (q de treat)")
    ax.set_title("CPTAC, 94 pares, valores observados")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_scatter(path: Path) -> None:
    df = pd.read_csv(CERT / "cptac_pride_overlap.tsv", sep="\t")
    both = df.loc[df["_merge"].astype(str) == "both"].copy()
    x = pd.to_numeric(both["log2fc_cptac"], errors="coerce")
    y = pd.to_numeric(both["log2fc_pride"], errors="coerce")
    ok = x.notna() & y.notna()
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    ax.scatter(x[ok], y[ok], s=10, c="#5c6b7a", linewidths=0, alpha=0.45)
    ax.axhline(0, color="#bbbbbb", lw=0.6)
    ax.axvline(0, color="#bbbbbb", lw=0.6)
    lim = 4
    ax.plot([-lim, lim], [-lim, lim], color="#c47b2b", lw=0.8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2FC CPTAC")
    ax.set_ylabel("log2FC PRIDE")
    ax.set_title("Misma proteína, dos cohortes")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_counts(path: Path, audit: dict, secondary: dict) -> None:
    labels = ["CPTAC\nprimario", "PRIDE\nréplica", "Mujer", "Hombre", "Estadio", "Plasma"]
    values = [
        audit["cptac_n_significant"],
        audit["pride_n_significant"],
        secondary["contrasts"]["female"]["n_significant"],
        secondary["contrasts"]["male"]["n_significant"],
        secondary["contrasts"]["stage"]["n_significant"],
        secondary["contrasts"]["plasma"]["n_significant"],
    ]
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    colors_bar = ["#0f6e6e", "#1f3a5f", "#3d6b8a", "#3d6b8a", "#8a8f98", "#c47b2b"]
    ax.bar(labels, values, color=colors_bar)
    for i, v in enumerate(values):
        ax.text(i, v + 8, str(v), ha="center", fontsize=8)
    ax.set_ylabel("Proteínas con q < 0,05 en treat")
    ax.set_title("El mismo umbral en cada contraste que sí tenía datos")
    ax.set_ylim(0, max(values) * 1.15)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_stroma(path: Path) -> None:
    df = pd.read_csv(CERT / "cptac_unshared_treat_stroma_adjusted.tsv", sep="\t")
    base = pd.read_csv(CERT / "cptac_unshared_treat.tsv", sep="\t")[["protein_id", "log2fc"]]
    base = base.rename(columns={"log2fc": "lfc_crudo"})
    m = df.merge(base, on="protein_id", how="inner")
    m = m.loc[m["certified_unadjusted"].astype(bool)]
    keep = m["significant"].astype(bool)
    fig, ax = plt.subplots(figsize=(6.6, 6.0))
    ax.scatter(m.loc[~keep, "lfc_crudo"], m.loc[~keep, "log2fc"], s=14, c="#c47b2b", linewidths=0, alpha=0.65, label="se cae al igualar composición")
    ax.scatter(m.loc[keep, "lfc_crudo"], m.loc[keep, "log2fc"], s=14, c="#0f6e6e", linewidths=0, alpha=0.7, label="sobrevive")
    lim = 4
    ax.plot([-lim, lim], [-lim, lim], color="#999999", lw=0.8)
    for v in (0.585, -0.585):
        ax.axhline(v, color="#bbbbbb", lw=0.6, ls="--")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2FC sin ajustar")
    ax.set_ylabel("log2FC con el par igualado en estroma")
    ax.set_title("El efecto encoge hacia cero al quitar la composición")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_adenoma(path: Path) -> None:
    tumor = pd.read_csv(AXES / "pride_paired_tumor_vs_normal_treat.tsv", sep="\t")
    polyp = pd.read_csv(AXES / "pride_paired_polyp_vs_normal_treat.tsv", sep="\t")
    t = tumor.loc[tumor["significant"].astype(bool), ["protein_id", "log2fc"]].rename(columns={"log2fc": "tumor"})
    p = polyp[["protein_id", "log2fc"]].rename(columns={"log2fc": "polipo"})
    m = t.merge(p, on="protein_id", how="inner").dropna()
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    ax.scatter(m["tumor"], m["polipo"], s=16, c="#1f3a5f", linewidths=0, alpha=0.6)
    lim = float(np.nanmax(np.abs(np.concatenate([m["tumor"], m["polipo"]])))) * 1.05
    ax.plot([-lim, lim], [-lim, lim], color="#c47b2b", lw=0.9, label="mismo tamaño")
    ax.axhline(0, color="#bbbbbb", lw=0.6)
    ax.axvline(0, color="#bbbbbb", lw=0.6)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2FC tumor − normal (mismo paciente)")
    ax.set_ylabel("log2FC pólipo − normal (mismo paciente)")
    ax.set_title("El adenoma ya apunta a donde apunta el tumor, con menos fuerza")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_gap(path: Path) -> None:
    m = pd.read_csv(AXES / "tejido_a_plasma.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    same = np.sign(m["log2fc"]) == np.sign(m["lfc_plasma"])
    ax.scatter(m.loc[same, "log2fc"], m.loc[same, "lfc_plasma"], s=22, c="#0f6e6e", linewidths=0, alpha=0.7, label="mismo signo")
    ax.scatter(m.loc[~same, "log2fc"], m.loc[~same, "lfc_plasma"], s=22, c="#c47b2b", linewidths=0, alpha=0.7, label="signo opuesto")
    ax.axhline(0, color="#999999", lw=0.7)
    ax.axvline(0, color="#999999", lw=0.7)
    ax.set_xlabel("log2FC en tejido (tumor − adyacente)")
    ax.set_ylabel("log2FC en suero (cáncer − sano)")
    ax.set_title("La dirección en tejido no predice la dirección en sangre")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_models(path: Path, axes_audit: dict) -> None:
    models = axes_audit["model_selection"]["models"]
    names = list(models.keys())
    aucs = [models[n]["auc_external"] for n in names]
    lows = [models[n]["ci_low"] for n in names]
    highs = [models[n]["ci_high"] for n in names]
    y = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    for i, n in enumerate(names):
        ax.plot([lows[i], highs[i]], [y[i], y[i]], color="#1f3a5f", lw=2.2)
        ax.plot(aucs[i], y[i], "o", color="#0f6e6e", ms=7)
        ax.text(highs[i] + 0.01, y[i], f"{aucs[i]:.3f}".replace(".", ","), va="center", fontsize=8.5)
    ax.axvline(0.5, color="#c47b2b", lw=0.9, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("_", " ") for n in names], fontsize=9)
    ax.set_xlim(0.35, 1.08)
    ax.set_xlabel("AUC en la cohorte externa, con intervalo del 95 % por bootstrap")
    ax.set_title("Los intervalos se solapan: no hay modelo ganador")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_subindices(path: Path, cert: dict) -> None:
    labels, values, states = [], [], []
    for sub in cert["subindices"]:
        labels.append(sub["code"].replace("_", "\n"))
        values.append(0.0 if sub["value"] is None else float(sub["value"]))
        states.append(sub["status"])
    fig, ax = plt.subplots(figsize=(7.6, 3.6))
    colors_bar = ["#0f6e6e" if st == "calculado" else "#cfd4da" for st in states]
    ax.bar(labels, values, color=colors_bar)
    for i, (v, st) in enumerate(zip(values, states)):
        ax.text(i, v + 0.02, "no calculable" if st != "calculado" else f"{v:.3f}".replace(".", ","), ha="center", fontsize=8)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("valor del subíndice")
    ax.set_title("Certificado R: cuatro de cinco subíndices con número, sin colapsarlos")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def styles() -> dict:
    return {
        "h1": ParagraphStyle("h1", fontName="Arial-Bold", fontSize=16, leading=20, textColor=colors.HexColor("#1f3a5f"), spaceAfter=8, spaceBefore=12),
        "h2": ParagraphStyle("h2", fontName="Arial-Bold", fontSize=12.5, leading=16, textColor=colors.HexColor("#0f6e6e"), spaceAfter=6, spaceBefore=10),
        "body": ParagraphStyle("body", fontName="Arial", fontSize=10.5, leading=14.5, alignment=TA_JUSTIFY, spaceAfter=8),
        "cover": ParagraphStyle("cover", fontName="Arial-Bold", fontSize=18, leading=23, alignment=TA_CENTER, textColor=colors.HexColor("#1f3a5f"), spaceAfter=10),
        "sub": ParagraphStyle("sub", fontName="Arial", fontSize=11, leading=15, alignment=TA_CENTER, textColor=colors.HexColor("#333333"), spaceAfter=6),
        "cap": ParagraphStyle("cap", fontName="Arial", fontSize=8.5, leading=11, alignment=TA_LEFT, textColor=colors.HexColor("#444444"), spaceBefore=2, spaceAfter=8),
        "small": ParagraphStyle("small", fontName="Arial", fontSize=9, leading=12, alignment=TA_LEFT),
        "smallw": ParagraphStyle("smallw", fontName="Arial", fontSize=9, leading=12, alignment=TA_LEFT, textColor=colors.white),
    }


def bullet(text: str, style: ParagraphStyle) -> ListItem:
    return ListItem(Paragraph(esc(text), style), leftIndent=12)


def build() -> None:
    setup_font()
    FIG.mkdir(parents=True, exist_ok=True)
    audit = load_json(CERT / "inference_audit.json")
    cert = load_json(ROOT / "results" / "certification" / "certificado_resumen.json")
    secondary = load_json(SEC / "secondary_audit.json")
    axes = load_json(AXES / "axes_audit.json")
    inventory = pd.read_csv(SEC / "variables_inventario.tsv", sep="\t")
    panel = pd.read_csv(AXES / "panel_externo_plasma.tsv", sep="\t")
    fig_architecture(FIG / "arquitectura.png")
    fig_volcano(FIG / "volcano_cptac.png")
    fig_scatter(FIG / "scatter_cptac_pride.png")
    fig_counts(FIG / "conteos.png", audit, secondary)
    fig_stroma(FIG / "estroma.png")
    fig_adenoma(FIG / "adenoma.png")
    fig_gap(FIG / "tejido_plasma.png")
    fig_models(FIG / "modelos.png", axes)
    fig_subindices(FIG / "subindices.png", cert)
    s = styles()
    story: list = []

    def h1(text: str) -> None:
        story.append(Paragraph(esc(text), s["h1"]))

    def h2(text: str) -> None:
        story.append(Paragraph(esc(text), s["h2"]))

    def p(text: str) -> None:
        story.append(Paragraph(esc(text), s["body"]))

    def img(name: str, width: float) -> None:
        story.append(Image(str(FIG / name), width=width, height=width * 0.62))

    story.append(Spacer(1, 2.2 * cm))
    story.append(Paragraph("OncoGO-CRC", s["cover"]))
    story.append(Paragraph("Qué se busca, cómo se analiza y qué se puede afirmar", s["cover"]))
    story.append(Spacer(1, 0.4 * cm))
    story.append(Paragraph(
        "Documento explicativo de la corrida del 1 de octubre de 2026, escrito para quien no conoce el tema.",
        s["sub"],
    ))
    story.append(Paragraph(
        "Los números salen de las matrices públicas ya descargadas en el disco JANCODATA. Esas matrices no se editaron.",
        s["sub"],
    ))
    story.append(PageBreak())

    h1("1. Qué se busca")
    p("El proyecto no busca un fármaco ni un examen listo para el hospital. Busca saber, con datos públicos de proteómica de cáncer colorrectal, si una proteína que cambia entre tumor y tejido adyacente en un estudio cambia en la misma dirección en otro estudio hecho con otra tecnología. El producto es una metodología y una tabla de hipótesis que un laboratorio puede auditar.")
    p("La hipótesis central, llamada H1, decía dos cosas a la vez: que las proteínas coincidirían poco (menos del 70 % en el signo del cambio) y que las vías biológicas coincidirían mucho (más del 80 %). Las dos cláusulas van juntas. Si la de proteína no se cumple, H1 no se sostiene, aunque alguna vía coincida.")
    p("Hay otras preguntas declaradas de antemano. H0 pregunta si dos estudios públicos comparten ficheros sin decirlo. H2 pregunta si una red aporta más que un modelo simple. H4 pregunta si exigir que la proteína sea medible en más de una plataforma elimina candidatos o si casi todos ya lo eran. Esas preguntas solo se responden con el dato que de verdad existe. Si el dato no está, la pregunta queda abierta.")

    h1("2. Conceptos, en el orden en que se usan")
    h2("Proteína, péptido y por qué el péptido importa")
    p("Una proteína es la molécula que la célula usa para hacer el trabajo. Un péptido es un fragmento de esa proteína, que es lo que el espectrómetro de masas mide. Si un péptido aparece en dos proteínas distintas, verlo subir no dice cuál de las dos cambió. Un péptido no compartido pertenece a una sola proteína. En la matriz de CPTAC esa medición tiene un nombre propio: Unshared Log Ratio. Esa es la columna que decide el contraste. La otra columna, Log Ratio, mezcla péptidos compartidos y no compartidos; se guardó como comprobación y no se usó para aceptar o rechazar H1.")
    p("El archivo de PRIDE, proteinGroups de MaxQuant, trae una columna llamada Unique peptides: cuántos péptidos de esa proteína no están en otra. Eso permite preguntar si la proteína se vio con al menos un péptido propio en la segunda cohorte. No permite afirmar que la secuencia medida fue la misma en las dos plataformas, porque los ficheros públicos de esta corrida no traen la secuencia alineada lado a lado.")
    h2("log2FC, q-valor y treat")
    p("El log2FC es el cambio de abundancia en logaritmo base 2. Un log2FC de 1 significa el doble. Un log2FC de 0,585 significa un 50 % más o un 50 % menos (1,5 veces). El diseño fijó 0,585 como efecto mínimo antes de mirar estas tablas.")
    p("Con miles de proteínas, un p-valor pequeño aparece por azar. El q-valor de Benjamini–Hochberg controla la proporción esperada de llamadas falsas entre las que se declaran significativas. Aquí el umbral es 0,05.")
    p("Un cambio diminuto puede tener un q-valor pequeño si hay muchas muestras. treat no pregunta si el cambio es distinto de cero. Pregunta si es más grande que 0,585. El p-valor ya incluye ese umbral. No se filtra la tabla después. Por eso un contraste que antes marcaba miles de proteínas ahora marca cientos, o ninguna.")
    h2("Por qué no se rellenan los huecos y por qué no se juntan los estudios")
    p("Cuando una proteína no se detecta en una muestra, la celda está vacía. Rellenarla con el valor de las muestras vecinas inventa una medición. El contraste certificado no rellena. Si un par tumor–normal tiene un lado vacío, ese par no entra en el promedio de esa proteína.")
    p("CPTAC midió con TMT: cada muestra se compara con una mezcla de referencia y el número ya es un log-ratio. PRIDE midió con LFQ: el número es una intensidad, y solo se toma log2 cuando es mayor que cero. El cero de MaxQuant significa que no se cuantificó, no que la abundancia sea cero. Juntar las dos matrices en un solo modelo mezclaría la tecnología con la biología. Cada cohorte se contrasta por dentro. Después se comparan los resultados, proteína por proteína, mediante el identificador UniProt.")
    h2("Tejido adyacente, estroma y pureza")
    p("El normal de CPTAC es tejido de al lado del tumor, no sangre de una persona sana. Ese tejido trae músculo, matriz y células que no son el epitelio del tumor. Un cambio puede ser del epitelio del tumor o puede ser que la mezcla de tejidos sea distinta a los dos lados. La tabla clínica pública de PDC no trae el porcentaje de pureza del patólogo.")
    p("La literatura reciente de proteómica de colon dice que esa mezcla no es un detalle. Los subtipos proteómicos de CCR se separan por contenido estromal e inmune, y proteínas como CAVIN1 y PRELP reflejan células mesenquimales más que biología del tumor. En esta corrida las dos aparecen muy disminuidas en tumor: CAVIN1 con log2FC −1,75 y PRELP con −1,55, las dos significativas. Eso no es un hallazgo del tumor: es la firma de que el tejido de al lado tiene más estroma.")
    p("Por eso se construyó un sustituto de composición con los datos que hay: para cada espécimen, el promedio de las proteínas anotadas en matriz extracelular (término GO:0031012), estandarizadas proteína por proteína. Son 282 proteínas de esta matriz. Después, para cada par tumor–normal se calcula el contraste estromal del par y se usa como covariable. El efecto que se reporta es el que quedaría si el par no tuviera diferencia de composición. Esto no es una pureza medida, y las proteínas del índice salen del mismo proteoma que se está contrastando: es un análisis de sensibilidad, no un patrón de oro.")

    h1("3. Datos que sí entraron, y los que no")
    p("Descubrimiento: PDC000116, proteoma de colon de CPTAC, tecnología TMT. Se usaron 94 pares de la misma persona, columna de péptido no compartido, valores tal como vienen en el fichero. Las filas Mean, Median y StdDev no son proteínas y se excluyeron del contraste. No se volvió a aplicar log2: el fichero ya está en esa escala.")
    p("Réplica: PXD034575, cohorte 1 de Tanaka 2024, LFQ de MaxQuant, 16 tumores y 16 normales. Los pólipos de ese estudio no entran en el contraste tumor frente a normal. En el diseño del proyecto esta cohorte es un control de independencia del mismo grupo de laboratorio, no la réplica primaria PXD019504. El cruce dice lo que esta cohorte muestra. No sustituye el diseño de cuatro plataformas.")
    p("Plasma: PXD055218, suero humano medido por DIA, 20 muestras de cáncer y 20 de personas sanas. Es otro líquido y otro proteoma. Se analiza aparte y no se usa para H1.")
    bullets = [
        "Edad: las 198 filas de péptido no compartido tienen la edad vacía. No hay contraste de edad ni de cáncer de inicio temprano.",
        "MSI, subtipo CMS, metástasis, respuesta al tratamiento, país y centro de reclutamiento: todos desconocidos o vacíos.",
        "Sitio del tumor: todo lo informado es colon. No hay recto con qué comparar.",
        "Lote: un solo lote TMT. No hay un segundo lote que corregir.",
        "Etnicidad hispana en tumor: 7 frente a 87. Por debajo del mínimo de 15 por grupo. No se contrasta.",
        "Raza del censo: no es ancestría genética, y los grupos pequeños tienen 15 y 7 tumores. No se contrasta.",
        "Pureza del patólogo: no está en la tabla. Se usa un sustituto estromal del propio proteoma, declarado como tal.",
        "Pólipos de PXD034575: 16, emparejados con el normal y el tumor del mismo paciente. Sí se contrastan, en su propia sección.",
    ]
    story.append(ListFlowable([bullet(t, s["body"]) for t in bullets], bulletType="bullet", start="•"))
    p("El mínimo de 15 observaciones por lado se fijó para los contrastes secundarios antes de leer sus q-valores. Una variable por debajo de ese mínimo se describe y no se modela.")

    h1("4. Cómo está armado el análisis")
    p("El código vive en el repositorio. Los ficheros pesados viven en /Volumes/JANCODATA/OncoGO-CRC. La figura siguiente es el camino que produce el certificado. El contraste antiguo, el que aplicaba un segundo logaritmo y trataba las dos columnas del mismo espécimen como si fueran dos personas, queda fuera de este camino.")
    story.append(Image(str(FIG / "arquitectura.png"), width=16.5 * cm, height=9.0 * cm))
    story.append(Paragraph("Figura 1. Camino del contraste certificado. Las matrices de origen no se modifican.", s["cap"]))
    p("La t moderada estima la varianza de cada proteína y la acerca un poco a la varianza típica del experimento. En CPTAC ese ajuste fue pequeño (grados de libertad a priori 3,22). Si el ajuste no se puede estimar, el método vuelve a una t ordinaria con el mismo umbral de efecto. No es una red neuronal ni un modelo que se elige por tener el AUC más alto.")
    p("Se usaron estas técnicas, y no otras, por una razón concreta. Una t contra cero, como la de la primera corrida, declara significativo un cambio minúsculo. Un bosque aleatorio dentro de CPTAC separa tumor y normal con AUC 1 porque las dos clases ya están muy separadas en la misma matriz; eso no elige un modelo ni demuestra que la firma viaje. Juntar cohortes o corregir lote con ComBat para inferir mezclaría plataforma y biología. msqrob2, el modelo pensado para TMT con el plex como efecto aleatorio, necesita las intensidades de cada canal; el fichero público trae el log-ratio de proteína ya resumido, y fabricar los canales sería inventar datos. Una red de grafos solo se interpretaría si ganara a un modelo lineal, a un bosque aleatorio y a un nulo que conserva el grado de cada nodo. Ese nulo no está calculado, así que no se interpreta una red.")

    h1("5. Resultado primario")
    p(f"En CPTAC, treat marca {audit['cptac_n_significant']} proteínas de {audit['cptac_n_tested']} contrastadas, con {audit['cptac_pairs']} pares. De esas, la mayoría está más baja en tumor que en el tejido de al lado. Las primeras por q-valor son SCGN, NCAM1, CHGA, CADM3 y ADH1B, todas disminuidas. En PRIDE, el mismo treat marca {audit['pride_n_significant']} proteínas de {audit['pride_n_tested']}. Solo 28 proteínas pasan el umbral en las dos cohortes a la vez.")
    story.append(Image(str(FIG / "volcano_cptac.png"), width=14.5 * cm, height=10.2 * cm))
    story.append(Paragraph("Figura 2. Cada punto es una proteína en CPTAC. La línea punteada es el efecto mínimo de 0,585.", s["cap"]))
    p(f"Entre las proteínas significativas en CPTAC que además tienen signo en PRIDE, {cert['n_same_sign']} de {cert['n_sign_evaluable']} conservan la dirección. La concordancia es {comma(cert['sign_concordance'])}, con intervalo exacto de 95 % de {comma(cert['sign_ci_low'])} a {comma(cert['sign_ci_high'])}. H1 necesitaba una concordancia por debajo de 0,70. {comma(cert['sign_concordance'])} está por encima. La cláusula de proteína se rechaza.")
    p("La correlación de la magnitud, en las 1.263 proteínas con log2FC en las dos cohortes, es Pearson 0,584 y Spearman 0,521. El signo coincide más de lo que H1 predecía. La lista de proteínas que pasan el umbral en las dos cohortes sigue siendo corta: 28.")
    story.append(Image(str(FIG / "scatter_cptac_pride.png"), width=12.2 * cm, height=11.8 * cm))
    story.append(Paragraph("Figura 3. Cada punto es una proteína presente en las dos cohortes. La recta es la igualdad de efecto.", s["cap"]))
    p("Cien permutaciones del signo de los 94 pares, con semilla 20260914, dan una mediana de 0 proteínas significativas. El recuento observado, 688, tiene p empírico 0,0099. El contraste no es un azar de etiquetas. La comprobación con la columna de todos los péptidos, que no decide la hipótesis, marca 679 proteínas. El recuento no depende de haber elegido la columna de péptido no compartido.")
    p("De vías, 34 términos de proceso biológico salen en CPTAC y solo 4 se pueden comparar con PRIDE. Esos 4 coinciden de signo. El intervalo va de 0,40 a 1,0. Cuatro términos no certifican la cláusula de vía. Como la cláusula de proteína falla, H1 en conjunto queda rechazada.")

    h1("6. Péptidos, con el dato que hay")
    p(f"De las {audit['peptide_n_discovery_significant']} proteínas significativas en CPTAC, {audit['peptide_n_measurable']} tienen al menos un péptido único en el proteinGroups de PXD034575. La fracción es {comma(audit['peptide_measurable_fraction'])}. 458 no están en esa tabla. De las que sí están, casi todas tienen péptido único: la columna no filtra más que la presencia en el fichero. H4 se falsificaría si al menos el 90 % de los candidatos ya fuera medible en la segunda cohorte. 0,33 no llega a 0,90. El filtro de presencia no es redundante. La identidad de secuencia entre plataformas sigue sin estar en los ficheros.")

    r4 = audit["r4_stroma"]
    h1("7. Cuánto del resultado es tumor y cuánto es composición del tejido")
    p(f"El contraste estromal medio de los pares es {comma(r4['stroma_delta_mean'])}, con desviación {comma(r4['stroma_delta_sd'])}. El signo negativo dice que el tumor tiene menos matriz extracelular que su propio tejido adyacente, igual que indican CAVIN1 y PRELP. Esa diferencia no es ruido: está en casi todos los pares.")
    p(f"Al añadir ese contraste como covariable y evaluar el efecto en diferencia de composición cero, de las {r4['n_certified']} proteínas certificadas siguen pasando treat {r4['n_certified_surviving']}. La fracción es {comma(r4['fraction_surviving'])}. La magnitud mediana baja de 1,08 a {comma(r4['median_abs_log2fc_adjusted'])}. Y en {r4['n_certified_stroma_dependent_slope']} de {r4['n_certified']} proteínas, es decir el {comma(r4['fraction_stroma_dependent_slope'])}, el efecto del par se mueve de forma clara con el contraste estromal de ese par.")
    story.append(Image(str(FIG / "estroma.png"), width=12.4 * cm, height=11.3 * cm))
    story.append(Paragraph("Figura 5. Cada punto es una de las 688 proteínas certificadas. Los puntos naranjas dejan de pasar el umbral cuando el par se iguala en composición.", s["cap"]))
    p("Esta es la corrección más importante de la corrida. Más de la mitad de la señal tumor frente a adyacente en esta matriz es composición de tejido, no biología del epitelio tumoral. Un trabajo que no haga este ajuste publicará proteínas de estroma y de músculo liso como si fueran del tumor. El subíndice R4 pasa de no calculable a 0,446, y lo que mide es exactamente eso.")
    p("Hay que decir la limitación con claridad: el índice estromal se construye con proteínas de la misma matriz, así que el ajuste quita varianza compartida y puede quitar también señal real del epitelio. El número correcto de R4 exige el porcentaje de pureza del patólogo o una deconvolución con referencias externas. Lo que hay aquí es una cota, y la cota dice que la mitad de la lista es composición.")

    seq = axes["adenoma_carcinoma_sequence"]
    h1("8. La secuencia adenoma–carcinoma, con datos que nadie había usado")
    p(f"PXD034575 trae tres tejidos del mismo paciente: tumor, normal y pólipo. Los 16 pacientes tienen los tres. El pólipo nunca se había contrastado en este proyecto. Con la misma regla y emparejando dentro del paciente, el brazo tumor marca {axes['pride_paired_tumor']['n_significant']} proteínas y el brazo pólipo {axes['pride_paired_polyp']['n_significant']}.")
    p(f"La pregunta declarada es si el adenoma ya apunta en la misma dirección que el carcinoma. De las {seq['direction']['n']} proteínas significativas en tumor que tienen efecto en pólipo, {seq['direction']['same']} coinciden en el signo: {comma(seq['direction']['concordance'])}, con intervalo de {comma(seq['direction']['ci_low'])} a {comma(seq['direction']['ci_high'])}. De esas que coinciden, {seq['intermediate_given_same_direction']['same']} de {seq['intermediate_given_same_direction']['n']} tienen además menor magnitud en el pólipo que en el tumor: {comma(seq['intermediate_given_same_direction']['concordance'])}. La magnitud mediana es {comma(seq['median_abs_lfc_polyp'])} en pólipo y {comma(seq['median_abs_lfc_tumor'])} en tumor.")
    story.append(Image(str(FIG / "adenoma.png"), width=12.0 * cm, height=11.3 * cm))
    story.append(Paragraph("Figura 6. El pólipo del mismo paciente ya se mueve hacia donde está el tumor, y se queda a medio camino.", s["cap"]))
    p("Esto es un gradiente ordenado: normal, adenoma, carcinoma. Es el resultado más limpio de la corrida porque el emparejamiento es dentro del mismo paciente, lo que elimina de golpe la variación entre personas. Tiene un límite serio: son 16 pacientes de un solo estudio, y el pólipo es del mismo paciente que ya tiene cáncer, no de una persona con solo un adenoma. No dice que un adenoma cualquiera vaya a progresar.")

    h1("9. Del tejido a la sangre: la brecha que nadie debería ignorar")
    gap = axes["tissue_to_plasma_gap"]
    p("Esta es la pregunta más útil para quien piense en un marcador en sangre. Si una proteína sube en el tumor, ¿sube en el suero? La literatura reciente dice que no se puede dar por hecho: la sobreexpresión en la superficie del tumor casi nunca se refleja en la concentración circulante, y CEACAM5 en particular muestra una disociación entre tejido y plasma.")
    p(f"Con los datos de esta corrida, de las {gap['n_tissue_significant_with_uniprot']} proteínas certificadas en tejido con identificador UniProt, solo {gap['n_measured_in_plasma']} aparecen medidas en el suero DIA: una fracción de {comma(gap['fraction_tissue_hits_seen_in_plasma'])}. De esas, {gap['same']} de {gap['n']} coinciden en el signo, o sea {comma(gap['concordance'])}, con intervalo de {comma(gap['ci_low'])} a {comma(gap['ci_high'])} y p de {comma(gap['p_vs_half'])} frente al azar. Es una moneda al aire. Y {gap['n_also_significant_in_plasma']} proteínas del tejido son además significativas en el suero: ninguna.")
    story.append(Image(str(FIG / "tejido_plasma.png"), width=12.0 * cm, height=10.5 * cm))
    story.append(Paragraph("Figura 7. Las 86 proteínas del tejido que el suero también mide. Si el tejido predijera la sangre, los puntos estarían en dos cuadrantes.", s["cap"]))
    p("CEACAM5, la proteína del antígeno CEA, lo ilustra bien. En el tejido de CPTAC sube en tumor con log2FC 0,88 y q de 0,009. En la matriz de suero medida sin enriquecimiento no aparece: la proteómica DIA no dirigida no la ve. Un marcador que la clínica mide a diario con un ensayo dirigido no es visible en un barrido no dirigido del mismo líquido. Eso no contradice la clínica: dice que cada plataforma tiene su ventana, y que una lista de tejido no se convierte en una prueba de sangre.")

    h2("Panel externo: lo que otros ya publicaron, probado en esta matriz")
    p("Para no pescar en los datos, se tomó una lista de proteínas con dirección ya publicada en plasma o suero de CCR y se comprobó esa dirección en la cohorte de esta corrida. La lista se fijó antes de mirar la tabla y no se ajustó después.")
    pan = axes["external_plasma_panel"]
    header2 = [Paragraph(esc(x), s["smallw"]) for x in ("Proteína", "Dirección esperada", "log2FC en suero", "¿Coincide?", "log2FC en tejido")]
    data2 = [header2]
    for _, row in panel.iterrows():
        measured = bool(row["plasma_measured"])
        data2.append([
            Paragraph(esc(row["gene_symbol"]), s["small"]),
            Paragraph("sube" if int(row["expected_direction_plasma"]) > 0 else "baja", s["small"]),
            Paragraph(comma(float(row["plasma_log2fc"])) if measured and pd.notna(row["plasma_log2fc"]) else "no medida", s["small"]),
            Paragraph("sí" if row["direction_matches"] is True or str(row["direction_matches"]) == "True" else ("no" if measured else "—"), s["small"]),
            Paragraph(comma(float(row["tissue_log2fc"])) if pd.notna(row.get("tissue_log2fc")) else "—", s["small"]),
        ])
    t2 = Table(data2, colWidths=[3.0 * cm, 3.4 * cm, 3.4 * cm, 2.4 * cm, 3.6 * cm])
    t2.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1f3a5f")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Arial-Bold"),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.2, colors.HexColor("#cccccc")),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(t2)
    story.append(Spacer(1, 0.3 * cm))
    p(f"De las {pan['n_panel']} proteínas del panel, el suero de esta cohorte mide {pan['n_measured_in_plasma']}. De esas, {pan['same']} de {pan['n']} van en la dirección publicada: {comma(pan['concordance'])}, con intervalo de {comma(pan['ci_low'])} a {comma(pan['ci_high'])} y p de {comma(pan['p_vs_half'])}. La dirección es compatible con lo publicado, y con siete proteínas el intervalo es demasiado ancho para llamarlo una réplica. Ninguna pasa el umbral de efecto en esta cohorte de 20 contra 20.")
    p("La tabla contiene el hallazgo más instructivo del documento. LRG1, SERPINA1 y ORM1 están claramente disminuidas en el tejido tumoral frente al adyacente, con log2FC de −1,14, −1,16 y −0,97, y las tres suben en el suero: +1,23, +0,55 y +0,76. Son proteínas de fase aguda que el hígado produce en respuesta a la inflamación sistémica. La sangre no está leyendo el tumor: está leyendo la respuesta del organismo. Por eso el tejido y el plasma pueden apuntar en direcciones opuestas sin que ninguno de los dos esté mal medido.")

    h1("10. Las otras variables que sí se pudieron contrastar")
    p("El contraste primario no se volvió a calcular. Lo que sigue usa la misma regla, sobre la misma matriz sin editar, y se guarda en una carpeta distinta.")
    sex = secondary["sex_sign_agreement_female_reference"]
    female = secondary["contrasts"]["female"]
    male = secondary["contrasts"]["male"]
    stage = secondary["contrasts"]["stage"]
    plasma = secondary["contrasts"]["plasma"]
    p(f"Sexo. Hay {female['n_pairs']} pares de mujer y {male['n_pairs']} de hombre. treat marca {female['n_significant']} proteínas en mujeres y {male['n_significant']} en hombres. De las significativas en mujeres, {sex['same']} de {sex['n']} tienen el mismo signo en hombres: concordancia {comma(sex['concordance'])} (intervalo {comma(sex['ci_low'])} a {comma(sex['ci_high'])}). La dirección tumor frente a adyacente no cambia con el sexo dentro de CPTAC. Esto es una réplica interna del mismo estudio, no una cohorte externa.")
    p(f"Estadio. Entre tumores, {stage['n_advanced']} son avanzados y {stage['n_early']} tempranos. Con el mismo umbral de 0,585, ninguna proteína pasa el q-valor. El estadio, en esta matriz y con esta regla, no produce una lista.")
    p(f"Plasma. Veinte muestras de cáncer y veinte de personas sanas, log2 solo donde la intensidad es positiva, sin imputar. treat marca {plasma['n_significant']} proteínas de {plasma['n_tested']} contrastadas: SAA1 (P0DJI8), LBP (P18428) y CHMP4A (Q9BY43). SAA1 y LBP son proteínas de fase aguda, esperables en suero, y no se trasladan al tejido. CHMP4A se apoya en menos muestras que las otras dos. Tres proteínas no son un panel diagnóstico. El plasma no se mezcla con H1.")
    story.append(Image(str(FIG / "conteos.png"), width=15.2 * cm, height=8.4 * cm))
    story.append(Paragraph("Figura 4. Proteínas que pasan treat. Estadio en cero y plasma en tres, con la misma regla que CPTAC.", s["cap"]))

    h2("Inventario de variables")
    header = [Paragraph(esc(x), s["smallw"]) for x in ("Variable", "¿Se contrastó?", "Por qué")]
    data = [header]
    for _, row in inventory.iterrows():
        flag = str(row["aplica"]).strip().lower() == "true"
        data.append([
            Paragraph(esc(row["variable"]), s["small"]),
            Paragraph("sí" if flag else "no", s["small"]),
            Paragraph(esc(row["detalle"]), s["small"]),
        ])
    table = Table(data, colWidths=[4.2 * cm, 2.4 * cm, 10.2 * cm])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1f3a5f")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Arial-Bold"),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.2, colors.HexColor("#cccccc")),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(table)
    story.append(Spacer(1, 0.3 * cm))

    h1("11. Modelos de aprendizaje y cómo se elige uno")
    ms = axes["model_selection"]
    p("No se elige un modelo por el AUC de dentro. Dentro de CPTAC, tumor y tejido adyacente se separan tanto que un clasificador llega a un AUC de 1. Ese número no distingue un buen modelo de uno malo: solo dice que las dos clases ya estaban separadas en la misma matriz.")
    p(f"La regla se escribió antes de mirar la tabla. La métrica primaria es el AUC en la cohorte que el modelo no vio. Cada modelo lleva un intervalo del 95 % por bootstrap de la muestra de prueba. Solo hay ganador si el intervalo del primero no se solapa con el del segundo. El entrenamiento usa {ms['n_train']} especímenes de CPTAC y la prueba {ms['n_test']} de PXD034575, sobre {ms['n_features_common_uniprot']} proteínas con el mismo UniProt en las dos.")
    story.append(Image(str(FIG / "modelos.png"), width=15.0 * cm, height=6.7 * cm))
    story.append(Paragraph("Figura 8. AUC externo con intervalo. El punto es la estimación; la barra, el bootstrap.", s["cap"]))
    p(f"El resultado es {ms['verdict']}. Elastic net da {comma(ms['models']['elastic_net']['auc_external'])} con intervalo de {comma(ms['models']['elastic_net']['ci_low'])} a {comma(ms['models']['elastic_net']['ci_high'])}; el bosque aleatorio {comma(ms['models']['random_forest']['auc_external'])} de {comma(ms['models']['random_forest']['ci_low'])} a {comma(ms['models']['random_forest']['ci_high'])}; la logística L1 {comma(ms['models']['logistica_l1']['auc_external'])} de {comma(ms['models']['logistica_l1']['ci_low'])} a {comma(ms['models']['logistica_l1']['ci_high'])}. Con 32 muestras de prueba los tres intervalos se pisan. Declarar ganador al mejor número sería leer ruido.")
    p("Las dos matrices están en escalas distintas, log-ratio frente a log2 de intensidad, así que cada cohorte se estandariza por proteína dentro de sí misma. Eso convierte la transferencia en una prueba de patrón relativo y no de calibración absoluta, y así hay que escribirlo. Una red o un modelo profundo no se ajusta aquí: hay del orden de cien muestras y miles de proteínas, y sin un nulo que conserve el grado de cada nodo un grafo que parece mejor puede estar premiando a los nodos más conectados. El producto de la tesis sigue siendo la tabla de hipótesis, no un clasificador.")

    h1("12. El certificado, subíndice por subíndice")
    story.append(Image(str(FIG / "subindices.png"), width=15.2 * cm, height=7.2 * cm))
    story.append(Paragraph("Figura 9. Los cinco subíndices no se promedian en un número único. R3 sigue sin dato.", s["cap"]))
    p(f"R1, dirección: {comma(cert['sign_concordance'])}. R2, se mide en la segunda cohorte: {comma(cert['measurable_fraction'])}. R4, sobrevive al igualar composición: {comma(cert['subindices'][3]['value'])}. R5, p empírico frente a 100 permutaciones: {comma(audit['r5_empirical_p'], 4)}. R3, homogeneidad entre estudios, necesita un meta-análisis de efectos aleatorios que esta corrida no hace. Colapsar los cinco en un índice único esconde justo lo que importa: que la dirección es firme, que la medibilidad es baja y que la mitad del efecto es composición.")

    h1("13. Conclusión")
    p("Con los datos públicos disponibles y sin modificarlos, hay cuatro cosas que se pueden escribir y defender.")
    p("Primera. En CPTAC, 688 proteínas cambian más de un 50 % entre tumor y tejido adyacente con un efecto mínimo declarado, el patrón no es un azar de etiquetas (p empírico 0,0099 sobre 100 permutaciones) y la dirección se repite entre mujeres y hombres del mismo estudio con una concordancia de 0,998. La medición es sólida.")
    p("Segunda. H1 queda rechazada, y se rechaza por el lado contrario al que el diseño esperaba. La hipótesis predecía que las proteínas coincidirían poco entre cohortes, menos del 70 %. Coinciden en el 87,4 %, con intervalo de 80,1 a 92,8. La dirección del efecto es más reproducible de lo que el proyecto suponía. Lo que no se reproduce es la lista: solo 28 proteínas pasan el umbral en las dos cohortes, porque la segunda tiene 16 pares frente a 94 y menos potencia.")
    p("Tercera, y es la que cambia cómo hay que leer todo lo demás. Al igualar cada par en contenido estromal, solo el 44,6 % de las proteínas certificadas sigue pasando el umbral, y el 84,3 % tiene un efecto que se mueve con la composición del par. El tumor tiene menos matriz extracelular que su propio tejido adyacente. Más de la mitad de lo que un análisis bulk llamaría biología del tumor es mezcla de tejidos.")
    p("Cuarta. El tejido no predice la sangre. Solo el 12,5 % de las proteínas certificadas en tejido aparece medida en el suero, la dirección coincide en el 48,8 % —azar puro, p 0,91— y ninguna es significativa en los dos líquidos. Proteínas de fase aguda como LRG1, SERPINA1 y ORM1 bajan en el tejido tumoral y suben en el suero: la sangre lee la respuesta del organismo, no el tumor. CEACAM5 sube en tejido y la proteómica no dirigida de suero no la ve.")
    p("A eso se suma un resultado nuevo con datos que estaban sin usar: en los 16 pacientes de PXD034575 con pólipo, normal y tumor del mismo individuo, el adenoma ya apunta en la dirección del carcinoma en el 94,3 % de las proteínas y con menor magnitud en el 91,9 % de ellas. Es un gradiente ordenado, medido dentro del mismo paciente.")
    p("Lo que sigue sin certificarse hay que decirlo igual. No hay meta-análisis de efectos aleatorios, así que R3 no existe. No hay pureza del patólogo: R4 es una cota construida con el propio proteoma. No hay intensidades de canal TMT, así que no se puede ajustar el modelo pensado para esa tecnología. La concordancia de vía se apoya en cuatro términos. PXD034575 es el control de independencia del diseño, no la réplica primaria PXD019504, y la auditoría de ficheros duplicados no se repitió. No queda certificado ningún marcador, ninguna diana y ninguna recomendación clínica. Nada de este documento sirve para decidir sobre un paciente.")

    h1("14. Literatura reciente que respalda estas decisiones")
    refs = [
        "Composite proteomic and metabolomic plasma biomarkers for detection of colorectal, lung and ovarian cancers. Molecular Cancer, 2026. doi:10.1186/s12943-026-02654-1. Da la dirección publicada de CEACAM5, ferritina, SPP1 y PECAM1 en plasma, con AUC de 0,66 a 0,69. Es la fuente del panel externo.",
        "Translational Gap in Biomarker Discovery: Tumor Surface Markers Rarely Mirror Circulating Levels. Muestra que CEACAM5 se disocia entre tejido y plasma y que CEACAM1 fue el único marcador concordante. Es la base de la sección 9.",
        "Proteomics of colorectal tumors identifies the role of CAVIN1 in tumor relapse. Molecular Systems Biology, 2025. doi:10.1038/s44320-025-00102-8. Los subtipos proteómicos de CCR se separan por contenido estromal e inmune. Justifica el ajuste por composición.",
        "ProteoMixture: a cell type deconvolution tool for bulk tissue proteomic data. iScience, 2024. Deconvolución de proteoma bulk con firmas de estroma e inmune. Es el método de referencia que aquí se aproxima con el índice de matriz extracelular.",
        "Evaluating and Optimizing Mass Spectrometry Proteomics Data to Deconvolve Cell-Type-Specific Protein Expression in Tumors. PMC12772123. Deconvolución de CCR en epitelio, estroma e inmune, y la advertencia sobre la transformación de datos TMT.",
        "Mass Spectrometry Proteomics Characterization of Plasma Biomarkers for Colorectal Cancer Associated With Inflammation. PMC11191626. Complemento C5 elevado en plasma de CCR y confirmado por ELISA; C8A y C4B separan estadio temprano de avanzado.",
        "Label-free quantitative proteomics reveals aberrant expression levels of LRG, C9, FN, A1AT and AGP1 in the plasma of patients with colorectal cancer. Clinical Proteomics, 2023. doi:10.1186/s12014-023-09407-y. Direcciones de LRG1, C9, SERPINA1, ORM1 y FN1 en plasma.",
        "A plasma metabolite-based test to detect minimal residual disease in post-surgery patients with colorectal cancer. 160 pacientes, sensibilidad 62 % y especificidad 80 %, con adelanto mediano de 471 días. Hoy la evidencia más fuerte en enfermedad residual es metabolómica, no proteómica.",
        "Insights into ctDNA assessment to detect minimal residual disease in localized colorectal cancer. Scientific Reports, 2026. doi:10.1038/s41598-026-58101-0. 47 pacientes de estadio II y III de alto riesgo: sensibilidad 37,5 % solo con sangre y 50 % con criterio informado por tejido; las recaídas pulmonares y peritoneales no se detectaron.",
        "Unexplained increases in serum carcinoembryonic antigen levels in colorectal cancer patients during the postoperative follow-up period. PMC10898814. El 20 % de los pacientes sin recaída tuvo al menos una elevación de CEA, el 98 % de ellas entre 5 y 15 ng/ml.",
        "The diagnostic accuracy of a single CEA blood test in detecting colorectal cancer recurrence: FACS trial. PLOS One. Con umbral de 5 µg/l, el CEA aislado pierde cerca de la mitad de las recaídas y obliga a combinarlo con imagen.",
    ]
    story.append(ListFlowable([bullet(t, s["small"]) for t in refs], bulletType="bullet", start="•"))
    p("Estas referencias no se usaron para ajustar ningún umbral. Se usaron para decidir qué preguntar y para fijar la dirección esperada del panel de plasma antes de mirar la matriz. Las dos últimas describen el comportamiento conocido del CEA en seguimiento y están citadas por contexto del marco teórico, no como resultado de esta corrida.")

    h1("15. Dónde está cada archivo")
    p("Matrices y tablas: /Volumes/JANCODATA/OncoGO-CRC. Contraste primario y ajuste estromal: results/certified/. Certificado R: results/certification/informe_brecha.md. Sexo, estadio y plasma: results/secondary/. Pólipos, brecha tejido–plasma, panel externo y modelos: results/axes/. Figuras: results/secondary/figures/. Este PDF: OncoGO-CRC_tesis_explicada.pdf, en la raíz de esa carpeta.")

    def footer(canvas, doc) -> None:
        canvas.saveState()
        canvas.setFont("Arial", 8)
        canvas.setFillColor(colors.HexColor("#666666"))
        canvas.drawString(2 * cm, 1.2 * cm, "OncoGO-CRC · corrida certificada · 1 de octubre de 2026")
        canvas.drawRightString(A4[0] - 2 * cm, 1.2 * cm, str(doc.page))
        canvas.restoreState()

    doc = SimpleDocTemplate(
        str(PDF_PATH),
        pagesize=A4,
        leftMargin=1.8 * cm,
        rightMargin=1.8 * cm,
        topMargin=1.6 * cm,
        bottomMargin=1.8 * cm,
        title="OncoGO-CRC. Qué se busca, cómo se analiza y qué se puede afirmar",
        author="OncoGO-CRC",
    )
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(PDF_PATH, PDF_PATH.stat().st_size)


def main() -> int:
    ap = argparse.ArgumentParser(description="Genera el PDF explicativo de la corrida certificada.")
    ap.add_argument(
        "--data-root",
        type=Path,
        default=None,
        help="Raíz de los resultados. Por omisión: la caché del repo si está completa, si no el disco de datos.",
    )
    ap.add_argument("--out", type=Path, default=PDF_DEFAULT, help="Ruta del PDF.")
    ap.add_argument("--figures", type=Path, default=REPO / "docs" / "tesis_figuras", help="Carpeta de figuras.")
    ap.add_argument(
        "--sync-cache",
        action="store_true",
        help="Copia los insumos del disco de datos a docs/tesis_datos antes de generar.",
    )
    args = ap.parse_args()

    if args.sync_cache:
        if missing_inputs(WORKSPACE_DEFAULT):
            print(f"No se puede sincronizar: faltan insumos en {WORKSPACE_DEFAULT}.")
            return 2
        print(f"Insumos copiados a la caché del repo: {sync_cache(WORKSPACE_DEFAULT, CACHE_DEFAULT)}")

    root = args.data_root
    if root is None:
        root = CACHE_DEFAULT if not missing_inputs(CACHE_DEFAULT) else WORKSPACE_DEFAULT

    faltan = missing_inputs(root)
    if faltan:
        print(f"No se genera el PDF. Faltan {len(faltan)} insumos en {root}:")
        for rel in faltan:
            print(f"  - {rel}")
        print(
            "\nEstos archivos los produce el pipeline (etapas certified_inference, secondary_contrasts,\n"
            "translational_axes, certification). Con el disco de datos montado:\n"
            f"  .venv/bin/python scripts/build_thesis_pdf.py --sync-cache\n"
            "Eso deja la copia en docs/tesis_datos y desde entonces el PDF se rehace sin el disco."
        )
        return 1

    set_paths(root, args.out, args.figures)
    build()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
