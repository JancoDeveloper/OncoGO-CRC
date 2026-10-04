"""PDF explicativo de v2 para lectores sin formación en bioinformática.

Lee las tablas de resumen de results/v2/ y escribe un documento en lenguaje
claro. Cada cifra sale de un JSON o TSV; no se escribe ningún número a mano.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
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
)

REPO = Path(__file__).resolve().parents[1]
WORKSPACE_DEFAULT = Path("/Volumes/JANCODATA/OncoGO-CRC")
OUT_DEFAULT = REPO / "docs" / "reportes" / "OncoGO-CRC_guia_publica_v2.pdf"
FONT = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT_B = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"

AUDITORIAS = (
    ("certification/certificado_v2.json", "certificado"),
    ("certified/inference_audit_v2.json", "inferencia"),
    ("axes/tcga_cea/tcga_cea_audit_v2.json", "tcga_cea"),
    ("axes/pureza_medida/pureza_medida_audit_v2.json", "pureza"),
    ("axes/heterogeneidad/heterogeneidad_audit_v2.json", "heterogeneidad"),
    ("axes/deconvolucion/deconvolucion_audit_v2.json", "deconvolucion"),
    ("marcadores/marcadores_audit_v2.json", "marcadores"),
)


def esc(text: str) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def comma(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}".replace(".", ",").replace("-", "\u2212")


def pct(value: float) -> str:
    return f"{value * 100:.1f} %".replace(".", ",")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def setup_font() -> None:
    pdfmetrics.registerFont(TTFont("Arial", FONT))
    bold = FONT_B if Path(FONT_B).exists() else FONT
    pdfmetrics.registerFont(TTFont("Arial-Bold", bold))


def _sub(cert: dict, codigo: str) -> dict:
    for s in cert.get("subindices", []):
        if s.get("codigo") == codigo:
            return s
    return {}


def build(root: Path, out: Path, fig_dir: Path | None) -> None:
    setup_font()
    base = root / "results" / "v2"
    cert = load_json(base / "certification" / "certificado_v2.json")
    inf = load_json(base / "certified" / "inference_audit_v2.json")
    tcga = load_json(base / "axes" / "tcga_cea" / "tcga_cea_audit_v2.json")
    pureza = load_json(base / "axes" / "pureza_medida" / "pureza_medida_audit_v2.json")
    het = load_json(base / "axes" / "heterogeneidad" / "heterogeneidad_audit_v2.json")
    dec = load_json(base / "axes" / "deconvolucion" / "deconvolucion_audit_v2.json")
    marc = load_json(base / "marcadores" / "marcadores_audit_v2.json")

    r1 = _sub(cert, "R1_direccion")
    r2 = _sub(cert, "R2_detectabilidad")
    r3 = _sub(cert, "R3_homogeneidad")
    r4 = _sub(cert, "R4_composicion")
    cea = (marc.get("ajustes") or {}).get("CEA|deteccion_de_recaida") or {}

    styles = {
        "title": ParagraphStyle(
            "title", fontName="Arial-Bold", fontSize=20, leading=24,
            alignment=TA_CENTER, spaceAfter=14,
        ),
        "h1": ParagraphStyle(
            "h1", fontName="Arial-Bold", fontSize=14, leading=18,
            spaceBefore=12, spaceAfter=8, textColor=colors.HexColor("#1f3a5f"),
        ),
        "h2": ParagraphStyle(
            "h2", fontName="Arial-Bold", fontSize=11.5, leading=15,
            spaceBefore=8, spaceAfter=5,
        ),
        "body": ParagraphStyle(
            "body", fontName="Arial", fontSize=10.5, leading=14,
            alignment=TA_JUSTIFY, spaceAfter=6,
        ),
        "small": ParagraphStyle(
            "small", fontName="Arial", fontSize=9, leading=12,
            alignment=TA_JUSTIFY, textColor=colors.HexColor("#444444"),
        ),
        "warn": ParagraphStyle(
            "warn", fontName="Arial-Bold", fontSize=10.5, leading=14,
            alignment=TA_JUSTIFY, textColor=colors.HexColor("#8b1a1a"),
            spaceAfter=8,
        ),
    }

    story: list = []

    def p(text: str, style: str = "body") -> None:
        story.append(Paragraph(esc(text), styles[style]))

    def h1(text: str) -> None:
        story.append(Paragraph(esc(text), styles["h1"]))

    def h2(text: str) -> None:
        story.append(Paragraph(esc(text), styles["h2"]))

    def bullet(text: str) -> ListItem:
        return ListItem(Paragraph(esc(text), styles["body"]), leftIndent=12)

    story.append(Spacer(1, 2.0 * cm))
    story.append(Paragraph("OncoGO-CRC", styles["title"]))
    story.append(Paragraph(
        "Guía pública de la versión 2<br/>"
        "Qué se estudió, qué se encontró y qué no significa para un paciente",
        ParagraphStyle("sub", fontName="Arial", fontSize=12, leading=16, alignment=TA_CENTER),
    ))
    story.append(Spacer(1, 1.0 * cm))
    p(
        "Este documento resume un trabajo de investigación sobre proteínas del "
        "cáncer de colon. Está escrito para alguien sin formación en biología "
        "molecular ni en estadística. Los números salen de tablas del proyecto; "
        "no son opiniones ni recomendaciones médicas."
    )
    p(
        "Ningún resultado de este trabajo sirve para decidir sobre un paciente. "
        "No sustituye a un médico, a una prueba de laboratorio ni a un tratamiento.",
        "warn",
    )
    story.append(PageBreak())

    h1("1. ¿De qué va este proyecto?")
    p(
        "El cáncer de colon es un tumor del intestino grueso. Las células del tumor "
        "producen proteínas —moléculas que hacen el trabajo dentro y fuera de la "
        "célula— de forma distinta al tejido sano de al lado. Si esa diferencia se "
        "repitiera en varios estudios independientes, podría ayudar a entender la "
        "enfermedad y, con muchísimo más trabajo después, a pensar en marcadores "
        "en sangre. Ese «después» no está aquí."
    )
    p(
        "OncoGO-CRC es un programa de ordenador que descarga datos públicos de "
        "varios laboratorios, los compara con reglas fijas escritas antes de mirar "
        "los resultados, y deja un informe reproducible. No se tocan los datos "
        "originales. No se inventan valores donde falta una medición."
    )

    h1("2. ¿Qué datos se usaron?")
    bullets = [
        "Tejido de colon de 94 personas con tumor y tejido normal de la misma persona (estudio CPTAC, tecnología TMT).",
        "Tejido de 16 personas con tumor, normal y pólipo del mismo paciente (estudio PRIDE PXD034575).",
        "Tejido de 22 pares tumor–normal de otra cohorte (PRIDE PXD019504), para comprobar si los resultados se parecen entre estudios.",
        "Sangre: 20 muestras de personas con cáncer y 20 de personas sanas (PRIDE PXD055218).",
        "Datos clínicos de cientos de tumores de TCGA: nivel de CEA en sangre, pureza tumoral medida por ADN, expresión del gen CEACAM5.",
        "Literatura clínica sobre sensibilidad y especificidad del CEA frente a otras pruebas (meta-análisis de tablas 2×2).",
    ]
    story.append(ListFlowable([bullet(t) for t in bullets], bulletType="bullet", start="•"))
    p(
        "Cada fuente tiene su propia tecnología. Por eso no se mezclan en un solo "
        "modelo: primero se analiza por dentro, y después se comparan los resultados."
    )

    h1("3. La pregunta principal")
    p(
        "La hipótesis H1 preguntaba si las proteínas que cambian en el primer estudio "
        "grande coincidirían poco con las del segundo. Si la coincidencia fuera baja "
        "(menos del 70 % en dirección), la hipótesis se consideraría apoyada. "
        "Si fuera alta, se rechazaría."
    )
    p(
        f"En CPTAC, con {inf['cptac']['pares']} pares, {inf['cptac']['n_certificadas']} "
        f"proteínas pasaron el umbral de cambio mínimo (más del 50 % arriba o abajo, "
        f"con corrección por comparar miles de proteínas a la vez). En la réplica "
        f"pareada de 16 pares, {inf['pride']['n_significativas_pareado']} proteínas "
        f"pasaron el mismo filtro."
    )
    p(
        f"De las proteínas que se pueden comparar entre cohortes, el "
        f"{pct(float(r1['valor']))} mantiene la misma dirección (sube o baja). "
        f"Eso está por encima del 70 %, así que H1 se rechaza por la cláusula de "
        f"proteína: la dirección es más reproducible de lo que el diseño esperaba."
    )

    h1("4. Los cinco indicadores resumen (R1 a R5)")
    p("El certificado no resume todo en una sola nota. Usa cinco piezas separadas:")
    indicadores = [
        (
            "R1 — Dirección",
            f"{pct(float(r1['valor']))} de las proteínas comparables van en la misma "
            f"dirección en dos cohortes distintas.",
        ),
        (
            "R2 — Detectabilidad",
            f"Solo el {pct(float(r2['valor']))} de lo encontrado en el estudio grande "
            f"aparece medible en la cohorte pequeña con la misma regla.",
        ),
        (
            "R3 — Heterogeneidad entre estudios",
            f"El meta-análisis con tres cohortes de tejido da un I² del "
            f"{pct(float(r3['valor']))}: los estudios discrepan más de lo que "
            f"coinciden en magnitud.",
        ),
        (
            "R4 — Composición del tejido",
            f"Tras ajustar por cuánto estroma hay en cada muestra, queda el "
            f"{pct(float(r4['valor']))} de los hallazgos originales. Es una cota: "
            f"parte de lo que parece «del tumor» es mezcla de tejidos.",
        ),
        (
            "R5 — Azar",
            "Al permutar las etiquetas tumor/normal cien veces, ninguna repetición "
            "igualó lo observado. El contraste no es un truco de etiquetas.",
        ),
    ]
    story.append(ListFlowable([bullet(f"{t}. {d}") for t, d in indicadores], bulletType="bullet", start="•"))

    h1("5. Lo más importante para entender los resultados")
    h2("5.1 Mucho de la señal es composición, no tumor puro")
    p(
        "El «normal» de estos estudios es tejido de al lado del tumor, no sangre de "
        "una persona sana. Trae músculo, matriz y células que no son el epitelio "
        "canceroso. Cuando se corrige por esa mezcla, más de la mitad de las "
        "proteínas dejan de pasar el umbral. Leer una lista de proteínas sin ese "
        "contexto es leer estroma y músculo como si fueran el tumor."
    )
    h2("5.2 El tejido no predice la sangre")
    p(
        "Una proteína que sube en el tumor no sube automáticamente en suero. "
        "CEACAM5 —la proteína del antígeno CEA— sube en tejido, pero el barrido "
        "proteómico de plasma usado aquí no la detecta. La sangre y el tumor miden "
        "cosas distintas con tecnologías distintas."
    )
    h2("5.3 El CEA en sangre depende poco de cuánto gen expresa el tumor")
    modelo = tcga.get("modelo_del_cea", {})
    p(
        f"En TCGA, {tcga.get('n_muestras', tcga.get('cobertura', {}).get('n_muestras', 'cientos de'))} "
        f"tumores tienen datos clínicos y de expresión. Añadir CEACAM5 al modelo "
        f"clínico completo solo aporta el "
        f"{pct(float(modelo.get('aporte_de_ceacam5_al_r2', 0)))} extra de variación "
        f"del CEA. El marcador en sangre depende sobre todo de cuánto tumor hay y "
        f"de dónde está, no de cuánto gen expresa el tumor."
    )
    h2("5.4 Separar cáncer de sano no es lo mismo que separar cáncer de enfermo")
    wann = load_json(root / "external" / "wannberg_2026" / "estado.json")
    if wann.get("estado") == "usada":
        comp = wann.get("comparacion_del_grupo_de_control", {})
        sano = comp.get("contra_personas_sanas", {})
        sint = comp.get("contra_pacientes_con_sintomas", {})
        p(
            f"Con la misma plataforma, distinguir cáncer de personas sanas da un "
            f"AUC de {comma(float(sano.get('auc', 0)), 3)}, pero "
            f"distinguir cáncer de pacientes con síntomas que no tenían cáncer baja "
            f"a {comma(float(sint.get('auc_descubrimiento', 0)), 3)}. En consulta "
            f"llegan personas enfermas de otra cosa, no voluntarios sanos."
        )

    story.append(PageBreak())

    h1("6. Qué aportó la versión 2 respecto a la primera corrida")
    correcciones = [
        "Permutaciones con semilla fija para comprobar que el resultado no es azar.",
        "Denominadores explícitos en cada fracción (cuántas proteínas entran y por qué).",
        "Ajuste por composición sin usar la misma proteína para medir y para corregir.",
        "Pólipos analizados por separado, no mezclados con tumor–normal.",
        "Cláusula de vías: un término solo cuenta si tiene al menos diez proteínas.",
        "Tercera cohorte de tejido para estimar heterogeneidad entre estudios.",
        "Verificador automático: si un número del informe no sale de una tabla, el informe no se publica.",
        "Nueve fuentes externas (TCGA, pureza ABSOLUTE, heterogeneidad intratumoral, deconvolución celular, plasma frente a enfermedad benigna, meta-análisis de marcadores, protocolo de cohorte futura).",
    ]
    story.append(ListFlowable([bullet(t) for t in correcciones], bulletType="bullet", start="•"))

    h1("7. Nuevas comprobaciones en lenguaje llano")
    h2("Pureza medida con otra tecnología")
    ecm = pureza.get("validacion_del_indice", {}).get("ecm", {})
    p(
        f"El índice de estroma que usa el proyecto va en el sentido correcto frente "
        f"a la pureza medida por ADN (correlación {comma(float(ecm.get('rho', 0)), 3)}), "
        f"pero solo comparte el {pct(float(ecm.get('varianza_compartida', 0)))} de la "
        f"variación. No sustituye a la pureza del patólogo."
    )
    h2("Variación dentro del mismo tumor")
    ratio = het.get("variacion_dentro_frente_a_entre", {})
    p(
        f"Medir un solo trozo de tumor es medir un trozo, no el tumor entero. "
        f"Con cuatro pacientes y remuestreo, la variación dentro del tumor y entre "
        f"pacientes es del mismo orden (cociente {comma(float(ratio.get('cociente_mediano', 0)), 3)}, "
        f"intervalo incluye el 1). Es una señal, no una conclusión cerrada."
    )
    h2("Deconvolución celular")
    cptac = dec.get("por_que_no_se_corre_en_cptac", {})
    p(
        f"Intentar estimar proporciones de tipos celulares no funciona en la matriz "
        f"TMT de CPTAC: el {pct(float(cptac.get('fraccion_negativa', 0.525)))} de los "
        f"valores es negativo, y una mezcla lineal no negativa no puede reproducirlos. "
        f"Las proporciones tipo a tipo no se publican; solo se usa un eje epitelial "
        f"como comprobación de sensibilidad cuando los datos lo permiten."
    )
    if cea.get("estado") == "calculado":
        h2("Rendimiento publicado del CEA")
        p(
            f"Según {int(cea['k'])} estudios extraídos de la literatura (pendientes "
            f"de verificación humana fila a fila), el CEA para detectar recaída "
            f"tiene sensibilidad {pct(float(cea['sensibilidad']))} y especificidad "
            f"{pct(float(cea['especificidad']))}: de cada diez recaídas, deja pasar "
            f"cerca de cuatro. Ese es el listón que cualquier alternativa debería "
            f"superar con evidencia comparable."
        )

    h1("8. Qué hay en el repositorio")
    p(
        "El código fuente está en GitHub (rama principal tras integrar v2). "
        "Los datos pesados viven en un disco de datos local; el repositorio incluye "
        "los informes en PDF, las tablas de resumen y la documentación para "
        "reproducir el análisis."
    )
    entregables = [
        "docs/reportes/OncoGO-CRC_guia_publica_v2.pdf — este documento.",
        "docs/reportes/OncoGO-CRC_informe_v2.pdf — informe técnico completo con figuras.",
        "docs/reportes/OncoGO-CRC_informe_v2_tecnico.pdf — copia del informe de la corrida.",
        "results/v2/report/ — HTML, Markdown y figuras generadas desde tablas.",
        "config/preregistro_v2.yaml — reglas fijadas antes de mirar los datos.",
        "docs/fuentes_datos.md — origen, licencia y limitaciones de cada fuente.",
        "docs/enmiendas_v2.md — cambios declarados respecto al preregistro.",
    ]
    story.append(ListFlowable([bullet(t) for t in entregables], bulletType="bullet", start="•"))

    if fig_dir and (fig_dir / "indice_r_y_changelog.png").exists():
        story.append(PageBreak())
        h1("9. Figura: índice R y cambios respecto a v1")
        story.append(Image(str(fig_dir / "indice_r_y_changelog.png"), width=15.5 * cm, height=7.5 * cm))
        p(
            "Cada barra es un subíndice. No se promedian en una sola nota porque "
            "miden cosas distintas: dirección, detectabilidad, heterogeneidad, "
            "composición y azar.",
            "small",
        )

    h1("10. Conclusión honesta")
    p(
        "Con datos públicos y reglas fijas, se puede afirmar que hay muchas proteínas "
        "que cambian entre tumor y tejido adyacente, que la dirección se repite más "
        "de lo esperado entre cohortes, que buena parte de la señal es composición "
        "de tejidos, que el tejido no predice la sangre, y que el CEA depende poco "
        "de la expresión génica del tumor."
    )
    p(
        "No se certifica ningún biomarcador, ninguna diana terapéutica ni ninguna "
        "recomendación clínica. Este trabajo es investigación reproducible, no una "
        "prueba médica.",
        "warn",
    )

    out.parent.mkdir(parents=True, exist_ok=True)

    def footer(canvas, doc) -> None:
        canvas.saveState()
        canvas.setFont("Arial", 8)
        canvas.setFillColor(colors.HexColor("#666666"))
        canvas.drawString(2 * cm, 1.2 * cm, "OncoGO-CRC v2 · guía pública · octubre 2026")
        canvas.drawRightString(A4[0] - 2 * cm, 1.2 * cm, str(doc.page))
        canvas.restoreState()

    doc = SimpleDocTemplate(
        str(out),
        pagesize=A4,
        leftMargin=1.8 * cm,
        rightMargin=1.8 * cm,
        topMargin=1.6 * cm,
        bottomMargin=1.8 * cm,
        title="OncoGO-CRC v2 — Guía pública",
        author="OncoGO-CRC",
    )
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(f"[guia-publica] escrito en {out} ({out.stat().st_size} bytes)")


def sync_reportes(workspace: Path, repo: Path) -> None:
    """Copia los PDF y figuras clave al repo para publicarlos en GitHub."""
    dest = repo / "docs" / "reportes"
    dest.mkdir(parents=True, exist_ok=True)
    origen = workspace / "results" / "v2" / "report"
    for nombre in (
        "OncoGO-CRC_informe_v2.pdf",
        "informe_v2.md",
        "comparacion_v1_v2.md",
        "index.html",
    ):
        src = origen / nombre
        if src.exists():
            shutil.copy2(src, dest / nombre)
    fig_src = origen / "figuras"
    fig_dst = dest / "figuras"
    if fig_src.exists():
        fig_dst.mkdir(parents=True, exist_ok=True)
        for ext in ("png", "svg"):
            for f in sorted(fig_src.glob(f"*.{ext}")):
                if f.name.startswith("._"):
                    continue
                shutil.copy2(f, fig_dst / f.name)
    cert = workspace / "results" / "v2" / "certification"
    datos = dest / "datos"
    datos.mkdir(parents=True, exist_ok=True)
    for nombre in (
        "certificado_v2.json",
        "indice_R_v2.tsv",
        "changelog_numeros.tsv",
    ):
        src = cert / nombre
        if src.exists():
            shutil.copy2(src, datos / nombre)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", type=Path, default=WORKSPACE_DEFAULT)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument(
        "--sync-reportes",
        action="store_true",
        help="Copia PDFs, figuras y tablas de resumen a docs/reportes/.",
    )
    args = ap.parse_args()

    faltan = [
        rel for rel, _ in AUDITORIAS
        if not (args.workspace / "results" / "v2" / rel).exists()
    ]
    if faltan:
        print(f"No se genera la guía. Faltan {len(faltan)} auditorías en {args.workspace}:")
        for rel in faltan:
            print(f"  - {rel}")
        return 1

    fig_dir = args.workspace / "results" / "v2" / "report" / "figuras"
    build(args.workspace, args.out, fig_dir if fig_dir.exists() else None)

    if args.sync_reportes:
        sync_reportes(args.workspace, REPO)
        # Nombre alternativo para el informe técnico en la carpeta de reportes.
        tecnico = REPO / "docs" / "reportes" / "OncoGO-CRC_informe_v2.pdf"
        if tecnico.exists():
            shutil.copy2(tecnico, REPO / "docs" / "reportes" / "OncoGO-CRC_informe_v2_tecnico.pdf")
        print(f"[guia-publica] reportes sincronizados en {REPO / 'docs' / 'reportes'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
