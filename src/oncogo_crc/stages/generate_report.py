from __future__ import annotations

from pathlib import Path

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv
from oncogo_crc.utils.report_utils import write_executive_summary
from oncogo_crc.utils.thesis_web_digest import build_quantitative_bullets


def _preview_md(path: Path, n: int = 8) -> str:
    if not path.exists():
        return "_No generado aún._"
    df = read_tsv(path).head(n)
    return "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "| " + " | ".join(["---"] * len(df.columns)) + " |"] + ["| " + " | ".join(map(str, row)) + " |" for row in df.values.tolist()])


def run(settings: Settings) -> None:
    """Genera el informe Markdown consolidado para la tesis / comité."""
    res = settings.resolve(settings.config["paths"]["results_dir"])
    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")

    sections = {
        "Contexto": (
            "OncoGO-CRC integra proteómica pública de CCR con Gene Ontology (ORA + GSEA pre-ranked), "
            "firma multi-cohorte CPTAC∩PRIDE y priorización mediante modelos interpretables."
        ),
        "Resumen cuantitativo (FDR)": build_quantitative_bullets(settings),
        "Diferencial (extracto)": _preview_md(res / "differential" / "cptac_tumor_vs_normal_diff.tsv"),
        "GSEA (extracto)": _preview_md(res / "gsea" / "cptac_tumor_vs_normal_gsea_go.tsv"),
        "GO ORA (extracto)": _preview_md(res / "go" / "cptac_tumor_vs_normal_go.tsv"),
        "Reproducibilidad CPTAC∩PRIDE": _preview_md(
            res / "reproducibility" / "cptac_pride_differential_overlap_summary.tsv"
        ),
        "Red (hubs)": _preview_md(res / "networks" / f"{primary}_hub_proteins.tsv"),
        "ML (ranking)": _preview_md(res / "ml" / f"{primary}_feature_ranking.tsv"),
        "Traslación": _preview_md(res / "clinical_translation" / f"{primary}_translation_ranked.tsv"),
    }
    out = settings.resolve(settings.config["paths"]["reports_dir"]) / "executive_summary" / "results_summary.md"
    write_executive_summary(out, sections=sections)
    print(f"[report] {settings.relative_to_project(out)}")
