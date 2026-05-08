from __future__ import annotations

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv, write_tsv
from oncogo_crc.utils.translation_utils import build_translation_table


def run(settings: Settings) -> None:
    """Fusión de capas (diferencial, GO implícito vía genes, red, ML) y niveles de evidencia."""
    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")
    res = settings.resolve(settings.config["paths"]["results_dir"])
    diff_path = res / "differential" / "cptac_tumor_vs_normal_diff.tsv"
    if not diff_path.exists():
        print("[clinical] sin resultado diferencial núcleo; omitiendo.")
        return
    diff = read_tsv(diff_path)
    hubs_path = res / "networks" / f"{primary}_hub_proteins.tsv"
    hubs = read_tsv(hubs_path) if hubs_path.exists() else None
    ml_path = res / "ml" / f"{primary}_feature_ranking.tsv"
    ml = read_tsv(ml_path) if ml_path.exists() else None
    if ml is not None:
        ml = ml.set_index("protein_id")
    tiers = list(settings.config.get("translation", {}).get("evidence_tiers", []))
    out = build_translation_table(diff, hubs, ml, tiers=tiers or ["validated", "hypothesis_generating", "future_validation"])
    out_dir = res / "clinical_translation"
    write_tsv(out, out_dir / f"{primary}_translation_ranked.tsv")
    print(f"[clinical] ranking -> results/clinical_translation/{primary}_translation_ranked.tsv")
