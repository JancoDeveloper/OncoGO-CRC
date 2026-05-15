from __future__ import annotations

from oncogo_crc.settings import Settings
from oncogo_crc.utils.export_bundle import write_protein_context_bundle


def run(settings: Settings) -> None:
    cfg = settings.config.get("export_bundle", {}) or {}
    if not bool(cfg.get("enabled", True)):
        print("[export_bundle] deshabilitado (export_bundle.enabled=false).")
        return
    res = settings.resolve(settings.config["paths"]["results_dir"])
    comp = str(cfg.get("primary_comparison", "cptac_tumor_vs_normal"))
    out = res / "export" / "protein_context_bundle.json"
    write_protein_context_bundle(
        settings,
        out,
        comparison=comp,
        top=int(cfg.get("diff_top", 250)),
        uniprot_top_n=int(cfg.get("uniprot_top_n", 0)),
    )
    print(f"[export_bundle] -> {settings.relative_to_project(out)}")
