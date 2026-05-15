from __future__ import annotations

from oncogo_crc.settings import Settings
from oncogo_crc.utils.ids import effective_protein_mapping, harmonize_protein_table
from oncogo_crc.utils.io import read_tsv, write_tsv
from oncogo_crc.utils.mapping_diagnostics import build_mapping_diagnostics_row, write_mapping_diagnostics_rows


def run(settings: Settings) -> None:
    """Armoniza identificadores proteicos con el mapa maestro y colapsa por gen."""
    mapping = effective_protein_mapping(settings)
    map_cfg = settings.config.get("mapping", {}) or {}
    n_sup = sum(
        1
        for rel in map_cfg.get("supplementary_map_paths", []) or []
        if settings.resolve(str(rel).strip()).exists()
    )
    if n_sup:
        print(f"[harmonize] mapa efectivo ({len(mapping)} filas; suplementos presentes: {n_sup})")

    interim = settings.resolve(settings.config["paths"]["interim_dir"])
    repro_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "reproducibility"
    diag_rows: list[dict] = []
    for key, ds in settings.datasets["datasets"].items():
        matrix_path = settings.resolve(ds["matrix_path"])
        if not matrix_path.exists():
            continue
        protein_col = ds.get("protein_column", "protein_id")
        expr = read_tsv(matrix_path)
        out = harmonize_protein_table(expr, mapping, protein_col=protein_col)
        out_path = interim / f"{key}_harmonized.tsv"
        write_tsv(out, out_path)
        print(f"[harmonize] {key} -> {settings.relative_to_project(out_path)}")
        diag_rows.append(
            build_mapping_diagnostics_row(dataset_key=key, harmonized_path=out_path, mapping=mapping)
        )
    if diag_rows:
        cov_path = repro_dir / "mapping_coverage.tsv"
        write_mapping_diagnostics_rows(diag_rows, cov_path)
        print(f"[harmonize] mapping_coverage -> {settings.relative_to_project(cov_path)}")
