from __future__ import annotations

from oncogo_crc.settings import Settings
from oncogo_crc.utils.ids import harmonize_protein_table
from oncogo_crc.utils.io import read_tsv, write_tsv


def run(settings: Settings) -> None:
    """Armoniza identificadores proteicos con el mapa maestro y colapsa por gen."""
    mapping = read_tsv(settings.protein_id_map_path())
    interim = settings.resolve(settings.config["paths"]["interim_dir"])
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
