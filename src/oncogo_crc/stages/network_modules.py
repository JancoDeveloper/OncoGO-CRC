from __future__ import annotations


from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import read_tsv, write_tsv
from oncogo_crc.utils.network_utils import correlation_network, hub_proteins


def run(settings: Settings) -> None:
    """Red proteína–proteína por correlación y detección de hubs."""
    primary = settings.config.get("project", {}).get("primary_dataset", "cptac_colorectal")
    mat_path = settings.processed_matrix_path(primary, "zscore")
    if not mat_path.exists():
        print("[network] sin matriz z-score primaria; omitiendo.")
        return
    df = read_tsv(mat_path)
    meta = {"protein_id", "gene_symbol", "uniprot_accession"}
    sample_cols = [c for c in df.columns if c not in meta]
    mat = df.set_index("protein_id")[sample_cols].astype(float)
    net_cfg = settings.config.get("network", {})
    G = correlation_network(mat, threshold=float(net_cfg.get("correlation_threshold", 0.6)))
    hubs = hub_proteins(G, min_degree=int(net_cfg.get("min_degree_for_hub", 8)))
    out_dir = settings.resolve(settings.config["paths"]["results_dir"]) / "networks"
    write_tsv(hubs, out_dir / f"{primary}_hub_proteins.tsv")
    print(f"[network] hubs -> results/networks/{primary}_hub_proteins.tsv")
