from __future__ import annotations

import json

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import write_tsv


def _write_manifest(settings: Settings, rows: list[dict]) -> None:
    out = settings.resolve("data/raw/manifest.tsv")
    write_tsv(pd.DataFrame(rows), out)


def _demo_synthetic(settings: Settings) -> None:
    rng = np.random.default_rng(settings.config.get("project", {}).get("seed", 42))
    dcfg = settings.config.get("demo", {})
    n_prot = int(dcfg.get("n_proteins", 120))
    nt = int(dcfg.get("n_samples_tumor", 24))
    nn = int(dcfg.get("n_samples_normal", 24))
    n_shift = int(dcfg.get("n_shifted_proteins", 25))

    proteins = [f"P{i:04d}" for i in range(1, n_prot + 1)]
    genes = [f"GENE{i:03d}" for i in range(1, n_prot + 1)]
    mapping = pd.DataFrame({"protein_id": proteins, "gene_symbol": genes, "uniprot_accession": [""] * n_prot})
    map_path = settings.resolve("data/interim/.demo_protein_mapping.tsv")
    write_tsv(mapping, map_path)

    tumor_cols = [f"CPTAC-T{i:02d}" for i in range(1, nt + 1)]
    normal_cols = [f"CPTAC-N{i:02d}" for i in range(1, nn + 1)]
    base = rng.lognormal(mean=1.0, sigma=0.45, size=(n_prot, nt + nn))
    mat = pd.DataFrame(base, index=proteins, columns=tumor_cols + normal_cols)
    for i in range(n_shift):
        mat.iloc[i, :nt] *= rng.uniform(1.4, 2.2)

    expr = mat.reset_index()
    expr = expr.rename(columns={expr.columns[0]: "protein_id"})
    cptac_path = settings.resolve(settings.datasets["datasets"]["cptac_colorectal"]["matrix_path"])
    write_tsv(expr, cptac_path)

    meta_rows = [{"sample_id": s, "condition": "tumor", "stage_bin": "advanced", "msi_status": "MSS", "batch": "b1"} for s in tumor_cols]
    meta_rows += [{"sample_id": s, "condition": "normal", "stage_bin": "early", "msi_status": "MSS", "batch": "b1"} for s in normal_cols]
    demo_meta = settings.resolve("data/interim/_demo_cptac_colorectal_samples.tsv")
    write_tsv(pd.DataFrame(meta_rows), demo_meta)

    go_pool = []
    for go_id, ns in [
        ("GO:0006955", "BP"),
        ("GO:0006281", "BP"),
        ("GO:0007165", "BP"),
        ("GO:0008284", "BP"),
        ("GO:0043065", "BP"),
        ("GO:0001525", "BP"),
        ("GO:0006355", "BP"),
        ("GO:0005975", "BP"),
        ("GO:0005506", "MF"),
        ("GO:0005524", "MF"),
        ("GO:0003723", "MF"),
        ("GO:0005886", "CC"),
        ("GO:0005737", "CC"),
        ("GO:0005634", "CC"),
    ]:
        go_pool.append((go_id, ns))
    go_rows = []
    for g in genes:
        k = int(rng.integers(2, 6))
        picks = rng.choice(len(go_pool), size=k, replace=False)
        for j in picks:
            gid, ns = go_pool[int(j)]
            go_rows.append({"gene_symbol": g, "go_id": gid, "namespace": ns})
    go_path = settings.resolve("data/interim/.demo_go_gene_long.tsv")
    write_tsv(pd.DataFrame(go_rows), go_path)

    n_pt, n_pn = 18, 18
    pt = [f"PLASMA-CRC-{i:02d}" for i in range(1, n_pt + 1)]
    pn = [f"PLASMA-HC-{i:02d}" for i in range(1, n_pn + 1)]
    sub_prot = proteins[: min(40, n_prot)]
    basep = rng.lognormal(mean=1.0, sigma=0.5, size=(len(sub_prot), n_pt + n_pn))
    mp = pd.DataFrame(basep, index=sub_prot, columns=pt + pn)
    for i in range(min(8, len(sub_prot))):
        mp.iloc[i, :n_pt] *= rng.uniform(1.2, 1.9)
    plasma_expr = mp.reset_index()
    plasma_expr = plasma_expr.rename(columns={plasma_expr.columns[0]: "protein_id"})
    plasma_path = settings.resolve(settings.datasets["datasets"]["plasma_crc_2024"]["matrix_path"])
    write_tsv(plasma_expr, plasma_path)
    plasma_meta = [{"sample_id": s, "condition": "crc", "batch": "b1"} for s in pt] + [
        {"sample_id": s, "condition": "healthy", "batch": "b1"} for s in pn
    ]
    write_tsv(pd.DataFrame(plasma_meta), settings.resolve("data/interim/_demo_plasma_crc_2024_samples.tsv"))


def run(settings: Settings) -> None:
    """Registra datasets: en producción valida rutas; en `--demo` genera datos sintéticos reproducibles."""
    rows = []
    if settings.demo:
        _demo_synthetic(settings)
        for key, ds in settings.datasets["datasets"].items():
            p = settings.resolve(ds["matrix_path"])
            rows.append(
                {
                    "dataset": key,
                    "matrix_path": settings.relative_to_project(p),
                    "status": "demo_synthetic",
                }
            )
        _write_manifest(settings, rows)
        (settings.resolve("data/raw/demo_readme.json")).write_text(
            json.dumps({"note": "Matrices sintéticas para validar el pipeline sin descargar CPTAC/PRIDE."}, indent=2),
            encoding="utf-8",
        )
        return

    for key, ds in settings.datasets["datasets"].items():
        p = settings.resolve(ds["matrix_path"])
        status = "missing"
        if p.exists():
            status = "present"
        rows.append({"dataset": key, "matrix_path": settings.relative_to_project(p), "status": status})
    _write_manifest(settings, rows)
    missing = [r for r in rows if r["status"] == "missing"]
    if missing:
        msg = (
            "Faltan matrices en data/raw/. Coloque los TSV según config/datasets.yaml "
            "o ejecute `oncogo-crc run --demo` para datos sintéticos.\n"
            + "\n".join(f"  - {m['matrix_path']}" for m in missing)
        )
        print(msg)
