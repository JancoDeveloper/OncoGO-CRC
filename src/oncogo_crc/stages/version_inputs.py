from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _matrix_columns(path: Path) -> list[str]:
    suf = path.suffix.lower()
    if suf in (".xlsx", ".xls"):
        return pd.read_excel(path, engine="openpyxl", nrows=1).columns.astype(str).tolist()
    head = path.read_bytes()[:2]
    if head == b"\x1f\x8b":
        return pd.read_csv(path, sep="\t", compression="gzip", nrows=0).columns.tolist()
    return pd.read_csv(path, sep="\t", nrows=0).columns.tolist()


def _read_metadata(path: Path) -> pd.DataFrame:
    suf = path.suffix.lower()
    if suf in (".xlsx", ".xls"):
        return pd.read_excel(path, engine="openpyxl")
    head = path.read_bytes()[:2]
    if head == b"\x1f\x8b":
        return pd.read_csv(path, sep="\t", compression="gzip")
    return pd.read_csv(path, sep="\t")


def _read_matrix_rowcount(path: Path) -> int:
    suf = path.suffix.lower()
    if suf in (".xlsx", ".xls"):
        return int(pd.read_excel(path, engine="openpyxl").shape[0])
    head = path.read_bytes()[:2]
    if head == b"\x1f\x8b":
        return int(pd.read_csv(path, sep="\t", compression="gzip", low_memory=False).shape[0])
    return int(pd.read_csv(path, sep="\t", low_memory=False).shape[0])


def run(settings: Settings) -> None:
    """
    Versiona snapshots de matrices/metadata actuales para trazabilidad.
    """
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    base = settings.resolve(f"data/versions/{ts}")
    base.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    docs: list[str] = [
        "# Inventario de matrices versionadas",
        "",
        f"- version_id: `{ts}`",
        "- objetivo: congelar matrices/metadata activas para definir la versión de trabajo reproducible.",
        "",
    ]

    for ds, cfg in settings.datasets.get("datasets", {}).items():
        if cfg.get("enabled") is False:
            continue
        matrix_path = settings.resolve(cfg["matrix_path"])
        metadata_path = settings.resolve(cfg["sample_metadata"])
        if not matrix_path.exists() or not metadata_path.exists():
            rows.append(
                {
                    "dataset": ds,
                    "status": "missing_inputs",
                    "matrix_path": str(matrix_path),
                    "metadata_path": str(metadata_path),
                }
            )
            docs.extend([f"## {ds}", "- estado: missing_inputs", ""])
            continue

        ds_dir = base / ds
        ds_dir.mkdir(parents=True, exist_ok=True)
        matrix_dst = ds_dir / "matrix.tsv"
        metadata_dst = ds_dir / "metadata.tsv"
        matrix_dst.write_bytes(matrix_path.read_bytes())
        metadata_dst.write_bytes(metadata_path.read_bytes())

        matrix_cols = _matrix_columns(matrix_dst)
        metadata_df = _read_metadata(metadata_dst)
        sample_cols = max(len(matrix_cols) - 1, 0)
        first_col = matrix_cols[0] if matrix_cols else ""

        rows.append(
            {
                "dataset": ds,
                "status": "versioned",
                "version_id": ts,
                "matrix_snapshot": settings.relative_to_project(matrix_dst),
                "metadata_snapshot": settings.relative_to_project(metadata_dst),
                "matrix_sha256": _sha256(matrix_dst),
                "metadata_sha256": _sha256(metadata_dst),
                "n_proteins": _read_matrix_rowcount(matrix_dst),
                "n_samples_matrix": int(sample_cols),
                "n_rows_metadata": int(metadata_df.shape[0]),
                "matrix_first_column": first_col,
                "matrix_function": "abundancia proteica por muestra",
                "metadata_function": "anotacion clinica y grupos analiticos",
            }
        )

        docs.extend(
            [
                f"## {ds}",
                "- estado: versioned",
                "- matrix_function: abundancia proteica por muestra",
                "- metadata_function: anotacion clinica y grupos analiticos",
                f"- matrix_first_column: `{first_col}`",
                f"- n_samples_matrix: {sample_cols}",
                f"- n_rows_metadata: {int(metadata_df.shape[0])}",
                "",
            ]
        )

    inv_dir = settings.resolve("data/versions")
    inv_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(inv_dir / "matrix_versions.tsv", sep="\t", index=False)
    (inv_dir / "latest_version.txt").write_text(ts, encoding="utf-8")
    report_path = settings.resolve("reports/executive_summary/matrix_inventory.md")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(docs), encoding="utf-8")
    print(f"[version_inputs] version_id={ts} -> data/versions/{ts}")
