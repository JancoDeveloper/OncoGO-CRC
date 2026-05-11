from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any

import pandas as pd

from oncogo_crc.settings import Settings


def _cptac_version() -> str:
    try:
        import importlib.metadata as im

        return im.version("cptac")
    except Exception:
        return "unknown"


def _flatten_proteomics_columns(df: pd.DataFrame) -> pd.DataFrame:
    if isinstance(df.columns, pd.MultiIndex):
        # Típico cptac: nivel 0 = gen símbolo
        df = df.copy()
        df.columns = [str(x[0]) if x[0] is not None else str(x) for x in df.columns]
    return df


def _proteomics_to_pipeline_matrix(prot: pd.DataFrame) -> pd.DataFrame:
    """Muestras en columnas → protein_id + columnas sample_id."""
    prot = _flatten_proteomics_columns(prot)
    m = prot.T.copy()
    m.insert(0, "protein_id", m.index.astype(str))
    m = m.reset_index(drop=True)
    return m


def _clinical_to_metadata(clin: pd.DataFrame) -> pd.DataFrame:
    out = clin.reset_index()
    if out.columns[0] != "sample_id":
        out = out.rename(columns={out.columns[0]: "sample_id"})
    out["sample_id"] = out["sample_id"].astype(str)
    if "condition" not in out.columns:
        candidates = [
            c
            for c in out.columns
            if any(
                k in c.lower()
                for k in ("tissue", "tumor", "sample_type", "histologic", "diagnosis", "type")
            )
        ]
        if candidates:
            col = candidates[0]
            ser = out[col].astype(str).str.lower()
            out["condition"] = ser.apply(
                lambda x: "tumor"
                if any(k in x for k in ("tumor", "cancer", "crc", "carcinoma", "primary"))
                else ("normal" if any(k in x for k in ("normal", "adjacent", "healthy", "control")) else "unknown")
            )
        else:
            out["condition"] = "unknown"
    return out


def run(settings: Settings) -> None:
    """
    Exporta CPTAC Colon desde el paquete `cptac` a las rutas declaradas en datasets.yaml (MVP rápido).

    Requiere: pip install 'oncogo-crc[ingest]' o pip install cptac
    """
    if settings.demo:
        print("[ingest_cptac] omitido en modo demo.")
        return

    try:
        import cptac
    except ImportError:
        print("[ingest_cptac] Falta dependencia: pip install 'oncogo-crc[ingest]'  (o pip install cptac)")
        return

    ds_cfg = settings.datasets.get("datasets", {}).get("cptac_colorectal")
    if not ds_cfg or ds_cfg.get("enabled") is False:
        print("[ingest_cptac] cptac_colorectal deshabilitado o ausente en datasets.yaml.")
        return

    co = cptac.Colon()
    prot = co.get_proteomics()
    clin = co.get_clinical()
    if prot is None or clin is None:
        print("[ingest_cptac] get_proteomics() o get_clinical() devolvió None.")
        return

    matrix = _proteomics_to_pipeline_matrix(prot)
    meta = _clinical_to_metadata(clin)

    m_path = settings.resolve(ds_cfg["matrix_path"])
    md_path = settings.resolve(ds_cfg["sample_metadata"])
    m_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(m_path, sep="\t", index=False)
    meta.to_csv(md_path, sep="\t", index=False)

    prov: dict[str, Any] = {
        "ingest": "cptac.Colon",
        "cptac_package_version": _cptac_version(),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "matrix_path": str(m_path.relative_to(settings.root)),
        "metadata_path": str(md_path.relative_to(settings.root)),
        "n_proteins": int(matrix.shape[0]),
        "n_matrix_columns_samples": int(matrix.shape[1] - 1),
        "n_metadata_rows": int(meta.shape[0]),
        "note": "Validar en tesis que la capa proteómica exportada corresponde al análisis deseado (global vs fosfo, etc.).",
    }
    prov_path = m_path.parent / "cptac_ingest_provenance.json"
    prov_path.write_text(json.dumps(prov, indent=2), encoding="utf-8")
    print(f"[ingest_cptac] exportado -> {m_path.name}, {md_path.name}; provenance -> {prov_path.name}")
