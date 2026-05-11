"""Une mapping PDC TMT10 + biospecimen/clínica y normaliza matrices PRIDE/plasma."""

from __future__ import annotations

import os

from oncogo_crc.settings import Settings
from oncogo_crc.utils.public_metadata_harmonize import (
    harmonize_cptac_colon_pdc,
    harmonize_plasma_dia_pq500,
    harmonize_pride_maxquant_lfq,
)


def run(settings: Settings) -> None:
    cfg = settings.config.get("harmonize_public_metadata") or {}
    if not bool(cfg.get("enabled", True)):
        print("[harmonize_public_metadata] deshabilitado en config.yaml.")
        return

    token = os.environ.get(settings.config.get("download_real", {}).get("auth", {}).get("pdc_token_env", "PDC_TOKEN") or "PDC_TOKEN")
    clinical_cfg = settings.config.get("clinical_filters", {}) or {}

    ds_all = settings.datasets.get("datasets", {})

    cpt = cfg.get("cptac") or {}
    if bool(cpt.get("enabled", True)) and "cptac_colorectal" in ds_all:
        d = ds_all["cptac_colorectal"]
        dr = d.get("download_real", {}) or {}
        study = str(cpt.get("pdc_study_submitter_id") or dr.get("pdc_study_submitter_id") or "PDC000116").strip()
        m_path = settings.resolve(d["matrix_path"])
        meta_path = settings.resolve(d["sample_metadata"])
        print(f"[harmonize_public_metadata] CPTAC PDC {study} -> {meta_path.relative_to(settings.root)}")
        info = harmonize_cptac_colon_pdc(
            matrix_path=m_path,
            metadata_out=meta_path,
            pdc_study_submitter_id=study,
            clinical_cfg=clinical_cfg,
            token=token or None,
        )
        print(f"  rows={info['n_metadata_rows']} tumor={info['n_tumor']} normal={info['n_normal']} unknown={info['n_unknown']}")

    pr = cfg.get("pride") or {}
    if bool(pr.get("enabled", True)) and "pride_crc_meta" in ds_all:
        d = ds_all["pride_crc_meta"]
        m_path = settings.resolve(d["matrix_path"])
        meta_path = settings.resolve(d["sample_metadata"])
        print(f"[harmonize_public_metadata] PRIDE MaxQuant LFQ -> {m_path.name} + {meta_path.relative_to(settings.root)}")
        info = harmonize_pride_maxquant_lfq(matrix_path=m_path, metadata_out=meta_path)
        print(f"  proteins={info['n_proteins']} lfq_columns={info['n_lfq_samples']}")

    pl = cfg.get("plasma") or {}
    if bool(pl.get("enabled", True)) and "plasma_crc_2024" in ds_all:
        d = ds_all["plasma_crc_2024"]
        m_path = settings.resolve(d["matrix_path"])
        meta_path = settings.resolve(d["sample_metadata"])
        print(f"[harmonize_public_metadata] plasma DIA -> {m_path.name} + {meta_path.relative_to(settings.root)}")
        info = harmonize_plasma_dia_pq500(matrix_path=m_path, metadata_out=meta_path)
        print(
            f"  peptides={info['n_peptides']} ratio_cols={info['n_ratio_columns']} "
            f"run_prefixes={info['n_run_prefixes']}"
        )
