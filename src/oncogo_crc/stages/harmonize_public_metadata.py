"""Une mapping PDC TMT10 + biospecimen/clínica y normaliza matrices PRIDE/plasma."""

from __future__ import annotations

import os

from oncogo_crc.settings import Settings
from oncogo_crc.utils.clinical_demographics import write_demographics_artifacts
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
    if settings.demo:
        print(
            "[harmonize_public_metadata] modo demo: se omite unión PDC/PRIDE (las muestras sintéticas "
            "vienen de `download` → `_demo_*_samples.tsv`)."
        )
        return

    token = os.environ.get(settings.config.get("download_real", {}).get("auth", {}).get("pdc_token_env", "PDC_TOKEN") or "PDC_TOKEN")
    clinical_cfg = dict(settings.config.get("clinical_filters", {}) or {})
    cpt_h = cfg.get("cptac") or {}
    cms_sup = str(cpt_h.get("cms_subtype_supplement_tsv", "") or "").strip()
    if cms_sup:
        clinical_cfg["cms_subtype_supplement_tsv"] = str(settings.resolve(cms_sup))

    ds_all = settings.datasets.get("datasets", {})

    cpt = cfg.get("cptac") or {}
    if bool(cpt.get("enabled", True)) and "cptac_colorectal" in ds_all:
        d = ds_all["cptac_colorectal"]
        dr = d.get("download_real", {}) or {}
        study = str(cpt.get("pdc_study_submitter_id") or dr.get("pdc_study_submitter_id") or "PDC000116").strip()
        m_path = settings.resolve(d["matrix_path"])
        meta_path = settings.resolve(d["sample_metadata"])
        print(f"[harmonize_public_metadata] CPTAC PDC {study} -> {settings.relative_to_project(meta_path)}")
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
        print(f"[harmonize_public_metadata] PRIDE MaxQuant LFQ -> {m_path.name} + {settings.relative_to_project(meta_path)}")
        info = harmonize_pride_maxquant_lfq(matrix_path=m_path, metadata_out=meta_path)
        print(f"  proteins={info['n_proteins']} lfq_columns={info['n_lfq_samples']}")

    pl = cfg.get("plasma") or {}
    if bool(pl.get("enabled", True)) and "plasma_crc_2024" in ds_all:
        d = ds_all["plasma_crc_2024"]
        m_path = settings.resolve(d["matrix_path"])
        meta_path = settings.resolve(d["sample_metadata"])
        accs = (d.get("download_real") or {}).get("accessions") or []
        accs_s = {str(a).strip() for a in accs}
        use_peptide_dia = "PXD036594" in accs_s
        if use_peptide_dia:
            print(f"[harmonize_public_metadata] plasma DIA (PXD036594) -> {m_path.name} + {settings.relative_to_project(meta_path)}")
            info = harmonize_plasma_dia_pq500(matrix_path=m_path, metadata_out=meta_path)
            print(
                f"  peptides={info['n_peptides']} ratio_cols={info['n_ratio_columns']} "
                f"run_prefixes={info['n_run_prefixes']}"
            )
        else:
            rel_meta = settings.relative_to_project(meta_path)
            if meta_path.exists() and meta_path.stat().st_size > 50:
                tail = (
                    f"Metadata ya presente ({rel_meta}); no se aplica harmonizador péptido/ratios "
                    "(típico tras DIA-NN + PXD015163 en config)."
                )
            else:
                tail = (
                    f"Curar {rel_meta} desde PRIDE SDRF o el artículo (condition: crc|healthy o equivalente)."
                )
            print(
                "[harmonize_public_metadata] plasma: harmonizador péptido/ratios sólo aplica a PXD036594 "
                "(Skyline ratios). "
                f"Accession(es) declarada(s) en datasets.yaml: {sorted(accs_s) or '—'}. "
                + tail
            )

    try:
        demo_path = write_demographics_artifacts(settings)
        print(
            f"[harmonize_public_metadata] resumen demografía/clínica -> "
            f"{settings.relative_to_project(demo_path)}"
        )
    except OSError as exc:
        print(f"[harmonize_public_metadata] resumen demografía omitido: {exc}")
