"""
Harmonización de metadata y matrices públicas (PDC CPTAC, PRIDE MaxQuant, plasma DIA).

Descarga auxiliares PDC vía GraphQL (URLs firmadas), cruza mapping TMT10 + biospecimen
y escribe TSV alineado a columnas de la matriz proteómica.
"""

from __future__ import annotations

import json
import re
import urllib.request
from io import BytesIO
from pathlib import Path
from typing import Any

import pandas as pd

PDC_GRAPHQL_ENDPOINT = "https://pdc.cancer.gov/graphql"


def _post_graphql(query: str, variables: dict[str, Any] | None = None, token: str | None = None) -> dict[str, Any]:
    payload = json.dumps({"query": query, "variables": variables or {}}).encode("utf-8")
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    req = urllib.request.Request(PDC_GRAPHQL_ENDPOINT, data=payload, headers=headers, method="POST")
    with urllib.request.urlopen(req, timeout=120) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    if "errors" in data:
        raise RuntimeError(str(data["errors"]))
    return data.get("data", {})


def _pdc_list_files(pdc_study_id: str, token: str | None) -> list[dict[str, Any]]:
    fq = (
        "query StudyFiles($pdc_study_id: String!) { "
        "filesPerStudy(pdc_study_id: $pdc_study_id, acceptDUA: true) { "
        "file_name file_size signedUrl { url } } }"
    )
    fdata = _post_graphql(fq, {"pdc_study_id": pdc_study_id}, token=token)
    return list(fdata.get("filesPerStudy") or [])


def _pdc_download_named_file(files: list[dict[str, Any]], substr: str, *, ends: str = ".xlsx") -> bytes:
    for f in files:
        name = str(f.get("file_name") or "")
        if substr in name and name.endswith(ends):
            surl = f.get("signedUrl") or {}
            url = surl.get("url") if isinstance(surl, dict) else None
            if url:
                with urllib.request.urlopen(url, timeout=180) as r:
                    return r.read()
    raise FileNotFoundError(f"PDC: no file matching {substr!r} ({ends})")


def _tmt_mapping_long(mapping_df: pd.DataFrame) -> pd.DataFrame:
    cols = list(mapping_df.columns)
    rows_out: list[dict[str, str]] = []
    i = 0
    while i < len(cols):
        c0 = str(cols[i])
        if "Participant ID" in c0 and i + 1 < len(cols):
            c1 = str(cols[i + 1])
            if "Specimen Label" in c1:
                pid_col, slab_col = cols[i], cols[i + 1]
                for ridx in range(len(mapping_df)):
                    pid = mapping_df.iloc[ridx][pid_col]
                    slab = mapping_df.iloc[ridx][slab_col]
                    if pd.isna(slab):
                        continue
                    s = str(slab).strip()
                    if not s or s.lower() == "colonref":
                        continue
                    rows_out.append({"participant_id": str(pid).strip(), "specimen_label": s})
                i += 2
                continue
        i += 1
    return pd.DataFrame(rows_out)


def _pathology_to_condition(status: str) -> str:
    s = (status or "").strip().lower()
    # "non-malignant" contiene la subcadena "malignant": comprobar primero el caso benigno.
    if "non-malignant" in s or "non malignant" in s or "nonmalignant" in s.replace("-", ""):
        return "normal"
    if "malignant" in s:
        return "tumor"
    return "unknown"


def _stage_to_bin(stage_val: Any, early_tokens: list[str], adv_tokens: list[str]) -> str:
    if pd.isna(stage_val):
        return ""
    s = str(stage_val).strip().lower()
    for t in sorted((str(x).lower() for x in adv_tokens if x), key=len, reverse=True):
        if t and t in s:
            return "advanced"
    for t in sorted((str(x).lower() for x in early_tokens if x), key=len, reverse=True):
        if t and t in s:
            return "early"
    if "stage" in s:
        if any(x in s for x in [" iii", " iv", "iiib", "iiic", "iva", "ivb"]):
            return "advanced"
        if any(x in s for x in ["stage iia", "stage iib", "stage iic", "stage ia", "stage ib"]):
            return "early"
    return ""


def harmonize_cptac_colon_pdc(
    *,
    matrix_path: Path,
    metadata_out: Path,
    pdc_study_submitter_id: str,
    clinical_cfg: dict[str, Any],
    token: str | None,
) -> dict[str, Any]:
    files = _pdc_list_files(pdc_study_submitter_id, token)
    map_bytes = _pdc_download_named_file(files, "TMT10_Label_to_Sample_Mapping_File_COprospective")
    bio_bytes = _pdc_download_named_file(files, "Biospecimens")
    clin_bytes = _pdc_download_named_file(files, "Clinical_Data")

    map_df = pd.read_excel(BytesIO(map_bytes), engine="openpyxl")
    bio_df = pd.read_excel(BytesIO(bio_bytes), engine="openpyxl")
    clin_df = pd.read_excel(BytesIO(clin_bytes), engine="openpyxl", sheet_name="Sheet1")

    long_map = _tmt_mapping_long(map_df)
    long_map = long_map.drop_duplicates(subset=["specimen_label"], keep="first")

    bio_sub = bio_df[["Participant ID", "Specimen Label", "Pathological Status"]].copy()
    bio_sub.columns = ["participant_id", "specimen_label", "pathological_status"]
    bio_sub["specimen_label"] = bio_sub["specimen_label"].astype(str).str.strip()
    bio_sub["condition"] = bio_sub["pathological_status"].map(_pathology_to_condition)

    sm = clinical_cfg.get("stage_mapping", {}) if isinstance(clinical_cfg, dict) else {}
    early_tokens = list(sm.get("early", []) or [])
    adv_tokens = list(sm.get("advanced", []) or [])

    clin_df = clin_df.copy()
    clin_df["_stage_bin"] = clin_df["Tumor Stage (Pathological)"].map(
        lambda x: _stage_to_bin(x, early_tokens, adv_tokens)
    )
    clin_agg = clin_df.sort_values("Participant ID").groupby("Participant ID", as_index=False).first()
    clin_small = clin_agg[["Participant ID", "_stage_bin"]].rename(
        columns={"Participant ID": "participant_id", "_stage_bin": "stage_bin"}
    )

    mdf = pd.read_csv(matrix_path, sep="\t", nrows=0, low_memory=False)
    sample_cols = [c for c in mdf.columns if c != "protein_id"]

    def _specimen_from_matrix_col(col: str) -> str | None:
        c = str(col)
        # Orden: el sufijo más largo primero (evita que "... Unshared Log Ratio" coincida con "... Log Ratio").
        for suf in (" Unshared Log Ratio", " Log Ratio"):
            if c.endswith(suf):
                return c[: -len(suf)]
        return None

    out_rows: list[dict[str, Any]] = []
    for col in sample_cols:
        slab = _specimen_from_matrix_col(col)
        if not slab:
            continue
        b = bio_sub[bio_sub["specimen_label"] == slab]
        if b.empty:
            cond = "unknown"
            pid = ""
        else:
            cond = str(b.iloc[0]["condition"])
            pid = str(b.iloc[0]["participant_id"]).strip()
        st_bin = ""
        if pid:
            cmatch = clin_small[clin_small["participant_id"] == pid]
            if not cmatch.empty:
                st_bin = str(cmatch.iloc[0]["stage_bin"] or "")
        mm = long_map[long_map["specimen_label"] == slab]
        map_pid = str(mm.iloc[0]["participant_id"]) if not mm.empty else pid
        out_rows.append(
            {
                "sample_id": col,
                "condition": cond,
                "participant_id": map_pid or pid,
                "stage_bin": st_bin,
                "msi_status": "",
                "batch": "pdc_tmt10",
            }
        )

    out = pd.DataFrame(out_rows)
    metadata_out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(metadata_out, sep="\t", index=False)

    # La matriz PDC a veces trae columnas de anotación al final (NCBIGeneID, …): quitarlas para alinear con metadata.
    sample_keep = [c for c in sample_cols if _specimen_from_matrix_col(c)]
    full_mat = pd.read_csv(matrix_path, sep="\t", low_memory=False)
    keep = ["protein_id"] + sample_keep
    keep = [c for c in keep if c in full_mat.columns]
    full_mat[keep].to_csv(matrix_path, sep="\t", index=False)

    return {
        "n_matrix_columns": len(sample_keep),
        "n_metadata_rows": len(out),
        "n_tumor": int((out["condition"] == "tumor").sum()),
        "n_normal": int((out["condition"] == "normal").sum()),
        "n_unknown": int((out["condition"] == "unknown").sum()),
    }


def harmonize_pride_maxquant_lfq(
    *,
    matrix_path: Path,
    metadata_out: Path,
) -> dict[str, Any]:
    df = pd.read_csv(matrix_path, sep="\t", low_memory=False)
    df.columns = [str(c).replace("\ufeff", "").strip() for c in df.columns]
    first_col = str(df.columns[0])
    lfq_re = re.compile(r"^LFQ intensity CRC\d+[MNP]$")

    def _first_protein_id(cell: Any) -> str:
        s = "" if pd.isna(cell) else str(cell).strip()
        if not s:
            return ""
        return s.split(";")[0].strip()

    if first_col == "Protein IDs":
        df = df.copy()
        df.insert(0, "protein_id", df[first_col].map(_first_protein_id))
        lfq_cols = [c for c in df.columns if isinstance(c, str) and lfq_re.match(c)]
        if not lfq_cols:
            raise ValueError("PRIDE harmonize: no se hallaron columnas LFQ intensity CRC##[MNP]")
        out_df = df[["protein_id"] + lfq_cols].copy()
        out_df.to_csv(matrix_path, sep="\t", index=False)
    elif first_col == "protein_id":
        lfq_cols = [c for c in df.columns if isinstance(c, str) and lfq_re.match(c)]
        if not lfq_cols:
            raise ValueError(
                "PRIDE harmonize: matriz ya tiene protein_id pero no columnas LFQ; "
                "vuelve a ejecutar download_real para PXD034575."
            )
        out_df = df[["protein_id"] + lfq_cols].copy()
        out_df.to_csv(matrix_path, sep="\t", index=False)
    else:
        raise ValueError(
            f"PRIDE harmonize: primera columna debe ser 'Protein IDs' o 'protein_id', hallada {first_col!r}"
        )

    meta_rows: list[dict[str, str]] = []
    for col in lfq_cols:
        if col.endswith("M"):
            cond = "tumor"
        elif col.endswith("N"):
            cond = "normal"
        elif col.endswith("P"):
            cond = "polyp"
        else:
            cond = "unknown"
        meta_rows.append({"sample_id": col, "condition": cond, "batch": "PXD034575"})

    meta = pd.DataFrame(meta_rows)
    metadata_out.parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(metadata_out, sep="\t", index=False)
    return {"n_proteins": len(out_df), "n_lfq_samples": len(lfq_cols)}


def harmonize_plasma_dia_pq500(
    *,
    matrix_path: Path,
    metadata_out: Path,
) -> dict[str, Any]:
    """
    Plasma PXD036594: matriz péptido + ratios Skyline exportados.
    Reparte CRC vs control en 10+10 según orden lexicográfico del prefijo de run
    (cohorte bioRxiv n=20); documentar en tesis y sustituir si el suplemento da otro mapeo.
    """
    df = pd.read_csv(matrix_path, sep="\t", low_memory=False)
    ratio_cols = [c for c in df.columns if isinstance(c, str) and c.endswith(" Ratio To Standard")]
    if not ratio_cols:
        raise ValueError("Plasma harmonize: no hay columnas terminadas en ' Ratio To Standard'")

    if "Protein Name" in df.columns and "Peptide Sequence" in df.columns:
        df = df.copy()
        df = df.dropna(subset=["Protein Name", "Peptide Sequence"], how="any")
        base = df["Protein Name"].astype(str).str.strip() + "|" + df["Peptide Sequence"].astype(str).str.strip()
        df.insert(0, "protein_id", base + "|r" + pd.Series(range(len(df)), index=df.index, dtype="int64").astype(str))
        out_df = df[["protein_id"] + ratio_cols].copy()
        out_df.to_csv(matrix_path, sep="\t", index=False)
    elif "protein_id" not in df.columns:
        raise ValueError(
            "Plasma harmonize: faltan 'Protein Name'/'Peptide Sequence' o columna 'protein_id' derivada."
        )
    else:
        dfix = df.copy()
        dfix = dfix[dfix["protein_id"].notna() & (dfix["protein_id"].astype(str).str.strip() != "")]
        dfix = dfix[~dfix["protein_id"].astype(str).str.lower().isin({"nan", "none"})]
        if dfix["protein_id"].duplicated().any():
            dfix = dfix.copy()
            dfix["protein_id"] = dfix["protein_id"].astype(str) + "|r" + dfix.index.astype(str)
        dfix[["protein_id"] + ratio_cols].to_csv(matrix_path, sep="\t", index=False)

    prefixes = sorted({c[: -len(" Ratio To Standard")] for c in ratio_cols})
    if len(prefixes) != 20:
        # Aún escribimos metadata por columna; el split 10/10 solo si son 20 runs.
        crc_set: set[str] = set()
    else:
        crc_set = set(prefixes[:10])

    meta_rows: list[dict[str, str]] = []
    for col in ratio_cols:
        pref = col[: -len(" Ratio To Standard")]
        if len(prefixes) == 20:
            cond = "crc" if pref in crc_set else "healthy"
        else:
            cond = "unknown"
        meta_rows.append(
            {
                "sample_id": col,
                "condition": cond,
                "batch": "PXD036594_decile_split" if len(prefixes) == 20 else "PXD036594_needs_curation",
            }
        )

    meta = pd.DataFrame(meta_rows)
    metadata_out.parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(metadata_out, sep="\t", index=False)
    n_pep = len(pd.read_csv(matrix_path, sep="\t", low_memory=False))
    return {"n_peptides": n_pep, "n_ratio_columns": len(ratio_cols), "n_run_prefixes": len(prefixes)}
