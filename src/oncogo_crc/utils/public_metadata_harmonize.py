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

from oncogo_crc.utils.clinical_metadata import (
    bin_age,
    harmonize_clinical_metadata,
    normalize_msi,
    normalize_sex,
)

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


def _first_clinical_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Match column con prioridad: exacto → tokenizado → contiene-todos-los-tokens.

    El Excel clínico de PDC para CPTAC varía ligeramente entre estudios y versiones
    (``S037-r2``, ``S037-r3``, ``CPTAC3``). Necesitamos un matcher tolerante a:

    - whitespace inconsistente (incluye non-breaking space U+00A0)
    - mayúsculas/minúsculas
    - guiones bajos vs espacios
    - paréntesis aclaratorios al final (p. ej. ``Diagnosis Age (Years)``)
    - sufijos descriptivos extra (p. ej. ``MSI Type per BAT26``)
    """
    def _norm(s: str) -> str:
        s = str(s).replace("\u00a0", " ").replace("_", " ").strip().lower()
        # quita paréntesis y su contenido al final, y comprime espacios
        s = re.sub(r"\([^)]*\)\s*$", "", s)
        s = re.sub(r"\s+", " ", s).strip()
        return s

    norm_to_orig: dict[str, str] = {_norm(c): str(c) for c in df.columns}

    # 1) Match exacto tras normalización
    for cand in candidates:
        k = _norm(cand)
        if k in norm_to_orig:
            return norm_to_orig[k]

    # 2) Match por tokens (cada token del candidato presente en la columna)
    for cand in candidates:
        toks = [t for t in _norm(cand).split() if t]
        if not toks:
            continue
        for norm_col, orig in norm_to_orig.items():
            if all(t in norm_col.split() for t in toks):
                return orig

    # 3) Fallback: la columna contiene todas las subcadenas
    for cand in candidates:
        toks = [t for t in _norm(cand).split() if t]
        if not toks:
            continue
        for norm_col, orig in norm_to_orig.items():
            if all(t in norm_col for t in toks):
                return orig

    return None


def _stage_to_bin(stage_val: Any, early_tokens: list[str], adv_tokens: list[str]) -> str:
    from oncogo_crc.utils.clinical_metadata import bin_stage_pathological

    return bin_stage_pathological(stage_val, {"early": early_tokens, "advanced": adv_tokens})


def harmonize_cptac_colon_pdc(
    *,
    matrix_path: Path,
    metadata_out: Path,
    pdc_study_submitter_id: str,
    clinical_cfg: dict[str, Any],
    token: str | None,
) -> dict[str, Any]:
    # Fallar antes de PDC si la matriz local no tiene columnas TMT (ahorra red y evita metadata vacía).
    mdf0 = pd.read_csv(matrix_path, sep="\t", nrows=0, low_memory=False)
    sample_cols0 = [c for c in mdf0.columns if c != "protein_id"]

    def _specimen_from_matrix_col_early(col: str) -> str | None:
        c = str(col)
        for suf in (" Unshared Log Ratio", " Log Ratio", "_Unshared_Log_Ratio", "_Log_Ratio", ".Unshared.Log.Ratio", ".Log.Ratio"):
            if c.endswith(suf):
                return c[: -len(suf)]
        return None

    tmt0 = [c for c in sample_cols0 if _specimen_from_matrix_col_early(c)]
    if not tmt0:
        rel = matrix_path
        try:
            rel = matrix_path.resolve()
        except OSError:
            pass
        raise ValueError(
            "CPTAC harmonize: la matriz no tiene columnas de muestra TMT reconocibles "
            "(sufijos `… Log Ratio` / `… Unshared Log Ratio`). "
            f"Columnas de muestra halladas: {len(sample_cols0)}. "
            "Si solo existe `protein_id`, reemplaza el archivo por la matriz PDC completa o configura "
            f"la descarga correcta. Archivo: `{rel}`."
        )

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
    # Detectar columna de stage (puede ser "Tumor Stage (Pathological)", "Pathological Stage", etc.)
    stage_col = _first_clinical_col(
        clin_df,
        [
            "Tumor Stage (Pathological)",
            "Pathological Stage",
            "tumor_stage_pathological",
            "ajcc_pathologic_tumor_stage",
            "Tumor Stage Pathological",
            "Tumor Stage",
        ],
    )
    if stage_col:
        clin_df["_stage_raw"] = clin_df[stage_col].astype(str)
        clin_df["_stage_bin"] = clin_df[stage_col].map(
            lambda x: _stage_to_bin(x, early_tokens, adv_tokens)
        )
    else:
        clin_df["_stage_raw"] = ""
        clin_df["_stage_bin"] = ""

    age_col = _first_clinical_col(
        clin_df,
        [
            "Diagnosis Age",
            "Age at Initial Pathologic Diagnosis",
            "Age At Diagnosis",
            "Age at Diagnosis",
            "Age (Years)",
            "age_at_diagnosis",
            "Age at Index",
            "Patient Age",
        ],
    )
    # Evitar "Age in Months…" como edad en años (EOCRC / age_bin).
    if age_col:
        ac_norm = str(age_col).replace("\u00a0", " ").strip().lower()
        if "month" in ac_norm and "year" not in ac_norm:
            age_col = _first_clinical_col(
                clin_df,
                [
                    "Age at Diagnosis (Years)",
                    "Age at Diagnosis in Years",
                    "Age at Initial Pathologic Diagnosis",
                    "Diagnosis Age",
                    "Age in Years",
                ],
            )
            if not age_col:
                print(
                    "  ⚠ La columna de edad detectada parece estar en **meses**; "
                    "no hay alternativa en años — age/age_bin pueden quedar vacíos."
                )
                age_col = None
    dtb_col = _first_clinical_col(clin_df, ["days_to_birth", "Days to Birth"])
    if age_col:
        clin_df["_age_raw"] = pd.to_numeric(clin_df[age_col], errors="coerce")
    elif dtb_col:
        clin_df["_age_raw"] = (-pd.to_numeric(clin_df[dtb_col], errors="coerce") / 365.25).round()
    else:
        clin_df["_age_raw"] = float("nan")

    msi_col = _first_clinical_col(
        clin_df,
        [
            # PDC / CPTAC2 colon prospective S037 usa "MSI Type" (no "MSI Status")
            "MSI Type",
            "MSI status",
            "MSI Status",
            "MSI Subtype",
            "MSI Classification",
            "Microsatellite Instability Status",
            "Microsatellite Instability",
            "MSS/MSI",
            "MSI per BAT26",
            "MSI",
        ],
    )
    if msi_col:
        clin_df["_msi_raw"] = clin_df[msi_col].astype(str)
    else:
        clin_df["_msi_raw"] = ""

    sex_col = _first_clinical_col(clin_df, ["Gender", "Sex"])
    if sex_col:
        clin_df["_sex_raw"] = clin_df[sex_col].astype(str)
    else:
        clin_df["_sex_raw"] = ""

    site_col = _first_clinical_col(
        clin_df,
        [
            "Tumor Primary Site",
            "Tumor Anatomic Site",
            "Anatomical Location",
            "Primary Site",
            "Neoplasm Disease Site",
            "Site of Resection or Biopsy",
        ],
    )
    if site_col:
        clin_df["_tumor_site"] = clin_df[site_col].astype(str)
    else:
        clin_df["_tumor_site"] = ""

    ethnicity_col = _first_clinical_col(
        clin_df,
        [
            "Ethnicity",
            "ethnicity",
            "Ethnic Group",
            "Ethnic Background",
            "demographics_ethnicity",
        ],
    )
    if ethnicity_col:
        clin_df["_ethnicity_raw"] = clin_df[ethnicity_col].astype(str)
    else:
        clin_df["_ethnicity_raw"] = ""

    race_col = _first_clinical_col(
        clin_df,
        ["Race", "race", "Race Category", "demographics_race"],
    )
    if race_col:
        clin_df["_race_raw"] = clin_df[race_col].astype(str)
    else:
        clin_df["_race_raw"] = ""

    country_col = _first_clinical_col(
        clin_df,
        [
            "Country",
            "Country of Origin",
            "Country of Residence",
            "country",
            "Country at Enrollment",
            "Nation",
        ],
    )
    if country_col:
        clin_df["_country_raw"] = clin_df[country_col].astype(str)
    else:
        clin_df["_country_raw"] = ""

    # No incluir "Site" suelto: coincide con "Primary Site of Disease" por token "site" y rompe el merge.
    enroll_col = _first_clinical_col(
        clin_df,
        [
            "Enrollment Site",
            "Treating Institution",
            "Institution",
            "Study Site",
            "Biospecimen Source Site",
            "Site of Enrollment",
        ],
    )
    if enroll_col:
        clin_df["_enrollment_site"] = clin_df[enroll_col].astype(str)
    else:
        clin_df["_enrollment_site"] = ""

    # Log explícito de qué columnas detectó/no detectó para trazabilidad del run.
    print(
        "[harmonize_public_metadata] CPTAC clinical columns detected: "
        f"stage={stage_col!r} age={age_col or dtb_col!r} msi={msi_col!r} "
        f"sex={sex_col!r} site={site_col!r} ethnicity={ethnicity_col!r} "
        f"race={race_col!r} country={country_col!r} enrollment={enroll_col!r}"
    )
    if not stage_col:
        print("  ⚠ stage no detectado — revisa columnas del Excel clínico.")
    if not (age_col or dtb_col):
        print("  ⚠ age no detectado — comparación EOCRC/LORC quedará omitida.")
    if not msi_col:
        print("  ⚠ MSI no detectado — comparación MSI vs MSS quedará omitida.")
    if not ethnicity_col:
        print(
            "  ⚠ ethnicity no detectado — estrato hispano/latino y comparación opcional quedarán en unknown."
        )
    if not country_col and not enroll_col:
        print(
            "  ⚠ país/sitio de enrolamiento no detectados — geographic_region_bin quedará unknown "
            "(típico en CPTAC US si el Excel no trae país)."
        )

    clin_agg = clin_df.sort_values("Participant ID").groupby("Participant ID", as_index=False).first()
    clin_small = clin_agg[
        [
            "Participant ID",
            "_stage_raw",
            "_stage_bin",
            "_age_raw",
            "_msi_raw",
            "_sex_raw",
            "_tumor_site",
            "_ethnicity_raw",
            "_race_raw",
            "_country_raw",
            "_enrollment_site",
        ]
    ].rename(
        columns={
            "Participant ID": "participant_id",
            "_stage_raw": "stage",
            "_stage_bin": "stage_bin",
            "_age_raw": "age",
            "_msi_raw": "msi_status",
            "_sex_raw": "gender",
            "_tumor_site": "tumor_site",
            "_ethnicity_raw": "ethnicity",
            "_race_raw": "race",
            "_country_raw": "country",
            "_enrollment_site": "enrollment_site",
        }
    )

    mdf = pd.read_csv(matrix_path, sep="\t", nrows=0, low_memory=False)
    sample_cols = [c for c in mdf.columns if c != "protein_id"]

    def _specimen_from_matrix_col(col: str) -> str | None:
        c = str(col)
        # Orden: el sufijo más largo primero (evita que "... Unshared Log Ratio" coincida con "... Log Ratio").
        for suf in (" Unshared Log Ratio", " Log Ratio"):
            if c.endswith(suf):
                return c[: -len(suf)]
        # Algunas exportaciones PDC usan sufijos compactos o guiones bajos.
        for suf in ("_Unshared_Log_Ratio", "_Log_Ratio", ".Unshared.Log.Ratio", ".Log.Ratio"):
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
        st_raw = ""
        st_bin = ""
        age_v: Any = float("nan")
        msi_v, gend, site_v = "", "", ""
        eth_v, race_v, country_v, enroll_v = "", "", "", ""
        mm = long_map[long_map["specimen_label"] == slab]
        map_pid = str(mm.iloc[0]["participant_id"]) if not mm.empty else pid
        lid = str(map_pid or pid).strip()
        if lid:
            cmatch = clin_small[clin_small["participant_id"] == lid]
            if not cmatch.empty:
                st_raw = str(cmatch.iloc[0]["stage"] or "")
                st_bin = str(cmatch.iloc[0]["stage_bin"] or "")
                age_v = cmatch.iloc[0]["age"]
                msi_v = str(cmatch.iloc[0]["msi_status"] or "")
                gend = str(cmatch.iloc[0]["gender"] or "")
                site_v = str(cmatch.iloc[0]["tumor_site"] or "")
                eth_v = str(cmatch.iloc[0].get("ethnicity", "") or "")
                race_v = str(cmatch.iloc[0].get("race", "") or "")
                country_v = str(cmatch.iloc[0].get("country", "") or "")
                enroll_v = str(cmatch.iloc[0].get("enrollment_site", "") or "")
        out_rows.append(
            {
                "sample_id": col,
                "condition": cond,
                "participant_id": map_pid or pid,
                "stage": st_raw,
                "stage_bin": st_bin,
                "age": age_v,
                "gender": gend,
                "msi_status": msi_v,
                "tumor_site": site_v,
                "ethnicity": eth_v,
                "race": race_v,
                "country": country_v,
                "enrollment_site": enroll_v,
                "batch": "pdc_tmt10",
            }
        )

    out = pd.DataFrame(out_rows)
    global_cfg = {"clinical_filters": clinical_cfg} if isinstance(clinical_cfg, dict) else {}
    out = harmonize_clinical_metadata(out, global_cfg)
    bins_cfg = clinical_cfg.get("age_bins", {}).get("bins", []) if isinstance(clinical_cfg, dict) else []
    if "sex" not in out.columns:
        out["sex"] = out["gender"].map(normalize_sex) if "gender" in out.columns else "unknown"
    if "msi_status" in out.columns:
        out["msi_status"] = out["msi_status"].map(normalize_msi)
    out["age_bin"] = out["age"].map(lambda x: bin_age(x, bins_cfg))
    sup = ""
    if isinstance(clinical_cfg, dict):
        sup = str(clinical_cfg.get("cms_subtype_supplement_tsv") or "").strip()
    if sup:
        sup_path = Path(sup)
        if sup_path.is_file():
            cms = pd.read_csv(sup_path, sep="\t")
            if "participant_id" in cms.columns and "cms_subtype" in cms.columns:
                out = out.merge(
                    cms[["participant_id", "cms_subtype"]].drop_duplicates("participant_id"),
                    on="participant_id",
                    how="left",
                )
    # Columnas estándar para validación / comparaciones opcionales (EOCRC, MSI, CMS).
    if "cms_subtype" in out.columns:
        def _cms_bin(v: Any) -> str:
            if v is None or (isinstance(v, float) and pd.isna(v)):
                return ""
            s = str(v).strip().upper()
            if s in {"CMS4", "4"}:
                return "CMS4"
            if s in {"CMS1", "CMS2", "CMS3", "1", "2", "3"}:
                return "nonCMS4"
            if s.startswith("CMS"):
                return "nonCMS4"
            return ""

        out["cms_subtype_bin"] = out["cms_subtype"].map(_cms_bin)
    else:
        out["cms_subtype_bin"] = ""
    for col in ("metastatic_status", "treatment_response", "stage"):
        if col not in out.columns:
            out[col] = "unknown"
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
    lfq_re = re.compile(r"^LFQ intensity CRC\d+[MNP]$", re.IGNORECASE)
    # Matrices ya anchas (demo, post-proceso download, o export interno) PRIDE-T-* / PRIDE-N-*.
    pride_wide_re = re.compile(r"^PRIDE-[TN]-\d+$", re.IGNORECASE)

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
            wide_cols = [c for c in df.columns if isinstance(c, str) and pride_wide_re.match(c)]
            if wide_cols:
                print(
                    "[harmonize_public_metadata] PRIDE: matriz ya ancha (PRIDE-T-/PRIDE-N-); "
                    "se genera metadata sin requerir MaxQuant LFQ."
                )
                lfq_cols = wide_cols
            else:
                raise ValueError(
                    "PRIDE harmonize: la matriz tiene `protein_id` pero no columnas "
                    "`LFQ intensity CRC##[MNP]` (MaxQuant) ni `PRIDE-T-##` / `PRIDE-N-##`. "
                    "Borra `data/raw/pride/pride_crc_matrix.tsv` y ejecuta de nuevo "
                    "`oncogo-crc run --stage download_real` para PXD034575, o restaura el TSV "
                    "proteinGroups.txt original."
                )
        if not lfq_re.match(lfq_cols[0]):
            # Ya está en formato final; no reescribir matriz.
            out_df = df[["protein_id"] + lfq_cols].copy()
        else:
            out_df = df[["protein_id"] + lfq_cols].copy()
            out_df.to_csv(matrix_path, sep="\t", index=False)
    else:
        raise ValueError(
            f"PRIDE harmonize: primera columna debe ser 'Protein IDs' o 'protein_id', hallada {first_col!r}"
        )

    meta_rows: list[dict[str, str]] = []
    for col in lfq_cols:
        if pride_wide_re.match(col):
            cond = "tumor" if str(col).upper().startswith("PRIDE-T") else "normal"
        elif col.endswith("M"):
            cond = "tumor"
        elif col.endswith("N"):
            cond = "normal"
        elif col.endswith("P"):
            # Tanaka et al., Cell Reports 2024: P es el tumor primario (pCRC)
            # y M la metástasis hepática. v2 guardó P como «polyp». El cierre
            # de tesis no usa esta etiqueta; usa etiqueta_tanaka.
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
    Plasma PXD036594: matriz péptido + ratios Skyline de un panel dirigido.
    El orden alfabético de las corridas no es una etiqueta clínica. Antes se
    partía en 10 y 10 y se llamaba crc/healthy. Eso fabricaba un contraste.
    La condición queda en sin_etiqueta_clinica hasta que el suplemento del
    artículo asigne cada corrida.
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
    meta_rows: list[dict[str, str]] = []
    for col in ratio_cols:
        meta_rows.append(
            {
                "sample_id": col,
                "condition": "sin_etiqueta_clinica",
                "batch": "PXD036594_sin_etiqueta_clinica",
            }
        )

    meta = pd.DataFrame(meta_rows)
    metadata_out.parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(metadata_out, sep="\t", index=False)
    n_pep = len(pd.read_csv(matrix_path, sep="\t", low_memory=False))
    return {"n_peptides": n_pep, "n_ratio_columns": len(ratio_cols), "n_run_prefixes": len(prefixes)}
