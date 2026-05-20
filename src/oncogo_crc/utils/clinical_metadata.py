from __future__ import annotations

import re
from typing import Any

import pandas as pd

SEX_MAP = {
    "male": "male",
    "m": "male",
    "female": "female",
    "f": "female",
    "unknown": "unknown",
    "not reported": "unknown",
    "na": "unknown",
    "none": "unknown",
}


def normalize_sex(value: Any) -> str:
    if value is None or pd.isna(value):
        return "unknown"
    return SEX_MAP.get(str(value).strip().lower(), "unknown")


def normalize_age(value: Any) -> int | None:
    if value is None or pd.isna(value):
        return None
    try:
        return int(float(value))
    except Exception:  # noqa: BLE001
        return None


def bin_age(age: Any, bins_cfg: list[dict[str, Any]] | None = None) -> str:
    a = normalize_age(age)
    if a is None:
        return "unknown"
    if bins_cfg:
        for item in bins_cfg:
            if a >= int(item.get("min", -999)) and a <= int(item.get("max", 999)):
                return str(item.get("name", "unknown"))
    return "older" if a >= 50 else "young"


def normalize_stage(value: Any) -> str:
    if value is None or pd.isna(value):
        return "unknown"
    return str(value).strip().lower()


def bin_stage(stage: Any, stage_mapping: dict[str, list[str]] | None) -> str:
    s = normalize_stage(stage)
    if s == "unknown":
        return "unknown"
    for label, values in (stage_mapping or {}).items():
        vals = {str(v).strip().lower() for v in values}
        if s in vals:
            return str(label)
    return "unknown"


def bin_stage_pathological(stage_val: Any, stage_mapping: dict[str, list[str]] | None) -> str:
    """Temprano vs avanzado por **subcadena** en el texto de estadio (compatible con PDC/TCGA)."""
    if stage_val is None or pd.isna(stage_val):
        return ""
    sm = stage_mapping or {}
    early_tokens = list(sm.get("early", []) or [])
    adv_tokens = list(sm.get("advanced", []) or [])
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


def normalize_ethnicity(value: Any) -> str:
    """Texto crudo normalizado (minúsculas) para etnicidad NIH/TCGA-style."""
    if value is None or pd.isna(value):
        return "unknown"
    s = str(value).strip().lower()
    if not s or s in {"unknown", "not reported", "na", "n/a", "none", "not available"}:
        return "unknown"
    return s


def bin_hispanic_latino(ethnicity: Any) -> str:
    """Estrato hispano/latino según categoría NIH (autorreportada en cohortes US).

    Valores: ``hispanic_latino`` | ``non_hispanic_latino`` | ``unknown``.

    **Limitación:** en CPTAC/PDC esto NO equivale a «paciente colombiano» ni a origen geográfico
    en Latinoamérica; es la variable demográfica estándar de estudios oncológicos en EE.UU.
    """
    s = normalize_ethnicity(ethnicity)
    if s == "unknown":
        return "unknown"
    if any(
        tok in s
        for tok in (
            "hispanic",
            "latino",
            "latina",
            "latin american",
            "latin american or hispanic",
        )
    ):
        if "not" in s and "hispanic" in s:
            return "non_hispanic_latino"
        return "hispanic_latino"
    if any(tok in s for tok in ("not hispanic", "non-hispanic", "non hispanic", "not latino")):
        return "non_hispanic_latino"
    return "unknown"


def normalize_race(value: Any) -> str:
    if value is None or pd.isna(value):
        return "unknown"
    s = str(value).strip().lower()
    if not s or s in {"unknown", "not reported", "na", "n/a", "none"}:
        return "unknown"
    return s


def normalize_country(value: Any) -> str:
    if value is None or pd.isna(value):
        return "unknown"
    s = str(value).strip()
    if not s or s.lower() in {"unknown", "not reported", "na", "n/a", "none"}:
        return "unknown"
    return s


# Países considerados América Latina / Caribe (ISO-style names en inglés, como en PDC/TCGA).
_LATIN_AMERICA_COUNTRIES = frozenset(
    {
        "argentina",
        "bolivia",
        "brazil",
        "chile",
        "colombia",
        "costa rica",
        "cuba",
        "dominican republic",
        "ecuador",
        "el salvador",
        "guatemala",
        "honduras",
        "mexico",
        "nicaragua",
        "panama",
        "paraguay",
        "peru",
        "puerto rico",
        "uruguay",
        "venezuela",
        "latin america",
    }
)

_US_COUNTRY_TOKENS = frozenset(
    {
        "united states",
        "united states of america",
        "usa",
        "u.s.",
        "u.s.a.",
        "us",
    }
)


def bin_geographic_region(country: Any, *, enrollment_site: Any = None) -> str:
    """Región geográfica aproximada a partir de país (y opcionalmente sitio de reclutamiento).

    Valores: ``latin_america`` | ``united_states`` | ``other`` | ``unknown``.

    CPTAC colon prospective recluta casi exclusivamente en centros de EE.UU.; si no hay país,
    se infiere ``united_states`` solo cuando el sitio de enrolamiento contiene tokens US explícitos.
    """
    c = normalize_country(country)
    if c != "unknown":
        cl = c.lower()
        if cl in _US_COUNTRY_TOKENS:
            return "united_states"
        if cl in _LATIN_AMERICA_COUNTRIES:
            return "latin_america"
        return "other"

    site = "" if enrollment_site is None or pd.isna(enrollment_site) else str(enrollment_site).strip().lower()
    if site:
        if any(tok in site for tok in ("colombia", "mexico", "brazil", "argentina", "chile", "peru", "latin")):
            return "latin_america"
        if any(
            tok in site
            for tok in (
                "united states",
                "usa",
                "u.s.",
                "washington",
                "baltimore",
                "boston",
                "stanford",
                "vanderbilt",
                "mayo",
            )
        ):
            return "united_states"
    return "unknown"


_STAGE_GROUP_RE = re.compile(
    r"\bstage\s*(?P<g>[ivx]+|[1-4])\b",
    re.IGNORECASE,
)


def parse_stage_group(stage_val: Any) -> str:
    """Extrae grupo AJCC romano (I–IV) del texto de estadio patológico."""
    if stage_val is None or pd.isna(stage_val):
        return ""
    s = str(stage_val).strip()
    if not s:
        return ""
    m = _STAGE_GROUP_RE.search(s)
    if not m:
        # Patrones compactos: "IIIA", "IVB", "Stage 3"
        m2 = re.search(r"\b([iv]{1,3}[abc]?|[1-4][abc]?)\b", s, re.IGNORECASE)
        if m2:
            g = m2.group(1).upper()
            if g[0].isdigit():
                return {"1": "I", "2": "II", "3": "III", "4": "IV"}.get(g[0], "")
            return g[0] if len(g) == 1 else g.rstrip("ABC")
        return ""
    g = m.group("g").upper()
    if g.isdigit():
        return {"1": "I", "2": "II", "3": "III", "4": "IV"}.get(g, "")
    return g[0] if len(g) >= 1 and g[0] in "IVX" else ""


def normalize_msi(value: Any) -> str:
    if value is None or pd.isna(value):
        return "unknown"
    s = str(value).strip().upper()
    if s in {"MSI-H", "MSIH", "HIGH"}:
        return "MSI-H"
    if s in {"MSS", "STABLE"}:
        return "MSS"
    if s in {"MSI-L", "MSIL", "LOW"}:
        return "MSI-L"
    return "unknown"


def harmonize_clinical_metadata(df: pd.DataFrame, global_cfg: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    bins_cfg = global_cfg.get("clinical_filters", {}).get("age_bins", {}).get("bins", [])
    stage_map = global_cfg.get("clinical_filters", {}).get("stage_mapping", {})

    if "sex" not in out.columns and "gender" in out.columns:
        out["sex"] = out["gender"]
    if "sex" not in out.columns and "gender_raw" in out.columns:
        out["sex"] = out["gender_raw"]

    out["sex"] = out.get("sex", pd.Series(index=out.index, dtype="object")).apply(normalize_sex)
    out["age"] = out.get("age", pd.Series(index=out.index, dtype="object")).apply(normalize_age)
    out["age_bin"] = out["age"].apply(lambda x: bin_age(x, bins_cfg))
    out["stage"] = out.get("stage", pd.Series(index=out.index, dtype="object")).apply(normalize_stage)
    # Preferir stage_bin ya asignado en harmonize PDC; si no, inferir por subcadena AJCC.
    if "stage_bin" in out.columns:
        empty_bin = out["stage_bin"].isna() | (out["stage_bin"].astype(str).str.strip() == "")
        out.loc[empty_bin, "stage_bin"] = out.loc[empty_bin, "stage"].apply(
            lambda x: bin_stage_pathological(x, stage_map)
        )
    else:
        out["stage_bin"] = out["stage"].apply(lambda x: bin_stage_pathological(x, stage_map))
    if "stage_group" not in out.columns:
        out["stage_group"] = out["stage"].apply(parse_stage_group)
    else:
        out["stage_group"] = out["stage_group"].fillna(out["stage"].apply(parse_stage_group))

    if "ethnicity" in out.columns:
        out["ethnicity"] = out["ethnicity"].apply(normalize_ethnicity)
    if "race" in out.columns:
        out["race"] = out["race"].apply(normalize_race)
    if "country" in out.columns:
        out["country"] = out["country"].apply(normalize_country)
    if "ethnicity" in out.columns:
        out["hispanic_latino_bin"] = out["ethnicity"].apply(bin_hispanic_latino)
    elif "hispanic_latino_bin" not in out.columns:
        out["hispanic_latino_bin"] = "unknown"

    enroll_col = "enrollment_site" if "enrollment_site" in out.columns else None
    if "country" in out.columns or enroll_col:
        out["geographic_region_bin"] = out.apply(
            lambda r: bin_geographic_region(
                r.get("country"),
                enrollment_site=r.get(enroll_col) if enroll_col else None,
            ),
            axis=1,
        )
    elif "geographic_region_bin" not in out.columns:
        out["geographic_region_bin"] = "unknown"

    out["msi_status"] = out.get("msi_status", pd.Series(index=out.index, dtype="object")).apply(normalize_msi)
    if "metastatic_status" not in out.columns:
        out["metastatic_status"] = "unknown"
    if "treatment_response" not in out.columns:
        out["treatment_response"] = "unknown"
    if "condition" not in out.columns:
        out["condition"] = "unknown"

    ordered = [
        "sample_id",
        "condition",
        "participant_id",
        "sex",
        "age",
        "age_bin",
        "stage",
        "stage_bin",
        "stage_group",
        "msi_status",
        "ethnicity",
        "race",
        "hispanic_latino_bin",
        "country",
        "geographic_region_bin",
        "enrollment_site",
        "tumor_site",
        "cms_subtype",
        "cms_subtype_bin",
        "metastatic_status",
        "treatment_response",
    ]
    existing = [c for c in ordered if c in out.columns]
    rest = [c for c in out.columns if c not in existing]
    return out[existing + rest]
