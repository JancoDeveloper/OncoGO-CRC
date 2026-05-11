from __future__ import annotations

from typing import Any


def _metadata_row_is_archive_bundle(row: dict[str, Any]) -> bool:
    """
    Excluye paquetes de instrumento (.zip de PRM/Skyline, etc.) como metadata tabular.
    Permite .tsv.gz / .txt.gz / SDRF comprimido.
    """
    nm = (row.get("name") or "").lower()
    if nm.endswith((".zip", ".tar", ".tgz", ".7z")):
        return True
    if nm.endswith(".gz") and not any(
        x in nm for x in (".tsv.gz", ".txt.gz", ".csv.gz", "sdrf", "metadata", "annotation")
    ):
        return True
    return False


def score_file(name_or_url: str, *, prefer_extensions: list[str], include_keywords: list[str], exclude_keywords: list[str]) -> float:
    s = (name_or_url or "").lower()
    score = 0.0

    # Extension priority
    for i, ext in enumerate(prefer_extensions):
        if s.endswith(ext.lower()):
            score += float(len(prefer_extensions) - i) * 4.0
            break

    # Positive/negative keywords
    for kw in include_keywords:
        if kw.lower() in s:
            score += 2.0
    for kw in exclude_keywords:
        if kw.lower() in s:
            score -= 4.0

    # Penalize tiny/non-informative names
    if any(x in s for x in ["readme", "description", "note"]):
        score -= 3.0
    return score


def pick_best_files(candidates: list[dict[str, Any]], file_selection: dict[str, Any]) -> dict[str, dict[str, Any] | None]:
    prefer_extensions = list(file_selection.get("prefer_extensions", [".tsv", ".txt", ".csv", ".xlsx"]))
    include_keywords = list(file_selection.get("include_keywords", []))
    exclude_keywords = list(file_selection.get("exclude_keywords", []))

    def _best(kind_keywords: list[str], *, extra_exclude: list[str] | None = None) -> dict[str, Any] | None:
        best = None
        best_score = float("-inf")
        ex = list(exclude_keywords) + list(extra_exclude or [])
        for row in candidates:
            label = f"{row.get('name','')} {row.get('url','')}"
            score = score_file(
                label,
                prefer_extensions=prefer_extensions,
                include_keywords=include_keywords + kind_keywords,
                exclude_keywords=ex,
            )
            # avoid tiny files when size known (<10KB)
            size = row.get("size")
            if isinstance(size, int) and size > 0 and size < 10_000:
                score -= 2.0
            if score > best_score:
                best_score = score
                best = row
        return best

    matrix = _best(
        ["matrix", "abundance", "protein", "intensity", "quant", "tmt10", "proteome", "colon", "cptac"],
        extra_exclude=[
            "sample.txt",
            "proteome.sample",
            "clinical_data",
            "label_to_sample",
            "mapping_file",
            "biospecimen",
            "mzid",
            "pep.xml",
            "spectra",
            "psm",
        ],
    )
    def _meta_best() -> dict[str, Any] | None:
        best = None
        best_score = float("-inf")
        ex = list(exclude_keywords) + [
            "proteingroups",
            "protein_groups",
            "peptides.tsv",
            "evidence.txt",
            "qc_protein",
            "maxquant",
            "mztab",
            "proteome.sample",
            "label_to_sample",
            ".zip",
            ".prm",
            "skyline",
            "dia.zip",
        ]
        kind_keywords = ["metadata", "sample", "clinical", "sdrf", "experimental", "design", "annotation"]
        for row in candidates:
            if _metadata_row_is_archive_bundle(row):
                continue
            label = f"{row.get('name','')} {row.get('url','')}"
            score = score_file(
                label,
                prefer_extensions=prefer_extensions,
                include_keywords=include_keywords + kind_keywords,
                exclude_keywords=ex,
            )
            size = row.get("size")
            if isinstance(size, int) and size > 0 and size < 10_000:
                score -= 2.0
            if score > best_score:
                best_score = score
                best = row
        return best

    metadata = _meta_best()
    if matrix and metadata and matrix.get("url") == metadata.get("url"):
        # Segundo mejor archivo distinto de la matriz; no re-aplicar exclusiones agresivas (p. ej. proteingroups)
        # o no quedaría ningún candidato en proyectos solo-MaxQuant.
        alt_candidates = [
            c
            for c in candidates
            if c.get("url") != matrix.get("url")
            and "qc" not in (c.get("name") or "").lower()
            and not _metadata_row_is_archive_bundle(c)
        ]
        metadata = _best_alt(
            ["metadata", "sample", "clinical", "sdrf", "experimental", "design", "annotation"],
            alt_candidates,
            prefer_extensions=prefer_extensions,
            base_exclude=list(exclude_keywords) + [".zip", ".prm", "skyline", "dia.zip"],
        )

    return {"matrix": matrix, "metadata": metadata}


def _best_alt(
    kind_keywords: list[str],
    candidates: list[dict[str, Any]],
    *,
    prefer_extensions: list[str],
    base_exclude: list[str],
) -> dict[str, Any] | None:
    if not candidates:
        return None
    best = None
    best_score = float("-inf")
    for row in candidates:
        if _metadata_row_is_archive_bundle(row):
            continue
        label = f"{row.get('name','')} {row.get('url','')}"
        score = score_file(
            label,
            prefer_extensions=prefer_extensions,
            include_keywords=kind_keywords,
            exclude_keywords=base_exclude,
        )
        size = row.get("size")
        if isinstance(size, int) and size > 0 and size < 10_000:
            score -= 2.0
        if score > best_score:
            best_score = score
            best = row
    return best
