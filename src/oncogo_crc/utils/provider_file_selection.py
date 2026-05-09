from __future__ import annotations

from typing import Any


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

    def _best(kind_keywords: list[str]) -> dict[str, Any] | None:
        best = None
        best_score = float("-inf")
        for row in candidates:
            label = f"{row.get('name','')} {row.get('url','')}"
            score = score_file(
                label,
                prefer_extensions=prefer_extensions,
                include_keywords=include_keywords + kind_keywords,
                exclude_keywords=exclude_keywords,
            )
            # avoid tiny files when size known (<10KB)
            size = row.get("size")
            if isinstance(size, int) and size > 0 and size < 10_000:
                score -= 2.0
            if score > best_score:
                best_score = score
                best = row
        return best

    return {"matrix": _best(["matrix", "abundance", "protein", "intensity"]), "metadata": _best(["metadata", "sample", "clinical"])}
