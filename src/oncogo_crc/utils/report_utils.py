from __future__ import annotations

from pathlib import Path


def write_executive_summary(
    path: Path,
    *,
    sections: dict[str, str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    parts = ["# OncoGO-CRC — resumen ejecutivo\n"]
    for title, body in sections.items():
        parts.append(f"\n## {title}\n\n{body}\n")
    path.write_text("".join(parts), encoding="utf-8")
