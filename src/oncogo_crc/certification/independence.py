"""Nivel A del detector de independencia: colisión de ficheros entre accesiones.

Un checksum idéntico es redepósito literal. El mismo tamaño y extensión en un
fichero crudo es sospecha fuerte. El mismo nombre normalizado es sospecha débil.
La ausencia de colisión no demuestra que los pacientes sean distintos.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

RAW_EXTENSIONS = {".raw", ".d", ".wiff", ".wiff2", ".mzml", ".mzxml", ".baf", ".tdf", ".lcd"}

GENERIC_NAME_PATTERNS = [
    re.compile(p, re.IGNORECASE)
    for p in (
        r"^(readme|licen[cs]e|checksum|md5sums?|file[_-]?list)",
        r"^(parameters?|settings?|mqpar|sdrf|experimental[_-]?design)",
        r"^(protein|peptide)groups?\.txt$",
        r"^evidence\.txt$",
        r"^msms\.txt$",
    )
]


@dataclass(frozen=True)
class FileRecord:
    accession: str
    file_name: str
    file_size: int | None = None
    checksum: str | None = None

    @property
    def extension(self) -> str:
        name = self.file_name.lower()
        for ext in sorted(RAW_EXTENSIONS, key=len, reverse=True):
            if name.endswith(ext):
                return ext
        return Path(name).suffix.lower()

    @property
    def is_raw(self) -> bool:
        return self.extension in RAW_EXTENSIONS

    @property
    def normalized_name(self) -> str:
        stem = self.file_name
        for ext in sorted(RAW_EXTENSIONS, key=len, reverse=True):
            if stem.lower().endswith(ext):
                stem = stem[: -len(ext)]
                break
        else:
            stem = Path(stem).stem
        stem = stem.split("/")[-1].split("\\")[-1]
        stem = re.sub(r"^\d{6,8}[_-]", "", stem)
        return re.sub(r"[^0-9a-z]+", "", stem.lower())

    @property
    def is_generic(self) -> bool:
        return any(p.search(Path(self.file_name).name) for p in GENERIC_NAME_PATTERNS)


@dataclass(frozen=True)
class Collision:
    kind: str
    accession_a: str
    accession_b: str
    file_a: str
    file_b: str
    evidence: str
    severity: str


def detect_file_collisions(records: list[FileRecord]) -> list[Collision]:
    """Colisiones entre accesiones distintas. No compara un fichero consigo mismo."""
    by_checksum: dict[str, list[FileRecord]] = defaultdict(list)
    by_size_ext: dict[tuple[int, str], list[FileRecord]] = defaultdict(list)
    by_name: dict[str, list[FileRecord]] = defaultdict(list)

    for rec in records:
        if not rec.file_name or rec.is_generic:
            continue
        if rec.checksum:
            by_checksum[rec.checksum.strip().lower()].append(rec)
        if rec.is_raw and rec.file_size:
            by_size_ext[(int(rec.file_size), rec.extension)].append(rec)
        if len(rec.normalized_name) >= 6:
            by_name[rec.normalized_name].append(rec)

    found: list[Collision] = []

    def pairs(bucket: list[FileRecord]):
        for a, b in combinations(bucket, 2):
            if a.accession != b.accession:
                yield a, b

    for checksum, bucket in by_checksum.items():
        for a, b in pairs(bucket):
            found.append(
                Collision(
                    "checksum",
                    a.accession,
                    b.accession,
                    a.file_name,
                    b.file_name,
                    f"checksum={checksum}",
                    "prueba",
                )
            )
    for (size, ext), bucket in by_size_ext.items():
        for a, b in pairs(bucket):
            found.append(
                Collision(
                    "size_ext",
                    a.accession,
                    b.accession,
                    a.file_name,
                    b.file_name,
                    f"bytes={size} ext={ext}",
                    "fuerte",
                )
            )
    for name, bucket in by_name.items():
        for a, b in pairs(bucket):
            found.append(
                Collision(
                    "normalized_name",
                    a.accession,
                    b.accession,
                    a.file_name,
                    b.file_name,
                    f"nombre_normalizado={name}",
                    "debil",
                )
            )
    return found
