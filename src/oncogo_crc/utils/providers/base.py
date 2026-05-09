from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class ProviderFile:
    name: str
    url: str | None = None
    size: int | None = None
    file_type: str | None = None
    checksum: str | None = None
    source: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class ProviderDatasetBundle:
    dataset_name: str
    provider: str
    raw_dir: Path
    files: list[ProviderFile] = field(default_factory=list)
    metadata_records: list[dict[str, Any]] = field(default_factory=list)
    manifest: dict[str, Any] = field(default_factory=dict)
    status: str = "pending"
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)


class BaseProviderAdapter:
    provider_name = "base"

    def __init__(self, dataset_name: str, dataset_cfg: dict[str, Any], global_cfg: dict[str, Any]):
        self.dataset_name = dataset_name
        self.dataset_cfg = dataset_cfg
        self.global_cfg = global_cfg

    def materialize_dataset(self, output_dir: Path) -> ProviderDatasetBundle:  # pragma: no cover - interface
        raise NotImplementedError
