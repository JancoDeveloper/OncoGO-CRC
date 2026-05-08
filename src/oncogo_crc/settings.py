from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


@dataclass
class Settings:
    root: Path
    config: dict[str, Any]
    datasets: dict[str, Any]
    comparisons: dict[str, Any]
    demo: bool

    @property
    def paths(self) -> dict[str, Any]:
        return self.config.get("paths", {})

    def resolve(self, rel: str) -> Path:
        return (self.root / rel).resolve()

    def protein_id_map_path(self) -> Path:
        if self.demo:
            p = self.resolve("data/interim/.demo_protein_mapping.tsv")
            if p.exists():
                return p
        rel = self.datasets.get("protein_id_map", "metadata/protein_id_maps/protein_mapping_master.tsv")
        return self.resolve(rel)

    def go_gene_long_path(self) -> Path:
        if self.demo:
            return self.resolve("data/interim/.demo_go_gene_long.tsv")
        return self.resolve("metadata/gene_ontology/go_gene_long.tsv")

    def sample_metadata_path(self, dataset_key: str) -> Path:
        """En modo demo, si existe un TSV generado en interim, tiene prioridad sobre metadata versionada."""
        interim = self.resolve(f"data/interim/_demo_{dataset_key}_samples.tsv")
        if self.demo and interim.exists():
            return interim
        ds = self.datasets["datasets"][dataset_key]
        return self.resolve(ds["sample_metadata"])


def load_settings(root: Path | None = None, *, demo: bool = False) -> Settings:
    base = root or _repo_root()
    cfg = yaml.safe_load((base / "config/config.yaml").read_text())
    dsets = yaml.safe_load((base / "config/datasets.yaml").read_text())
    comps = yaml.safe_load((base / "config/comparisons.yaml").read_text())
    demo_flag = bool(demo or cfg.get("demo", {}).get("enabled", False))
    return Settings(root=base.resolve(), config=cfg, datasets=dsets, comparisons=comps, demo=demo_flag)
