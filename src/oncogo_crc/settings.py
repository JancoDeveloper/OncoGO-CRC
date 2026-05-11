from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _normalize_optional_dir(root: Path, value: str | Path | None) -> Path | None:
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    p = Path(s).expanduser()
    if not p.is_absolute():
        p = (root / p).resolve()
    else:
        p = p.resolve()
    return p


def _coalesce_data_workspace(
    root: Path,
    yaml_paths: dict[str, Any],
    *,
    cli_value: str | Path | None,
) -> Path | None:
    """Prioridad: CLI > ONCOGO_CRC_DATA_WORKSPACE > config paths.data_workspace."""
    for candidate in (
        cli_value,
        os.environ.get("ONCOGO_CRC_DATA_WORKSPACE"),
        yaml_paths.get("data_workspace"),
    ):
        out = _normalize_optional_dir(root, candidate)
        if out is not None:
            return out
    return None


def _is_workspace_relative(rel: str) -> bool:
    """Rutas pesadas (datasets, procesados, informes) bajo data/, results/ o reports/."""
    norm = rel.replace("\\", "/").lstrip("/")
    return norm.startswith("data/") or norm.startswith("results/") or norm.startswith("reports/")


@dataclass
class Settings:
    root: Path
    config: dict[str, Any]
    datasets: dict[str, Any]
    comparisons: dict[str, Any]
    demo: bool
    data_workspace: Path | None = None

    @property
    def paths(self) -> dict[str, Any]:
        return self.config.get("paths", {})

    def resolve(self, rel: str) -> Path:
        p = Path(rel)
        if p.is_absolute():
            return p.resolve()
        base = self.root
        if self.data_workspace is not None and _is_workspace_relative(rel):
            base = self.data_workspace
        return (base / rel).resolve()

    def relative_to_project(self, path: Path | str) -> str:
        """Ruta legible respecto al workspace de datos o a la raíz del repo (logs, manifiestos)."""
        resolved = Path(path).resolve()
        if self.data_workspace is not None:
            dw = self.data_workspace.resolve()
            try:
                return str(resolved.relative_to(dw))
            except ValueError:
                pass
        try:
            return str(resolved.relative_to(self.root.resolve()))
        except ValueError:
            return str(resolved)

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


def load_settings(
    root: Path | None = None,
    *,
    demo: bool = False,
    data_workspace: str | Path | None = None,
) -> Settings:
    base = (root or _repo_root()).resolve()
    cfg = yaml.safe_load((base / "config/config.yaml").read_text())
    dsets = yaml.safe_load((base / "config/datasets.yaml").read_text())
    comps = yaml.safe_load((base / "config/comparisons.yaml").read_text())
    demo_flag = bool(demo or cfg.get("demo", {}).get("enabled", False))
    yaml_paths = cfg.get("paths", {}) if isinstance(cfg.get("paths"), dict) else {}
    dw = _coalesce_data_workspace(base, yaml_paths, cli_value=data_workspace)
    return Settings(root=base, config=cfg, datasets=dsets, comparisons=comps, demo=demo_flag, data_workspace=dw)
