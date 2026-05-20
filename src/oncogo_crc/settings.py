from __future__ import annotations

import os
import shutil
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


def _copy_metadata_file_if_needed(src: Path, dst: Path, *, refresh_if_repo_newer: bool) -> str:
    """Copia un archivo de metadata al workspace. Devuelve ``new``, ``updated`` o vacío.

    Si el destino no es escribible (volumen externo no montado, sin permisos), retorna sin error.
    """
    if not src.is_file():
        return ""
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if not dst.exists():
            shutil.copy2(src, dst)
            return "new"
        if refresh_if_repo_newer:
            if src.stat().st_mtime > dst.stat().st_mtime:
                shutil.copy2(src, dst)
                return "updated"
    except (OSError, PermissionError):
        return ""
    return ""


def _workspace_writable(workspace: Path) -> bool:
    """¿Podemos escribir en el workspace ahora mismo? (volumen montado y con permisos)."""
    try:
        workspace.mkdir(parents=True, exist_ok=True)
        probe = workspace / ".oncogo_write_probe"
        probe.write_text("ok", encoding="utf-8")
        probe.unlink(missing_ok=True)
        return True
    except (OSError, PermissionError):
        return False


def _seed_workspace_metadata_from_repo(repo_root: Path, workspace: Path) -> None:
    """Prepara ``<workspace>/metadata/`` para una corrida completa.

    - **``protein_id_maps``**, **``gene_ontology``** y **``tcga_auxiliary``** (plantillas): se copian si faltan
      o si la copia del **repositorio es más reciente** (p. ej. tras ``git pull`` o regenerar el bulk UniProt).
    - **``sample_metadata``** y **``clinical_layers``**: solo si **faltan** en el workspace (no se pisan
      curaciones ni harmonizaciones ya escritas allí).
    """
    src = repo_root / "metadata"
    dst = workspace / "metadata"
    if not src.is_dir():
        return
    n_new = 0
    n_updated = 0
    refresh_subs = {"protein_id_maps", "gene_ontology", "tcga_auxiliary"}
    for sub in refresh_subs:
        sdir = src / sub
        if not sdir.is_dir():
            continue
        for p in sdir.rglob("*"):
            if not p.is_file():
                continue
            rel = p.relative_to(sdir)
            out = dst / sub / rel
            tag = _copy_metadata_file_if_needed(p, out, refresh_if_repo_newer=True)
            if tag == "new":
                n_new += 1
            elif tag == "updated":
                n_updated += 1
    for sub in ("sample_metadata", "clinical_layers"):
        sdir = src / sub
        if not sdir.is_dir():
            continue
        for p in sdir.rglob("*"):
            if not p.is_file():
                continue
            rel = p.relative_to(sdir)
            out = dst / sub / rel
            tag = _copy_metadata_file_if_needed(p, out, refresh_if_repo_newer=False)
            if tag == "new":
                n_new += 1
    if n_new or n_updated:
        parts = [f"{n_new} nuevo(s)"] if n_new else []
        if n_updated:
            parts.append(f"{n_updated} actualizado(s) desde repo (mtime)")
        print(f"[settings] metadata → workspace ({dst}): {', '.join(parts)}.")


def _is_workspace_relative(rel: str) -> bool:
    """Rutas pesadas y metadata: data/, results/, reports/, metadata/, y docs/web (MkDocs)."""
    norm = rel.replace("\\", "/").lstrip("/")
    return (
        norm.startswith("data/")
        or norm.startswith("results/")
        or norm.startswith("reports/")
        or norm.startswith("metadata/")
        or norm.startswith("docs/web")
    )


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

    def processed_matrix_path(self, dataset_key: str, kind: str) -> Path:
        """``kind`` ∈ ``log2norm`` (diferencial, GO, GSEA, overlap) | ``zscore`` (ML, red, clustering)."""
        proc = self.resolve(self.config["paths"]["processed_dir"])
        if kind == "log2norm":
            return proc / f"{dataset_key}_log2norm.tsv"
        if kind == "zscore":
            return proc / f"{dataset_key}_zscore.tsv"
        raise ValueError(f"processed_matrix_path: kind desconocido {kind!r} (use log2norm o zscore)")


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
    if dw is not None and demo_flag:
        # Salvaguarda crítica: el modo demo escribe matrices sintéticas (GENE001…) y sobrescribiría
        # resultados reales acumulados en el workspace. Forzar fallback al repo evita pisar el sitio
        # publicado y los TSV de corridas reales.
        print(
            f"[settings] data_workspace={dw} ignorado porque --demo está activo. "
            "El modo demo siempre usa rutas relativas al repositorio para no contaminar el workspace "
            "con artefactos sintéticos."
        )
        dw = None
    if dw is not None:
        if _workspace_writable(dw):
            _seed_workspace_metadata_from_repo(base, dw)
            print(
                f"[settings] data_workspace={dw} — descargas y artefactos pesados bajo "
                "`data/`, `results/`, `reports/`, `metadata/`, `docs/web` (config y código siguen en el repo)."
            )
        else:
            print(
                f"[settings] data_workspace={dw} configurado pero NO accesible (volumen sin montar o sin permisos). "
                "Fallback: rutas relativas al repositorio. Conecta el volumen o ajusta `paths.data_workspace` para usarlo."
            )
            dw = None
    return Settings(root=base, config=cfg, datasets=dsets, comparisons=comps, demo=demo_flag, data_workspace=dw)
