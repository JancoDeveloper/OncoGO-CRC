from __future__ import annotations

import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.clinical_metadata import harmonize_clinical_metadata
from oncogo_crc.utils.providers.iprox_adapter import IProXAdapter
from oncogo_crc.utils.providers.manual_adapter import ManualAdapter
from oncogo_crc.utils.providers.pdc_adapter import PDCAdapter
from oncogo_crc.utils.providers.pride_adapter import PRIDEAdapter
from oncogo_crc.utils.providers.px_adapter import PXAdapter


def _build_url(base_url: str, query_params: dict[str, Any] | None) -> str:
    if not query_params:
        return base_url
    q = urllib.parse.urlencode({k: v for k, v in query_params.items() if v is not None}, doseq=True)
    sep = "&" if "?" in base_url else "?"
    return f"{base_url}{sep}{q}"


def _resolve_env_placeholders(value: str) -> str:
    # Permite usar "${TOKEN}" dentro del YAML de headers/query.
    def repl(match: re.Match[str]) -> str:
        key = match.group(1)
        return os.environ.get(key, "")

    return re.sub(r"\$\{([A-Z0-9_]+)\}", repl, value)


def _resolve_env_recursive(obj: Any) -> Any:
    if isinstance(obj, str):
        return _resolve_env_placeholders(obj)
    if isinstance(obj, dict):
        return {k: _resolve_env_recursive(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_resolve_env_recursive(v) for v in obj]
    return obj


def _fetch_json(url: str, *, timeout_s: int, headers: dict[str, str] | None = None) -> Any:
    req = urllib.request.Request(url=url, headers=headers or {})
    with urllib.request.urlopen(req, timeout=timeout_s) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _fetch_text(url: str, *, timeout_s: int, headers: dict[str, str] | None = None) -> str:
    req = urllib.request.Request(url=url, headers=headers or {})
    with urllib.request.urlopen(req, timeout=timeout_s) as resp:
        return resp.read().decode("utf-8", errors="ignore")


def _extract_links_from_html(html: str, base_url: str) -> list[str]:
    hrefs = re.findall(r"""href=["']([^"']+)["']""", html, flags=re.IGNORECASE)
    out = []
    for h in hrefs:
        if h.startswith("javascript:") or h.startswith("#"):
            continue
        out.append(urllib.parse.urljoin(base_url, h))
    return out


def _score_url(url: str, *, kind: str) -> float:
    u = url.lower()
    score = 0.0
    if any(x in u for x in [".css", ".js", ".svg", ".png", ".jpg", ".jpeg", ".gif", ".html", ".htm"]):
        return -100.0
    if any(ext in u for ext in [".tsv", ".txt", ".csv"]):
        score += 4
    if ".xlsx" in u or ".xls" in u or ".zip" in u or ".gz" in u:
        score -= 2
    if kind == "matrix":
        for kw, w in [
            ("matrix", 5),
            ("abundance", 4),
            ("protein", 3),
            ("quant", 3),
            ("expression", 3),
        ]:
            if kw in u:
                score += w
    else:
        for kw, w in [
            ("metadata", 5),
            ("sample", 4),
            ("clinical", 3),
            ("annotation", 3),
            ("design", 2),
        ]:
            if kw in u:
                score += w
    return score


def _autodiscover_pride_urls(accession: str, *, timeout_s: int) -> tuple[str | None, str | None]:
    # API pública PRIDE v2: lista de archivos del proyecto.
    api = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}/files?page=0&pageSize=1000"
    data = _fetch_json(api, timeout_s=timeout_s)
    files: list[dict[str, Any]] = []
    if isinstance(data, dict):
        emb = data.get("_embedded")
        if isinstance(emb, dict):
            files = list(emb.get("files") or [])
    elif isinstance(data, list):
        files = data

    candidate_urls: list[str] = []
    for f in files:
        # Campos usuales en PRIDE responses:
        direct = f.get("downloadLink")
        if isinstance(direct, str) and direct:
            candidate_urls.append(direct)
        for loc in f.get("publicFileLocations", []) or []:
            if isinstance(loc, str):
                candidate_urls.append(loc)
            elif isinstance(loc, dict):
                # Algunos payloads usan value/name/url.
                for k in ("value", "url", "name"):
                    v = loc.get(k)
                    if isinstance(v, str) and v.startswith(("http://", "https://", "ftp://")):
                        candidate_urls.append(v)

    candidate_urls = list(dict.fromkeys(candidate_urls))
    if not candidate_urls:
        return None, None

    matrix = max(candidate_urls, key=lambda u: _score_url(u, kind="matrix"), default=None)
    metadata = max(candidate_urls, key=lambda u: _score_url(u, kind="metadata"), default=None)
    return matrix, metadata


def _autodiscover_from_catalog_page(catalog_url: str, *, timeout_s: int) -> tuple[str | None, str | None]:
    html = _fetch_text(catalog_url, timeout_s=timeout_s)
    links = _extract_links_from_html(html, catalog_url)
    if not links:
        return None, None
    matrix = max(links, key=lambda u: _score_url(u, kind="matrix"), default=None)
    metadata = max(links, key=lambda u: _score_url(u, kind="metadata"), default=None)
    return matrix, metadata


def _normalize_table_if_needed(out_path: Path) -> str:
    """
    Si el archivo descargado está en CSV, lo normaliza a TSV para el pipeline.
    Retorna detalle de normalización.
    """
    try:
        import pandas as pd

        # Heurística simple: extensión csv o cabecera con comas y sin tabs.
        is_csv = out_path.suffix.lower() == ".csv"
        if not is_csv:
            sample = out_path.read_text(encoding="utf-8", errors="ignore").splitlines()[:3]
            joined = "\n".join(sample)
            is_csv = ("," in joined) and ("\t" not in joined)
        if is_csv:
            df = pd.read_csv(out_path)
            df.to_csv(out_path, sep="\t", index=False)
            return "normalized_csv_to_tsv"
    except Exception:  # noqa: BLE001
        return "normalization_skipped"
    return "already_tsv_like"


def _download_with_retries(
    *,
    url: str,
    out_path: Path,
    timeout_s: int,
    retries: int,
    headers: dict[str, str] | None,
    overwrite: bool,
) -> tuple[bool, str]:
    if out_path.exists() and not overwrite:
        return True, "already_exists"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url=url, headers=headers or {})
    last_err = ""
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout_s) as resp:
                payload = resp.read()
            out_path.write_bytes(payload)
            norm = _normalize_table_if_needed(out_path)
            return True, f"downloaded_attempt_{attempt};{norm}"
        except Exception as exc:  # noqa: BLE001
            last_err = str(exc)
            if attempt < retries:
                time.sleep(min(2**attempt, 8))
    return False, last_err


PROVIDER_REGISTRY = {
    "pdc": PDCAdapter,
    "pride": PRIDEAdapter,
    "px": PXAdapter,
    "iprox": IProXAdapter,
    "manual": ManualAdapter,
}
PROVIDER_ALIASES = {
    "publication_supplement": "px",
    "metastatic_repository": "px",
    "crcdb": "manual",
    "pharmacoproteomics": "manual",
}


def _run_provider_inventory(
    settings: Settings,
    dataset_name: str,
    ds_cfg: dict[str, Any],
    global_cfg: dict[str, Any],
) -> dict[str, Any]:
    provider_name = str(ds_cfg.get("provider") or ds_cfg.get("download_real", {}).get("provider") or "").strip().lower()
    provider_name = PROVIDER_ALIASES.get(provider_name, provider_name)
    if provider_name not in PROVIDER_REGISTRY:
        return {"status": "skipped_provider_unknown", "provider": provider_name}
    adapter_cls = PROVIDER_REGISTRY[provider_name]
    output_dir = settings.resolve(ds_cfg["matrix_path"]).parent
    try:
        adapter = adapter_cls(dataset_name, ds_cfg, global_cfg)
        bundle = adapter.materialize_dataset(output_dir)
        if bundle.metadata_records:
            md = pd.DataFrame(bundle.metadata_records)
            if "sample_id" in md.columns and md["sample_id"].notna().any():
                md = harmonize_clinical_metadata(md, settings.config)
                meta_path = settings.resolve(ds_cfg["sample_metadata"])
                meta_path.parent.mkdir(parents=True, exist_ok=True)
                md.to_csv(meta_path, sep="\t", index=False)
        return {
            "status": bundle.status,
            "provider": bundle.provider,
            "n_files": len(bundle.files),
            "n_metadata_records": len(bundle.metadata_records),
            "warnings": bundle.warnings,
            "errors": bundle.errors,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "status": "failed",
            "provider": provider_name,
            "n_files": 0,
            "n_metadata_records": 0,
            "warnings": [],
            "errors": [str(exc)],
        }


def run(settings: Settings) -> None:
    """
    Descarga matrices/metadatos reales desde URLs o APIs definidas en datasets.yaml.

    Convención esperada por dataset (opcional):
      download_real:
        matrix_url: "https://.../matrix.tsv"
        metadata_url: "https://.../metadata.tsv"
        headers: {"Authorization": "Bearer ..."}
        query_params: {"study": "xyz"}
        overwrite: false
        timeout_s: 60
        retries: 3
    """
    if settings.demo:
        print("[download_real] modo demo activo: omitiendo descarga real.")
        return

    dl_cfg = settings.config.get("download_real", {})
    enabled = bool(dl_cfg.get("enabled", False))
    if not enabled:
        print("[download_real] deshabilitado en config.yaml (download_real.enabled=false).")
        return

    default_timeout = int(dl_cfg.get("timeout_s", 60))
    default_retries = int(dl_cfg.get("retries", 3))
    default_overwrite = bool(dl_cfg.get("overwrite", False))

    rows: list[dict[str, Any]] = []
    summary: dict[str, Any] = {"datasets": {}}
    manifests_dir = settings.resolve(dl_cfg.get("manifests_dir", "data/manifests"))
    manifests_dir.mkdir(parents=True, exist_ok=True)
    for key, ds in settings.datasets.get("datasets", {}).items():
        if ds.get("enabled") is False:
            rows.append(
                {
                    "dataset": key,
                    "matrix_path": ds["matrix_path"],
                    "metadata_path": ds["sample_metadata"],
                    "matrix_status": "skipped_dataset_disabled",
                    "metadata_status": "skipped_dataset_disabled",
                }
            )
            continue
        spec = ds.get("download_real", {}) or {}
        if not bool(spec.get("enabled", True)):
            rows.append(
                {
                    "dataset": key,
                    "matrix_path": ds["matrix_path"],
                    "metadata_path": ds["sample_metadata"],
                    "matrix_status": "skipped_disabled",
                    "metadata_status": "skipped_disabled",
                }
            )
            print(f"[download_real] {key}: skipped (download_real.enabled=false)")
            continue
        matrix_url = spec.get("matrix_url")
        metadata_url = spec.get("metadata_url")
        provider = str(spec.get("provider", "")).strip().lower()
        source_accession = str(spec.get("source_accession", "")).strip()
        source_catalog_url = str(spec.get("source_catalog_url", "")).strip()
        headers = _resolve_env_recursive(spec.get("headers") or {})
        qparams = _resolve_env_recursive(spec.get("query_params") or {})
        timeout_s = int(spec.get("timeout_s", default_timeout))
        retries = int(spec.get("retries", default_retries))
        overwrite = bool(spec.get("overwrite", default_overwrite))

        m_path = settings.resolve(ds["matrix_path"])
        md_path = settings.resolve(ds["sample_metadata"])

        row: dict[str, Any] = {
            "dataset": key,
            "matrix_path": ds["matrix_path"],
            "metadata_path": ds["sample_metadata"],
            "matrix_status": "skipped_no_url",
            "metadata_status": "skipped_no_url",
        }

        # Autodiscovery cuando no hay URL directa definida.
        if not matrix_url or not metadata_url:
            try:
                auto_matrix = auto_metadata = None
                if provider == "pride" and source_accession:
                    auto_matrix, auto_metadata = _autodiscover_pride_urls(source_accession, timeout_s=timeout_s)
                elif source_catalog_url:
                    auto_matrix, auto_metadata = _autodiscover_from_catalog_page(
                        source_catalog_url, timeout_s=timeout_s
                    )
                if not matrix_url and auto_matrix:
                    matrix_url = auto_matrix
                    row["matrix_autodiscovered"] = True
                if not metadata_url and auto_metadata:
                    metadata_url = auto_metadata
                    row["metadata_autodiscovered"] = True
            except Exception as exc:  # noqa: BLE001
                row["autodiscovery_error"] = str(exc)

        if matrix_url:
            url = _build_url(str(matrix_url), qparams)
            ok, status = _download_with_retries(
                url=url,
                out_path=m_path,
                timeout_s=timeout_s,
                retries=retries,
                headers=headers,
                overwrite=overwrite,
            )
            row["matrix_status"] = "ok" if ok else "error"
            row["matrix_detail"] = status
            row["matrix_url"] = url

        if metadata_url:
            url = _build_url(str(metadata_url), qparams)
            ok, status = _download_with_retries(
                url=url,
                out_path=md_path,
                timeout_s=timeout_s,
                retries=retries,
                headers=headers,
                overwrite=overwrite,
            )
            row["metadata_status"] = "ok" if ok else "error"
            row["metadata_detail"] = status
            row["metadata_url"] = url

        rows.append(row)
        provider_summary = _run_provider_inventory(settings, key, ds, settings.config)
        summary["datasets"][key] = {
            "provider": provider_summary.get("provider", spec.get("provider")),
            "status": provider_summary.get("status", "unknown"),
            "n_files": provider_summary.get("n_files", 0),
            "n_metadata_records": provider_summary.get("n_metadata_records", 0),
            "warnings": provider_summary.get("warnings", []),
            "errors": provider_summary.get("errors", []),
            "matrix_status": row.get("matrix_status"),
            "metadata_status": row.get("metadata_status"),
            "matrix_url": row.get("matrix_url"),
            "metadata_url": row.get("metadata_url"),
        }
        (manifests_dir / f"{key}.json").write_text(
            json.dumps(summary["datasets"][key], indent=2),
            encoding="utf-8",
        )
        print(
            f"[download_real] {key}: matrix={row['matrix_status']} metadata={row['metadata_status']}"
        )

    out = settings.resolve("data/raw/download_real_manifest.tsv")
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    settings.resolve("data/raw/download_real_manifest.json").write_text(
        json.dumps(rows, indent=2), encoding="utf-8"
    )
    (manifests_dir / "_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[download_real] manifest -> {out.relative_to(settings.root)}")
