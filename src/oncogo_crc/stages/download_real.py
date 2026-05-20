from __future__ import annotations

import hashlib
import json
import os
import re
import socket
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeoutError
from pathlib import Path
from typing import Any

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.clinical_metadata import harmonize_clinical_metadata
from oncogo_crc.utils.providers.iprox_adapter import IProXAdapter
from oncogo_crc.utils.providers.manual_adapter import ManualAdapter
from oncogo_crc.utils.providers.base import ProviderDatasetBundle
from oncogo_crc.utils.providers.pdc_adapter import PDCAdapter
from oncogo_crc.utils.providers.pride_adapter import PRIDEAdapter
from oncogo_crc.utils.providers.px_adapter import PXAdapter
from oncogo_crc.utils.provider_file_selection import pick_best_files


def _sha256_file(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _matrix_has_tmt_style_samples(path: Path, *, protein_col: str = "protein_id") -> bool:
    """Heurística: matriz CPTAC/PDC TMT10 con muchas columnas `… Log Ratio` / `… Unshared Log Ratio`.

    Evita que ``overwrite: false`` deje un **stub** (solo ``protein_id`` o pocas columnas) y
    ``download_real`` marque ``matrix=ok`` sin haber descargado nada útil.
    """
    if not path.exists() or path.stat().st_size < 10_000:
        return False
    try:
        hdr = pd.read_csv(path, sep="\t", nrows=0, low_memory=False)
    except Exception:  # noqa: BLE001
        return False
    cols = [
        c
        for c in hdr.columns
        if str(c) != protein_col and str(c).strip().lower() not in {"gene", "gene_name"}
    ]
    if len(cols) < 16:
        return False
    suf = (
        " Log Ratio",
        " Unshared Log Ratio",
        "_Log_Ratio",
        "_Unshared_Log_Ratio",
        ".Log.Ratio",
        ".Unshared.Log.Ratio",
    )
    hits = sum(1 for c in cols if str(c).endswith(suf))
    return hits >= 8


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


def _is_useful_data_url(url: str | None) -> bool:
    if not url:
        return False
    u = str(url).strip().lower()
    if not u.startswith(("http://", "https://", "ftp://")):
        return False
    bad_tokens = [
        ".css",
        ".js",
        ".svg",
        ".png",
        ".jpg",
        ".jpeg",
        ".gif",
        ".html",
        ".htm",
        "fonts.gstatic",
        "/cgi/",
        "view=datasets",
    ]
    if any(tok in u for tok in bad_tokens):
        return False
    parsed = urllib.parse.urlparse(u)
    path = parsed.path or ""
    if path in {"", "/"}:
        return False
    good_ext = [".tsv", ".txt", ".csv", ".xlsx", ".xls", ".tsv.gz", ".zip", ".gz", ".mztab"]
    if any(path.endswith(ext) for ext in good_ext):
        return True
    if "amazonaws.com" in (parsed.netloc or "").lower():
        return any(ext in u for ext in [".tsv", ".txt", ".csv", ".gz", ".zip"])
    return any(k in path for k in ["/download/", "/files/"])


def _file_row_download_url(f: dict[str, Any]) -> str | None:
    direct = f.get("downloadLink") or f.get("publicLink")
    if isinstance(direct, str) and direct.startswith(("http://", "https://", "ftp://")):
        return direct
    for loc in f.get("publicFileLocations", []) or []:
        if isinstance(loc, str) and loc.startswith(("http://", "https://", "ftp://")):
            return loc
        if isinstance(loc, dict):
            for k in ("value", "url", "name"):
                v = loc.get(k)
                if isinstance(v, str) and v.startswith(("http://", "https://", "ftp://")):
                    return v
    return None


def _autodiscover_pride_urls(
    accession: str,
    *,
    timeout_s: int,
    file_selection: dict[str, Any],
) -> tuple[str | None, str | None]:
    """Elige matriz y metadata (p. ej. SDRF) usando nombre de archivo + URL, no solo la URL."""
    page = 0
    page_size = 500
    files: list[dict[str, Any]] = []
    while True:
        api = (
            f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}/files"
            f"?page={page}&pageSize={page_size}"
        )
        data = _fetch_json(api, timeout_s=timeout_s)
        chunk: list[dict[str, Any]] = []
        if isinstance(data, dict):
            emb = data.get("_embedded")
            if isinstance(emb, dict):
                chunk = list(emb.get("files") or [])
        elif isinstance(data, list):
            chunk = data
        if not chunk:
            break
        files.extend(chunk)
        if len(chunk) < page_size:
            break
        page += 1

    candidates: list[dict[str, Any]] = []
    for f in files:
        name = str(f.get("fileName") or f.get("name") or "")
        url = _file_row_download_url(f)
        if not url or not _is_useful_data_url(url):
            continue
        size = f.get("fileSize") or f.get("fileSizeBytes")
        try:
            isize = int(size) if size is not None else None
        except (TypeError, ValueError):
            isize = None
        candidates.append(
            {
                "name": name,
                "url": url,
                "size": isize,
                "file_type": f.get("fileCategory"),
            }
        )

    if not candidates:
        return None, None

    selected = pick_best_files(candidates, file_selection)
    mat = selected.get("matrix") or {}
    meta = selected.get("metadata") or {}
    matrix_url = mat.get("url") if mat else None
    metadata_url = meta.get("url") if meta else None
    if metadata_url and matrix_url and metadata_url == matrix_url:
        metadata_url = None
    return matrix_url, metadata_url


def _autodiscover_from_catalog_page(catalog_url: str, *, timeout_s: int) -> tuple[str | None, str | None]:
    html = _fetch_text(catalog_url, timeout_s=timeout_s)
    links = _extract_links_from_html(html, catalog_url)
    if not links:
        return None, None
    filtered = [u for u in links if _is_useful_data_url(u)]
    if not filtered:
        return None, None
    matrix = max(filtered, key=lambda u: _score_url(u, kind="matrix"), default=None)
    metadata = max(filtered, key=lambda u: _score_url(u, kind="metadata"), default=None)
    return matrix, metadata


def _maybe_gunzip_file_inplace(path: Path) -> None:
    """Si el archivo es gzip (aunque la extensión sea .tsv), descomprimir in-place."""
    if not path.exists() or path.stat().st_size < 3:
        return
    head = path.read_bytes()[:3]
    if head[:2] != b"\x1f\x8b":
        return
    import gzip

    data = gzip.decompress(path.read_bytes())
    path.write_bytes(data)


def _ensure_tsv_from_excel_if_needed(path: Path) -> None:
    """Si el archivo descargado es Excel pero la ruta destino es análisis TSV, convertir."""
    if not path.exists() or path.stat().st_size == 0:
        return
    try:
        pd.read_csv(path, sep="\t", nrows=2)
        return
    except Exception:  # noqa: BLE001
        pass
    try:
        df = pd.read_excel(path, engine="openpyxl")
    except Exception:  # noqa: BLE001
        return
    df.to_csv(path, sep="\t", index=False)


def _harmonize_protein_id_column(matrix_path: Path, *, protein_column: str) -> bool:
    """
    CPTAC/PDC exporta la primera columna como `Gene`; el pipeline espera `protein_id`.
    """
    if protein_column != "protein_id" or not matrix_path.exists():
        return False
    suf = matrix_path.suffix.lower()
    if suf not in {".tsv", ".txt", ""}:
        return False
    try:
        peek = pd.read_csv(matrix_path, sep="\t", nrows=0)
    except Exception:  # noqa: BLE001
        return False
    if len(peek.columns) == 0:
        return False
    first = str(peek.columns[0])
    if first == "protein_id":
        return False
    aliases = {"Gene", "gene", "Gene Symbol", "gene_symbol"}
    if first not in aliases:
        return False
    df = pd.read_csv(matrix_path, sep="\t", low_memory=False)
    df.rename(columns={first: "protein_id"}, inplace=True)
    df.to_csv(matrix_path, sep="\t", index=False)
    return True


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


def _fetch_url_bytes(url: str, headers: dict[str, str] | None, timeout_s: int) -> bytes:
    """Lectura completa con timeout fuerte; usado dentro de un thread vigilado."""
    req = urllib.request.Request(url=url, headers=headers or {})
    old_default = socket.getdefaulttimeout()
    # Para FTP, urllib NO siempre pasa el timeout al socket de control; default global asegura
    # que cualquier read/recv pendiente respete el wall-clock.
    socket.setdefaulttimeout(timeout_s)
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            return resp.read()
    finally:
        socket.setdefaulttimeout(old_default)


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
    # Wall-clock máximo por intento: 2× timeout efectivo del socket (cubre FTP control hangs).
    wallclock_s = max(int(timeout_s) * 2, int(timeout_s) + 30)
    last_err = ""
    for attempt in range(1, retries + 1):
        try:
            with ThreadPoolExecutor(max_workers=1, thread_name_prefix="oncogo-dl") as ex:
                fut = ex.submit(_fetch_url_bytes, url, headers, timeout_s)
                try:
                    payload = fut.result(timeout=wallclock_s)
                except FuturesTimeoutError:
                    fut.cancel()
                    last_err = f"wallclock_timeout_{wallclock_s}s"
                    if attempt < retries:
                        time.sleep(min(2**attempt, 8))
                    continue
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
    *,
    precomputed: ProviderDatasetBundle | None = None,
) -> dict[str, Any]:
    provider_name = str(ds_cfg.get("provider") or ds_cfg.get("download_real", {}).get("provider") or "").strip().lower()
    provider_name = PROVIDER_ALIASES.get(provider_name, provider_name)
    if provider_name not in PROVIDER_REGISTRY:
        return {"status": "skipped_provider_unknown", "provider": provider_name}
    adapter_cls = PROVIDER_REGISTRY[provider_name]
    output_dir = settings.resolve(ds_cfg["matrix_path"]).parent
    try:
        if precomputed is not None and precomputed.provider == provider_name:
            bundle = precomputed
        else:
            adapter = adapter_cls(dataset_name, ds_cfg, global_cfg)
            bundle = adapter.materialize_dataset(output_dir)
        # PDC: la metadata analítica (p. ej. tumor vs normal) debe venir del archivo clinical
        # descargado o curada manualmente; casesPerStudy queda solo en provider_manifest.json.
        if bundle.metadata_records and bundle.provider != "pdc":
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
    # Wall-clock por dataset (override por dataset.download_real.dataset_wallclock_s);
    # cubre todo el bloque matrix+metadata+provider para evitar que un solo accession bloquee la corrida.
    default_dataset_wallclock = int(dl_cfg.get("dataset_wallclock_s", 900))

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
        if isinstance(matrix_url, str) and not matrix_url.strip():
            matrix_url = None
        if isinstance(metadata_url, str) and not metadata_url.strip():
            metadata_url = None

        provider_top = str(ds.get("provider") or "").strip().lower()
        provider_raw = str(spec.get("provider") or provider_top or "").strip().lower()
        provider_resolved = PROVIDER_ALIASES.get(provider_raw, provider_raw)

        source_accession = str(spec.get("source_accession", "")).strip()
        if not source_accession:
            accs = spec.get("accessions") or []
            if isinstance(accs, list) and accs:
                source_accession = str(accs[0]).strip()
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

        precomputed_bundle: ProviderDatasetBundle | None = None
        if provider_resolved == "pdc" and (not matrix_url or not metadata_url):
            try:
                adapter = PDCAdapter(key, ds, settings.config)
                precomputed_bundle = adapter.materialize_dataset(m_path.parent)
                man = precomputed_bundle.manifest
                if not matrix_url and man.get("resolved_matrix_url"):
                    matrix_url = man["resolved_matrix_url"]
                    row["matrix_autodiscovered"] = True
                if not metadata_url and man.get("resolved_metadata_url"):
                    metadata_url = man["resolved_metadata_url"]
                    row["metadata_autodiscovered"] = True
            except Exception as exc:  # noqa: BLE001
                row["pdc_resolution_error"] = str(exc)

        # Autodiscovery cuando no hay URL directa definida.
        meta_auto = bool(spec.get("metadata_autodiscovery", True))
        need_matrix_url = not matrix_url
        need_metadata_url = not metadata_url and meta_auto
        if need_matrix_url or need_metadata_url:
            try:
                auto_matrix = auto_metadata = None
                if provider_resolved == "pride":
                    acc_list = list(spec.get("accessions") or [])
                    if not acc_list and source_accession:
                        acc_list = [source_accession]
                    for acc in acc_list:
                        if not acc:
                            continue
                        am, ad = _autodiscover_pride_urls(
                            str(acc).strip(),
                            timeout_s=timeout_s,
                            file_selection=spec.get("file_selection", {}) or {},
                        )
                        if am or ad:
                            auto_matrix, auto_metadata = am, ad
                            row["pride_accession_used"] = str(acc).strip()
                            break
                elif source_catalog_url and provider_resolved not in {"pdc", "pride"}:
                    auto_matrix, auto_metadata = _autodiscover_from_catalog_page(
                        source_catalog_url, timeout_s=timeout_s
                    )
                if need_matrix_url and auto_matrix:
                    matrix_url = auto_matrix
                    row["matrix_autodiscovered"] = True
                if need_metadata_url and auto_metadata:
                    metadata_url = auto_metadata
                    row["metadata_autodiscovered"] = True
            except Exception as exc:  # noqa: BLE001
                row["autodiscovery_error"] = str(exc)

        if matrix_url:
            url = _build_url(str(matrix_url), qparams)
            if not _is_useful_data_url(url):
                row["matrix_status"] = "error_invalid_url"
                row["matrix_detail"] = "url_not_useful_for_data"
            else:
                print(f"[download_real] {key}: matrix GET {url}")
                matrix_overwrite = overwrite
                if (
                    provider_resolved == "pdc"
                    and m_path.exists()
                    and not overwrite
                    and not _matrix_has_tmt_style_samples(m_path)
                ):
                    matrix_overwrite = True
                    row["matrix_redownload_forced"] = True
                    print(
                        f"[download_real] {key}: matriz local no parece TMT10 PDC válida "
                        f"({m_path.name}, {m_path.stat().st_size} bytes) — forzando re-descarga "
                        "(stub o archivo incompleto; overwrite seguía en false)."
                    )
                ok, status = _download_with_retries(
                    url=url,
                    out_path=m_path,
                    timeout_s=timeout_s,
                    retries=retries,
                    headers=headers,
                    overwrite=matrix_overwrite,
                )
                row["matrix_status"] = "ok" if ok else "error"
                row["matrix_detail"] = status
                if ok:
                    _maybe_gunzip_file_inplace(m_path)
                    _ensure_tsv_from_excel_if_needed(m_path)
                    if bool(spec.get("dia_nn_postprocess_to_wide", False)):
                        try:
                            from oncogo_crc.utils.dia_nn_plasma import postprocess_dia_nn_plasma_download

                            postprocess_dia_nn_plasma_download(m_path, md_path)
                            row["matrix_detail"] = (row.get("matrix_detail") or "") + ";dia_nn_pivoted_wide"
                            if md_path.exists():
                                row["metadata_status"] = "ok_generated"
                        except Exception as exc:  # noqa: BLE001
                            row["matrix_status"] = "error"
                            row["matrix_detail"] = f"dia_nn_postprocess_failed:{exc}"
            row["matrix_url"] = url
        if m_path.exists():
            _harmonize_protein_id_column(m_path, protein_column=str(ds.get("protein_column", "protein_id")))
            if row.get("matrix_status") == "ok":
                row["matrix_sha256"] = _sha256_file(m_path)

        if metadata_url:
            url = _build_url(str(metadata_url), qparams)
            if not _is_useful_data_url(url):
                row["metadata_status"] = "error_invalid_url"
                row["metadata_detail"] = "url_not_useful_for_data"
            else:
                print(f"[download_real] {key}: metadata GET {url}")
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
                if ok:
                    _maybe_gunzip_file_inplace(md_path)
                    _ensure_tsv_from_excel_if_needed(md_path)
                    row["metadata_sha256"] = _sha256_file(md_path)
            row["metadata_url"] = url

        if (
            not metadata_url
            and md_path.exists()
            and row.get("metadata_status") == "skipped_no_url"
            and not bool(spec.get("metadata_autodiscovery", True))
        ):
            row["metadata_status"] = "ok_local"

        rows.append(row)
        provider_summary = _run_provider_inventory(
            settings,
            key,
            ds,
            settings.config,
            precomputed=precomputed_bundle if provider_resolved == "pdc" else None,
        )
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
            "matrix_sha256": row.get("matrix_sha256"),
            "metadata_sha256": row.get("metadata_sha256"),
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
    print(f"[download_real] manifest -> {settings.relative_to_project(out)}")
