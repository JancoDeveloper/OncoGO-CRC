from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from oncogo_crc.utils.provider_file_selection import pick_best_files
from oncogo_crc.utils.providers.base import BaseProviderAdapter, ProviderDatasetBundle, ProviderFile

PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v2"


def _get_json(url: str, timeout: int = 60) -> dict[str, Any]:
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


class PRIDEAdapter(BaseProviderAdapter):
    provider_name = "pride"

    def materialize_dataset(self, output_dir: Path) -> ProviderDatasetBundle:
        output_dir.mkdir(parents=True, exist_ok=True)
        dr = self.dataset_cfg.get("download_real", {})
        accessions = dr.get("accessions") or ([dr.get("source_accession")] if dr.get("source_accession") else [])

        bundle = ProviderDatasetBundle(
            dataset_name=self.dataset_name,
            provider=self.provider_name,
            raw_dir=output_dir,
        )
        for acc in accessions:
            if not acc:
                continue
            project = _get_json(f"{PRIDE_API_BASE}/projects/{urllib.parse.quote(str(acc))}")
            bundle.metadata_records.append(
                {
                    "sample_id": None,
                    "condition": None,
                    "dataset_origin": self.dataset_name,
                    "provider": self.provider_name,
                    "project_accession": acc,
                    "title": project.get("title"),
                }
            )
            page = 0
            while True:
                files_json = _get_json(f"{PRIDE_API_BASE}/projects/{acc}/files?page={page}&pageSize=200")
                if isinstance(files_json, list):
                    files = files_json
                else:
                    files = (files_json.get("_embedded") or {}).get("files") or files_json.get("list") or []
                if not files:
                    break
                for f in files:
                    # Compatibilidad con distintas formas del payload de PRIDE
                    file_name = str(f.get("fileName") or f.get("name") or "")
                    dl = f.get("downloadLink") or f.get("publicLink")
                    if not dl:
                        for loc in f.get("publicFileLocations", []) or []:
                            if isinstance(loc, str) and loc.startswith(("http://", "https://", "ftp://")):
                                dl = loc
                                break
                            if isinstance(loc, dict):
                                for k in ("value", "url", "name"):
                                    v = loc.get(k)
                                    if isinstance(v, str) and v.startswith(("http://", "https://", "ftp://")):
                                        dl = v
                                        break
                                if dl:
                                    break
                    bundle.files.append(
                        ProviderFile(
                            name=file_name,
                            url=dl,
                            size=f.get("fileSize") or f.get("fileSizeBytes"),
                            file_type=f.get("fileCategory"),
                            checksum=f.get("checksum"),
                            source="PRIDE",
                            extra={"project_accession": acc},
                        )
                    )
                if len(files) < 200:
                    break
                page += 1
        selected = pick_best_files(
            [
                {"name": f.name, "url": f.url, "size": f.size, "file_type": f.file_type}
                for f in bundle.files
            ],
            dr.get("file_selection", {}),
        )
        bundle.manifest = {
            "dataset_name": self.dataset_name,
            "provider": self.provider_name,
            "n_files": len(bundle.files),
            "n_metadata_records": len(bundle.metadata_records),
            "selected_matrix_candidate": selected.get("matrix"),
            "selected_metadata_candidate": selected.get("metadata"),
        }
        (output_dir / "provider_manifest.json").write_text(json.dumps(bundle.manifest, indent=2), encoding="utf-8")
        bundle.status = "partial"
        if not bundle.files:
            bundle.warnings.append("No se detectaron archivos en PRIDE para las accesiones configuradas.")
        return bundle
