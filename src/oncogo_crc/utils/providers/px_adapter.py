from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

from oncogo_crc.utils.providers.base import BaseProviderAdapter, ProviderDatasetBundle, ProviderFile

PX_BASE = "https://proteomecentral.proteomexchange.org/cgi/GetDataset"


class PXAdapter(BaseProviderAdapter):
    provider_name = "px"

    def materialize_dataset(self, output_dir: Path) -> ProviderDatasetBundle:
        output_dir.mkdir(parents=True, exist_ok=True)
        dr = self.dataset_cfg.get("download_real", {})
        accession = str(dr.get("accession") or dr.get("source_accession") or "").strip()
        fallback = str(dr.get("fallback_project_uri") or dr.get("source_project_url") or "").strip()
        bundle = ProviderDatasetBundle(dataset_name=self.dataset_name, provider=self.provider_name, raw_dir=output_dir)
        if accession:
            url = f"{PX_BASE}?ID={urllib.parse.quote(accession)}&outputMode=XML"
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=60) as resp:
                xml = resp.read().decode("utf-8", errors="ignore")
            (output_dir / "px_dataset.xml").write_text(xml, encoding="utf-8")
            bundle.manifest["accession"] = accession
            bundle.manifest["xml_retrieved"] = True
        if fallback:
            bundle.files.append(ProviderFile(name=f"{accession or self.dataset_name}_fallback_project", url=fallback, source="PX"))
        bundle.manifest.update({"dataset_name": self.dataset_name, "provider": self.provider_name, "n_files": len(bundle.files)})
        (output_dir / "provider_manifest.json").write_text(json.dumps(bundle.manifest, indent=2), encoding="utf-8")
        bundle.status = "partial"
        return bundle
