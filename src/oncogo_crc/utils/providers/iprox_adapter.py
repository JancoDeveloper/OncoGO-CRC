from __future__ import annotations

import json
from pathlib import Path

from oncogo_crc.utils.providers.base import BaseProviderAdapter, ProviderDatasetBundle, ProviderFile


class IProXAdapter(BaseProviderAdapter):
    provider_name = "iprox"

    def materialize_dataset(self, output_dir: Path) -> ProviderDatasetBundle:
        output_dir.mkdir(parents=True, exist_ok=True)
        dr = self.dataset_cfg.get("download_real", {})
        accession = str(dr.get("accession") or dr.get("source_accession") or "").strip()
        uri = str(dr.get("fallback_project_uri") or dr.get("source_project_url") or "").strip()
        bundle = ProviderDatasetBundle(dataset_name=self.dataset_name, provider=self.provider_name, raw_dir=output_dir)
        if uri:
            bundle.files.append(ProviderFile(name=f"{accession or self.dataset_name}_iprox_project", url=uri, source="iProX"))
        bundle.manifest = {"dataset_name": self.dataset_name, "provider": self.provider_name, "accession": accession, "n_files": len(bundle.files)}
        (output_dir / "provider_manifest.json").write_text(json.dumps(bundle.manifest, indent=2), encoding="utf-8")
        bundle.status = "partial"
        return bundle
