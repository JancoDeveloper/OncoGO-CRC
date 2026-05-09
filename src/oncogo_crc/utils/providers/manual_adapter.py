from __future__ import annotations

import json
from pathlib import Path

from oncogo_crc.utils.providers.base import BaseProviderAdapter, ProviderDatasetBundle


class ManualAdapter(BaseProviderAdapter):
    provider_name = "manual"

    def materialize_dataset(self, output_dir: Path) -> ProviderDatasetBundle:
        output_dir.mkdir(parents=True, exist_ok=True)
        bundle = ProviderDatasetBundle(dataset_name=self.dataset_name, provider=self.provider_name, raw_dir=output_dir)
        bundle.manifest = {
            "dataset_name": self.dataset_name,
            "provider": self.provider_name,
            "note": "Dataset configurado en modo manual; completar matrix_url/metadata_url o carga local.",
        }
        (output_dir / "provider_manifest.json").write_text(json.dumps(bundle.manifest, indent=2), encoding="utf-8")
        bundle.status = "manual_pending"
        bundle.warnings.append("Provider manual: sin descarga automática configurada.")
        return bundle
