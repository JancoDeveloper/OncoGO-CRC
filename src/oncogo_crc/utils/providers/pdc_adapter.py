from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path
from typing import Any

from oncogo_crc.utils.providers.base import BaseProviderAdapter, ProviderDatasetBundle

PDC_GRAPHQL_ENDPOINT = "https://pdc.cancer.gov/graphql"


def _post_graphql(query: str, variables: dict[str, Any] | None = None, token: str | None = None) -> dict[str, Any]:
    payload = json.dumps({"query": query, "variables": variables or {}}).encode("utf-8")
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    req = urllib.request.Request(PDC_GRAPHQL_ENDPOINT, data=payload, headers=headers, method="POST")
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    if "errors" in data:
        raise RuntimeError(str(data["errors"]))
    return data.get("data", {})


class PDCAdapter(BaseProviderAdapter):
    provider_name = "pdc"

    def materialize_dataset(self, output_dir: Path) -> ProviderDatasetBundle:
        output_dir.mkdir(parents=True, exist_ok=True)
        dr = self.dataset_cfg.get("download_real", {})
        token_env = self.global_cfg.get("download_real", {}).get("auth", {}).get("pdc_token_env", "PDC_TOKEN")
        token = os.environ.get(token_env)
        study_submitter_id = dr.get("study_submitter_id") or dr.get("source_accession") or ""
        project_submitter_id = dr.get("project_submitter_id") or ""

        bundle = ProviderDatasetBundle(dataset_name=self.dataset_name, provider=self.provider_name, raw_dir=output_dir)
        query = (
            "query StudyCatalog { studyCatalog { study_id study_submitter_id project_submitter_id project_name "
            "study_name disease_type primary_site cases_count aliquots_count } }"
        )
        data = _post_graphql(query, token=token)
        studies = data.get("studyCatalog", []) or []
        matches = []
        for row in studies:
            if study_submitter_id and study_submitter_id.lower() not in str(row.get("study_submitter_id", "")).lower():
                continue
            if project_submitter_id and project_submitter_id.lower() not in str(row.get("project_submitter_id", "")).lower():
                continue
            matches.append(row)
        bundle.manifest = {
            "dataset_name": self.dataset_name,
            "provider": self.provider_name,
            "token_present": bool(token),
            "n_study_matches": len(matches),
            "study_matches": matches[:10],
        }
        (output_dir / "provider_manifest.json").write_text(json.dumps(bundle.manifest, indent=2), encoding="utf-8")
        bundle.status = "partial"
        if not matches:
            bundle.warnings.append("Sin coincidencias de estudio en PDC; revisar study_submitter_id/project_submitter_id.")
        return bundle
