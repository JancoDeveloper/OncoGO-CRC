from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path
from typing import Any

from oncogo_crc.utils.provider_file_selection import pick_best_files
from oncogo_crc.utils.providers.base import BaseProviderAdapter, ProviderDatasetBundle, ProviderFile

PDC_GRAPHQL_ENDPOINT = "https://pdc.cancer.gov/graphql"


def _post_graphql(query: str, variables: dict[str, Any] | None = None, token: str | None = None) -> dict[str, Any]:
    payload = json.dumps({"query": query, "variables": variables or {}}).encode("utf-8")
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    req = urllib.request.Request(PDC_GRAPHQL_ENDPOINT, data=payload, headers=headers, method="POST")
    with urllib.request.urlopen(req, timeout=120) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    if "errors" in data:
        raise RuntimeError(str(data["errors"]))
    return data.get("data", {})


def _list_pdc_study_ids(token: str | None) -> list[str]:
    data = _post_graphql("query { studyCatalog { pdc_study_id } }", token=token)
    rows = data.get("studyCatalog", []) or []
    return [str(r.get("pdc_study_id")) for r in rows if r.get("pdc_study_id")]


def _resolve_pdc_study_id(dr: dict[str, Any], token: str | None) -> str:
    """
    Identificador de estudio en el API actual: `pdc_study_id` (ej. PDC000116).
    Acepta claves `pdc_study_submitter_id` o `pdc_study_id` en YAML.
    """
    direct = str(dr.get("pdc_study_submitter_id") or dr.get("pdc_study_id") or "").strip()
    if direct:
        return direct
    # Legado: intentar emparejar por texto en el catálogo mínimo (solo IDs)
    study_sub = str(dr.get("study_submitter_id") or "").strip().lower()
    proj_sub = str(dr.get("project_submitter_id") or "").strip().lower()
    if not study_sub and not proj_sub:
        return ""
    # Sin metadatos por estudio en studyCatalog actual, no podemos resolver; el usuario debe fijar pdc_study_id.
    return ""


class PDCAdapter(BaseProviderAdapter):
    provider_name = "pdc"

    def materialize_dataset(self, output_dir: Path) -> ProviderDatasetBundle:
        output_dir.mkdir(parents=True, exist_ok=True)
        dr = self.dataset_cfg.get("download_real", {})
        token_env = self.global_cfg.get("download_real", {}).get("auth", {}).get("pdc_token_env", "PDC_TOKEN")
        token = os.environ.get(token_env)

        bundle = ProviderDatasetBundle(dataset_name=self.dataset_name, provider=self.provider_name, raw_dir=output_dir)

        pdc_study_id = _resolve_pdc_study_id(dr, token)
        bundle.manifest = {
            "dataset_name": self.dataset_name,
            "provider": self.provider_name,
            "token_present": bool(token),
            "pdc_study_id": pdc_study_id,
        }

        if not pdc_study_id:
            bundle.status = "partial"
            bundle.warnings.append(
                "Falta pdc_study_submitter_id o pdc_study_id (ej. PDC000116) en download_real para resolver archivos PDC."
            )
            (output_dir / "provider_manifest.json").write_text(
                json.dumps(bundle.manifest, indent=2), encoding="utf-8"
            )
            return bundle

        catalog_ids = _list_pdc_study_ids(token)
        if pdc_study_id not in catalog_ids:
            bundle.warnings.append(
                f"pdc_study_id={pdc_study_id!r} no aparece en studyCatalog ({len(catalog_ids)} estudios); "
                "verificar ID en el portal PDC."
            )

        fq = (
            "query StudyFiles($pdc_study_id: String!) { "
            "filesPerStudy(pdc_study_id: $pdc_study_id, acceptDUA: true) { "
            "file_name file_size file_type md5sum signedUrl { url } } }"
        )
        try:
            fdata = _post_graphql(fq, {"pdc_study_id": pdc_study_id}, token=token)
            raw_files = fdata.get("filesPerStudy", []) or []
        except Exception as exc:  # noqa: BLE001
            bundle.warnings.append(f"filesPerStudy falló: {exc}")
            raw_files = []

        file_selection = dr.get("file_selection", {})
        for f in raw_files:
            surl = f.get("signedUrl") or {}
            url = surl.get("url") if isinstance(surl, dict) else None
            bundle.files.append(
                ProviderFile(
                    name=str(f.get("file_name") or ""),
                    url=url,
                    size=f.get("file_size"),
                    file_type=f.get("file_type"),
                    checksum=f.get("md5sum"),
                    source="PDC",
                    extra={"pdc_study_id": pdc_study_id},
                )
            )

        candidates = [
            {"name": pf.name, "url": pf.url or "", "size": pf.size, "file_type": pf.file_type}
            for pf in bundle.files
            if pf.url and isinstance(pf.url, str) and pf.url.startswith("http")
        ]
        selected = pick_best_files(candidates, file_selection) if candidates else {"matrix": None, "metadata": None}
        mat = selected.get("matrix") or {}
        meta = selected.get("metadata") or {}
        # Preferir manifestación clínica explícita frente a solo biospecimen.
        if meta.get("name") and "biospecimen" in meta["name"].lower():
            for c in candidates:
                nm = (c.get("name") or "").lower()
                if "clinical_data" in nm or ("clinical" in nm and "data" in nm):
                    meta = c
                    selected["metadata"] = c
                    break
        bundle.manifest["selected_matrix_candidate"] = mat
        bundle.manifest["selected_metadata_candidate"] = meta
        bundle.manifest["resolved_matrix_url"] = mat.get("url") if mat else None
        bundle.manifest["resolved_metadata_url"] = meta.get("url") if meta else None
        bundle.manifest["n_files_in_study"] = len(raw_files)

        bundle.status = "partial"
        if not bundle.manifest.get("resolved_matrix_url"):
            bundle.warnings.append(
                "No se seleccionó matriz procesada; refinar file_selection (excluir .raw) o revisar archivos del estudio."
            )
        bundle.warnings.append(
            "Las URLs firmadas PDC caducan; re-ejecutar download_real antes del análisis final y documentar pdc_study_id."
        )

        (output_dir / "provider_manifest.json").write_text(
            json.dumps(bundle.manifest, indent=2, default=str), encoding="utf-8"
        )
        return bundle
