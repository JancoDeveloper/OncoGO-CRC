from __future__ import annotations

from collections.abc import Callable

from oncogo_crc.settings import Settings

from oncogo_crc.stages import (
    clinical_translation,
    differential_analysis,
    download_data,
    download_real,
    generate_report,
    go_enrichment,
    harmonize_ids,
    harmonize_public_metadata,
    ingest_cptac,
    ml_prioritization,
    network_modules,
    publish_results,
    qc_normalize,
    stability_analysis,
    validate_inputs,
    version_inputs,
)

StageFn = Callable[[Settings], None]

STAGES: dict[str, StageFn] = {
    "download_real": download_real.run,
    "harmonize_public_metadata": harmonize_public_metadata.run,
    "download": download_data.run,
    "ingest_cptac": ingest_cptac.run,
    "validate_inputs": validate_inputs.run,
    "version_inputs": version_inputs.run,
    "harmonize": harmonize_ids.run,
    "qc": qc_normalize.run,
    "differential": differential_analysis.run,
    "go": go_enrichment.run,
    "network": network_modules.run,
    "ml": ml_prioritization.run,
    "clinical": clinical_translation.run,
    "report": generate_report.run,
    "stability": stability_analysis.run,
    "publish_results": publish_results.run,
}

STAGE_ORDER = list(STAGES.keys())

# Etapas opcionales que no deben ejecutarse en `oncogo-crc run --stage all`
EXCLUDED_FROM_ALL: frozenset[str] = frozenset({"ingest_cptac"})
