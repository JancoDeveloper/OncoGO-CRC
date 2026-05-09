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
    ml_prioritization,
    network_modules,
    publish_results,
    qc_normalize,
    stability_analysis,
    validate_inputs,
)

StageFn = Callable[[Settings], None]

STAGES: dict[str, StageFn] = {
    "download_real": download_real.run,
    "download": download_data.run,
    "validate_inputs": validate_inputs.run,
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
