from __future__ import annotations

from collections.abc import Callable

from oncogo_crc.settings import Settings

from oncogo_crc.stages import (
    clinical_translation,
    differential_analysis,
    download_data,
    generate_report,
    go_enrichment,
    harmonize_ids,
    ml_prioritization,
    network_modules,
    qc_normalize,
)

StageFn = Callable[[Settings], None]

STAGES: dict[str, StageFn] = {
    "download": download_data.run,
    "harmonize": harmonize_ids.run,
    "qc": qc_normalize.run,
    "differential": differential_analysis.run,
    "go": go_enrichment.run,
    "network": network_modules.run,
    "ml": ml_prioritization.run,
    "clinical": clinical_translation.run,
    "report": generate_report.run,
}

STAGE_ORDER = list(STAGES.keys())
