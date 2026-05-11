from __future__ import annotations

from oncogo_crc.settings import Settings
from oncogo_crc.utils.web_publish import build_web_pages


def run(settings: Settings) -> None:
    web_cfg = settings.config.get("web_publish", {})
    if not bool(web_cfg.get("enabled", True)):
        print("[publish_results] deshabilitado en config.yaml (web_publish.enabled=false).")
        return
    out = build_web_pages(settings)
    print(f"[publish_results] web -> {settings.relative_to_project(out)}")
