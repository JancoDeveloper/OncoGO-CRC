"""Manifiesto de corrida del pipeline: trazabilidad para publicación y auditoría.

Se escribe en ``results/run_manifest.json`` (workspace). Las etapas actualizan secciones
para que ``publish_results`` / digest no mezclen artefactos de corridas anteriores.
"""

from __future__ import annotations

import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from oncogo_crc.settings import Settings


def manifest_path(settings: Settings) -> Path:
    return settings.resolve(settings.config.get("paths", {}).get("results_dir", "results")) / "run_manifest.json"


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def init_run(settings: Settings, *, stage_order: list[str]) -> None:
    """Inicializa manifiesto al comenzar ``oncogo-crc run`` (sobrescribe el anterior)."""
    run_id = f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}_{uuid.uuid4().hex[:10]}"
    doc: dict[str, Any] = {
        "run_id": run_id,
        "started_at": _utc_now(),
        "finished_at": None,
        "demo": bool(settings.demo),
        "stage_order": list(stage_order),
        "stages": {},
        "notes": [],
    }
    path = manifest_path(settings)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=2), encoding="utf-8")


def load_doc(settings: Settings) -> dict[str, Any]:
    path = manifest_path(settings)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def merge_doc(settings: Settings, patch: dict[str, Any]) -> None:
    doc = load_doc(settings)
    if not doc:
        return
    for k, v in patch.items():
        if isinstance(v, dict) and isinstance(doc.get(k), dict):
            doc[k] = {**doc[k], **v}
        else:
            doc[k] = v
    manifest_path(settings).write_text(json.dumps(doc, indent=2), encoding="utf-8")


def finalize_run(settings: Settings) -> None:
    merge_doc(settings, {"finished_at": _utc_now()})


def record_stage_meta(settings: Settings, stage: str, **meta: Any) -> None:
    path = manifest_path(settings)
    if not path.exists():
        doc: dict[str, Any] = {
            "run_id": f"adhoc_{uuid.uuid4().hex[:10]}",
            "started_at": _utc_now(),
            "finished_at": None,
            "demo": bool(settings.demo),
            "stage_order": [],
            "stages": {},
            "notes": ["run_manifest creado en caliente (etapa suelta sin init_run en CLI)."],
        }
    else:
        doc = load_doc(settings)
        if not doc:
            return
    st = doc.setdefault("stages", {})
    cur = st.get(stage, {})
    if isinstance(cur, dict):
        cur = {**cur, **meta}
    else:
        cur = dict(meta)
    st[stage] = cur
    doc["stages"] = st
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=2), encoding="utf-8")


def gsea_executed_this_run(settings: Settings) -> bool:
    doc = load_doc(settings)
    g = doc.get("stages", {}).get("gsea", {})
    if isinstance(g, dict):
        return bool(g.get("executed", False))
    return False


def gsea_comparison_trustworthy(settings: Settings, comp_key: str) -> bool:
    """True solo si el manifiesto indica que GSEA generó salida para esta comparación en esta corrida."""
    if not gsea_executed_this_run(settings):
        return False
    doc = load_doc(settings)
    comps = doc.get("stages", {}).get("gsea", {}).get("comparisons", {})
    if not isinstance(comps, dict):
        return False
    ent = comps.get(comp_key)
    return isinstance(ent, dict) and ent.get("status") == "ok"


def ml_transfer_executed_this_run(settings: Settings) -> bool:
    doc = load_doc(settings)
    m = doc.get("stages", {}).get("ml_transfer", {})
    return isinstance(m, dict) and bool(m.get("executed", False))


def ml_transfer_trustworthy(settings: Settings) -> bool:
    return ml_transfer_executed_this_run(settings)


def note(settings: Settings, text: str) -> None:
    doc = load_doc(settings)
    notes = doc.setdefault("notes", [])
    if isinstance(notes, list):
        notes.append(f"{_utc_now()}: {text}")
        doc["notes"] = notes
        manifest_path(settings).write_text(json.dumps(doc, indent=2), encoding="utf-8")
