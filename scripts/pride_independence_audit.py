#!/usr/bin/env python3
"""Auditoría de independencia nivel A sobre la API de PRIDE.

Compara checksum, tamaño y nombre normalizado entre accesiones. Un checksum
idéntico es redepósito literal. Cero colisiones no demuestran independencia
de pacientes: hace falta la huella a nivel de proteína (aún no implementada).

    python scripts/pride_independence_audit.py --accessions PXD019504 PXD031705 PXD034575

La salida queda en results/independence/. Repetir con red y TLS del sistema
antes de citar el hallazgo en la tesis.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from oncogo_crc.certification.independence import FileRecord, detect_file_collisions  # noqa: E402

PRIDE_V2 = "https://www.ebi.ac.uk/pride/ws/archive/v2"
DELAY_S = 0.35


def _get(url: str) -> dict | list | None:
    time.sleep(DELAY_S)
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "oncogo-crc/independence-audit"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        print(f"HTTP {exc.code} en {url}", file=sys.stderr)
        return None
    except urllib.error.URLError as exc:
        print(f"fallo de red en {url}: {exc}", file=sys.stderr)
        return None


def fetch_files(accession: str) -> list[FileRecord]:
    records: list[FileRecord] = []
    page = 0
    while page < 200:
        url = (
            f"{PRIDE_V2}/files/byProject?accession={urllib.parse.quote(accession)}"
            f"&pageSize=100&page={page}"
        )
        raw = _get(url)
        if raw is None:
            break
        if isinstance(raw, dict):
            items = raw.get("_embedded", {}).get("files") or raw.get("files") or raw.get("list") or []
        else:
            items = raw
        if not items:
            break
        for it in items:
            checksum = (it.get("checksum") or it.get("md5") or "").strip().lower() or None
            records.append(
                FileRecord(
                    accession=accession,
                    file_name=(it.get("fileName") or it.get("name") or "").strip(),
                    file_size=it.get("fileSizeBytes") or it.get("fileSize"),
                    checksum=checksum,
                )
            )
        if len(items) < 100:
            break
        page += 1
    print(f"{accession}: {len(records)} ficheros")
    return records


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accessions", nargs="+", required=True)
    ap.add_argument("--outdir", type=Path, default=Path("results/independence"))
    args = ap.parse_args(argv)

    records: list[FileRecord] = []
    for acc in args.accessions:
        records.extend(fetch_files(acc))
    if len({r.accession for r in records}) < 2:
        print("se necesitan ficheros de al menos dos accesiones", file=sys.stderr)
        return 1

    collisions = detect_file_collisions(records)
    args.outdir.mkdir(parents=True, exist_ok=True)
    tsv = args.outdir / "file_collisions.tsv"
    with tsv.open("w", newline="", encoding="utf-8") as fh:
        fields = ["kind", "accession_a", "accession_b", "file_a", "file_b", "evidence", "severity"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for col in collisions:
            writer.writerow(
                {
                    "kind": col.kind,
                    "accession_a": col.accession_a,
                    "accession_b": col.accession_b,
                    "file_a": col.file_a,
                    "file_b": col.file_b,
                    "evidence": col.evidence,
                    "severity": col.severity,
                }
            )
    n_checksum = sum(1 for c in collisions if c.kind == "checksum")
    summary = {
        "accessions": list(args.accessions),
        "n_files": len(records),
        "n_checksum_collisions": n_checksum,
        "n_collisions_total": len(collisions),
        "veredicto_nivel_A": "REDEPOSITO_LITERAL_DETECTADO" if n_checksum else "SIN_COLISION_EXACTA",
        "nota": (
            "El nivel A es concluyente en positivo y no concluyente en negativo. "
            "Declarar si esta pasada usó TLS verificado."
        ),
    }
    (args.outdir / "audit_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
