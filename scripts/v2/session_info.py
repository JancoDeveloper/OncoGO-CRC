"""Escribe results/v2/session_info.txt con las versiones exactas de la corrida.

Si R no está instalado, lo dice. Eso importa porque el diseño v2 nombra
limma, metafor, camera y mada, que son paquetes de R.
"""

from __future__ import annotations

import importlib.metadata as md
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

WORKSPACE = Path("/Volumes/JANCODATA/OncoGO-CRC")
REPO = Path(__file__).resolve().parents[2]

PAQUETES = (
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "statsmodels",
    "matplotlib",
    "seaborn",
    "networkx",
    "reportlab",
    "pypdf",
    "pypdfium2",
    "requests",
    "PyYAML",
    "gseapy",
    "shap",
    "pytest",
    "ruff",
)


def git(*args: str) -> str:
    try:
        return subprocess.run(
            ["git", *args], cwd=REPO, capture_output=True, text=True, check=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "no_disponible"


def main() -> int:
    lines = [
        "# OncoGO-CRC v2 — entorno de la corrida",
        f"fecha_utc: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"python: {sys.version.split()[0]}",
        f"ejecutable: {sys.executable}",
        f"plataforma: {platform.platform()}",
        f"maquina: {platform.machine()}",
        f"repo_commit: {git('rev-parse', 'HEAD')}",
        f"repo_rama: {git('rev-parse', '--abbrev-ref', 'HEAD')}",
        f"repo_sucio: {'si' if git('status', '--porcelain') else 'no'}",
        f"workspace: {WORKSPACE}",
        "",
        "# Paquetes de Python",
    ]
    for nombre in PAQUETES:
        try:
            lines.append(f"{nombre}: {md.version(nombre)}")
        except md.PackageNotFoundError:
            lines.append(f"{nombre}: no_instalado")

    lines += ["", "# R y paquetes de R"]
    rscript = shutil.which("Rscript")
    if rscript is None:
        lines += [
            "R: no_instalado",
            "limma: no_disponible (treat está reimplementado en Python: src/oncogo_crc/utils/moderated.py)",
            "metafor: no_disponible (REML reimplementado en Python: src/oncogo_crc/utils/meta.py)",
            "camera: no_disponible (prueba competitiva reimplementada en Python)",
            "mada: no_disponible (meta-análisis diagnóstico reimplementado en Python)",
            "renv.lock: no_aplica sin R",
        ]
    else:
        ver = subprocess.run(
            [rscript, "-e", "cat(R.version.string)"], capture_output=True, text=True
        ).stdout.strip()
        lines.append(f"R: {ver}")
        for pkg in ("limma", "metafor", "mada"):
            out = subprocess.run(
                [rscript, "-e", f'cat(tryCatch(as.character(packageVersion("{pkg}")), error=function(e) "no_instalado"))'],
                capture_output=True,
                text=True,
            ).stdout.strip()
            lines.append(f"{pkg}: {out or 'no_instalado'}")

    out_path = WORKSPACE / "results/v2/session_info.txt"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[session_info] -> {out_path}")

    lock = REPO / "requirements-lock.txt"
    frozen = subprocess.run(
        [sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True
    ).stdout
    keep = [
        linea
        for linea in frozen.splitlines()
        if linea and not linea.startswith("-e ") and "oncogo" not in linea
    ]
    lock.write_text(
        "# Generado por scripts/v2/session_info.py. Versiones exactas de la corrida v2.\n"
        + "\n".join(sorted(keep))
        + "\n",
        encoding="utf-8",
    )
    print(f"[session_info] -> {lock}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
