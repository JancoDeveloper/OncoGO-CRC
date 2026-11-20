#!/usr/bin/env python3
"""Verificación local: hashes, pruebas, corrida y comparación.

Escribe reports/verificacion_local_<fecha>.md. No descarga los crudos ni
escribe en el disco de datos. Si el árbol de git tiene cambios de contenido,
la etapa sellada no arranca: el propio pipeline lo rechaza.
"""

from __future__ import annotations

import argparse
import csv
import os
import platform
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]


def _correr(comando: list[str], cwd: Path, env: dict[str, str]) -> tuple[int, str, float]:
    inicio = time.perf_counter()
    resultado = subprocess.run(comando, cwd=cwd, env=env, capture_output=True, text=True)
    duracion = time.perf_counter() - inicio
    texto = (resultado.stdout or "") + (resultado.stderr or "")
    return resultado.returncode, texto[-8000:], duracion


def _git(args: list[str], env: dict[str, str]) -> str:
    resultado = subprocess.run(["git", *args], cwd=RAIZ, env=env, capture_output=True, text=True)
    return (resultado.stdout or "").strip()


def _entorno(py: Path) -> list[str]:
    lineas = [
        f"- SO: {platform.platform()}",
        f"- CPU: {platform.processor() or platform.machine()}",
        f"- máquina: {platform.machine()}",
        f"- Python: {sys.version.split()[0]} ({py})",
    ]
    try:
        memoria = subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, check=False)
        if memoria.returncode == 0 and memoria.stdout.strip().isdigit():
            gib = int(memoria.stdout.strip()) / (1024 ** 3)
            lineas.append(f"- RAM: {gib:.1f} GiB")
    except OSError:
        pass
    freeze = subprocess.run([str(py), "-m", "pip", "freeze"], capture_output=True, text=True)
    paquetes = [
        linea for linea in (freeze.stdout or "").splitlines()
        if linea.startswith(("numpy==", "pandas==", "scipy==", "scikit-learn==", "statsmodels==", "PyYAML==", "directlfq==", "gseapy==", "shap==", "pytest==", "ruff=="))
    ]
    lineas.append("- paquetes: " + ", ".join(paquetes))
    return lineas


def _tabla_hashes(path: Path) -> list[str]:
    if not path.is_file():
        return ["No se escribió la tabla de hashes."]
    lineas = ["| fuente | dataset | lado | estado | bytes |", "| --- | --- | --- | --- | --- |"]
    with path.open(encoding="utf-8") as fichero:
        for fila in csv.DictReader(fichero, delimiter="\t"):
            lineas.append(
                f"| {fila['fuente']} | {fila['dataset']} | {fila['lado']} | {fila['estado']} | {fila['bytes']} |"
            )
    return lineas


def _tabla_salidas(path: Path) -> list[str]:
    if not path.is_file():
        return ["No hubo comparación de salidas."]
    lineas = ["| fichero | clase | diferencia máxima |", "| --- | --- | --- |"]
    with path.open(encoding="utf-8") as fichero:
        for fila in csv.DictReader(fichero, delimiter="\t"):
            if fila["clase"] == "identico":
                continue
            lineas.append(f"| {fila['fichero']} | {fila['clase']} | {fila['diferencia_maxima']} |")
    if len(lineas) == 2:
        lineas.append("| — | todas idénticas en bytes | |")
    return lineas


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datos", type=Path, required=True)
    parser.add_argument("--salidas", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=RAIZ / ".venv" / "bin" / "python")
    parser.add_argument("--fecha", default=date.today().isoformat())
    args = parser.parse_args()
    py = args.python
    env = os.environ.copy()
    env.update({
        "PYTHONHASHSEED": "0",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "GIT_CONFIG_COUNT": "1",
        "GIT_CONFIG_KEY_0": "core.filemode",
        "GIT_CONFIG_VALUE_0": "false",
        "ONCOGO_DATA_DIR": str(args.datos),
        "ONCOGO_CRC_DATA_WORKSPACE": str(args.salidas),
        "TMPDIR": env.get("TMPDIR", "/Volumes/JANCODATA/tmp"),
    })
    informe = RAIZ / "reports" / f"verificacion_local_{args.fecha}.md"
    hashes = args.salidas / "reports" / "hashes_entrada.tsv"
    comparacion = args.salidas / "reports" / "comparacion_salidas.tsv"
    pasos: list[tuple[str, int, float, str]] = []

    def paso(nombre: str, comando: list[str], cwd: Path = RAIZ) -> int:
        codigo, texto, duracion = _correr(comando, cwd, env)
        pasos.append((nombre, codigo, duracion, texto))
        print(f"[{'ok' if codigo == 0 else 'falla'}] {nombre} ({duracion:.1f} s)", flush=True)
        return codigo

    paso("hashes", [str(py), "scripts/verificar_hashes_entrada.py", "--datos", str(args.datos), "--salida", str(hashes)])
    env_pruebas = env.copy()
    env_pruebas.pop("ONCOGO_DATA_DIR", None)
    env_pruebas["ONCOGO_CRC_DATA_WORKSPACE"] = str(args.salidas / "pytest-workspace")
    env_pruebas["ONCOGO_E2E"] = "1"

    def paso_prueba(nombre: str, comando: list[str]) -> int:
        codigo, texto, duracion = _correr(comando, RAIZ, env_pruebas)
        pasos.append((nombre, codigo, duracion, texto))
        print(f"[{'ok' if codigo == 0 else 'falla'}] {nombre} ({duracion:.1f} s)", flush=True)
        return codigo

    paso_prueba("ruff", [str(RAIZ / ".venv/bin/ruff"), "check", "src", "tests", "scripts"])
    paso_prueba("pytest", [str(py), "-m", "pytest", "-q", "--cov=oncogo_crc", "--cov-report=term-missing"])
    paso_prueba("e2e", [str(py), "-m", "pytest", "-q", "-m", "e2e", "--override-ini=addopts=", "tests/test_e2e_v3_7_sintetico.py"])
    sucio = _git(["status", "--porcelain"], env)
    if sucio:
        pasos.append(("reproducir", 1, 0.0, "árbol con cambios de contenido; la etapa no firma"))
        print("[falla] reproducir (árbol sucio)", flush=True)
    else:
        paso("regresion-v1", [str(py), "-m", "pytest", "tests/test_regresion_v1.py", "-q", "--override-ini=addopts="])
        crudos = args.datos / "data/raw/plasma/pxd055218_raw"
        paso(
            "crudos-solo-lectura",
            [str(py), "scripts/crudos_solo_lectura.py", "--crudos", str(crudos), "--salida", str(args.salidas / "reports" / "crudos_sha1.tsv")],
        )
        paso("v3-7", [str(RAIZ / ".venv/bin/oncogo-crc"), "run", "--stage", "v3_7", "--data-dir", str(args.datos), "--data-workspace", str(args.salidas)])
    nueva = args.salidas / "results" / "v3_7"
    if nueva.is_dir():
        paso(
            "comparar",
            [str(py), "scripts/comparar_salidas.py", "--nueva", str(nueva), "--referencia", str(RAIZ / "results" / "v3_7"), "--salida", str(comparacion)],
        )
    else:
        pasos.append(("comparar", 1, 0.0, "no está la corrida nueva"))

    commit = _git(["rev-parse", "HEAD"], env)
    rama = _git(["branch", "--show-current"], env) or "(sin rama)"
    lineas = [
        f"# Verificación local {args.fecha}",
        "",
        "## Entorno",
        *_entorno(py),
        f"- commit: `{commit}`",
        f"- rama: `{rama}`",
        f"- datos (solo lectura): `{args.datos}`",
        f"- salidas: `{args.salidas}`",
        "- R: no hace falta. El único uso es `scripts/v2/session_info.py`, que anota la ausencia y no entra en v3.7.",
        "",
        "## Pasos",
        "",
        "| paso | resultado | segundos |",
        "| --- | --- | --- |",
    ]
    for nombre, codigo, duracion, _texto in pasos:
        lineas.append(f"| {nombre} | {'pasa' if codigo == 0 else 'falla'} | {duracion:.1f} |")
    lineas.extend(["", "## Hashes de entrada", "", *_tabla_hashes(hashes), "", "## Salidas que no son idénticas en bytes", ""])
    lineas.extend(_tabla_salidas(comparacion))
    lineas.extend(["", "## Notas de los comandos", ""])
    for nombre, codigo, duracion, texto in pasos:
        if codigo == 0 and nombre not in {"pytest", "e2e", "comparar", "hashes"}:
            continue
        lineas.append(f"### {nombre}")
        lineas.append("")
        lineas.append("```")
        lineas.append(texto.strip() or "(sin salida)")
        lineas.append("```")
        lineas.append("")
    informe.parent.mkdir(parents=True, exist_ok=True)
    informe.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    print(informe)
    return 0 if all(codigo == 0 for _n, codigo, _d, _t in pasos) else 1


if __name__ == "__main__":
    raise SystemExit(main())
