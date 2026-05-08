from __future__ import annotations

import argparse
import sys
from pathlib import Path

from oncogo_crc.pipeline.registry import STAGE_ORDER, STAGES
from oncogo_crc.settings import load_settings


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="oncogo-crc",
        description="OncoGO-CRC — pipeline proteómico protein-centric con GO e IA explicable.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    runp = sub.add_parser("run", help="Ejecuta una o todas las etapas del pipeline.")
    runp.add_argument(
        "--stage",
        choices=STAGE_ORDER + ["all"],
        default="all",
        help="Etapa concreta o `all` en orden (download→…→report).",
    )
    runp.add_argument(
        "--demo",
        action="store_true",
        help="Genera datos sintéticos y metadatos runtime para validar el flujo sin CPTAC/PRIDE.",
    )
    runp.add_argument("--root", type=Path, default=None, help="Raíz del repositorio (por defecto: auto-detectada).")

    args = parser.parse_args(argv)
    if args.command != "run":
        parser.print_help()
        return 2

    settings = load_settings(args.root, demo=args.demo)
    order = STAGE_ORDER if args.stage == "all" else [args.stage]
    for name in order:
        print(f"\n=== [{name}] ===")
        STAGES[name](settings)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
