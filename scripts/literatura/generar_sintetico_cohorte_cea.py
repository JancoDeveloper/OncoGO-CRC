#!/usr/bin/env python3
"""Genera la cohorte sintética de D9: ``tests/SINTETICO_cohorte_cea.csv``.

Estos datos no son de nadie. Se generan con una semilla fija para que el fichero
sea reproducible, y el prefijo ``SINTETICO_`` y la carpeta ``tests/`` están ahí
para que nunca se confundan con datos reales. Son lo único que permite comprobar
que el análisis del protocolo funciona antes de que exista el estudio.

El mecanismo que se simula es el de la hipótesis que el protocolo quiere probar:
quien va a recaer empieza a subir el CEA de forma multiplicativa unos meses antes
de que la recaída sea visible, mientras que quien no va a recaer fluctúa alrededor
de su propio valor. Se añaden a propósito tres cosas que estropean el análisis en
la vida real y que el código tiene que saber tragar:

- Un fumador con el CEA basal alto y sin recaída, que es el falso positivo
  clásico del CEA.
- Pacientes censurados, que dejan el seguimiento antes de la ventana de doce
  meses.
- Visitas que faltan, porque la gente no acude a todas.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parents[2]
DESTINO = RAIZ / "tests" / "SINTETICO_cohorte_cea.csv"

# Derivada de la semilla base 20260914 con el desplazamiento de D9.
SEMILLA = 20260914 + 909

VISITAS = (3.0, 6.0, 9.0, 12.0, 18.0, 24.0, 30.0, 36.0)


def generar(n_pacientes: int = 300, semilla: int = SEMILLA) -> pd.DataFrame:
    rng = np.random.default_rng(semilla)
    filas = []
    for i in range(n_pacientes):
        paciente = f"SINT{i:04d}"
        estadio = rng.choice(["II", "III"], p=[0.45, 0.55])
        fumador = rng.random() < 0.2
        # Cada paciente tiene su propio valor de partida, que es la razón de que
        # un umbral común funcione mal.
        basal = rng.lognormal(mean=np.log(2.0) + (0.6 if fumador else 0.0), sigma=0.45)

        recae = rng.random() < (0.18 if estadio == "II" else 0.32)
        mes_recaida = float(rng.uniform(9, 34)) if recae else np.nan
        # Censura: una parte deja el seguimiento antes del final.
        ultimo = 36.0 if rng.random() > 0.18 else float(rng.uniform(12, 33))
        if recae:
            ultimo = min(ultimo, mes_recaida + 1.0)

        for v in VISITAS:
            if v > ultimo:
                break
            if rng.random() < 0.12:  # visita perdida
                continue
            nivel = np.log(basal)
            if recae:
                # Sube de forma multiplicativa en los seis meses previos.
                adelanto = mes_recaida - v
                if adelanto < 6.0:
                    nivel += 0.35 * (6.0 - max(adelanto, -2.0))
            cea = float(np.exp(nivel + rng.normal(0, 0.22)))
            filas.append(
                {
                    "paciente": paciente,
                    "visita_meses": v,
                    "cea": round(cea, 2),
                    "estadio": estadio,
                    "fumador": int(fumador),
                    "recaida": int(recae),
                    "mes_de_la_recaida": (
                        round(mes_recaida, 1) if recae else ""
                    ),
                    "ultimo_seguimiento_meses": round(ultimo, 1),
                }
            )
    return pd.DataFrame(filas)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--destino", default=str(DESTINO))
    ap.add_argument("--n", type=int, default=300)
    args = ap.parse_args()

    d = generar(args.n)
    destino = Path(args.destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(destino, index=False)
    print(
        f"[sintetico] {len(d)} visitas de {d['paciente'].nunique()} pacientes "
        f"sintéticos en {destino}"
    )
    print(
        f"[sintetico] con recaída: "
        f"{d.groupby('paciente')['recaida'].max().sum()} pacientes"
    )
    print("[sintetico] estos datos no son de nadie y solo valen para los tests")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
