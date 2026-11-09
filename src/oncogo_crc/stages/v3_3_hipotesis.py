"""H1 a H10. Lo que no tiene fichero queda no_evaluable o requiere_acceso."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.stages.v3_3_controles import _abrir
from oncogo_crc.v3.h5 import ejecutar as ejecutar_h5
from oncogo_crc.v3.h7 import ejecutar as ejecutar_h7
from oncogo_crc.v3.h10 import ejecutar as ejecutar_h10
from oncogo_crc.v3.salida import volcar_tsv


def _fila(codigo: str, nivel: str, estado: str, metrica: str, detalle: str, fuente: str, fraccion: str = "") -> dict:
    return {
        "hipotesis": codigo,
        "nivel": nivel,
        "estado": estado,
        "metrica": metrica,
        "detalle": detalle,
        "umbral": "",
        "fuente": fuente,
        "fraccion_multiverso": fraccion,
    }


def _heredadas(path: Path) -> list[dict]:
    if not path.is_file():
        return [
            _fila(codigo, "heredada", "no_evaluable", "", "no está la tabla de v3", "results/v3/hipotesis_v3.tsv")
            for codigo in ("H2", "H3", "H4")
        ]
    vieja = pd.read_csv(path, sep="\t")
    filas = []
    for codigo in ("H2", "H3", "H4"):
        bloque = vieja.loc[vieja["hipotesis"].eq(codigo)]
        if bloque.empty:
            filas.append(_fila(codigo, "heredada", "no_evaluable", "", "sin fila en v3", str(path)))
            continue
        estado = str(bloque["estado"].iloc[0])
        if estado == "no_evaluada":
            estado = "no_evaluable"
        filas.append(_fila(
            codigo,
            "heredada",
            estado,
            str(bloque["metrica"].iloc[0]),
            f"No se recalcula en v3.2. v3 decía: {bloque['detalle'].iloc[0]}",
            str(bloque["fuente"].iloc[0]),
        ))
    return filas


def _buscar(settings: Settings, nombre: str) -> Path | None:
    raiz = settings.resolve("external")
    hallados = list(raiz.rglob(nombre)) if raiz.is_dir() else []
    return hallados[0] if hallados else None


def run(settings: Settings) -> None:
    _prereg, salida, sello = _abrir(settings, "v3_3_hipotesis")
    k1 = salida / "k1_etiquetas.json"
    if not k1.is_file():
        raise FileNotFoundError("falta k1_etiquetas.json; corre v3-2-tejido antes")
    tejido = json.loads(k1.read_text(encoding="utf-8"))
    r1 = tejido["contrastes"]["primario_frente_a_mucosa"]
    filas = [_fila(
        "H1",
        "heredada",
        tejido["h1_veredicto"],
        "concordancia de signo, primario frente a mucosa",
        f"{r1['numerador']}/{r1['denominador']}",
        "results/v3_2/k1_etiquetas.json",
    )]
    filas.extend(_heredadas(settings.resolve("results/v3/hipotesis_v3.tsv")))
    print("[v3.2] H5: descarga, si falta, y treat por carcinoma")
    h5 = ejecutar_h5(settings, salida, sello, _prereg)
    filas.append(_fila(
        "H5",
        "secundaria",
        h5["estado"],
        "fraccion pan_carcinoma",
        str(h5.get("fraccion_pan_carcinoma", h5.get("motivo"))),
        "results/v3_2/h5_especificidad_organo.json",
    ))
    referencia = settings.resolve("external/referencia_unicelular")
    matrices = list(referencia.glob("*.tsv")) + list(referencia.glob("*.txt")) if referencia.is_dir() else []
    filas.append(_fila(
        "H6",
        "exploratoria",
        "no_disponible" if not matrices else "no_evaluable",
        "origen celular",
        "No hay una matriz de célula única con leyenda verificada en external/referencia_unicelular." if not matrices else "La matriz está, pero este corte no ajusta la interacción.",
        str(referencia),
    ))
    print("[v3.3] H7: modelo por proteína, fuera de recurso")
    h7 = ejecutar_h7(settings, salida, sello, _prereg)
    filas.append(_fila(
        "H7",
        "confirmatoria",
        h7["estado"],
        "AUC fuera de recurso y de proteína",
        str(h7.get("por_recurso", h7.get("motivo"))),
        "results/v3_3/h7_visibilidad_en_sangre.json",
    ))
    pxd = _buscar(settings, "*PXD064655*")
    filas.append(_fila(
        "H8",
        "secundaria",
        "no_disponible",
        "treat estadio IV frente a I-II",
        "PXD064655 no está en el disco." if pxd is None else f"El fichero {pxd} está, pero la leyenda no se reabrió en esta etapa.",
        "external",
    ))
    ukb = settings.resolve("external/ukb_ppp/estado.json")
    motivo_h9 = "requiere cuenta de Synapse"
    if ukb.is_file():
        motivo_h9 = str(json.loads(ukb.read_text(encoding="utf-8")).get("motivo"))
    filas.append(_fila("H9", "exploratoria", "requiere_acceso", "Wald o IVW con instrumentos cis", motivo_h9, str(ukb)))
    print("[v3.2] H10: atlas de enfermedades del HPA")
    h10 = ejecutar_h10(settings, salida, sello, _prereg)
    filas.append(_fila(
        "H10",
        "secundaria",
        "calculado" if h10["estado"] == "calculado" else h10["estado"],
        "proteínas específicas frente a sanos, otros cánceres y otras enfermedades",
        str(h10.get("n_especificas_crc", h10.get("motivo"))),
        "results/v3_2/h10_especificidad.json",
    ))
    # El estado de H10 en la tabla de hipótesis usa el vocabulario del preregistro.
    if filas[-1]["estado"] == "calculado":
        filas[-1]["estado"] = "no_evaluable"
        filas[-1]["detalle"] = (
            "El preregistro no fija un umbral para sostener o refutar H10. "
            f"El atlas da {h10.get('n_especificas')} específicas, "
            f"{h10.get('n_no_especificas')} no específicas y {h10.get('n_no_medidas')} no medidas. "
            "Fase de Pepe 1."
        )
    volcar_tsv(pd.DataFrame(filas), salida / "hipotesis_v3_2.tsv", sello)
    print("[v3.2] hipotesis_v3_2.tsv escrita")
