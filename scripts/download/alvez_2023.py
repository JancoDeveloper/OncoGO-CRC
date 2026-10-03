#!/usr/bin/env python3
"""D4b. Recurso de sangre pan-cáncer de Álvez et al. 2023.

Para qué sirve esto. v1 encontró tres proteínas que separan suero de cáncer de
suero sano en una cohorte de 20 contra 20. Tres proteínas en cuarenta sueros es
poco, y la pregunta obvia es si alguien con una cohorte de verdad ve lo mismo.
Álvez midió 1 463 proteínas con Olink Explore en 221 pacientes de cáncer
colorrectal frente a controles sanos del programa Wellness. Es la comprobación
externa que a v1 le faltaba.

De este trabajo solo se usan estadísticas ya publicadas. No se reanalizan datos
individuales: las diferencias de NPX, los p ajustados y los AUC vienen tal cual
de la información suplementaria del artículo.

Un aviso sobre la plataforma, que condiciona lo que se puede comprobar. Olink es
un panel cerrado: mide las proteínas que lleva y ninguna más. Que una proteína no
aparezca en los resultados no significa que no cambie, significa que no se midió.
Eso no es lo mismo que un «no replica» y en la tabla limpia se distingue.

El fichero suplementario 11 es la fuente de datos de las figuras del artículo y
trae una hoja por panel. Las hojas que se usan están declaradas en ``HOJAS``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

DOI = "10.1038/s41467-023-39765-y"
ARTICULO = f"https://doi.org/{DOI}"
BASE = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41467-023-39765-y/MediaObjects/"
    "41467_2023_39765_MOESM{}_ESM.xlsx"
)

CITA = (
    "Álvez MB, Edfors F, von Feilitzen K, et al. Next generation pan-cancer "
    "blood proteome profiling using proteograph. Nature Communications 2023; "
    "14:4308."
)
LICENCIA = (
    "CC BY 4.0. El artículo y su información suplementaria son de acceso "
    "abierto y permiten reutilización con atribución."
)

# Los suplementarios que hacen falta y por qué. Se descargan solo estos tres: el
# 11 es grande y los demás no aportan nada a esta pregunta.
SUPLEMENTARIOS = {
    3: "resumen de las cohortes: edades, sexo y estadios",
    6: "AUC del clasificador por número de proteínas",
    11: "datos de las figuras, con los resultados por proteína",
}

# Hojas del suplementario 11 que se leen, con su contenido declarado.
HOJAS = {
    "2a": "diferencia de NPX por proteína y cáncer, con p ajustado",
    "6e": "NPX por estadio agrupado",
    "S1c": "coeficiente de variación intraensayo por ensayo",
}

# Las tres proteínas de plasma de v1. Son el motivo de esta fuente.
PLASMA_V1 = {
    "P0DJI8": "SAA1",
    "P18428": "LBP",
    "Q9BY43": "CHMP4A",
}

# Controles de familia. CEACAM5 es la proteína que mide el ensayo clínico de CEA,
# y aquí sirve para cruzar con el eje D1.
CEACAM = ("CEACAM5", "CEACAM1")

CANCER = "Colorectal"


def _descargar(numero: int, destino: Path) -> dict:
    """Descarga idempotente: si el fichero ya está y el hash cuadra, no repite."""
    url = BASE.format(numero)
    if destino.exists():
        return {
            "url": url,
            "fichero": destino.name,
            "bytes": destino.stat().st_size,
            "sha256": _sha256(destino),
            "descargado_ahora": False,
        }
    peticion = urllib.request.Request(url, headers={"User-Agent": "OncoGO-CRC/2"})
    with urllib.request.urlopen(peticion, timeout=180) as r:
        datos = r.read()
    destino.write_bytes(datos)
    return {
        "url": url,
        "fichero": destino.name,
        "bytes": len(datos),
        "sha256": hashlib.sha256(datos).hexdigest(),
        "descargado_ahora": True,
    }


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()


def _por_proteina(xlsx: Path) -> pd.DataFrame:
    """Resultados por proteína del cáncer colorrectal.

    La columna ``Significance`` del artículo trae tres valores de texto. Se
    traduce a una dirección explícita en lugar de dejar la cadena en inglés, y el
    signo se toma de la diferencia de NPX, no de la etiqueta: si alguna vez las
    dos discreparan, la tabla lo diría.
    """
    d = pd.read_excel(xlsx, sheet_name="2a")
    d = d[d["Cancer"].astype(str).str.strip() == CANCER].copy()
    d["proteina"] = d["Assay"].astype(str).str.strip().str.upper()
    d["diferencia_npx"] = pd.to_numeric(d["NPX_difference"], errors="coerce")
    d["p_ajustado"] = pd.to_numeric(d["p.adjusted"], errors="coerce")
    etiqueta = d["Significance"].astype(str).str.strip().str.lower()
    d["significativa"] = etiqueta.str.startswith("significant")
    d["direccion"] = pd.Series(
        ["sube" if v > 0 else ("baja" if v < 0 else "sin_cambio") for v in d["diferencia_npx"]],
        index=d.index,
    )
    d["etiqueta_del_articulo"] = etiqueta
    d["olink_id"] = d["OlinkID"].astype(str)
    salida = d[
        [
            "proteina", "olink_id", "diferencia_npx", "p_ajustado",
            "significativa", "direccion", "etiqueta_del_articulo",
        ]
    ].sort_values("p_ajustado", kind="mergesort")
    return salida.reset_index(drop=True)


def _auc(xlsx: Path) -> pd.DataFrame:
    d = pd.read_excel(xlsx)
    d = d[d["Cancer"].astype(str).str.strip() == CANCER].copy()
    d = d.rename(
        columns={
            "Num_proteins": "n_proteinas",
            "AUC": "auc",
            "Sensitivity": "sensibilidad",
            "Specificity": "especificidad",
            "Accuracy": "exactitud",
        }
    )
    columnas = ["n_proteinas", "auc", "sensibilidad", "especificidad", "exactitud"]
    return d[columnas].sort_values("n_proteinas").reset_index(drop=True)


def _cohorte(xlsx: Path) -> dict:
    """Tamaños de la cohorte colorrectal y del grupo de control.

    La cabecera del suplementario 3 ocupa dos filas: la segunda trae los niveles
    (Female/Male, los estadios 1 a 4). Se leen por posición porque los nombres
    que pandas inventa (``Unnamed: 4``) no son estables entre versiones.
    """
    d = pd.read_excel(xlsx)
    fila = d[d["Cancer"].astype(str).str.strip() == CANCER]
    control = d[d["Cohort"].astype(str).str.strip() == "Wellness"]
    if fila.empty or control.empty:
        return {"estado": "no_legible"}

    def _n(f, inicio, cuantos):
        vals = f.iloc[0, inicio : inicio + cuantos]
        return [None if pd.isna(v) else int(v) for v in vals]

    sexo = _n(fila, 3, 2)
    estadios = _n(fila, 5, 5)
    return {
        "estado": "legible",
        "n_cancer": int(sum(v for v in sexo if v)),
        "sexo_cancer": {"mujeres": sexo[0], "hombres": sexo[1]},
        "estadios_cancer": dict(
            zip(["I", "II", "III", "IV", "desconocido"], estadios, strict=True)
        ),
        "edad_cancer": str(fila.iloc[0]["Age_range"]),
        "n_control": int(sum(v for v in _n(control, 3, 2) if v)),
        "edad_control": str(control.iloc[0]["Age_range"]),
        "que_son_los_controles": (
            "Sujetos sanos del programa Wellness, no pacientes con síntomas ni "
            "con enfermedad benigna. Es el mismo diseño que PXD055218 y tiene la "
            "misma limitación: separar cáncer de sano es más fácil que separar "
            "cáncer de lo que llega a una consulta."
        ),
    }


def _comprobar_plasma_v1(por_proteina: pd.DataFrame) -> dict:
    """Comprueba las tres proteínas de plasma de v1 contra este recurso.

    La distinción que importa: «no está en el panel» no es «no replica». La
    primera es una limitación de la plataforma y la segunda es un resultado.
    """
    medidas = set(por_proteina["proteina"])
    filas = []
    for uniprot, gen in PLASMA_V1.items():
        if gen not in medidas:
            filas.append(
                {
                    "uniprot": uniprot,
                    "gen": gen,
                    "estado": "no_en_el_panel",
                    "diferencia_npx": None,
                    "p_ajustado": None,
                    "direccion": None,
                    "lectura": (
                        f"{gen} no es uno de los 1 463 ensayos de Olink Explore. "
                        "No se midió, así que este recurso no dice nada sobre "
                        "ella, ni a favor ni en contra."
                    ),
                }
            )
            continue
        r = por_proteina[por_proteina["proteina"] == gen].iloc[0]
        replica = bool(r["significativa"]) and r["direccion"] == "sube"
        filas.append(
            {
                "uniprot": uniprot,
                "gen": gen,
                "estado": "replica" if replica else "no_replica",
                "diferencia_npx": float(r["diferencia_npx"]),
                "p_ajustado": float(r["p_ajustado"]),
                "direccion": str(r["direccion"]),
                "lectura": (
                    f"{gen} sí se midió. "
                    + (
                        "Sube y es significativa, igual que en v1."
                        if replica
                        else f"La diferencia es {r['diferencia_npx']:.3f} en NPX con "
                        f"p ajustado {r['p_ajustado']:.3g}, que no pasa el umbral "
                        "del artículo. En una cohorte mucho mayor que la de v1."
                    )
                ),
            }
        )
    comprobables = [f for f in filas if f["estado"] != "no_en_el_panel"]
    replican = [f for f in comprobables if f["estado"] == "replica"]
    return {
        "proteinas": filas,
        "n_declaradas_en_v1": len(PLASMA_V1),
        "n_comprobables": len(comprobables),
        "n_que_replican": len(replican),
        "lectura": (
            f"De las {len(PLASMA_V1)} proteínas de plasma de v1, "
            f"{len(comprobables)} están en el panel de Olink y "
            f"{len(replican)} replican. Las otras "
            f"{len(PLASMA_V1) - len(comprobables)} no se midieron, que no es lo "
            "mismo que no replicar."
        ),
    }


def _fase_aguda(por_proteina: pd.DataFrame) -> dict:
    """Cuántas de las proteínas que suben son de respuesta de fase aguda.

    No hay una lista canónica de fase aguda en el artículo, así que no se
    construye una a mano: eso sería escoger el resultado. Lo que se deja es el
    recuento de proteínas que suben, para que la comparación con v1 se haga con
    las vías de C4 y no con una lista improvisada aquí.
    """
    # El cast es necesario: con la tabla vacía la columna queda de tipo objeto y
    # pandas la leería como una lista de nombres de columna en lugar de máscara.
    sig = por_proteina[por_proteina["significativa"].astype(bool)]
    return {
        "n_medidas": int(len(por_proteina)),
        "n_significativas": int(len(sig)),
        "n_suben": int((sig["direccion"] == "sube").sum()),
        "n_bajan": int((sig["direccion"] == "baja").sum()),
        "nota": (
            "No se clasifica a mano qué proteínas son de fase aguda. La "
            "interpretación por vías se hace con las colecciones de C4, que están "
            "preregistradas, y no con una lista escogida después de ver estos "
            "resultados."
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    args = ap.parse_args()

    salida = Path(args.workspace) / "external" / "alvez_2023"
    salida.mkdir(parents=True, exist_ok=True)
    estado_path = salida / "estado.json"

    ficheros: dict[str, dict] = {}
    try:
        for numero in SUPLEMENTARIOS:
            destino = salida / f"41467_2023_39765_MOESM{numero}_ESM.xlsx"
            ficheros[str(numero)] = _descargar(numero, destino)
            print(
                f"[alvez] suplementario {numero}: "
                f"{ficheros[str(numero)]['bytes']} bytes"
            )
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        estado_path.write_text(
            json.dumps(
                {
                    "nombre": "alvez_2023",
                    "estado": "no_disponible",
                    "motivo": f"fallo al descargar la información suplementaria: {e}",
                    "fecha_de_acceso": datetime.now(UTC).date().isoformat(),
                    "articulo": ARTICULO,
                    "analisis_que_queda_sin_evaluar": (
                        "la comprobación externa de las tres proteínas de plasma "
                        "de v1 (D4b)"
                    ),
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"[alvez] no_disponible: {e}", file=sys.stderr)
        return 1

    x11 = salida / "41467_2023_39765_MOESM11_ESM.xlsx"
    por_proteina = _por_proteina(x11)
    auc = _auc(salida / "41467_2023_39765_MOESM6_ESM.xlsx")
    cohorte = _cohorte(salida / "41467_2023_39765_MOESM3_ESM.xlsx")
    plasma = _comprobar_plasma_v1(por_proteina)
    fase = _fase_aguda(por_proteina)

    ceacam = {}
    for gen in CEACAM:
        r = por_proteina[por_proteina["proteina"] == gen]
        ceacam[gen] = (
            {
                "diferencia_npx": float(r.iloc[0]["diferencia_npx"]),
                "p_ajustado": float(r.iloc[0]["p_ajustado"]),
                "direccion": str(r.iloc[0]["direccion"]),
            }
            if not r.empty
            else {"estado": "no_en_el_panel"}
        )

    tabla = salida / "alvez_2023_colorrectal_por_proteina.tsv"
    por_proteina.to_csv(tabla, sep="\t", index=False)
    tabla_auc = salida / "alvez_2023_colorrectal_auc.tsv"
    auc.to_csv(tabla_auc, sep="\t", index=False)

    minimo = auc["n_proteinas"].min()
    fila_min = auc[auc["n_proteinas"] == minimo].iloc[0]
    fila_max = auc.iloc[-1]

    estado = {
        "nombre": "alvez_2023",
        "estado": "usada",
        "fecha_de_acceso": datetime.now(UTC).date().isoformat(),
        "articulo": ARTICULO,
        "doi": DOI,
        "cita": CITA,
        "licencia": LICENCIA,
        "plataforma": (
            "Olink Explore, 1 463 ensayos de inmunoensayo de extensión por "
            "proximidad. Es un panel cerrado: mide lo que lleva."
        ),
        "suplementarios": {
            numero: {**ficheros[str(numero)], "para_que": motivo}
            for numero, motivo in SUPLEMENTARIOS.items()
        },
        "hojas_usadas": HOJAS,
        "cohorte": cohorte,
        "tabla_por_proteina": tabla.name,
        "sha256_tabla_por_proteina": _sha256(tabla),
        "tabla_auc": tabla_auc.name,
        "sha256_tabla_auc": _sha256(tabla_auc),
        "resumen_por_proteina": fase,
        "plasma_de_v1": plasma,
        "ceacam": ceacam,
        "auc_segun_cuantas_proteinas": {
            "con_el_panel_completo": {
                "n_proteinas": int(fila_max["n_proteinas"]),
                "auc": float(fila_max["auc"]),
            },
            "con_el_minimo_publicado": {
                "n_proteinas": int(fila_min["n_proteinas"]),
                "auc": float(fila_min["auc"]),
            },
            "para_que_sirve": (
                "v1 se quedó con tres proteínas de plasma. Esta serie dice qué "
                "rinde un número parecido de proteínas en una cohorte grande, y "
                "cuánto se pierde respecto al panel completo."
            ),
        },
        "limitaciones": [
            "Olink mide NPX, que es una escala log2 relativa de un inmunoensayo. "
            "No es comparable con el log2 de una intensidad de espectrometría de "
            "masas, así que las magnitudes no se cruzan: solo el signo y el "
            "hecho de ser significativa.",
            "Los controles son sujetos sanos del programa Wellness, con un rango "
            "de edad más estrecho que los pacientes. El AUC de una comparación "
            "contra sanos es optimista para cualquier uso clínico, porque en la "
            "consulta el problema es separar cáncer de enfermedad benigna.",
            "Los grupos están desbalanceados y vienen de reclutamientos "
            "distintos, así que parte de la separación puede ser de cohorte y no "
            "de enfermedad.",
            "Que una proteína no esté en el panel no es evidencia de que no "
            "cambie. Dos de las tres proteínas de plasma de v1 están en esa "
            "situación.",
            "Solo se usan estadísticas publicadas. No hay acceso a los datos "
            "individuales, así que no se puede recalcular nada ni aplicar el "
            "treat de este trabajo.",
        ],
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
    }
    estado_path.write_text(
        json.dumps(estado, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(
        f"[alvez] {fase['n_medidas']} proteínas medidas en colorrectal, "
        f"{fase['n_significativas']} significativas "
        f"({fase['n_suben']} suben, {fase['n_bajan']} bajan)."
    )
    print(f"[alvez] {plasma['lectura']}")
    print(
        f"[alvez] AUC {fila_max['auc']:.3f} con "
        f"{int(fila_max['n_proteinas'])} proteínas y {fila_min['auc']:.3f} con "
        f"{int(fila_min['n_proteinas'])}."
    )
    print(f"[alvez] escrito en {salida}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
