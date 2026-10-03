#!/usr/bin/env python3
"""D7. UK Biobank Pharma Proteomics Project: pQTL del plasma.

Esta fuente **requiere acceso** y el script no intenta obtenerlo. Las
estadísticas de resumen viven en Synapse (syn51364943) bajo los términos del
UKB-PPP, que exigen cuenta y aceptación de condiciones. La regla de este trabajo
es que una fuente con cuenta se deja preparada y marcada ``requiere_acceso``, no
se fuerza. Así que aquí solo se comprueba qué páginas son públicas, se deja el
análisis escrito y probado, y se sigue.

Para qué serviría. C8 encontró que la dirección del cambio en tejido no se
traslada al plasma, y hasta ahora hay dos explicaciones sobre la mesa: la ventana
de detección, que es lo que mide C8, y la elección del grupo de control, que es
lo que mide D6. Un pQTL aporta una tercera, y es la más incómoda. Si el nivel de
una proteína en plasma lo determina en buena parte la variante genética que cada
persona hereda, entonces las diferencias entre dos individuos no hablan de la
enfermedad: hablan del genotipo. Una proteína con un cis-pQTL fuerte es una mala
candidata a marcador por razones que no tienen nada que ver con el cáncer.

Las dos proteínas de plasma de v1 que son de fase aguda, SAA1 y LBP, son
precisamente de las que se sabe que tienen variantes de efecto grande. Si el día
que haya acceso se confirma, el resultado de v1 en plasma queda explicado por
tres cosas distintas y ninguna es el tumor.

El análisis está escrito y tiene tests sobre datos sintéticos. El día que alguien
con acceso descargue el fichero de Synapse, esto corre sin tocar nada.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

ARTICULO = "https://www.nature.com/articles/s41586-023-06592-6"
PORTAL = "http://ukb-ppp.gwas.eu/"
SYNAPSE = "https://www.synapse.org/Synapse:syn51364943"
SYNAPSE_ID = "syn51364943"
SYNAPSE_DOI = "10.7303/syn51364943"

CITA = (
    "Sun BB, Chiou J, Traylor M, et al. Plasma proteomic associations with "
    "genetics and health in the UK Biobank. Nature 2023; 622:329-338. "
    "doi:10.1038/s41586-023-06592-6"
)
LICENCIA = (
    "Las estadísticas de resumen están sujetas a los términos del UKB-PPP y de "
    "UK Biobank. Hace falta cuenta de Synapse y aceptar las condiciones. Este "
    "script no las acepta en nombre de nadie."
)

# Proteínas que se consultarían, y por qué cada una. La lista se fija aquí, antes
# de tener los datos, para que no se elija después de verlos.
OBJETIVO = {
    "SAA1": "proteína de plasma de v1; fase aguda; se le conocen variantes fuertes",
    "LBP": "proteína de plasma de v1; fase aguda; está en el panel de Olink de D4b",
    "CHMP4A": "proteína de plasma de v1; no está en el panel de Olink",
    "CEACAM5": "el marcador clínico del eje D1",
    "S100A12": "la única inflamatoria de las 22 de D6",
    "KRT19": "citoqueratina de D6; control de una proteína de recambio celular",
    "PARP1": "nuclear de D6; control de una proteína que no debería estar en plasma",
}

# Umbral de significación del artículo para un pQTL primario. Se declara en lugar
# de reinventarse: el objetivo es leer su resultado, no rehacer su GWAS.
P_UMBRAL = 1.7e-11

# A partir de qué fracción de varianza explicada se considera que el nivel en
# plasma está determinado en buena parte por la genética. No hay un valor
# canónico, así que se fija aquí y se justifica: con un 5 % de la varianza, dos
# personas de genotipos opuestos difieren lo bastante para mover un marcador
# usado con un umbral fijo.
UMBRAL_VARIANZA = 0.05

# Columnas que se esperan del fichero de Synapse. Se declaran para que, si el
# formato cambia, el script lo diga en lugar de leer una columna equivocada.
COLUMNAS = {
    "protein": "nombre del ensayo de Olink",
    "variant": "identificador de la variante",
    "cis_trans": "si la variante está en el locus del gen (cis) o fuera (trans)",
    "beta": "efecto sobre el nivel en NPX por copia del alelo",
    "se": "error estándar del efecto",
    "log10p": "menos el logaritmo en base 10 del p",
    "maf": "frecuencia del alelo menos común",
}


def varianza_explicada(beta: float, maf: float) -> float:
    """Fracción de la varianza del nivel en plasma que explica una variante.

    Para un nivel estandarizado, la varianza que aporta una variante bialélica en
    equilibrio es ``2·p·(1−p)·beta²``. El NPX del UKB-PPP viene estandarizado por
    rango inverso normal, así que la fórmula aplica directamente y no hace falta
    la varianza del fenotipo.

    La fórmula es una aproximación y lo es por una razón concreta: supone que el
    efecto es aditivo y que la variante está en equilibrio. Para un cis-pQTL de
    efecto grande las dos cosas se cumplen razonablemente; para una variante rara
    el resultado se infla y por eso se declara la frecuencia al lado.
    """
    if not np.isfinite(beta) or not np.isfinite(maf) or not 0 < maf < 1:
        return float("nan")
    return float(2.0 * maf * (1.0 - maf) * beta**2)


def resumir_por_proteina(d: pd.DataFrame) -> pd.DataFrame:
    """Un resumen por proteína: su cis-pQTL más fuerte y lo que explica.

    Se queda con el cis-pQTL de mayor varianza explicada y no con el de menor p.
    No son lo mismo: una variante muy común con efecto pequeño puede dar un p
    minúsculo en 54 000 personas y explicar casi nada. Lo que importa aquí es
    cuánto mueve el nivel, no cuán seguro está el artículo de que lo mueve.
    """
    faltan = sorted(set(COLUMNAS) - set(d.columns))
    if faltan:
        raise KeyError(
            f"el fichero de Synapse no trae las columnas {faltan}; el formato "
            f"ha cambiado y hay que revisar el parser antes de usar los números"
        )
    t = d.copy()
    t["proteina"] = t["protein"].astype(str).str.upper()
    t["varianza_explicada"] = [
        varianza_explicada(b, m)
        for b, m in zip(t["beta"], t["maf"], strict=True)
    ]
    t = t[t["log10p"] >= -np.log10(P_UMBRAL)]
    cis = t[t["cis_trans"].astype(str).str.lower() == "cis"]

    filas = []
    for proteina, g in cis.groupby("proteina"):
        mejor = g.loc[g["varianza_explicada"].idxmax()]
        filas.append(
            {
                "proteina": proteina,
                "variante": str(mejor["variant"]),
                "beta": float(mejor["beta"]),
                "se": float(mejor["se"]),
                "maf": float(mejor["maf"]),
                "log10p": float(mejor["log10p"]),
                "varianza_explicada": float(mejor["varianza_explicada"]),
                "n_cis_pqtl": int(len(g)),
                "n_trans_pqtl": int(
                    (
                        (t["proteina"] == proteina)
                        & (t["cis_trans"].astype(str).str.lower() == "trans")
                    ).sum()
                ),
                "determinada_por_genetica": bool(
                    mejor["varianza_explicada"] >= UMBRAL_VARIANZA
                ),
            }
        )
    return pd.DataFrame(
        filas,
        columns=[
            "proteina", "variante", "beta", "se", "maf", "log10p",
            "varianza_explicada", "n_cis_pqtl", "n_trans_pqtl",
            "determinada_por_genetica",
        ],
    ).sort_values("varianza_explicada", ascending=False).reset_index(drop=True)


def lectura(resumen: pd.DataFrame) -> dict:
    """Qué dice el resumen sobre las proteínas que importan a este trabajo.

    La ausencia de una proteína en el resumen no es que no tenga pQTL: puede ser
    que no esté en el panel de 2 923 del UKB-PPP. Igual que en D4b, no medido y
    no encontrado no son lo mismo, y aquí también se distinguen.
    """
    presentes = set(resumen["proteina"])
    por_proteina = {}
    for gen, motivo in OBJETIVO.items():
        if gen not in presentes:
            por_proteina[gen] = {
                "estado": "no_en_el_panel_o_sin_cis_pqtl",
                "motivo_de_interes": motivo,
                "aviso": (
                    "No aparece con cis-pQTL significativo. Puede ser que no esté "
                    "entre las 2 923 proteínas del UKB-PPP, o que no tenga uno. "
                    "Las dos cosas se distinguen con la lista de ensayos, no con "
                    "este fichero."
                ),
            }
            continue
        r = resumen[resumen["proteina"] == gen].iloc[0]
        por_proteina[gen] = {
            "estado": "con_cis_pqtl",
            "motivo_de_interes": motivo,
            "variante": r["variante"],
            "varianza_explicada": round(float(r["varianza_explicada"]), 4),
            "determinada_por_genetica": bool(r["determinada_por_genetica"]),
        }
    determinadas = [
        g for g, v in por_proteina.items() if v.get("determinada_por_genetica")
    ]
    return {
        "umbral_de_varianza": UMBRAL_VARIANZA,
        "por_proteina": por_proteina,
        "n_determinadas_por_genetica": len(determinadas),
        "determinadas_por_genetica": sorted(determinadas),
        "que_significa": (
            "Una proteína cuyo nivel en plasma depende en buena parte de la "
            "variante que cada persona hereda es una mala candidata a marcador, "
            "y lo es por una razón que no tiene nada que ver con el cáncer: dos "
            "personas sanas de genotipos opuestos ya difieren entre sí."
        ),
        "lo_que_no_dice": (
            "Un cis-pQTL no invalida un marcador. Dice que el umbral fijo no "
            "vale igual para todo el mundo, que es otra cosa."
        ),
    }


def _publico(url: str) -> dict:
    """Comprueba si una página es pública. Sin credenciales y sin intentarlas.

    El estado tiene tres valores, no dos, y la distinción importa. Que el
    servidor conteste 401 o 403 es una restricción de acceso. Que falle el
    certificado o se agote el tiempo es un problema de esta máquina, y decir que
    la página «no es accesible sin cuenta» por eso sería falso.
    """
    try:
        peticion = urllib.request.Request(
            url, headers={"User-Agent": "OncoGO-CRC/2"}
        )
        with urllib.request.urlopen(peticion, timeout=60) as r:
            return {"url": url, "codigo": r.status, "acceso": "publico"}
    except urllib.error.HTTPError as e:
        restringido = e.code in {401, 402, 403}
        return {
            "url": url,
            "codigo": e.code,
            "acceso": "restringido" if restringido else "error_del_servidor",
        }
    except (urllib.error.URLError, TimeoutError) as e:
        return {
            "url": url,
            "codigo": None,
            "acceso": "no_comprobable",
            "motivo": (
                f"no se pudo comprobar desde esta máquina: {e}. No dice nada "
                "sobre si la página es pública o no."
            ),
        }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    ap.add_argument(
        "--estadisticas",
        default=None,
        help=(
            "ruta a un fichero de estadísticas de resumen ya descargado de "
            "Synapse por alguien con acceso. Sin esto, la fuente queda "
            "requiere_acceso."
        ),
    )
    args = ap.parse_args()

    ws = Path(args.workspace)
    salida = ws / "external" / "ukb_ppp"
    salida.mkdir(parents=True, exist_ok=True)

    paginas = {
        "articulo": _publico(ARTICULO),
        "portal": _publico(PORTAL),
        "synapse": _publico(SYNAPSE),
    }
    for nombre, p in paginas.items():
        print(f"[ukb-ppp] {nombre}: {p['acceso']} (código {p['codigo']})")

    comun = {
        "nombre": "ukb_ppp",
        "fecha_de_acceso": datetime.now(UTC).date().isoformat(),
        "articulo": ARTICULO,
        "portal": PORTAL,
        "synapse": {"url": SYNAPSE, "id": SYNAPSE_ID, "doi": SYNAPSE_DOI},
        "cita": CITA,
        "licencia": LICENCIA,
        "contenido_declarado_por_la_fuente": {
            "n_participantes": 54219,
            "n_proteinas": 2923,
            "n_asociaciones_primarias": 14287,
            "de_donde_sale": (
                "del resumen del portal público. No se ha descargado ningún dato."
            ),
        },
        "paginas_publicas": paginas,
        "columnas_esperadas": COLUMNAS,
        "proteinas_que_se_consultarian": OBJETIVO,
        "umbral_de_p_del_articulo": P_UMBRAL,
        "umbral_de_varianza_explicada": UMBRAL_VARIANZA,
        "para_que_serviria": (
            "Aportar una tercera explicación a la brecha tejido-plasma de C8, "
            "junto a la ventana de detección (C8) y la elección del grupo de "
            "control (D6): que el nivel en plasma de una proteína dependa de la "
            "variante que cada persona hereda."
        ),
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
    }

    if args.estadisticas is None:
        estado = {
            **comun,
            "estado": "requiere_acceso",
            "motivo": (
                "las estadísticas de resumen están en Synapse "
                f"({SYNAPSE_ID}) bajo los términos del UKB-PPP, que exigen "
                "cuenta y aceptación de condiciones. Este script no intenta "
                "autenticarse."
            ),
            "como_completarlo": (
                "Quien tenga cuenta de Synapse y haya aceptado los términos "
                "descarga el fichero y vuelve a ejecutar este script con "
                "--estadisticas RUTA. El análisis está escrito y probado, así "
                "que no hay nada más que programar."
            ),
            "analisis_que_queda_sin_evaluar": (
                "el peso de la genética en el nivel en plasma de las proteínas "
                "de v1 (D7)"
            ),
            "lo_que_no_se_hace": (
                "No se sustituye por ninguna otra fuente ni por datos "
                "simulados. Los datos sintéticos de este análisis existen solo "
                "en tests/, para comprobar que el código funciona."
            ),
        }
        (salida / "estado.json").write_text(
            json.dumps(estado, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        print(
            "[ukb-ppp] requiere_acceso: las estadísticas están en Synapse "
            f"{SYNAPSE_ID} y hacen falta cuenta y condiciones aceptadas."
        )
        print("[ukb-ppp] el análisis queda escrito y con tests; no se simula nada.")
        print(f"[ukb-ppp] ficha en {salida / 'estado.json'}")
        return 0

    ruta = Path(args.estadisticas)
    if not ruta.exists():
        print(f"[ukb-ppp] no existe {ruta}", file=sys.stderr)
        return 1
    sep = "\t" if ruta.suffix in {".tsv", ".txt", ".gz"} else ","
    d = pd.read_csv(ruta, sep=sep, low_memory=False)
    try:
        resumen = resumir_por_proteina(d)
    except KeyError as e:
        print(f"[ukb-ppp] {e}", file=sys.stderr)
        return 1

    tabla = salida / "ukb_ppp_cis_pqtl_resumen.tsv"
    resumen.to_csv(tabla, sep="\t", index=False)
    estado = {
        **comun,
        "estado": "usada",
        "fichero_de_entrada": {
            "ruta_declarada": str(ruta),
            "n_filas": int(len(d)),
            "proporcionado_por": (
                "una persona con acceso al UKB-PPP. El script no descargó nada."
            ),
        },
        "tabla": tabla.name,
        "lectura": lectura(resumen),
    }
    (salida / "estado.json").write_text(
        json.dumps(estado, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[ukb-ppp] {len(resumen)} proteínas con cis-pQTL significativo.")
    print(f"[ukb-ppp] escrito en {salida}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
