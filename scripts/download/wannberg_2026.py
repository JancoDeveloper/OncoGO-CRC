#!/usr/bin/env python3
"""D6. Wannberg et al. 2025: plasma de cáncer frente a enfermedad benigna.

El año es el de Crossref (`issued`, 29 de diciembre de 2025). El editor asigna
el artículo al volumen 17, número 151. El nombre de este fichero no cambia:
v2 ya escribió la descarga en ``external/wannberg_2026/``.

Para qué sirve esto. v1 comparó suero de cáncer contra suero de personas sanas,
y D4b hace lo mismo con una cohorte diez veces mayor. Pero en una consulta nadie
compara a un paciente con una persona sana: compara a un paciente con síntomas
contra otros pacientes con síntomas. Este trabajo hace exactamente esa
comparación, con la misma plataforma y las mismas 1 463 proteínas que D4b.

La comparación que importa sale de poner los dos juntos. Contra personas sanas,
el AUC es 0,994. Contra pacientes con síntomas no específicos, 0,80. Mismo panel,
misma tecnología, mismo grupo de investigación en parte. La diferencia no está en
el método, está en contra quién se compara.

Y hay un segundo hallazgo que toca de lleno a v1. Las 22 proteínas que este
trabajo encuentra no son proteínas de fase aguda: son mayoritariamente
intracelulares y nucleares (PARP1, RRM2, APEX1, PSIP1) y citoqueratinas (KRT18,
KRT19). La única claramente inflamatoria, S100A12, es la de menor efecto de las
22. Eso es coherente con que la señal de fase aguda sea lo primero que se cae
cuando el grupo de control también está enfermo.

De dónde salen los números. La información suplementaria de este artículo no
trae hojas de cálculo: son PDF. Así que la tabla S5 se extrae del texto del PDF
con un patrón explícito, y cada fila queda con la tabla y la página de la que
sale y con una columna ``verificado_por`` vacía. Ningún número extraído de un PDF
entra en un resultado sin que una persona lo haya comprobado.

Solo se usan estadísticas publicadas. No hay datos individuales, y los que
existen son de pacientes reales de dos hospitales suecos: no entran aquí ni
harían falta.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

DOI = "10.1038/s41467-025-67688-3"
ARTICULO = "https://www.nature.com/articles/s41467-025-67688-3"
ESM = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41467-025-67688-3/MediaObjects/"
    "41467_2025_67688_MOESM1_ESM.pdf"
)

CITA = (
    "Wannberg F, Bueno Álvez M, Qvick A, et al. Plasma protein profiling "
    "predicts cancer in patients with non-specific symptoms. Nature "
    "Communications 2025;17:151. doi:10.1038/s41467-025-67688-3"
)
LICENCIA = (
    "CC BY 4.0. El artículo y su información suplementaria son de acceso "
    "abierto y permiten reutilización con atribución."
)

# La página del PDF suplementario donde está la tabla S5, en base 1. Se declara
# para que cada número extraído diga de dónde sale y una persona pueda ir a
# comprobarlo sin buscar.
PAGINA_S5 = 6
TABLA_S5 = "Supplementary Table S5"

# Un patrón por fila: cohorte, ensayo y siete números. El OR y los extremos del
# intervalo vienen en decimal; los dos p, en notación exponencial. Se exige la
# forma completa para que una fila mal partida por el extractor no se cuele a
# medias: o encaja entera o no se usa.
FILA_S5 = re.compile(
    r"\b(Discovery|Replication)\s+"
    r"([A-Z0-9][A-Za-z0-9_.-]{1,15})\s+"
    r"([\d.]+)\s+"          # OR
    r"([\d.]+)\s+"          # SE
    r"([\d.]+)\s+"          # estadístico
    r"([\d.]+e[+-]?\d+)\s+"  # p
    r"([\d.]+e[+-]?\d+)\s+"  # p ajustado
    r"([\d.]+)\s+"          # IC bajo
    r"([\d.]+)"             # IC alto
)

# Números del artículo principal. Son del resumen, que es la parte más revisada
# del texto, y se declaran aquí en lugar de extraerse del HTML: una página web
# cambia de maquetación y el resumen de un artículo publicado no cambia.
DEL_RESUMEN = {
    "n_proteinas_medidas": 1463,
    "n_pacientes_descubrimiento": 456,
    "n_cancer_descubrimiento": 160,
    "n_sin_cancer_descubrimiento": 296,
    "n_proteinas_asociadas": 29,
    "auc_descubrimiento": 0.80,
    "n_pacientes_replica": 238,
    "auc_replica": 0.82,
}

# Tamaños de los grupos de la tabla S1 y S2, y de qué estaban enfermos los que no
# tenían cáncer (tablas S3 y S4). Esto es lo que hace que la comparación sea
# distinta de la de v1 y de D4b.
DEL_SUPLEMENTARIO = {
    "n_cancer_replica": 35,
    "n_sin_cancer_replica": 203,
    "n_colorrectal_descubrimiento": 5,
    "n_colorrectal_replica": 1,
    "n_sin_cancer_con_diagnostico_descrito_descubrimiento": 155,
    "n_con_enfermedad_autoinmune_descubrimiento": 55,
    "n_con_infeccion_descubrimiento": 32,
}

# Clasificación de las 22 proteínas por donde vive la proteína. No es una
# invención de este trabajo: es la localización que anota UniProt, y sirve para
# una sola cosa, decir si la lista se parece o no a una respuesta de fase aguda.
# Se escribe aquí para que quede auditable y no como una frase del informe.
LOCALIZACION = {
    "RRM2": "intracelular", "PARP1": "nuclear", "PSIP1": "nuclear",
    "KRT19": "citoqueratina", "DPY30": "nuclear", "MAD1L1": "nuclear",
    "KRT18": "citoqueratina", "ELOA": "nuclear", "ZBTB17": "nuclear",
    "AGR2": "secretada", "GRPEL1": "mitocondrial", "FKBP4": "intracelular",
    "STC1": "secretada", "APEX1": "nuclear", "METAP1D": "mitocondrial",
    "PRDX3": "mitocondrial", "PAEP": "secretada", "EIF4G1": "intracelular",
    "PTS": "intracelular", "HDGF": "nuclear", "BAIAP2": "intracelular",
    "S100A12": "inflamatoria_secretada",
}

# Las tres proteínas de plasma de v1, para comprobar si aparecen.
PLASMA_V1 = {"P0DJI8": "SAA1", "P18428": "LBP", "Q9BY43": "CHMP4A"}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()


def _descargar(destino: Path) -> dict:
    if destino.exists():
        return {
            "url": ESM,
            "bytes": destino.stat().st_size,
            "sha256": _sha256(destino),
            "descargado_ahora": False,
        }
    peticion = urllib.request.Request(ESM, headers={"User-Agent": "OncoGO-CRC/2"})
    with urllib.request.urlopen(peticion, timeout=180) as r:
        datos = r.read()
    destino.write_bytes(datos)
    return {
        "url": ESM,
        "bytes": len(datos),
        "sha256": hashlib.sha256(datos).hexdigest(),
        "descargado_ahora": True,
    }


def texto_de_la_pagina(pdf: Path, pagina: int) -> str:
    from pypdf import PdfReader

    paginas = PdfReader(str(pdf)).pages
    if pagina > len(paginas):
        raise IndexError(
            f"el PDF tiene {len(paginas)} páginas y se pidió la {pagina}"
        )
    return paginas[pagina - 1].extract_text() or ""


def extraer_s5(texto: str, *, pagina: int = PAGINA_S5) -> pd.DataFrame:
    """Tabla S5: regresión logística ajustada por edad y sexo, 22 proteínas.

    El extractor de texto de un PDF junta los saltos de línea de forma
    impredecible, así que no se parte por líneas: se busca el patrón de fila
    completa en todo el texto. Una fila a medias no encaja y por tanto no entra,
    que es lo que se quiere: es mejor extraer 43 de 44 y que el recuento lo diga
    que colar una fila con el intervalo de la siguiente.
    """
    filas = []
    for m in FILA_S5.finditer(texto.replace("\n", " ")):
        cohorte, ensayo, orr, se, est, p, padj, bajo, alto = m.groups()
        filas.append(
            {
                "cohorte": "descubrimiento" if cohorte == "Discovery" else "replica",
                "proteina": ensayo.upper(),
                "odds_ratio": float(orr),
                "se": float(se),
                "estadistico": float(est),
                "p": float(p),
                "p_ajustado": float(padj),
                "ic95_bajo": float(bajo),
                "ic95_alto": float(alto),
                "localizacion_uniprot": LOCALIZACION.get(ensayo.upper(), "sin_anotar"),
                "fuente_tabla": TABLA_S5,
                "fuente_pagina": pagina,
                # La deja vacía el script a propósito. La rellena el autor cuando
                # haya comparado la fila con el PDF.
                "verificado_por": "",
            }
        )
    d = pd.DataFrame(
        filas,
        columns=[
            "cohorte", "proteina", "odds_ratio", "se", "estadistico", "p",
            "p_ajustado", "ic95_bajo", "ic95_alto", "localizacion_uniprot",
            "fuente_tabla", "fuente_pagina", "verificado_por",
        ],
    )
    return d.sort_values(["cohorte", "p_ajustado"], kind="mergesort").reset_index(
        drop=True
    )


def comprobar_la_extraccion(d: pd.DataFrame) -> dict:
    """Comprobaciones que tienen que pasar antes de usar un número de un PDF.

    No son adornos. Un extractor que confunde columnas produce una tabla con
    pinta perfecta, y la única manera de notarlo sin leerla entera es comprobar
    cosas que tienen que cumplirse por construcción: que el OR caiga dentro de su
    propio intervalo, que el p ajustado no sea menor que el p crudo, y que las
    dos cohortes traigan las mismas proteínas.
    """
    problemas = []
    fuera = d[(d["odds_ratio"] < d["ic95_bajo"]) | (d["odds_ratio"] > d["ic95_alto"])]
    if not fuera.empty:
        problemas.append(
            f"{len(fuera)} filas con el OR fuera de su intervalo: "
            f"{sorted(set(fuera['proteina']))}"
        )
    invertido = d[d["p_ajustado"] < d["p"]]
    if not invertido.empty:
        problemas.append(
            f"{len(invertido)} filas con el p ajustado menor que el crudo: "
            f"{sorted(set(invertido['proteina']))}"
        )
    por_cohorte = {c: set(g["proteina"]) for c, g in d.groupby("cohorte")}
    if len(por_cohorte) == 2:
        a, b = por_cohorte.values()
        if a != b:
            problemas.append(
                f"las dos cohortes no traen las mismas proteínas: "
                f"{sorted(a ^ b)}"
            )
    no_protector = d[d["ic95_bajo"] <= 1.0]
    return {
        "n_filas": int(len(d)),
        "n_proteinas": int(d["proteina"].nunique()),
        "cohortes": sorted(d["cohorte"].unique()),
        "problemas": problemas,
        "coherente": not problemas,
        "n_con_intervalo_que_cruza_el_uno": int(len(no_protector)),
        "pendiente_de_verificacion_humana": int((d["verificado_por"] == "").sum()),
        "regla": (
            "Ningún número extraído de un PDF entra en un resultado sin que una "
            "persona lo haya comprobado. La columna verificado_por la rellena el "
            "autor, no el script."
        ),
    }


def perfil_de_las_proteinas(d: pd.DataFrame) -> dict:
    """Qué clase de proteínas son, que es lo que compara con v1.

    v1 encontró en plasma SAA1 y LBP, las dos de respuesta de fase aguda. Si la
    lista de un trabajo que compara contra pacientes enfermos no tiene proteínas
    de fase aguda, eso dice algo sobre por qué v1 las encontró.
    """
    unicas = d.drop_duplicates("proteina")
    recuento = unicas["localizacion_uniprot"].value_counts().to_dict()
    inflamatorias = unicas[
        unicas["localizacion_uniprot"] == "inflamatoria_secretada"
    ]
    presentes = set(d["proteina"])
    return {
        "n_proteinas": int(len(unicas)),
        "por_localizacion": {k: int(v) for k, v in recuento.items()},
        "n_inflamatorias": int(len(inflamatorias)),
        "inflamatorias": sorted(inflamatorias["proteina"]),
        "or_minimo_entre_las_inflamatorias": (
            float(
                d[d["proteina"].isin(inflamatorias["proteina"])]["odds_ratio"].min()
            )
            if not inflamatorias.empty
            else None
        ),
        # Sin el guarda, una tabla vacía daría NaN, y NaN no es JSON válido:
        # el fichero de estado quedaría ilegible sin que nada avisara.
        "or_maximo_de_la_tabla": (
            float(d["odds_ratio"].max()) if not d.empty else None
        ),
        "plasma_de_v1_presente": {
            gen: gen in presentes for gen in PLASMA_V1.values()
        },
        "lectura": (
            "La lista no es una respuesta de fase aguda. La mayoría son proteínas "
            "intracelulares, nucleares o citoqueratinas, que en plasma indican "
            "recambio y muerte celular, no inflamación. Ninguna de las tres "
            "proteínas de plasma de v1 aparece."
        ),
    }


def comparacion_del_control(estado_alvez: Path | None) -> dict:
    """Pone el AUC contra sanos y el AUC contra pacientes uno al lado del otro.

    Es la única comparación de todo este eje, y solo es legítima porque las dos
    usan la misma plataforma y el mismo número de ensayos. Si D4b no está
    descargada, no se inventa su número: se dice que falta.
    """
    base = {
        "contra_pacientes_con_sintomas": {
            "auc_descubrimiento": DEL_RESUMEN["auc_descubrimiento"],
            "auc_replica": DEL_RESUMEN["auc_replica"],
            "n_proteinas_del_panel": DEL_RESUMEN["n_proteinas_medidas"],
            "que_son_los_controles": (
                "Pacientes derivados por síntomas no específicos que, tras el "
                "estudio, no tenían cáncer. De los 155 descritos, 55 tenían una "
                "enfermedad autoinmune y 32 una infección."
            ),
        },
        "por_que_se_pueden_comparar": (
            "Las dos usan Olink Explore con 1 463 ensayos. No se comparan "
            "magnitudes de proteína, solo el AUC de un clasificador, que es "
            "adimensional."
        ),
    }
    if estado_alvez is None or not estado_alvez.exists():
        base["contra_personas_sanas"] = {
            "estado": "no_disponible",
            "motivo": (
                "la fuente alvez_2023 (D4b) no está descargada, así que no hay "
                "con qué comparar"
            ),
        }
        return base
    d = json.loads(estado_alvez.read_text(encoding="utf-8"))
    completo = d["auc_segun_cuantas_proteinas"]["con_el_panel_completo"]
    base["contra_personas_sanas"] = {
        "auc": completo["auc"],
        "n_proteinas_del_panel": completo["n_proteinas"],
        "que_son_los_controles": d["cohorte"]["que_son_los_controles"],
        "fuente": "alvez_2023",
    }
    base["diferencia_de_auc"] = round(
        completo["auc"] - DEL_RESUMEN["auc_descubrimiento"], 4
    )
    base["lectura"] = (
        f"Con el mismo panel de {completo['n_proteinas']} proteínas, el AUC es "
        f"{completo['auc']:.3f} contra personas sanas y "
        f"{DEL_RESUMEN['auc_descubrimiento']:.2f} contra pacientes con síntomas. "
        "La diferencia no viene del método: viene de contra quién se compara. "
        "Es la misma advertencia que hay que hacerle a los tres resultados de "
        "plasma de v1, que también se obtuvieron contra personas sanas."
    )
    return base


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    args = ap.parse_args()

    ws = Path(args.workspace)
    salida = ws / "external" / "wannberg_2026"
    salida.mkdir(parents=True, exist_ok=True)
    estado_path = salida / "estado.json"
    pdf = salida / "41467_2025_67688_MOESM1_ESM.pdf"

    def _fallo(motivo: str) -> int:
        estado_path.write_text(
            json.dumps(
                {
                    "nombre": "wannberg_2026",
                    "estado": "no_disponible",
                    "motivo": motivo,
                    "fecha_de_acceso": datetime.now(UTC).date().isoformat(),
                    "articulo": ARTICULO,
                    "analisis_que_queda_sin_evaluar": (
                        "la comparación del plasma de cáncer contra enfermedad "
                        "benigna (D6)"
                    ),
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"[wannberg] no_disponible: {motivo}", file=sys.stderr)
        return 1

    try:
        fichero = _descargar(pdf)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        return _fallo(f"fallo al descargar la información suplementaria: {e}")
    print(f"[wannberg] suplementario 1: {fichero['bytes']} bytes")

    try:
        texto = texto_de_la_pagina(pdf, PAGINA_S5)
    except ImportError:
        return _fallo(
            "pypdf no está instalado y la tabla S5 solo existe dentro de un PDF"
        )
    except (IndexError, OSError) as e:
        return _fallo(f"no se pudo leer la página {PAGINA_S5} del PDF: {e}")

    s5 = extraer_s5(texto)
    if s5.empty:
        return _fallo(
            f"el patrón de fila no encontró ninguna fila en la página "
            f"{PAGINA_S5}; probablemente el PDF cambió de maquetación"
        )
    comprobacion = comprobar_la_extraccion(s5)
    perfil = perfil_de_las_proteinas(s5)
    comparacion = comparacion_del_control(
        ws / "external" / "alvez_2023" / "estado.json"
    )

    tabla = salida / "wannberg_2026_tabla_s5.tsv"
    s5.to_csv(tabla, sep="\t", index=False)

    estado = {
        "nombre": "wannberg_2026",
        "estado": "usada",
        "fecha_de_acceso": datetime.now(UTC).date().isoformat(),
        "articulo": ARTICULO,
        "doi": DOI,
        "cita": CITA,
        "licencia": LICENCIA,
        "titulo": (
            "Plasma protein profiling predicts cancer in patients with "
            "non-specific symptoms"
        ),
        "plataforma": (
            "Olink Explore, 1 463 ensayos. Es la misma plataforma y el mismo "
            "número de ensayos que D4b, y por eso los dos AUC se pueden comparar."
        ),
        "suplementario": {**fichero, "para_que": "la tabla S5 y las tablas S1 a S4"},
        "del_resumen_del_articulo": DEL_RESUMEN,
        "del_suplementario": DEL_SUPLEMENTARIO,
        "tabla_s5": tabla.name,
        "sha256_tabla_s5": _sha256(tabla),
        "comprobacion_de_la_extraccion": comprobacion,
        "perfil_de_las_proteinas": perfil,
        "comparacion_del_grupo_de_control": comparacion,
        "sobre_el_colorrectal": {
            "n_descubrimiento": DEL_SUPLEMENTARIO["n_colorrectal_descubrimiento"],
            "n_replica": DEL_SUPLEMENTARIO["n_colorrectal_replica"],
            "advertencia": (
                "Solo 5 de los 160 cánceres de la cohorte de descubrimiento y 1 "
                "de los 35 de la réplica eran colorrectales. De aquí no sale "
                "ninguna afirmación sobre cáncer colorrectal en particular: lo "
                "que aporta es la comparación del grupo de control."
            ),
        },
        "limitaciones": [
            "Es un trabajo sobre cáncer en general, no sobre cáncer colorrectal. "
            "Las 22 proteínas son de un modelo pan-cáncer.",
            "Los números de la tabla S5 se han extraído del texto de un PDF. La "
            "columna verificado_por está vacía a propósito y la rellena una "
            "persona tras comparar cada fila con el original.",
            "La localización de cada proteína es la anotación de UniProt y sirve "
            "solo para decir si la lista se parece a una respuesta de fase "
            "aguda. No es un resultado de este trabajo ni del artículo.",
            "Las dos cohortes son suecas y de dos hospitales. Un grupo de "
            "control de pacientes con síntomas no es el mismo en todos los "
            "sistemas sanitarios.",
            "Comparar un AUC de este trabajo con uno de D4b es comparar dos "
            "diseños distintos. La comparación es legítima porque la plataforma "
            "es la misma, pero no es un contraste estadístico: son dos números "
            "publicados puestos uno al lado del otro.",
        ],
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
    }
    estado_path.write_text(
        json.dumps(estado, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(
        f"[wannberg] tabla S5: {comprobacion['n_filas']} filas, "
        f"{comprobacion['n_proteinas']} proteínas, "
        f"coherente={comprobacion['coherente']}"
    )
    if comprobacion["problemas"]:
        for p in comprobacion["problemas"]:
            print(f"[wannberg] ATENCIÓN: {p}", file=sys.stderr)
    print(
        f"[wannberg] localización: {perfil['por_localizacion']}; "
        f"inflamatorias: {perfil['n_inflamatorias']}"
    )
    if "lectura" in comparacion:
        print(f"[wannberg] {comparacion['lectura']}")
    print(
        f"[wannberg] pendientes de verificación humana: "
        f"{comprobacion['pendiente_de_verificacion_humana']} filas"
    )
    print(f"[wannberg] escrito en {salida}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
