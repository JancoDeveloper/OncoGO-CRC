#!/usr/bin/env python3
"""D3. Búsqueda de cohortes de tejido adicionales, y lo que apareció en su lugar.

Qué se buscaba. El preregistro fija los criterios de inclusión de D3 en
``d3_inclusion_cohortes``: tejido humano, tumor y normal, matriz descargable y
**al menos 10 pacientes**. Se revisaron los 571 depósitos de PRIDE que responden
a «colorectal» y los estudios colorrectales del PDC.

Qué se encontró. Ningún depósito nuevo cumple los cuatro criterios. Los que
tienen matriz publicada no llegan a 10 pacientes o no traen tejido normal; los
que tienen el diseño adecuado no publican matriz de proteínas, solo espectros. La
lista de candidatos y el motivo de cada descarte están en ``BUSQUEDA``, porque
«no encontramos más» y «buscamos y esto es lo que hay» no son lo mismo.

Qué se usa. ``PXD036434`` no cumple el mínimo de pacientes —son cuatro— así que
**no entra en el meta-análisis**, y el umbral no se dobla para que entre. Pero
trae algo que ninguna cohorte del certificado tiene: ocho regiones distintas del
mismo tumor de un mismo paciente. Eso permite preguntar algo que toca al diseño
de toda la tesis: si un trozo de tumor representa al tumor.

La pregunta importa porque las tres cohortes del certificado miden una biopsia
por paciente. Si las regiones de un tumor difieren entre sí tanto como difieren
pacientes distintos, entonces lo que mide una cohorte de biopsias únicas es más
ruidoso de lo que su tamaño sugiere, y eso no se arregla con más pacientes.

Nada de este script identifica a nadie: son cuatro pacientes anónimos de un
depósito público, sin edad, sexo ni ninguna otra característica.
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

import numpy as np
import pandas as pd

ACCESION = "PXD036434"
API = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{}/files?pageSize=500"
PORTAL = f"https://www.ebi.ac.uk/pride/archive/projects/{ACCESION}"
FICHERO = "Protein_quantity_matrix.xlsx"
DOI = "10.1186/s12967-022-03855-0"
CITA = (
    "Praus F, et al. Panomics reveals patient individuality as the major "
    "driver of colorectal cancer progression. Journal of Translational "
    "Medicine 2023;21:41. doi:10.1186/s12967-022-03855-0"
)
LICENCIA = (
    "PRIDE Archive, datos públicos. La reutilización está permitida citando el "
    "depósito y el artículo."
)

# Lo que se revisó y por qué no entró. Cada línea es un candidato real.
BUSQUEDA: list[dict] = [
    {
        "fuente": "PXD019504",
        "titulo": "Global proteome analysis of paired colorectal cancer sample",
        "veredicto": "incluida",
        "motivo": "cumple los cuatro criterios; es la cohorte que C10 ya usa",
    },
    {
        "fuente": "PDC000278",
        "titulo": "VU Normal Colon Epithelium - Proteome",
        "veredicto": "no_disponible",
        "motivo": (
            "60 alícuotas de epitelio de colon normal de 30 personas, sin solape "
            "con PDC000111, pero el depósito no publica ninguna matriz de "
            "proteínas: solo espectros crudos (393 GB), procesados (77 GB) y "
            "ficheros de PSM. Construir la matriz desde los PSM sería un "
            "reanálisis con decisiones distintas de las del CDAP, así que "
            "compararla con PDC000111 compararía dos pipelines y no dos tejidos."
        ),
    },
    {
        "fuente": "PDC000109",
        "titulo": "Prospective COAD Proteome S037-1",
        "veredicto": "excluida",
        "motivo": (
            "97 de sus 100 pacientes son los mismos de PDC000116, la cohorte de "
            "descubrimiento. No es independiente."
        ),
    },
    {
        "fuente": "PXD002137",
        "titulo": "Proteome Landscapes of the Human Colorectal Mucosa, and its "
                  "Adenoma and Cancer",
        "veredicto": "no_disponible",
        "motivo": (
            "el diseño es el ideal —mucosa, adenoma y cáncer— pero el depósito "
            "solo tiene 99 ficheros crudos y uno de búsqueda: no hay matriz de "
            "proteínas publicada"
        ),
    },
    {
        "fuente": "PXD001676",
        "titulo": "Quantitative proteomics analysis of paired colorectal cancer "
                  "and non-tumorigenic tissues",
        "veredicto": "no_disponible",
        "motivo": "sin matriz de proteínas publicada; solo crudos, búsqueda y picos",
    },
    {
        "fuente": "PXD060933",
        "titulo": "In-Depth Proteomic Analysis of Stage II Colorectal Cancer",
        "veredicto": "excluida",
        "motivo": (
            "tiene proteinGroups publicado, pero sus 12 muestras son todas de "
            "cáncer colorrectal y no hay tejido normal"
        ),
    },
    {
        "fuente": "PXD059755",
        "titulo": "Glycoproteomics and proteomics analysis of Colorectal Cancer "
                  "and adjacent normal FFPE tissues",
        "veredicto": "excluida",
        "motivo": (
            "solo 594 grupos de proteínas, y la anotación publicada deja el "
            "campo de enfermedad como «not available»: asignar qué canal de TMT "
            "es tumor y cuál normal sería adivinarlo"
        ),
    },
    {
        "fuente": "PXD051882",
        "titulo": "Glycoproteomics and Proteomics analysis of Colorectal Cancer "
                  "tissues",
        "veredicto": "excluida",
        "motivo": (
            "el diseño publicado da el canal de TMT de cada muestra (A21, B21...) "
            "pero no dice cuál de esas etiquetas es tumor y cuál mucosa normal"
        ),
    },
    {
        "fuente": ACCESION,
        "titulo": "Panomics reveals patient individuality as the major driver of "
                  "colorectal cancer progression",
        "veredicto": "no_cumple_el_minimo_de_pacientes",
        "motivo": (
            "son 4 pacientes y el preregistro exige 10, así que no entra en el "
            "meta-análisis. Se usa para una pregunta distinta, que no es un "
            "contraste entre cohortes: la variación dentro de un mismo tumor."
        ),
    },
]

# Las ocho regiones del tumor del paciente 4, y el resto del diseño.
REGIONES_DEL_MISMO_TUMOR = [f"T4_{i}" for i in range(1, 9)]
TUMORES_DE_UN_PACIENTE = ["T1", "T2", "T3"]
MUCOSAS = ["NM1", "NM2", "NM3", "NM4"]
METASTASIS = ["LM1", "LM2", "LM3", "LM4"]
# El paciente 4 aporta su mucosa y su metástasis pero su tumor son ocho columnas,
# así que para los contrastes pareados hace falta elegir una y declararlo.
TUMOR_DEL_PACIENTE_4 = "T4_1"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def localizar() -> str:
    """Pide a PRIDE la URL del fichero, en vez de construirla a mano.

    La ruta del archivo lleva el año y el mes del depósito, así que una URL
    escrita a mano se rompe en silencio cuando PRIDE reorganiza.
    """
    req = urllib.request.Request(
        API.format(ACCESION), headers={"Accept": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=180) as r:
        datos = json.loads(r.read().decode("utf-8"))
    ficheros = (
        datos.get("_embedded", {}).get("files", [])
        if isinstance(datos, dict)
        else datos
    )
    for f in ficheros:
        if f.get("fileName") == FICHERO:
            for loc in f.get("publicFileLocations", []):
                url = str(loc.get("value", ""))
                if url.startswith("ftp://"):
                    return "https://" + url[len("ftp://") :]
            raise RuntimeError(f"{FICHERO} existe pero sin ubicación pública FTP")
    raise RuntimeError(
        f"{FICHERO} no está en {ACCESION}; el depósito cambió de contenido y hay "
        "que volver a mirar qué fichero trae la matriz"
    )


def _descargar(url: str, destino: Path) -> dict:
    if destino.exists() and destino.stat().st_size > 0:
        return {
            "fichero": destino.name,
            "bytes": destino.stat().st_size,
            "sha256": _sha256(destino),
            "ya_estaba": True,
        }
    with urllib.request.urlopen(url, timeout=600) as r:
        datos = r.read()
    destino.write_bytes(datos)
    return {
        "fichero": destino.name,
        "bytes": len(datos),
        "sha256": hashlib.sha256(datos).hexdigest(),
        "ya_estaba": False,
    }


def parsear(xlsx: Path) -> tuple[pd.DataFrame, dict]:
    """Matriz de genes por muestra, en log2.

    Los valores del depósito son cocientes contra una abundancia de referencia
    común, así que el log2 los convierte en lo mismo que maneja el resto del
    trabajo: cambios multiplicativos. La referencia es la misma para todas las
    columnas, de modo que las comparaciones entre columnas no la arrastran.

    El identificador se toma del símbolo de gen que viene en ``DESCRIPTION``
    (``GN=PLEC``) y no de la columna ``PROTEIN ID``, porque esa trae nombres de
    entrada de UniProt (``PLEC_HUMAN``) y el resto del trabajo cruza por símbolo
    de gen. Cruzar un nombre de entrada con un símbolo no falla: devuelve cero
    coincidencias y parece que las cohortes no comparten proteínas.
    """
    hojas = pd.ExcelFile(xlsx)
    d = hojas.parse(hojas.sheet_names[0])
    for c in ("DESCRIPTION", "PROTEIN ID"):
        if c not in d.columns:
            raise RuntimeError(
                f"la matriz no tiene la columna «{c}»; el depósito cambió de "
                "formato y seguir daría genes mal emparejados"
            )
    d["gen"] = d["DESCRIPTION"].astype(str).str.extract(r"GN=([A-Za-z0-9_.-]+)")
    sin_gen = int(d["gen"].isna().sum())
    d = d.dropna(subset=["gen"])
    duplicados = int(d["gen"].duplicated().sum())
    d = d.drop_duplicates("gen", keep="first").set_index("gen")

    columnas = [
        c
        for c in d.columns
        if re.fullmatch(r"(LM|NM|T)\d+(_\d+)?", str(c))
    ]
    esperadas = set(
        REGIONES_DEL_MISMO_TUMOR + TUMORES_DE_UN_PACIENTE + MUCOSAS + METASTASIS
    )
    faltan = sorted(esperadas - set(columnas))
    if faltan:
        raise RuntimeError(
            f"faltan columnas del diseño: {faltan}. El análisis depende de que "
            "T, NM y LM signifiquen tumor, mucosa normal y metástasis hepática, "
            "y si el depósito renombró algo hay que volver a leerlo"
        )
    m = d[columnas].apply(pd.to_numeric, errors="coerce")
    ceros = int((m <= 0).sum().sum())
    m = m.mask(m <= 0.0)
    m = pd.DataFrame(np.log2(m.to_numpy(dtype=float)), index=m.index, columns=columnas)

    resumen = {
        "n_genes": int(m.shape[0]),
        "n_columnas": int(m.shape[1]),
        "sin_simbolo_de_gen": sin_gen,
        "genes_duplicados_descartados": duplicados,
        "valores_no_positivos_tratados_como_ausencia": ceros,
        "cuantificacion": (
            "Cociente contra una abundancia de referencia común, en log2. La "
            "referencia es la misma para todas las columnas."
        ),
        "identificador": (
            "Símbolo de gen extraído de GN= en DESCRIPTION, no el nombre de "
            "entrada de UniProt de la columna PROTEIN ID."
        ),
        "diseno": {
            "regiones_del_mismo_tumor": REGIONES_DEL_MISMO_TUMOR,
            "tumores_de_un_paciente": TUMORES_DE_UN_PACIENTE,
            "mucosas_normales": MUCOSAS,
            "metastasis_hepaticas": METASTASIS,
            "tumor_elegido_para_el_paciente_4": TUMOR_DEL_PACIENTE_4,
            "por_que_hay_que_elegir": (
                "El paciente 4 tiene ocho regiones de tumor y una sola mucosa. "
                "Para un contraste pareado hay que escoger una región, y la "
                "elección cambia el resultado: por eso el análisis la promedia "
                "sobre todas las elecciones posibles en lugar de fijar una."
            ),
        },
    }
    return m, resumen


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    args = ap.parse_args()

    ws = Path(args.workspace)
    salida = ws / "external" / "heterogeneidad_intratumoral"
    (salida / "crudos").mkdir(parents=True, exist_ok=True)
    estado_path = salida / "estado.json"

    def _fallo(motivo: str) -> int:
        estado_path.write_text(
            json.dumps(
                {
                    "nombre": "heterogeneidad_intratumoral",
                    "estado": "no_disponible",
                    "motivo": motivo,
                    "fecha_de_acceso": datetime.now(UTC).date().isoformat(),
                    "accesion": ACCESION,
                    "portal": PORTAL,
                    "busqueda_de_cohortes_d3": BUSQUEDA,
                    "analisis_que_queda_sin_evaluar": (
                        "si un trozo de tumor representa al tumor (D3)"
                    ),
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"[heterogeneidad] no_disponible: {motivo}", file=sys.stderr)
        return 1

    try:
        url = localizar()
        fichero = _descargar(url, salida / "crudos" / FICHERO)
    except (RuntimeError, urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        return _fallo(f"fallo al localizar o descargar la matriz: {e}")
    print(f"[heterogeneidad] matriz: {fichero['bytes']} bytes")

    try:
        matriz, resumen = parsear(salida / "crudos" / FICHERO)
    except (RuntimeError, ValueError) as e:
        return _fallo(f"la matriz no se pudo parsear: {e}")

    limpia = salida / "heterogeneidad_intratumoral_proteoma.tsv"
    matriz.to_csv(limpia, sep="\t")
    print(
        f"[heterogeneidad] {resumen['n_genes']} genes x {resumen['n_columnas']} "
        "muestras de 4 pacientes."
    )

    estado = {
        "nombre": "heterogeneidad_intratumoral",
        "estado": "usada",
        "fecha_de_acceso": datetime.now(UTC).date().isoformat(),
        "accesion": ACCESION,
        "portal": PORTAL,
        "doi": DOI,
        "cita": CITA,
        "licencia": LICENCIA,
        "fichero_crudo": fichero,
        "matriz_limpia": limpia.name,
        "sha256_matriz_limpia": _sha256(limpia),
        "parseo": resumen,
        "n_pacientes": 4,
        "no_entra_en_el_meta_analisis": (
            "El preregistro exige al menos 10 pacientes por cohorte "
            "(d3_inclusion_cohortes.min_pacientes). Estos son 4, así que esta "
            "cohorte no entra en el meta-análisis de C10 y el umbral no se dobla "
            "para que entre. Se usa para una pregunta que no es un contraste "
            "entre cohortes."
        ),
        "para_que_sirve": (
            "Ninguna cohorte del certificado tiene más de un trozo del mismo "
            "tumor. Esta tiene ocho del mismo tumor de un paciente, y eso permite "
            "preguntar si una biopsia representa al tumor."
        ),
        "busqueda_de_cohortes_d3": BUSQUEDA,
        "criterios_del_preregistro": (
            "tejido humano, tumor y normal, matriz descargable y al menos 10 "
            "pacientes; el contraste se hace dentro de cada cohorte y las "
            "cohortes solo se combinan a nivel de efectos"
        ),
        "limitaciones": [
            "Son cuatro pacientes. De aquí no sale ninguna estimación de efecto "
            "que se pueda comparar con las cohortes del certificado.",
            "Las ocho regiones son de UN tumor de UN paciente. Lo que se puede "
            "decir es que existe al menos un tumor así, no cuántos tumores son "
            "así. Generalizarlo sería inventar.",
            "La cuantificación es un cociente contra una referencia común y la "
            "plataforma no es la de ninguna cohorte del certificado, así que las "
            "magnitudes no son comparables entre fuentes. Lo que se compara aquí "
            "es variación contra variación dentro de la misma matriz.",
            "Las regiones del mismo tumor se midieron probablemente más juntas en "
            "el tiempo que las muestras de pacientes distintos. Si eso sesga "
            "algo, las haría parecer más parecidas entre sí, no menos: el sesgo "
            "va en contra del resultado y no a favor.",
        ],
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
    }
    estado_path.write_text(
        json.dumps(estado, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[heterogeneidad] usada: {limpia}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
