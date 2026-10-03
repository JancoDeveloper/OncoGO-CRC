#!/usr/bin/env python3
"""D5. Referencia unicelular para deconvolver la composición del tejido.

Para qué sirve. R4 ajusta por un índice estromal construido del propio proteoma,
y el preregistro lista ``deconvolucion_referencia_externa`` entre los métodos de
sensibilidad de C7. Una referencia externa rompería la circularidad: los perfiles
de cada tipo celular vendrían de fuera del proteoma que se está contrastando.

De dónde sale la referencia. Del Human Protein Atlas, fichero
``rna_single_cell_type.tsv``, que publica la expresión media de cada gen en cada
tipo celular a partir de datos unicelulares. Es público y no requiere cuenta.

Lo que no se pudo usar, y por qué. El atlas unicelular de cáncer colorrectal de
Pelka et al. (``GSE178341``) sería la referencia específica de este tejido. No se
ha podido descargar desde esta máquina: el servidor ``ftp.ncbi.nlm.nih.gov`` no
resuelve en el DNS de esta red. **Eso no es una restricción de acceso**: el
fichero es público y la ruta queda escrita abajo para completarlo con
``--pelka RUTA`` cuando se descargue desde otra red. Se marca
``no_comprobable_desde_esta_maquina`` y no «no disponible», porque decir que no
está disponible sería falso.

Las dos limitaciones de fondo, que no se arreglan con otra referencia.

La primera es que la referencia es de ARN y lo que se deconvuelve es proteína. El
acuerdo entre el transcrito y la proteína del mismo gen es flojo: D2 lo midió en
CEACAM5 y salió una correlación de rangos de 0,24. Una firma de ARN aplicada a un
proteoma es una aproximación, y hay que decirlo cada vez.

La segunda es que los tipos celulares del atlas son de todos los tejidos juntos.
Un fibroblasto de piel y uno de colon no tienen el mismo perfil, así que
«fibroblastos» aquí es un promedio entre tejidos.

Este script no descarga ningún dato de ninguna persona: son medias por tipo
celular ya agregadas y publicadas.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import socket
import sys
import urllib.error
import urllib.request
import zipfile
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

HPA_URL = "https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip"
HPA_INTERNO = "rna_single_cell_type.tsv"
HPA_VERSION_INFO = "https://www.proteinatlas.org/about/download"
LICENCIA_HPA = (
    "Human Protein Atlas, licencia Creative Commons Attribution-ShareAlike 4.0 "
    "International. Permite reutilización citando la fuente."
)
CITA_HPA = (
    "Karlsson M, Zhang C, Méar L, et al. A single-cell type transcriptomics map "
    "of human tissues. Science Advances 2021; 7(31):eabh2169."
)

# La referencia específica de colon, que no se ha podido traer desde esta red.
PELKA = {
    "acceso": "GSE178341",
    "portal": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178341",
    "titulo": "A Single Cell Atlas of MMRd and MMRp Colorectal Cancer",
    "cita": (
        "Pelka K, Hofree M, Chen JH, et al. Spatially organized multicellular "
        "immune hubs in human colorectal cancer. Cell 2021; 184(18):4734-4752."
    ),
    "ficheros": [
        "GSE178341_crc10x_full_c295v4_submit.h5",
        "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz",
        "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz",
    ],
    "ruta": (
        "ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178341/suppl/"
    ),
    "intentos": (
        "2026-10-03: el DNS de esta red resuelve el servidor de forma "
        "intermitente. En un momento en que resolvía, el fichero respondió con "
        "código 200 y 1203550558 bytes anunciados, y la descarga se cortó a los "
        "100 MB. El parcial se borró para que ninguna etapa leyera un fichero "
        "truncado que parecería válido."
    ),
    "como_completarlo": (
        "Descargar los tres ficheros desde una red que resuelva "
        "ftp.ncbi.nlm.nih.gov y volver a ejecutar este script con "
        "--pelka CARPETA. El h5 trae la matriz de células por genes y el "
        "cluster.csv la etiqueta de tipo celular de cada célula; la firma se "
        "construye promediando por etiqueta."
    ),
}

# Los tipos celulares del colon. Se fijan aquí, antes de mirar ningún resultado,
# porque elegirlos después de ver qué sale es elegir el resultado.
TIPOS_DEL_COLON = [
    "colonocytes",
    "enterocytes",
    "goblet cells",
    "paneth cells",
    "fibroblasts",
    "macrophages",
    "t-cells",
    "b-cells",
    "plasma cells",
    "vascular endothelial cells",
    "smooth muscle cells",
]

# Los cuatro que son epitelio. El eje epitelio frente a no epitelio es lo único
# que esta deconvolución sostiene, y se declara aquí y no al final.
TIPOS_EPITELIALES = ["colonocytes", "enterocytes", "goblet cells", "paneth cells"]

# La regla de selección de marcadores. No estaba en el preregistro y queda
# declarada como enmienda E9: un gen entra si su tipo con más expresión llega a
# NCPM_MINIMO y es al menos ESPECIFICIDAD_MINIMA veces la media de los demás.
NCPM_MINIMO = 10.0
ESPECIFICIDAD_MINIMA = 5.0
ENMIENDA_DE_LA_REGLA = "E9"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _descargar(url: str, destino: Path) -> dict:
    if destino.exists() and destino.stat().st_size > 0:
        return {
            "fichero": destino.name,
            "bytes": destino.stat().st_size,
            "sha256": _sha256(destino),
            "ya_estaba": True,
        }
    with urllib.request.urlopen(url, timeout=900) as r:
        datos = r.read()
    destino.write_bytes(datos)
    return {
        "fichero": destino.name,
        "bytes": len(datos),
        "sha256": hashlib.sha256(datos).hexdigest(),
        "ya_estaba": False,
    }


def _estado_de_pelka() -> dict:
    """Distingue «no resuelve el DNS» de «hace falta una cuenta».

    Son dos cosas distintas y la diferencia importa: con la primera el dato es
    público y solo hace falta otra red; con la segunda haría falta un permiso.
    Confundirlas pondría en el informe una barrera que no existe.
    """
    host = "ftp.ncbi.nlm.nih.gov"
    try:
        socket.gethostbyname(host)
    except OSError as e:
        return {
            **PELKA,
            "estado": "no_comprobable_desde_esta_maquina",
            "motivo": (
                f"{host} no resuelve en el DNS de esta red ({e}). El fichero es "
                "público y no requiere cuenta: lo que falla es la red, no el "
                "permiso."
            ),
            "no_es_requiere_acceso": (
                "No se marca «requiere_acceso» porque sería falso. GEO es "
                "público."
            ),
        }
    return {
        **PELKA,
        "estado": "pendiente_de_descarga",
        "motivo": (
            "El servidor resuelve desde esta máquina y el fichero responde con "
            "código 200, así que la descarga es posible. Este script no la hace "
            "automáticamente: el h5 de la matriz completa son 1,2 GB y un intento "
            "interrumpido deja un fichero truncado que parece válido. Se descarga "
            "aparte, se comprueba el tamaño y se pasa con --pelka."
        ),
    }


def construir_firma(tsv: Path) -> tuple[pd.DataFrame, dict]:
    """Media de expresión por tipo celular, restringida a los tipos del colon."""
    d = pd.read_csv(tsv, sep="\t")
    for c in ("Gene name", "Cell type", "nCPM"):
        if c not in d.columns:
            raise RuntimeError(
                f"el fichero del HPA no tiene la columna «{c}»; cambió de formato "
                "y seguir construiría una firma de otra cosa"
            )
    presentes = sorted(set(d["Cell type"]) & set(TIPOS_DEL_COLON))
    faltan = sorted(set(TIPOS_DEL_COLON) - set(presentes))
    if len(presentes) < 5:
        raise RuntimeError(
            f"solo {len(presentes)} de los {len(TIPOS_DEL_COLON)} tipos "
            "celulares del colon están en esta versión del HPA; con menos de "
            "cinco la deconvolución no separa epitelio de estroma"
        )
    sub = d[d["Cell type"].isin(presentes)]
    firma = sub.pivot_table(
        index="Gene name", columns="Cell type", values="nCPM", aggfunc="mean"
    )
    firma.index = [str(g).strip().upper() for g in firma.index]
    firma = firma[~firma.index.duplicated(keep="first")]
    return firma, {
        "tipos_pedidos": TIPOS_DEL_COLON,
        "tipos_presentes": presentes,
        "tipos_ausentes": faltan,
        "n_genes": int(firma.shape[0]),
        "tipos_epiteliales": [t for t in TIPOS_EPITELIALES if t in presentes],
    }


def seleccionar_marcadores(firma: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Se queda con los genes que distinguen un tipo celular de los demás.

    Un gen expresado igual en todos los tipos no aporta nada a la deconvolución y
    sí aporta ruido, porque engorda el sistema sin añadir información. La regla
    es explícita y queda declarada como enmienda, porque el preregistro no la
    traía.
    """
    maximo = firma.max(axis=1)
    n = firma.shape[1]
    media_de_los_demas = (firma.sum(axis=1) - maximo) / max(n - 1, 1)
    # El uno del denominador evita dividir por cero en genes casi no expresados,
    # que es donde el cociente se dispararía sin significar nada.
    especificidad = maximo / (media_de_los_demas + 1.0)
    entra = (maximo >= NCPM_MINIMO) & (especificidad >= ESPECIFICIDAD_MINIMA)
    marcadores = firma[entra]
    if len(marcadores) < 100:
        raise RuntimeError(
            f"solo {len(marcadores)} genes pasan la regla de marcadores; con "
            "menos de cien la deconvolución no tiene con qué separar los tipos"
        )
    return marcadores, {
        "regla": (
            f"el tipo con más expresión llega a {NCPM_MINIMO} nCPM y es al menos "
            f"{ESPECIFICIDAD_MINIMA} veces la media de los demás tipos"
        ),
        "no_estaba_preregistrada": (
            f"Esta regla no está en config/preregistro_v2.yaml. Queda declarada "
            f"como enmienda {ENMIENDA_DE_LA_REGLA} en docs/enmiendas_v2.md."
        ),
        "ncpm_minimo": NCPM_MINIMO,
        "especificidad_minima": ESPECIFICIDAD_MINIMA,
        "n_marcadores": int(len(marcadores)),
        "marcadores_por_tipo": {
            str(k): int(v)
            for k, v in marcadores.idxmax(axis=1).value_counts().items()
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    ap.add_argument(
        "--pelka",
        default=None,
        help=(
            "carpeta con los ficheros de GSE178341 ya descargados desde otra "
            "red; sin ella la referencia específica de colon queda declarada y "
            "sin usar"
        ),
    )
    args = ap.parse_args()

    ws = Path(args.workspace)
    salida = ws / "external" / "referencia_unicelular"
    (salida / "crudos").mkdir(parents=True, exist_ok=True)
    estado_path = salida / "estado.json"
    pelka = _estado_de_pelka()
    if args.pelka:
        pelka["carpeta_indicada"] = str(args.pelka)
        if not Path(args.pelka).exists():
            pelka["estado"] = "ruta_indicada_no_existe"

    def _fallo(motivo: str) -> int:
        estado_path.write_text(
            json.dumps(
                {
                    "nombre": "referencia_unicelular",
                    "estado": "no_disponible",
                    "motivo": motivo,
                    "fecha_de_acceso": datetime.now(UTC).date().isoformat(),
                    "hpa_url": HPA_URL,
                    "referencia_especifica_de_colon": pelka,
                    "analisis_que_queda_sin_evaluar": (
                        "la deconvolución con referencia externa como R4 "
                        "alternativo (D5)"
                    ),
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"[unicelular] no_disponible: {motivo}", file=sys.stderr)
        return 1

    zip_path = salida / "crudos" / "rna_single_cell_type.tsv.zip"
    try:
        fichero = _descargar(HPA_URL, zip_path)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
        return _fallo(f"fallo al descargar la referencia del HPA: {e}")
    print(f"[unicelular] HPA: {fichero['bytes']} bytes")

    tsv = salida / "crudos" / HPA_INTERNO
    if not tsv.exists():
        try:
            with zipfile.ZipFile(zip_path) as z:
                if HPA_INTERNO not in z.namelist():
                    return _fallo(
                        f"el zip del HPA no contiene {HPA_INTERNO}; trae "
                        f"{z.namelist()}"
                    )
                tsv.write_bytes(z.read(HPA_INTERNO))
        except zipfile.BadZipFile as e:
            return _fallo(f"el fichero del HPA no es un zip válido: {e}")

    try:
        firma, info = construir_firma(tsv)
        marcadores, regla = seleccionar_marcadores(firma)
    except (RuntimeError, ValueError) as e:
        return _fallo(f"la firma no se pudo construir: {e}")

    destino = salida / "firma_tipos_celulares_colon.tsv"
    marcadores.to_csv(destino, sep="\t")
    print(
        f"[unicelular] {regla['n_marcadores']} marcadores para "
        f"{len(info['tipos_presentes'])} tipos celulares."
    )

    estado = {
        "nombre": "referencia_unicelular",
        "estado": "usada",
        "fecha_de_acceso": datetime.now(UTC).date().isoformat(),
        "fuente": "Human Protein Atlas, rna_single_cell_type.tsv",
        "url": HPA_URL,
        "version": HPA_VERSION_INFO,
        "cita": CITA_HPA,
        "licencia": LICENCIA_HPA,
        "fichero_crudo": fichero,
        "firma": destino.name,
        "sha256_firma": _sha256(destino),
        "tipos_celulares": info,
        "seleccion_de_marcadores": regla,
        "referencia_especifica_de_colon": pelka,
        "para_que_sirve": (
            "Es la referencia externa que el preregistro lista como método de "
            "sensibilidad de C7 (deconvolucion_referencia_externa). Al venir de "
            "fuera del proteoma que se contrasta, rompe la circularidad que C7 "
            "corrigió en el índice estromal."
        ),
        "limitaciones": [
            "La referencia es de ARN y lo que se deconvuelve es proteína. D2 "
            "midió el acuerdo entre transcrito y proteína del mismo gen en "
            "CEACAM5 y salió una correlación de rangos de 0,24. Una firma de ARN "
            "aplicada a un proteoma es una aproximación.",
            "Los tipos celulares del atlas son de todos los tejidos juntos. Un "
            "fibroblasto de piel y uno de colon no tienen el mismo perfil, así "
            "que «fibroblastos» aquí es un promedio entre tejidos.",
            "La referencia específica de colon (GSE178341) no se ha podido "
            "descargar desde esta red y queda declarada, no usada.",
            "Varios tipos celulares comparten marcadores (fibroblastos y músculo "
            "liso, entre otros). Eso no es un defecto de los datos sino del "
            "problema: con perfiles parecidos, cualquier método reparte la señal "
            "entre ellos de forma inestable.",
            "La regla de selección de marcadores no estaba preregistrada y se "
            "declara como enmienda E9.",
        ],
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
    }
    estado_path.write_text(
        json.dumps(estado, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[unicelular] referencia de colon: {pelka['estado']}")
    print(f"[unicelular] usada: {destino}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
