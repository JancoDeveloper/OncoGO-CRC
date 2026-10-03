#!/usr/bin/env python3
"""D2. Proteoma retrospectivo de colon de TCGA/CPTAC (Zhang et al. 2014).

Para qué sirve esto. R4 mide cuánto del efecto tumor-normal sobrevive al igualar
la composición del tejido, y lo hace con un índice estromal construido del propio
proteoma. El preregistro obliga a llamarlo «una cota del efecto de composición, no
una pureza medida, salvo cuando la covariable es la pureza ABSOLUTE de TCGA». Esta
cohorte es la que permite quitar ese «salvo»: son pacientes de TCGA, así que
tienen proteoma y pureza medida por ADN en la misma persona. Con eso se puede
preguntar lo que hasta ahora no se había preguntado: ¿el índice estromal mide
pureza, o mide otra cosa?

Y de paso cierra lo que D1 dejó abierto. D1 encontró que el CEACAM5 del tejido
explica solo el 2,4 % de la varianza del CEA en sangre, pero lo midió con mRNA
porque es lo único que TCGA tiene. Aquí hay la proteína, en los mismos pacientes.

Tres discrepancias con el enunciado, encontradas al consultar la API del PDC en
vez de suponer el identificador. Quedan escritas porque cada una habría producido
un resultado de aspecto normal y equivocado.

1. **El identificador no era adivinable.** El estudio es ``PDC000111``,
   ``TCGA COAD Proteome S016-1``, de CPTAC2 Retrospective y label free. El
   candidato que más se parece por el nombre, ``PDC000289``
   («TCGA Colorectal Cancer CompRef Samples»), son muestras de referencia para
   control de calidad y no pacientes: usarlo habría dado una cohorte de
   duplicados técnicos presentada como cohorte de personas.
2. **Esta cohorte no tiene normales.** Son 95 alícuotas de 90 pacientes y las 95
   son tumor primario. Así que D2 **no puede recalcular R4**: sin pares no hay
   contraste tumor-normal. Lo que puede hacer es validar la covariable de R4,
   que es una pregunta distinta y se declara como tal.
3. **``PDC000109`` no es una cohorte independiente.** Es «Prospective COAD
   Proteome S037-1», label free, y 97 de sus 100 pacientes son los mismos
   pacientes de ``PDC000116``, la cohorte de descubrimiento del certificado. Es
   la misma gente recuantificada por otro laboratorio con otro método. Si entrara
   como réplica, contaría a esas personas dos veces y el acuerdo entre cohortes
   saldría altísimo por construcción.

Qué se descarga. Solo las matrices derivadas del CDAP, que son tres ficheros de
pocos megas. Los espectros crudos del estudio son más de un terabyte y no hacen
falta para nada de esto.

Ningún dato de este script identifica a una persona: son los códigos públicos de
TCGA y sus valores agregados por muestra, tal como los publica el consorcio.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

GRAPHQL = "https://proteomic.datacommons.cancer.gov/graphql"
PDC_STUDY_ID = "PDC000111"
PORTAL = f"https://pdc.cancer.gov/pdc/study/{PDC_STUDY_ID}"
DOI = "10.1038/nature13438"
CITA = (
    "Zhang B, Wang J, Wang X, et al. Proteogenomic characterization of human "
    "colon and rectal cancer. Nature 2014; 513(7518):382-387. "
    "doi:10.1038/nature13438"
)
LICENCIA = (
    "Datos abiertos del Proteomic Data Commons del NCI. Las matrices derivadas "
    "del CDAP son de acceso público y no requieren acuerdo de uso de datos."
)

# El estudio hermano que no se puede usar como réplica, y por qué.
ESTUDIO_NO_INDEPENDIENTE = "PDC000109"

# La cuantificación que se usa. «Unshared» son los péptidos que solo pueden
# venir de ese gen; es la misma decisión que v1 tomó en CPTAC con la columna
# «Unshared Log Ratio», y mantenerla evita comparar dos criterios distintos.
FICHERO = "TCGA_Colon_VU_Proteome_CDAP.r2.precursor_area.tsv"
SUFIJO = " Unshared Area"


def _consulta(query: str, *, timeout: int = 180) -> dict:
    """GraphQL por GET. El POST de este servidor rechaza algunas consultas."""
    url = f"{GRAPHQL}?{urllib.parse.urlencode({'query': query})}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        cuerpo = json.loads(r.read().decode("utf-8"))
    if "errors" in cuerpo:
        raise RuntimeError(f"la API del PDC rechazó la consulta: {cuerpo['errors']}")
    return cuerpo["data"]


def identificar_el_estudio() -> dict:
    """Comprueba el identificador contra la API en vez de darlo por bueno."""
    d = _consulta(
        "{study(pdc_study_id:\"%s\" acceptDUA:true){study_id study_submitter_id "
        "project_name disease_type analytical_fraction experiment_type "
        "cases_count aliquots_count}}" % PDC_STUDY_ID
    )
    estudios = d.get("study") or []
    if len(estudios) != 1:
        raise RuntimeError(
            f"{PDC_STUDY_ID} devolvió {len(estudios)} estudios y se esperaba uno"
        )
    return estudios[0]


def tipos_de_muestra(study_id: str) -> dict:
    """Cuenta tumores y normales. Sin normales no hay contraste tumor-normal."""
    d = _consulta(
        '{biospecimenPerStudy(study_id:"%s"){aliquot_submitter_id '
        "case_submitter_id sample_type disease_type}}" % study_id
    )
    b = d.get("biospecimenPerStudy") or []
    if not b:
        raise RuntimeError("la API no devolvió biospecímenes para este estudio")
    tipos: dict[str, int] = {}
    enfermedades: dict[str, int] = {}
    for x in b:
        tipos[x["sample_type"]] = tipos.get(x["sample_type"], 0) + 1
        enfermedades[x["disease_type"]] = enfermedades.get(x["disease_type"], 0) + 1
    casos = sorted({x["case_submitter_id"] for x in b})
    normales = sum(v for k, v in tipos.items() if "Normal" in k)
    return {
        "n_alicuotas": len(b),
        "n_casos": len(casos),
        "por_tipo_de_muestra": tipos,
        "por_enfermedad": enfermedades,
        "n_normales": normales,
        "hay_pares_tumor_normal": normales > 0,
        "consecuencia_de_no_tener_normales": (
            "Sin normales no hay contraste tumor-normal, así que esta cohorte no "
            "recalcula R4. Valida la covariable de R4 contra una pureza medida, "
            "que es una pregunta distinta."
        )
        if normales == 0
        else "",
        "replicas_tecnicas": len(b) - len(casos),
        "nota_de_replicas": (
            "Las alícuotas que sobran son segundas medidas del mismo tumor. "
            "Sirven para estimar cuánto ruido tiene el ensayo, que es lo que "
            "decide si una correlación débil es biología o es el método."
        ),
    }


def comprobar_la_independencia() -> dict:
    """Comprueba si el estudio hermano comparte pacientes con el descubrimiento.

    Es la comprobación que convierte una sospecha en un dato. Dos estudios con
    nombres distintos y laboratorios distintos pueden ser la misma gente, y la
    única forma de saberlo es cruzar los identificadores de caso.
    """
    try:
        hermano = _consulta(
            '{study(pdc_study_id:"%s" acceptDUA:true){study_id '
            "study_submitter_id experiment_type}}" % ESTUDIO_NO_INDEPENDIENTE
        )["study"][0]
        descubrimiento = _consulta(
            '{study(pdc_study_id:"PDC000116" acceptDUA:true){study_id '
            "study_submitter_id}}"
        )["study"][0]
        casos_h = {
            x["case_submitter_id"]
            for x in _consulta(
                '{biospecimenPerStudy(study_id:"%s"){case_submitter_id}}'
                % hermano["study_id"]
            )["biospecimenPerStudy"]
        }
        casos_d = {
            x["case_submitter_id"]
            for x in _consulta(
                '{biospecimenPerStudy(study_id:"%s"){case_submitter_id}}'
                % descubrimiento["study_id"]
            )["biospecimenPerStudy"]
        }
    except (RuntimeError, urllib.error.URLError, KeyError, IndexError, TimeoutError) as e:
        return {"estado": "no_comprobable", "motivo": str(e)}

    compartidos = casos_h & casos_d
    return {
        "estado": "comprobado",
        "estudio": ESTUDIO_NO_INDEPENDIENTE,
        "nombre": hermano["study_submitter_id"],
        "metodo": hermano["experiment_type"],
        "n_casos": len(casos_h),
        "n_compartidos_con_pdc000116": len(compartidos),
        "es_independiente": len(compartidos) == 0,
        "veredicto": (
            f"{len(compartidos)} de los {len(casos_h)} pacientes de "
            f"{ESTUDIO_NO_INDEPENDIENTE} son los mismos de PDC000116. No es una "
            "cohorte independiente: es la misma gente medida por otro "
            "laboratorio con otro método. No entra como réplica en ningún sitio."
        )
        if compartidos
        else (
            f"{ESTUDIO_NO_INDEPENDIENTE} no comparte pacientes con PDC000116."
        ),
    }


def url_firmada(study_id: str) -> str:
    """Las matrices del CDAP necesitan una URL firmada; la base devuelve 403."""
    d = _consulta(
        '{filesPerStudy(study_id:"%s" acceptDUA:true data_category:"Protein '
        'Assembly"){file_name signedUrl{url}}}' % study_id
    )
    for f in d.get("filesPerStudy") or []:
        if f["file_name"] == FICHERO:
            url = (f.get("signedUrl") or {}).get("url")
            if url:
                return url
            raise RuntimeError(f"{FICHERO} existe pero la API no dio URL firmada")
    raise RuntimeError(
        f"{FICHERO} no está entre los ficheros del estudio; el CDAP cambió de "
        "nombre de fichero y hay que volver a mirar cuál es la matriz derivada"
    )


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for bloque in iter(lambda: fh.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()


def _descargar(url: str, destino: Path) -> dict:
    """Idempotente: si el fichero ya está, no se vuelve a bajar."""
    if destino.exists() and destino.stat().st_size > 0:
        return {
            "fichero": destino.name,
            "bytes": destino.stat().st_size,
            "sha256": _sha256(destino),
            "ya_estaba": True,
        }
    req = urllib.request.Request(url, headers={"Accept": "*/*"})
    with urllib.request.urlopen(req, timeout=600) as r:
        datos = r.read()
    destino.write_bytes(datos)
    return {
        "fichero": destino.name,
        "bytes": len(datos),
        "sha256": hashlib.sha256(datos).hexdigest(),
        "ya_estaba": False,
    }


def parsear(crudo: Path) -> tuple[pd.DataFrame, dict]:
    """Matriz de genes por alícuota, en log2, con el cero tratado como ausencia.

    El fichero del CDAP está en área de precursor, que es una intensidad lineal.
    El log2 va aquí porque la abundancia de proteína se compara multiplicando, y
    el cero se convierte en ausente porque un área de cero significa que el
    péptido no se detectó, no que haya cero proteína. No se imputa nada.
    """
    d = pd.read_csv(crudo, sep="\t")
    if "Gene" not in d.columns:
        raise RuntimeError(
            "la matriz no tiene columna Gene; el CDAP cambió de formato y el "
            "resto del análisis daría genes mal emparejados sin avisar"
        )
    cols = [c for c in d.columns if c.endswith(SUFIJO)]
    if not cols:
        raise RuntimeError(
            f"ninguna columna acaba en «{SUFIJO}»; el criterio de cuantificación "
            "cambió y hay que decidirlo otra vez en vez de coger la primera que "
            "se parezca"
        )
    m = d.set_index("Gene")[cols].apply(pd.to_numeric, errors="coerce")
    ceros = int((m == 0).sum().sum())
    m = m.mask(m <= 0.0)
    m = pd.DataFrame(np.log2(m.to_numpy(dtype=float)), index=m.index, columns=m.columns)
    m.columns = [c[: -len(SUFIJO)] for c in m.columns]
    m = m[~m.index.duplicated(keep="first")]

    resumen = {
        "n_genes": int(m.shape[0]),
        "n_alicuotas": int(m.shape[1]),
        "cuantificacion": (
            "Área de precursor de péptidos no compartidos, en log2. Es el mismo "
            "criterio de «no compartido» que v1 usó en CPTAC."
        ),
        "ceros_tratados_como_ausencia": ceros,
        "fraccion_de_ausencias": round(ceros / float(m.size), 4),
        "no_se_imputa": (
            "Un área de cero es que el péptido no se detectó. Rellenarla con un "
            "número inventaría una medida que nadie hizo."
        ),
    }
    return m, resumen


def barcodes(columnas: list[str]) -> dict[str, str]:
    """Del nombre de alícuota al paciente de TCGA.

    Las columnas vienen como ``A6-3807-01A-22``, sin el prefijo ``TCGA-`` que sí
    llevan las tablas clínicas. Olvidarlo deja el cruce vacío y parece que no hay
    solape cuando lo que no hay es la misma forma de escribir el código.
    """
    out = {}
    for c in columnas:
        partes = str(c).split("-")
        if len(partes) < 2:
            continue
        out[c] = "TCGA-" + "-".join(partes[:2])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    args = ap.parse_args()

    ws = Path(args.workspace)
    salida = ws / "external" / "tcga_cptac_retro"
    (salida / "crudos").mkdir(parents=True, exist_ok=True)
    estado_path = salida / "estado.json"

    def _fallo(motivo: str) -> int:
        estado_path.write_text(
            json.dumps(
                {
                    "nombre": "tcga_cptac_retro",
                    "estado": "no_disponible",
                    "motivo": motivo,
                    "fecha_de_acceso": datetime.now(UTC).date().isoformat(),
                    "pdc_study_id": PDC_STUDY_ID,
                    "portal": PORTAL,
                    "analisis_que_queda_sin_evaluar": (
                        "la validación del índice estromal de R4 contra la "
                        "pureza ABSOLUTE medida, y el CEACAM5 de proteína "
                        "frente al CEA en sangre (D2)"
                    ),
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"[tcga-retro] no_disponible: {motivo}", file=sys.stderr)
        return 1

    try:
        estudio = identificar_el_estudio()
    except (RuntimeError, urllib.error.URLError, TimeoutError, KeyError) as e:
        return _fallo(f"no se pudo identificar {PDC_STUDY_ID} en la API del PDC: {e}")
    print(f"[tcga-retro] {PDC_STUDY_ID} = {estudio['study_submitter_id']}")

    try:
        muestras = tipos_de_muestra(estudio["study_id"])
    except (RuntimeError, urllib.error.URLError, TimeoutError) as e:
        return _fallo(f"no se pudieron leer los tipos de muestra: {e}")
    print(
        f"[tcga-retro] {muestras['n_alicuotas']} alícuotas de "
        f"{muestras['n_casos']} pacientes; normales: {muestras['n_normales']}"
    )

    try:
        url = url_firmada(estudio["study_id"])
        fichero = _descargar(url, salida / "crudos" / FICHERO)
    except (RuntimeError, urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        return _fallo(f"fallo al descargar la matriz derivada: {e}")
    print(f"[tcga-retro] matriz: {fichero['bytes']} bytes")

    try:
        matriz, resumen = parsear(salida / "crudos" / FICHERO)
    except (RuntimeError, ValueError) as e:
        return _fallo(f"la matriz no se pudo parsear: {e}")

    limpia = salida / "tcga_cptac_retro_proteoma.tsv"
    matriz.to_csv(limpia, sep="\t")
    mapa = barcodes(list(matriz.columns))
    pd.DataFrame(
        {"alicuota": list(mapa), "paciente": [mapa[k] for k in mapa]}
    ).to_csv(salida / "tcga_cptac_retro_alicuotas.tsv", sep="\t", index=False)

    independencia = comprobar_la_independencia()
    if independencia.get("veredicto"):
        print(f"[tcga-retro] {independencia['veredicto']}")

    estado = {
        "nombre": "tcga_cptac_retro",
        "estado": "usada",
        "fecha_de_acceso": datetime.now(UTC).date().isoformat(),
        "pdc_study_id": PDC_STUDY_ID,
        "study_id": estudio["study_id"],
        "study_submitter_id": estudio["study_submitter_id"],
        "portal": PORTAL,
        "doi": DOI,
        "cita": CITA,
        "licencia": LICENCIA,
        "proyecto": estudio["project_name"],
        "metodo": estudio["experiment_type"],
        "fraccion": estudio["analytical_fraction"],
        "fichero_crudo": fichero,
        "matriz_limpia": limpia.name,
        "sha256_matriz_limpia": _sha256(limpia),
        "parseo": resumen,
        "muestras": muestras,
        "por_que_estos_pacientes": (
            "Son pacientes de TCGA, así que tienen pureza medida por ADN "
            "(ABSOLUTE) y CEA sérico prequirúrgico en la misma persona. Los "
            "pacientes de PDC000116 no son de TCGA y no tienen ninguna de las "
            "dos cosas."
        ),
        "independencia_del_estudio_hermano": independencia,
        "discrepancias_con_el_enunciado": [
            "El identificador no era adivinable. El estudio es PDC000111 "
            "(TCGA COAD Proteome S016-1). PDC000289, que por nombre es el que "
            "más se parece, son muestras CompRef de control de calidad y no "
            "pacientes.",
            "Esta cohorte no tiene normales: las 95 alícuotas son tumor "
            "primario. D2 no puede recalcular R4, solo validar su covariable.",
            "PDC000109 no es una cohorte independiente: comparte pacientes con "
            "PDC000116, el descubrimiento del certificado.",
        ],
        "limitaciones": [
            "Es label free y PDC000116 es TMT. El índice estromal que se valida "
            "aquí es el mismo constructo, pero no son los mismos valores: lo que "
            "se comprueba es si la idea del índice mide pureza, no si los "
            "números concretos de R4 están bien.",
            "Son tumores sin normales, así que el índice se calcula entre "
            "tumores y no como diferencia dentro de un par. Es una validación "
            "del constructo, no del contraste.",
            "El 36 % de la matriz son ausencias. En label free eso es normal y "
            "se trata como ausencia, no se imputa.",
            "El CEA sérico solo existe en una parte de estos pacientes, así que "
            "el análisis del CEACAM5 corre sobre menos gente que el de la "
            "pureza. Los dos tamaños se publican por separado.",
        ],
        "no_entra_en_el_certificado": (
            "Nada de esta fuente cambia un subíndice. Lo que hace es decir si el "
            "lenguaje obligatorio del preregistro sobre R4 está justificado."
        ),
    }
    estado_path.write_text(
        json.dumps(estado, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[tcga-retro] usada: {limpia}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
