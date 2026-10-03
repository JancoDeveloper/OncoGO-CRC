#!/usr/bin/env python3
"""D1. TCGA COAD/READ: CEA sérico, pureza medida, MSI, subtipo y CEACAM.

Para qué sirve esto. R4 de v1 y de v2 mide cuánto del efecto sobrevive al igualar
la composición del tejido, pero lo hace con un índice estromal construido del
propio proteoma. Eso es una cota, no una pureza. TCGA sí tiene pureza medida por
ADN (ABSOLUTE), así que permite comparar la cota con una pureza de verdad. Y
además tiene CEA sérico prequirúrgico, que es lo único de todo este trabajo que
conecta con un marcador que se usa en clínica.

Tres discrepancias con el enunciado, encontradas al probar los endpoints en vez
de suponerlos. Quedan aquí escritas porque el nombre equivocado habría dado una
columna vacía y un «no disponible» falso.

1. El campo del CEA no se llama ``preoperative_pretreatment_cea_level``. En el
   BCR Biotab actual del GDC se llama ``cea_level_pretreatment``.
2. ``HEAD`` contra ``api.gdc.cancer.gov/data/<id>`` devuelve 400 aunque el
   fichero exista. Hay que comprobarlo con ``GET``.
3. El tarball de cBioPortal en ``cbioportal-datahub.s3.amazonaws.com`` responde
   403. La API REST de cBioPortal sí funciona y es la que se usa.

Y una limitación que no es un fallo de la descarga. El subtipo CMS no está en
TCGA. Lo que hay es el subtipo molecular de PanCanAtlas (CIN, MSI, GS, POLE),
que no es lo mismo. Las etiquetas CRCSC de CMS viven en Synapse y requieren
cuenta, así que CMS queda ``requiere_acceso`` y se usa el subtipo de PanCanAtlas
diciendo explícitamente que es otra variable.

Nada de lo que escribe este script es un dato de un paciente concreto
identificable: son los identificadores públicos de TCGA y sus valores agregados
por muestra, tal como los publica el consorcio.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data/{}"
CBIO = "https://www.cbioportal.org/api"
ESTUDIO = "coadread_tcga_pan_can_atlas_2018"
LISTA = f"{ESTUDIO}_all"

# Pureza ABSOLUTE de PanCanAtlas. El identificador es el de la página de
# publicaciones del GDC; si cambia, el script lo dirá en lugar de inventarlo.
PUREZA_ABSOLUTE = "4f277128-f793-4354-a13d-30cc7fe9f6b5"

# Entrez de las CEACAM que importan. CEACAM5 es la proteína que el ensayo de CEA
# mide en sangre; las otras dos entran como control de familia.
GENES = {1048: "CEACAM5", 4680: "CEACAM6", 634: "CEACAM1"}

ATRIBUTOS_PACIENTE = [
    "AGE", "SEX", "AJCC_PATHOLOGIC_TUMOR_STAGE", "PATH_M_STAGE",
    "PATH_N_STAGE", "PATH_T_STAGE", "SUBTYPE", "ICD_O_3_SITE",
]
ATRIBUTOS_MUESTRA = [
    "MSI_SCORE_MANTIS", "MSI_SENSOR_SCORE", "ANEUPLOIDY_SCORE",
    "TMB_NONSYNONYMOUS", "TUMOR_TISSUE_SITE",
]

CITA = (
    "The Cancer Genome Atlas Network. Comprehensive molecular characterization of "
    "human colon and rectal cancer. Nature 2012;487:330-7. "
    "Pureza: Carter SL et al. Nat Biotechnol 2012 (ABSOLUTE), vía PanCanAtlas. "
    "Acceso vía API del GDC y API de cBioPortal."
)
LICENCIA = (
    "Datos abiertos de TCGA (NIH Genomic Data Commons, nivel 3 agregado). "
    "cBioPortal redistribuye bajo ODC-ODbL."
)


def _get(url: str, *, cuerpo: dict | None = None, timeout: int = 180) -> bytes:
    datos = json.dumps(cuerpo).encode() if cuerpo is not None else None
    req = urllib.request.Request(
        url,
        data=datos,
        headers={
            "Accept": "application/json",
            "User-Agent": "OncoGO-CRC/2 (investigacion academica)",
            **({"Content-Type": "application/json"} if datos else {}),
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _clinico_biotab(destino: Path) -> tuple[pd.DataFrame | None, list[dict], str]:
    """El BCR Biotab de paciente de COAD y de READ, que es donde vive el CEA."""
    consulta = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "in", "content": {
                    "field": "cases.project.project_id",
                    "value": ["TCGA-COAD", "TCGA-READ"]}},
                {"op": "in", "content": {
                    "field": "data_type", "value": ["Clinical Supplement"]}},
                {"op": "in", "content": {
                    "field": "data_format", "value": ["BCR Biotab"]}},
            ],
        },
        "fields": "file_id,file_name,md5sum,file_size",
        "size": "60",
    }
    try:
        listado = json.loads(_get(GDC_FILES, cuerpo=consulta))
    except (urllib.error.URLError, json.JSONDecodeError, TimeoutError) as e:
        return None, [], f"el listado del GDC falló: {e}"

    hits = listado.get("data", {}).get("hits", [])
    queridos = [
        h for h in hits
        if h["file_name"].startswith("nationwidechildrens.org_clinical_patient_")
    ]
    if not queridos:
        return None, [], "el GDC no devolvió ningún clinical_patient en BCR Biotab"

    tablas, fichas = [], []
    for h in queridos:
        try:
            crudo = _get(GDC_DATA.format(h["file_id"]))
        except (urllib.error.URLError, TimeoutError) as e:
            fichas.append({"fichero": h["file_name"], "estado": f"falló: {e}"})
            continue
        (destino / h["file_name"]).write_bytes(crudo)
        # Las dos primeras filas tras la cabecera son descripciones, no pacientes.
        tabla = pd.read_csv(io.BytesIO(crudo), sep="\t", low_memory=False, skiprows=[1, 2])
        tabla["proyecto"] = (
            "TCGA-READ" if h["file_name"].endswith("_read.txt") else "TCGA-COAD"
        )
        tablas.append(tabla)
        fichas.append(
            {
                "fichero": h["file_name"],
                "file_id": h["file_id"],
                "md5_declarado_por_gdc": h["md5sum"],
                "sha256_descargado": _sha256(crudo),
                "n_filas": int(len(tabla)),
                "estado": "usada",
            }
        )
    if not tablas:
        return None, fichas, "ninguna descarga del clínico terminó"
    return pd.concat(tablas, ignore_index=True), fichas, ""


def _pureza(destino: Path) -> tuple[pd.DataFrame | None, dict, str]:
    """Pureza ABSOLUTE de PanCanAtlas, filtrada a COAD y READ más adelante."""
    try:
        crudo = _get(GDC_DATA.format(PUREZA_ABSOLUTE))
    except (urllib.error.URLError, TimeoutError) as e:
        return None, {}, f"la pureza ABSOLUTE no se pudo descargar: {e}"
    cabecera = crudo[:200].decode("utf-8", "replace")
    if "purity" not in cabecera:
        return None, {}, (
            "el fichero de pureza no tiene columna 'purity'; el identificador de "
            "PanCanAtlas puede haber cambiado"
        )
    (destino / "TCGA_mastercalls.abs_purity.txt").write_bytes(crudo)
    tabla = pd.read_csv(io.BytesIO(crudo), sep="\t", low_memory=False)
    ficha = {
        "fichero": "TCGA_mastercalls.abs_purity.txt",
        "file_id": PUREZA_ABSOLUTE,
        "sha256": _sha256(crudo),
        "n_filas_pancan": int(len(tabla)),
        "estado": "usada",
    }
    return tabla, ficha, ""


def _clinico_cbioportal() -> tuple[pd.DataFrame, list[str]]:
    """MSI, subtipo molecular y estadio desde la API de cBioPortal."""
    filas, fallos = [], []
    for tipo, atributos in (("PATIENT", ATRIBUTOS_PACIENTE), ("SAMPLE", ATRIBUTOS_MUESTRA)):
        for attr in atributos:
            url = (
                f"{CBIO}/studies/{ESTUDIO}/clinical-data"
                f"?clinicalDataType={tipo}&attributeId={attr}&pageSize=5000"
            )
            try:
                datos = json.loads(_get(url))
            except (urllib.error.URLError, json.JSONDecodeError, TimeoutError) as e:
                fallos.append(f"{attr}: {e}")
                continue
            for d in datos:
                filas.append(
                    {
                        "paciente": d["patientId"],
                        "muestra": d.get("sampleId"),
                        "atributo": attr,
                        "valor": d["value"],
                    }
                )
    if not filas:
        return pd.DataFrame(), fallos
    largo = pd.DataFrame(filas)
    # Los atributos de paciente se expanden a todas sus muestras al unir después.
    ancho = (
        largo.pivot_table(
            index=["paciente", "muestra"], columns="atributo", values="valor",
            aggfunc="first",
        )
        .reset_index()
    )
    de_paciente = (
        largo.loc[largo["muestra"].isna()]
        .pivot_table(index="paciente", columns="atributo", values="valor", aggfunc="first")
        .reset_index()
    )
    de_muestra = ancho.loc[ancho["muestra"].notna()].drop(
        columns=[c for c in ATRIBUTOS_PACIENTE if c in ancho.columns], errors="ignore"
    )
    return de_muestra.merge(de_paciente, on="paciente", how="outer"), fallos


def _expresion() -> tuple[pd.DataFrame, str]:
    """RNA-seq de las CEACAM, por muestra."""
    url = f"{CBIO}/molecular-profiles/{ESTUDIO}_rna_seq_v2_mrna/molecular-data/fetch"
    try:
        datos = json.loads(
            _get(url, cuerpo={"entrezGeneIds": list(GENES), "sampleListId": LISTA})
        )
    except (urllib.error.URLError, json.JSONDecodeError, TimeoutError) as e:
        return pd.DataFrame(), f"la expresión no se pudo descargar: {e}"
    if not isinstance(datos, list) or not datos:
        return pd.DataFrame(), "la API de expresión devolvió una respuesta vacía"
    tabla = pd.DataFrame(
        [
            {
                "muestra": d["sampleId"],
                "paciente": d["patientId"],
                "gen": GENES.get(d["entrezGeneId"], str(d["entrezGeneId"])),
                "rnaseq_rsem": d["value"],
            }
            for d in datos
        ]
    )
    return (
        tabla.pivot_table(
            index=["paciente", "muestra"], columns="gen", values="rnaseq_rsem",
            aggfunc="first",
        )
        .add_prefix("rnaseq_")
        .reset_index(),
        "",
    )


def _tabla_limpia(
    clinico: pd.DataFrame | None,
    pureza: pd.DataFrame | None,
    cbio: pd.DataFrame,
    expresion: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    """Una fila por muestra tumoral, con lo que haya de cada fuente."""
    base = expresion.copy() if not expresion.empty else cbio[["paciente", "muestra"]].copy()
    if not cbio.empty:
        base = base.merge(cbio, on=["paciente", "muestra"], how="outer")

    if clinico is not None and "bcr_patient_barcode" in clinico.columns:
        columnas = {"bcr_patient_barcode": "paciente"}
        for origen, nuevo in (
            ("cea_level_pretreatment", "cea_pretratamiento_ng_ml"),
            ("ajcc_pathologic_tumor_stage", "estadio_bcr"),
            ("tumor_tissue_site", "sitio_bcr"),
            ("proyecto", "proyecto"),
        ):
            if origen in clinico.columns:
                columnas[origen] = nuevo
        recorte = clinico[list(columnas)].rename(columns=columnas)
        recorte["cea_pretratamiento_ng_ml"] = pd.to_numeric(
            recorte.get("cea_pretratamiento_ng_ml"), errors="coerce"
        )
        base = base.merge(recorte.drop_duplicates("paciente"), on="paciente", how="left")

    if pureza is not None and "purity" in pureza.columns:
        p = pureza.copy()
        p["muestra"] = p["array"].astype(str).str.strip()
        p = p[["muestra", "purity", "ploidy", "Cancer DNA fraction"]].rename(
            columns={
                "purity": "pureza_absolute",
                "ploidy": "ploidia_absolute",
                "Cancer DNA fraction": "fraccion_adn_tumoral",
            }
        )
        for col in ("pureza_absolute", "ploidia_absolute", "fraccion_adn_tumoral"):
            p[col] = pd.to_numeric(p[col], errors="coerce")
        base = base.merge(p.drop_duplicates("muestra"), on="muestra", how="left")

    for col in ("MSI_SCORE_MANTIS", "MSI_SENSOR_SCORE", "AGE", "ANEUPLOIDY_SCORE",
                "TMB_NONSYNONYMOUS"):
        if col in base.columns:
            base[col] = pd.to_numeric(base[col], errors="coerce")

    base = base.rename(
        columns={
            "SUBTYPE": "subtipo_pancanatlas",
            "MSI_SCORE_MANTIS": "msi_mantis",
            "MSI_SENSOR_SCORE": "msi_sensor",
            "AJCC_PATHOLOGIC_TUMOR_STAGE": "estadio_ajcc",
            "PATH_M_STAGE": "estadio_m",
            "PATH_N_STAGE": "estadio_n",
            "PATH_T_STAGE": "estadio_t",
            "SEX": "sexo",
            "AGE": "edad",
            "TUMOR_TISSUE_SITE": "sitio_del_tumor",
            "ICD_O_3_SITE": "icd_o_3_sitio",
            "ANEUPLOIDY_SCORE": "aneuploidia",
            "TMB_NONSYNONYMOUS": "tmb_no_sinonima",
        }
    )
    # Solo muestras de tumor primario: el sufijo -01 del código de TCGA.
    base = base.loc[base["muestra"].astype(str).str.endswith("-01")]
    base = base.sort_values("muestra").reset_index(drop=True)

    def cuenta(col: str) -> int:
        return int(base[col].notna().sum()) if col in base.columns else 0

    cobertura = {
        "n_muestras": int(len(base)),
        "n_con_cea": cuenta("cea_pretratamiento_ng_ml"),
        "n_con_pureza_absolute": cuenta("pureza_absolute"),
        "n_con_cea_y_pureza": int(
            (
                base.get("cea_pretratamiento_ng_ml", pd.Series(dtype=float)).notna()
                & base.get("pureza_absolute", pd.Series(dtype=float)).notna()
            ).sum()
        ),
        "n_con_rnaseq_ceacam5": cuenta("rnaseq_CEACAM5"),
        "n_con_msi_mantis": cuenta("msi_mantis"),
        "n_con_subtipo": cuenta("subtipo_pancanatlas"),
    }
    if cobertura["n_muestras"]:
        cobertura["fraccion_con_cea"] = cobertura["n_con_cea"] / cobertura["n_muestras"]
    return base, cobertura


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    args = ap.parse_args()

    salida = Path(args.workspace) / "external" / "tcga_coadread"
    crudos = salida / "crudos"
    crudos.mkdir(parents=True, exist_ok=True)

    estado: dict = {
        "nombre": "tcga_coadread",
        "fecha_de_acceso": datetime.now(UTC).strftime("%Y-%m-%d"),
        "cita": CITA,
        "licencia": LICENCIA,
        "fuentes": {},
        "discrepancias_con_el_enunciado": [
            "El campo del CEA es cea_level_pretreatment, no "
            "preoperative_pretreatment_cea_level.",
            "HEAD contra api.gdc.cancer.gov/data/<id> devuelve 400 aunque el fichero "
            "exista; hay que usar GET.",
            "El tarball de cbioportal-datahub responde 403; se usa la API REST.",
        ],
        "cms": {
            "estado": "requiere_acceso",
            "motivo": (
                "Las etiquetas CMS del CRC Subtyping Consortium están en Synapse "
                "(syn2623706) y requieren cuenta. No se intenta autenticar."
            ),
            "sustituto_usado": (
                "Subtipo molecular de PanCanAtlas (CIN, MSI, GS, POLE). No es CMS y "
                "no se interpreta como si lo fuera."
            ),
        },
    }

    print("[tcga] clínico BCR Biotab (el CEA vive aquí)…")
    clinico, fichas, err = _clinico_biotab(crudos)
    estado["fuentes"]["clinico_bcr_biotab"] = {
        "estado": "usada" if clinico is not None else "no_disponible",
        "motivo": err,
        "ficheros": fichas,
    }

    print("[tcga] pureza ABSOLUTE de PanCanAtlas…")
    pureza, ficha_p, err_p = _pureza(crudos)
    estado["fuentes"]["pureza_absolute"] = {
        "estado": "usada" if pureza is not None else "no_disponible",
        "motivo": err_p,
        **({"fichero": ficha_p} if ficha_p else {}),
    }

    print("[tcga] clínico y subtipos de cBioPortal…")
    cbio, fallos = _clinico_cbioportal()
    estado["fuentes"]["cbioportal_clinico"] = {
        "estado": "usada" if not cbio.empty else "no_disponible",
        "estudio": ESTUDIO,
        "atributos_pedidos": ATRIBUTOS_PACIENTE + ATRIBUTOS_MUESTRA,
        "atributos_que_fallaron": fallos,
    }

    print("[tcga] expresión de las CEACAM…")
    expresion, err_e = _expresion()
    estado["fuentes"]["cbioportal_rnaseq"] = {
        "estado": "usada" if not expresion.empty else "no_disponible",
        "motivo": err_e,
        "genes": list(GENES.values()),
        "perfil": f"{ESTUDIO}_rna_seq_v2_mrna",
    }

    if clinico is None and cbio.empty and expresion.empty:
        estado["estado"] = "no_disponible"
        estado["motivo"] = "ninguna de las fuentes de D1 respondió"
        (salida / "estado.json").write_text(
            json.dumps(estado, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        print("[tcga] ninguna fuente respondió: queda no_disponible", file=sys.stderr)
        return 1

    tabla, cobertura = _tabla_limpia(clinico, pureza, cbio, expresion)
    destino = salida / "tcga_coadread_cea_pureza.tsv"
    tabla.to_csv(destino, sep="\t", index=False)

    estado["estado"] = "usada"
    estado["tabla_limpia"] = destino.name
    estado["sha256_de_la_tabla_limpia"] = _sha256(destino.read_bytes())
    estado["columnas"] = list(tabla.columns)
    estado["cobertura"] = cobertura
    estado["limitaciones"] = [
        "El CEA de TCGA es prequirúrgico y de una sola extracción. No es una serie "
        "de seguimiento, así que no dice nada sobre la cinética del marcador.",
        "El valor del CEA puede estar censurado por el límite inferior del ensayo y "
        "el límite no viene en la tabla. Cualquier modelo tiene que tratarlo como "
        "censura, no como un número exacto.",
        "La pureza ABSOLUTE es de ADN de la misma pieza, no de la sección exacta de "
        "la que sale el proteoma. Es una pureza medida, pero de otra alícuota.",
        "El subtipo disponible es el de PanCanAtlas, no CMS.",
    ]
    (salida / "estado.json").write_text(
        json.dumps(estado, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(f"[tcga] {cobertura['n_muestras']} muestras de tumor primario.")
    print(
        f"[tcga] CEA en {cobertura['n_con_cea']} "
        f"({cobertura.get('fraccion_con_cea', 0):.1%}), "
        f"pureza en {cobertura['n_con_pureza_absolute']}, "
        f"las dos en {cobertura['n_con_cea_y_pureza']}."
    )
    print(f"[tcga] escrito en {salida}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
