"""Ventana de detección tejido→plasma con concentraciones del HPA (C8).

v1 concluyó que la dirección del tejido no se transfiere al suero. La conclusión
se midió sobre todas las proteínas presentes en las dos matrices, sin separar
las que el método alcanza de las que no. Aquí se estima la curva de detección
frente a la concentración en sangre, se localiza el umbral empírico y la
concordancia se mide dentro y fuera de la ventana, por separado.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.quant_columns import is_summary_protein
from oncogo_crc.utils.reproducibility_overlap import (
    load_gene_uniprot_from_mapping,
    uniprot_base,
)
from oncogo_crc.utils.ids import effective_protein_mapping
from oncogo_crc.v2 import deteccion
from oncogo_crc.v2.prereg import cargar

# Mínimos propios de C8. No estaban en el preregistro con un número, así que se
# reutiliza el mínimo general de términos y se declara aquí de dónde sale.
N_MIN_CURVA = 50
N_MIN_VENTANA = 10


def _concentracion_hpa(path, gene_to_uni: dict[str, str]) -> tuple[pd.Series, dict]:
    """Concentración en sangre por UniProt, prefiriendo espectrometría de masas.

    Se prefiere la medida por espectrometría porque es la tecnología del propio
    experimento de suero. El inmunoensayo se usa solo cuando la otra falta.
    """
    hpa = pd.read_csv(path, sep="\t", low_memory=False)
    ms = pd.to_numeric(hpa.get("conc_sangre_ms_pg_l"), errors="coerce")
    im = pd.to_numeric(hpa.get("conc_sangre_inmunoensayo_pg_l"), errors="coerce")
    conc = ms.where(ms.notna(), im)
    origen = np.where(ms.notna(), "espectrometria", np.where(im.notna(), "inmunoensayo", ""))
    hpa = hpa.assign(_conc=conc, _origen=origen)

    # El HPA trae su propio UniProt; si falta, se traduce desde el símbolo.
    uni = hpa.get("uniprot", pd.Series("", index=hpa.index)).fillna("").astype(str)
    uni = uni.map(lambda v: uniprot_base(str(v).split(",")[0].strip()) or "")
    desde_gen = (
        hpa.get("gene", pd.Series("", index=hpa.index))
        .fillna("")
        .astype(str)
        .str.upper()
        .map(lambda g: gene_to_uni.get(g, ""))
    )
    hpa["_uni"] = np.where(uni.str.len() > 0, uni, desde_gen)

    util = hpa.loc[(hpa["_uni"].astype(str).str.len() > 0) & hpa["_conc"].notna()]
    util = util.sort_values("_conc", ascending=False, kind="stable").drop_duplicates("_uni")
    serie = util.set_index("_uni")["_conc"].astype(float)
    info = {
        "n_genes_en_el_hpa": int(len(hpa)),
        "n_con_concentracion": int(hpa["_conc"].notna().sum()),
        "n_con_concentracion_y_uniprot": int(len(serie)),
        "n_por_espectrometria": int((util["_origen"] == "espectrometria").sum()),
        "n_por_inmunoensayo": int((util["_origen"] == "inmunoensayo").sum()),
        "preferencia": "espectrometría cuando existe; inmunoensayo solo si falta",
    }
    return serie, info


def _medidas_en_plasma(path) -> set[str]:
    """UniProt de todas las proteínas que la matriz de suero llegó a cuantificar."""
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    raw = raw.loc[~raw["protein_id"].map(is_summary_protein)]
    columnas = [c for c in raw.columns if c != "protein_id"]
    valores = raw[columnas].apply(pd.to_numeric, errors="coerce")
    medida = (valores.notna() & (valores > 0)).any(axis=1)
    return {
        uniprot_base(str(p)) or ""
        for p in raw.loc[medida, "protein_id"]
    } - {""}


def run(settings: Settings) -> None:
    prereg = cargar(settings.root)
    salida = settings.resolve("results") / "v2" / "deteccion"
    salida.mkdir(parents=True, exist_ok=True)

    ficha_path = settings.resolve("external") / "hpa" / "estado.json"
    ficha = json.loads(ficha_path.read_text(encoding="utf-8")) if ficha_path.exists() else {}
    estado_hpa = ficha.get("hpa", {}).get("estado", "no_disponible")
    if estado_hpa != "usada":
        auditoria = {
            "version": "v2",
            "estado": "no_disponible",
            "motivo": ficha.get("hpa", {}).get(
                "motivo", "no hay tabla del Human Protein Atlas en external/hpa"
            ),
            "consecuencia": (
                "C8 queda sin evaluar. Sin concentraciones en sangre no se puede "
                "distinguir «no se detecta» de «no está», así que la brecha "
                "tejido→plasma de v1 se queda como estaba, sin reinterpretar."
            ),
            "como_obtenerlo": "python scripts/download/human_protein_atlas.py",
        }
        (salida / "deteccion_audit_v2.json").write_text(
            json.dumps(auditoria, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        print(f"[v2-deteccion] no disponible: {auditoria['motivo']}")
        return

    mapping = effective_protein_mapping(settings)
    gene_to_uni = load_gene_uniprot_from_mapping(mapping, exclude_ambiguous=True)
    conc, info_hpa = _concentracion_hpa(
        settings.resolve("external") / "hpa" / "hpa_concentracion_sangre.tsv", gene_to_uni
    )
    print(
        f"[v2-deteccion] HPA: {info_hpa['n_con_concentracion_y_uniprot']} proteínas con "
        f"concentración y UniProt ({info_hpa['n_por_espectrometria']} por espectrometría)."
    )

    plasma_ds = settings.datasets["datasets"]["plasma_crc_2024"]
    medidas = _medidas_en_plasma(settings.resolve(plasma_ds["matrix_path"]))
    print(f"[v2-deteccion] la matriz de suero cuantifica {len(medidas)} proteínas.")

    # La curva de detección se estima sobre el universo del HPA: de las proteínas
    # con concentración conocida, ¿cuáles llegó a medir el experimento de suero?
    detectada = pd.Series(
        [u in medidas for u in conc.index], index=conc.index, dtype=bool
    )
    curva = deteccion.curva_de_deteccion(conc, detectada, n_min=N_MIN_CURVA)
    if curva["estado"] == "calculado":
        print(
            f"[v2-deteccion] curva: pendiente {curva['pendiente_por_orden_de_magnitud']:.3f} "
            f"por orden de magnitud (z={curva['z_pendiente']:.1f}), umbral "
            f"{curva['umbral_pg_l']:.3g} pg/L"
        )
    else:
        print(f"[v2-deteccion] curva no evaluable: {curva['motivo']}")

    cert = settings.resolve("results") / "v2" / "certified"
    cptac = pd.read_csv(cert / "cptac_unshared_treat.tsv", sep="\t")
    tejido = cptac.loc[cptac["significant"]].copy()
    tejido["_uni"] = tejido["protein_id"].astype(str).str.upper().map(gene_to_uni)
    tejido = tejido.dropna(subset=["_uni"]).drop_duplicates("_uni")
    lfc_tejido = tejido.set_index("_uni")["log2fc"].astype(float)

    plasma_tabla = pd.read_csv(
        settings.resolve("results") / "secondary" / "plasma_crc_vs_healthy_treat.tsv", sep="\t"
    )
    plasma_tabla["_uni"] = plasma_tabla["protein_id"].map(lambda v: uniprot_base(str(v)) or "")
    plasma_tabla = plasma_tabla.loc[plasma_tabla["_uni"].str.len() > 0].drop_duplicates("_uni")
    lfc_plasma = plasma_tabla.set_index("_uni")["log2fc"].astype(float)

    umbral = curva.get("umbral_log10_pg_l") if curva["estado"] == "calculado" else None
    ventana = deteccion.concordancia_dentro_de_la_ventana(
        lfc_tejido, lfc_plasma, conc, umbral_log10=umbral, n_min=N_MIN_VENTANA
    )
    print(f"[v2-deteccion] todas: {ventana['todas']['valor']} ({ventana['todas']['denominador']})")
    if "dentro_de_la_ventana" in ventana:
        d = ventana["dentro_de_la_ventana"]
        f = ventana["fuera_de_la_ventana"]
        print(
            f"[v2-deteccion] dentro de la ventana: {d['valor']} ({d['denominador']}); "
            f"fuera: {f['valor']} ({f['denominador']})"
        )

    tabla = deteccion.tabla_por_proteina(
        lfc_tejido, lfc_plasma, conc, detectada, umbral_log10=umbral
    )
    write_tsv(tabla, salida / "ventana_deteccion_por_proteina.tsv")

    auditoria = {
        "version": "v2",
        "estado": "calculado",
        "fuente_concentraciones": {
            "nombre": "human_protein_atlas",
            "estado": estado_hpa,
            "cita": ficha.get("hpa", {}).get("cita"),
            "licencia": ficha.get("hpa", {}).get("licencia"),
            "sha256_de_la_tabla_limpia": ficha.get("hpa", {}).get("sha256_de_la_tabla_limpia"),
            **info_hpa,
        },
        "n_proteinas_cuantificadas_en_suero": len(medidas),
        "curva_de_deteccion": curva,
        "concordancia": ventana,
        "umbral_deteccion": prereg.get("c8_ventana_deteccion.umbral_deteccion"),
        "lenguaje_obligatorio": prereg.get("c8_ventana_deteccion.lenguaje_obligatorio"),
        "n_min_curva": N_MIN_CURVA,
        "n_min_ventana": N_MIN_VENTANA,
        "limitaciones": [
            "La concentración del HPA es de sujetos sanos de referencia, no de los "
            "pacientes de PXD055218. Es una aproximación al orden de magnitud.",
            "Solo una fracción pequeña del proteoma tiene concentración publicada, "
            "así que la ventana se define sobre un subconjunto sesgado hacia las "
            "proteínas abundantes y ya estudiadas.",
            "El tejido es CPTAC y el suero es otra cohorte. No hay ningún paciente "
            "con las dos medidas, así que la comparación es entre cohortes.",
            "Ningún resultado de este fichero sirve para decidir sobre un paciente.",
        ],
    }
    (salida / "deteccion_audit_v2.json").write_text(
        json.dumps(auditoria, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    print(f"[v2-deteccion] escrito en {salida}")
