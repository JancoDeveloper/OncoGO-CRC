"""Meta-análisis de efectos aleatorios por proteína y R3 (C10, C11).

Qué desbloquea esto. R3 es la mediana de 1 − I², y con dos cohortes I² es tan
inestable que el preregistro exige tres. v1 tenía dos, así que R3 figuraba
``no_calculable``. Aquí entra una tercera cohorte de tejido, PXD019504, y R3 pasa
a ser calculable.

Antes de usarla hay que probar que es una cohorte distinta. El repositorio ya
documentaba que PXD019504 y PXD031705 comparten ficheros, así que la sospecha
estaba justificada. La comprobación tiene dos niveles y vive en
``scripts/download/independencia_cohortes.py``: las sumas de verificación de
PRIDE y, lo que de verdad decide, la correlación de cada muestra de una cohorte
contra cada muestra de la otra. Este stage se niega a incluir una cohorte que no
haya pasado las dos.

Las matrices no se fusionan. Cada cohorte produce su propio contraste pareado y
lo único que se combina son los efectos y sus errores estándar, proteína por
proteína, con REML.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.moderated import moderate_table, paired_observed_stats
from oncogo_crc.utils.reproducibility_overlap import (
    load_gene_uniprot_from_mapping,
    uniprot_base,
)
from oncogo_crc.utils.ids import effective_protein_mapping
from oncogo_crc.stages.certified_inference import MIN_COMPLETE
from oncogo_crc.stages.v2_certified import _pride_pares
from oncogo_crc.v2.meta import meta_efectos_aleatorios
from oncogo_crc.v2.prereg import cargar

# En PXD019504 el tumor es T; en PXD034575 es M. No son intercambiables y
# confundirlos invierte el signo de toda la cohorte.
PATRON_019504 = ("N", "T")


def _acceso(valor: str) -> str:
    """Accession de un identificador de MaxQuant.

    PXD019504 trae accesiones sueltas separadas por punto y coma. Otros depósitos
    usan el formato de FASTA ``sp|P12345|NOMBRE_HUMAN``, donde la accesión es el
    campo del medio y no el último: quedarse con el último daría el nombre de la
    proteína y no cruzaría con nada.
    """
    partes = str(valor).split("|")
    return partes[1] if len(partes) >= 3 else partes[0]


def _matriz_maxquant(path) -> tuple[pd.DataFrame, dict[str, dict[str, str]]]:
    """Log2 de las columnas LFQ de un proteinGroups, con sus pares por paciente."""
    columnas = pd.read_csv(path, sep="\t", nrows=0, low_memory=False).columns
    lfq = [c for c in columnas if str(c).startswith("LFQ intensity ")]
    clave = "Majority protein IDs" if "Majority protein IDs" in columnas else "Protein IDs"
    tabla = pd.read_csv(path, sep="\t", usecols=[clave, *lfq], low_memory=False)
    banderas = [
        b
        for b in ("Reverse", "Potential contaminant", "Only identified by site")
        if b in columnas
    ]
    if banderas:
        marcas = pd.read_csv(path, sep="\t", usecols=banderas, low_memory=False)
        tabla = tabla.loc[marcas.isna().all(axis=1).to_numpy()]
    ids = tabla[clave].astype(str).str.split(";").str[0].map(_acceso).str.upper()
    valores = tabla[lfq].apply(pd.to_numeric, errors="coerce")
    # El cero de MaxQuant es ausencia de identificación, no abundancia cero.
    with np.errstate(divide="ignore", invalid="ignore"):
        log2 = np.where(valores > 0, np.log2(valores), np.nan)
    matriz = pd.DataFrame(log2, index=ids, columns=lfq)
    matriz = matriz.loc[~matriz.index.duplicated(keep="first")]

    pacientes: dict[str, dict[str, str]] = {}
    for col in lfq:
        nombre = str(col).replace("LFQ intensity ", "").strip()
        sufijo = nombre[-1].upper()
        if sufijo in PATRON_019504:
            pacientes.setdefault(nombre[:-1], {})[sufijo] = col
    return matriz, pacientes


def _contraste_019504(path, *, lfc: float, alpha: float) -> tuple[pd.DataFrame, int]:
    matriz, pacientes = _matriz_maxquant(path)
    ids = [p for p, v in pacientes.items() if {"T", "N"} <= set(v)]
    if not ids:
        raise RuntimeError(
            "PXD019504: no hay ningún par T/N. En esta cohorte el tumor es T, no M."
        )
    diffs = (
        matriz[[pacientes[p]["T"] for p in ids]].to_numpy()
        - matriz[[pacientes[p]["N"] for p in ids]].to_numpy()
    )
    coef, s2, df, stdev, n_obs = paired_observed_stats(diffs, min_pairs=MIN_COMPLETE)
    tabla, _d0, _s02 = moderate_table(
        protein_ids=matriz.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method="moderated_t_treat_paired_lfq_pxd019504",
    )
    tabla["comparison_id"] = "pxd019504_tumor_vs_normal"
    return tabla, len(ids)


def _por_uniprot(tabla: pd.DataFrame, *, a_uniprot=None) -> pd.DataFrame:
    """Efecto y error estándar por UniProt, una fila por proteína."""
    t = tabla[["protein_id", "log2fc", "se"]].copy()
    if a_uniprot is None:
        t["uniprot"] = t["protein_id"].map(lambda v: uniprot_base(str(v)) or "")
    else:
        t["uniprot"] = t["protein_id"].astype(str).str.upper().map(a_uniprot).fillna("")
    t = t.loc[t["uniprot"].str.len() > 0]
    t = t.dropna(subset=["log2fc", "se"])
    t = t.loc[np.isfinite(t["se"]) & (t["se"] > 0)]
    # Si varias filas dan al mismo UniProt gana la de menor error estándar.
    return t.sort_values("se", kind="stable").drop_duplicates("uniprot")


def _cohortes_permitidas(settings: Settings) -> tuple[list[str], dict]:
    """Solo entran las cohortes que pasaron los dos niveles de independencia."""
    path = settings.resolve("external") / "independencia" / "estado.json"
    if not path.exists():
        return [], {
            "estado": "no_disponible",
            "motivo": (
                "no hay auditoría de independencia. Ejecute "
                "scripts/download/independencia_cohortes.py"
            ),
        }
    ficha = json.loads(path.read_text(encoding="utf-8"))
    nivel_a = {d["par"]: d for d in ficha.get("nivel_a_hashes_de_pride", [])}
    nivel_b = ficha.get("nivel_b_correlacion_entre_muestras", {})
    par = "PXD019504_vs_PXD034575"
    a_ok = nivel_a.get(par, {}).get("n_hashes_compartidos", 1) == 0
    b = nivel_b.get(par, {})
    b_ok = b.get("estado") == "evaluado" and b.get("n_pares_sospechosos") == 0
    permitidas = ["PXD019504"] if (a_ok and b_ok) else []
    resumen = {
        "estado": "evaluado",
        "nivel_a_sin_ficheros_compartidos": bool(a_ok),
        "nivel_b_sin_muestras_repetidas": bool(b_ok),
        "r_maximo_entre_muestras": b.get("r_maximo"),
        "r_mediano_entre_muestras": b.get("r_mediano"),
        "umbral_de_sospecha": b.get("umbral"),
        "pxd031705_excluida": (
            "Comparte ficheros con PXD019504, así que no es una cohorte distinta. "
            "Entra una de las dos, nunca las dos."
        ),
        "veredicto": (
            "PXD019504 entra como tercera cohorte de tejido"
            if permitidas
            else "PXD019504 no entra: no supera la comprobación de independencia"
        ),
    }
    return permitidas, resumen


def _diagnostico_tecnologia(
    con_tres: pd.DataFrame, tabla_lfq: pd.DataFrame, cohortes_lfq: list[str]
) -> dict:
    """¿La heterogeneidad viene de la tecnología o de las poblaciones?

    Se comparan las mismas proteínas con y sin CPTAC. Si I² cae al quitar la
    cohorte de TMT, lo que R3 está midiendo es sobre todo la diferencia entre
    marcaje isobárico y LFQ, no una diferencia entre pacientes.
    """
    if not len(con_tres) or not len(tabla_lfq):
        return {"estado": "no_calculable", "motivo": "no hay proteínas en los dos conjuntos"}
    juntas = con_tres.merge(tabla_lfq, on="uniprot", how="inner")
    if not len(juntas):
        return {"estado": "no_calculable", "motivo": "ninguna proteína en común"}
    i2_tres = juntas["i2"].dropna()
    i2_lfq = juntas["i2_solo_lfq"].dropna()
    if not len(i2_tres) or not len(i2_lfq):
        return {"estado": "no_calculable", "motivo": "I² no estimable en uno de los dos"}
    caida = float(i2_tres.median() - i2_lfq.median())
    return {
        "estado": "calculado",
        "n_proteinas_comparadas": int(len(juntas)),
        "cohortes_lfq": cohortes_lfq,
        "i2_mediano_con_las_tres": float(i2_tres.median()),
        "i2_mediano_solo_lfq": float(i2_lfq.median()),
        "caida_al_quitar_cptac": caida,
        "advertencia": (
            "Al quitar la cohorte de TMT la heterogeneidad baja mucho, así que una "
            "parte grande de I² es diferencia de tecnología y no de poblaciones. R3 "
            "no debe leerse como «los estudios discrepan biológicamente»."
            if caida > 0.2
            else "La heterogeneidad apenas cambia al quitar la cohorte de TMT, así que "
            "no se explica por la diferencia de tecnología."
        ),
        "nota": (
            "Con dos cohortes I² es inestable, así que este número es un diagnóstico "
            "y no un valor del certificado. R3 sigue siendo el de las tres."
        ),
    }


def run(settings: Settings) -> None:
    prereg = cargar(settings.root)
    alpha = float(prereg.get("inferencia.alpha"))
    lfc = float(prereg.get("inferencia.treat_lfc"))
    min_cohortes = int(prereg.get("c10_meta.min_cohortes_para_reportar_r3"))
    salida = settings.resolve("results") / "v2" / "meta"
    salida.mkdir(parents=True, exist_ok=True)

    permitidas, independencia = _cohortes_permitidas(settings)
    print(f"[v2-meta] independencia: {independencia.get('veredicto', independencia.get('motivo'))}")

    mapping = effective_protein_mapping(settings)
    gene_to_uni = load_gene_uniprot_from_mapping(mapping, exclude_ambiguous=True)

    cert = settings.resolve("results") / "v2" / "certified"
    cptac = pd.read_csv(cert / "cptac_unshared_treat.tsv", sep="\t")
    certificadas = set(
        cptac.loc[cptac["significant"], "protein_id"].astype(str).str.upper()
    )
    efectos: dict[str, pd.DataFrame] = {
        "cptac_pdc000116": _por_uniprot(cptac, a_uniprot=gene_to_uni)
    }

    tabla_pride, _diffs, n_pares_pride = _pride_pares(settings, lfc=lfc, alpha=alpha)
    efectos["pride_pxd034575"] = _por_uniprot(tabla_pride)

    detalle_cohortes = {
        "cptac_pdc000116": {"tecnologia": "TMT, razón logarítmica", "n_pares": 94},
        "pride_pxd034575": {"tecnologia": "LFQ MaxQuant", "n_pares": n_pares_pride},
    }

    if "PXD019504" in permitidas:
        path = settings.resolve("external") / "pride" / "PXD019504" / "proteinGroups.txt"
        if path.exists():
            tabla_019504, n_pares = _contraste_019504(path, lfc=lfc, alpha=alpha)
            write_tsv(tabla_019504, salida / "pxd019504_tumor_vs_normal_treat.tsv")
            efectos["pride_pxd019504"] = _por_uniprot(tabla_019504)
            detalle_cohortes["pride_pxd019504"] = {
                "tecnologia": "LFQ MaxQuant",
                "n_pares": n_pares,
                "etiqueta_tumor": "T",
            }
            print(
                f"[v2-meta] PXD019504: {n_pares} pares T/N, "
                f"{int(tabla_019504['significant'].sum())} significativas."
            )
        else:
            print(f"[v2-meta] falta {path}")

    for nombre, t in efectos.items():
        print(f"[v2-meta] {nombre}: {len(t)} proteínas con efecto y error estándar.")

    # El meta-análisis se limita a las proteínas certificadas en descubrimiento.
    # Hacerlo sobre todo el proteoma mezclaría la pregunta de reproducibilidad con
    # un cribado nuevo.
    uni_certificadas = {
        gene_to_uni[g] for g in certificadas if g in gene_to_uni
    }
    indexadas = {k: v.set_index("uniprot") for k, v in efectos.items()}
    filas = []
    for uni in sorted(uni_certificadas):
        y, se, presentes = [], [], []
        for nombre, t in indexadas.items():
            if uni in t.index:
                y.append(float(t.loc[uni, "log2fc"]))
                se.append(float(t.loc[uni, "se"]))
                presentes.append(nombre)
        if len(y) < 2:
            continue
        res = meta_efectos_aleatorios(np.array(y), np.array(se), metodo="reml")
        if res is None:
            continue
        filas.append(
            {
                "uniprot": uni,
                "n_cohortes": len(y),
                "cohortes": ";".join(presentes),
                "efecto_combinado": res.efecto,
                "se_combinado": res.se,
                "ci_low": res.ci_low,
                "ci_high": res.ci_high,
                "tau2": res.tau2,
                "q_cochran": res.q,
                "i2": res.i2,
                "uno_menos_i2": None if res.i2 is None else 1.0 - res.i2,
                "metodo": res.metodo,
            }
        )
    tabla_meta = pd.DataFrame(filas)
    write_tsv(tabla_meta, salida / "meta_por_proteina.tsv")

    # Diagnóstico de tecnología. Si la heterogeneidad baja mucho al quitar CPTAC,
    # lo que mide I² es la diferencia entre TMT y LFQ, no entre poblaciones. Se
    # calcula sobre las mismas proteínas para que la comparación sea limpia.
    solo_lfq = [k for k in indexadas if k != "cptac_pdc000116"]
    filas_lfq = []
    for uni in sorted(uni_certificadas):
        y, se = [], []
        for nombre in solo_lfq:
            t = indexadas[nombre]
            if uni in t.index:
                y.append(float(t.loc[uni, "log2fc"]))
                se.append(float(t.loc[uni, "se"]))
        if len(y) < 2:
            continue
        res = meta_efectos_aleatorios(np.array(y), np.array(se), metodo="reml")
        if res is not None and res.i2 is not None:
            filas_lfq.append({"uniprot": uni, "i2_solo_lfq": res.i2, "tau2_solo_lfq": res.tau2})
    tabla_lfq = pd.DataFrame(filas_lfq)
    if len(tabla_lfq):
        write_tsv(tabla_lfq, salida / "i2_solo_cohortes_lfq.tsv")

    con_tres = tabla_meta.loc[tabla_meta["n_cohortes"] >= min_cohortes]
    i2_ok = con_tres["uno_menos_i2"].dropna()
    if len(i2_ok) and len(indexadas) >= min_cohortes:
        r3 = float(i2_ok.median())
        estado_r3 = "calculado"
        motivo_r3 = ""
    else:
        r3 = None
        estado_r3 = "no_calculable"
        motivo_r3 = (
            f"hay {len(indexadas)} cohortes de tejido y "
            f"{len(con_tres)} proteínas medidas en al menos {min_cohortes}"
        )
    print(
        f"[v2-meta] R3: {estado_r3} {'' if r3 is None else round(r3, 4)} "
        f"({len(con_tres)} proteínas en ≥{min_cohortes} cohortes de {len(tabla_meta)} con ≥2)"
    )

    auditoria = {
        "version": "v2",
        "estado": estado_r3,
        "n_cohortes": len(indexadas),
        "cohortes": detalle_cohortes,
        "independencia_de_cohortes": independencia,
        "min_cohortes_para_r3": min_cohortes,
        "n_proteinas_certificadas": len(certificadas),
        "n_proteinas_certificadas_con_uniprot": len(uni_certificadas),
        "n_proteinas_con_al_menos_dos_cohortes": int(len(tabla_meta)),
        "n_proteinas_evaluables": int(len(con_tres)),
        "mediana_1_menos_i2": r3,
        "motivo": motivo_r3,
        "mediana_i2": float(con_tres["i2"].dropna().median()) if len(con_tres) else None,
        "mediana_tau2": float(con_tres["tau2"].dropna().median()) if len(con_tres) else None,
        "n_con_efecto_combinado_distinto_de_cero": (
            int(((con_tres["ci_low"] > 0) | (con_tres["ci_high"] < 0)).sum())
            if len(con_tres)
            else 0
        ),
        "metodo": "REML, Higgins y Thompson 2002 para I². Sustituye a metafor::rma (enmienda E2).",
        "diagnostico_de_tecnologia": _diagnostico_tecnologia(con_tres, tabla_lfq, solo_lfq),
        "limitaciones": [
            "CPTAC mide razones logarítmicas con TMT y las otras dos miden log2 de "
            "intensidad LFQ. El TMT comprime las razones, así que la magnitud no es "
            "estrictamente comparable aunque el signo y el orden sí lo sean. Un I² "
            "alto puede reflejar esa diferencia de tecnología y no de biología.",
            "El meta-análisis se limita a las proteínas certificadas en CPTAC, así "
            "que el conjunto está seleccionado por el descubrimiento. No es un "
            "cribado nuevo.",
            "Solo las proteínas medidas en al menos tres cohortes entran en R3. Las "
            "que solo están en dos se reportan en la tabla pero no en el resumen.",
            "Ningún resultado de este fichero sirve para decidir sobre un paciente.",
        ],
    }
    (salida / "meta_audit_v2.json").write_text(
        json.dumps(auditoria, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    print(f"[v2-meta] escrito en {salida}")
