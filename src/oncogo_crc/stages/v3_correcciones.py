"""Correcciones K1, K4, K5, K6, K7 y K10.

Escribe solo bajo results/v3/. No toca v1, v2 ni el cierre. Los umbrales salen
de config/preregistro_v3.yaml.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal

from oncogo_crc.settings import Settings
from oncogo_crc.stages.certified_inference import _peptide_fraction, _unique_peptide_counts
from oncogo_crc.tesis.cierre import a_escala_log2
from oncogo_crc.utils.dia_nn_plasma import build_plasma_metadata_for_pxd055218_runs
from oncogo_crc.utils.ids import effective_protein_mapping
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.moderated import moderate_table, unpaired_observed_stats
from oncogo_crc.utils.reproducibility_overlap import load_gene_uniprot_from_mapping, uniprot_base
from oncogo_crc.v2.metricas import Fraccion
from oncogo_crc.v3.auc import bootstrap_auc, bootstrap_diferencia, delong_ic, nulo_del_maximo_plegado
from oncogo_crc.v3.concordancia import (
    alinear_por_clave,
    concordancia_de_lin,
    concordancia_de_signo,
    spearman,
)
from oncogo_crc.v3.etiquetas import (
    CONTRASTE_PRINCIPAL,
    contraste_no_pareado,
    contraste_pareado,
    filas_desde_columnas,
    matriz_log2_lfq,
)
from oncogo_crc.v3.faltantes import fisher_presencia
from oncogo_crc.v3.potencia import (
    alfa_efectivo,
    clasificar,
    efecto_minimo_detectable,
    p_tost,
)
from oncogo_crc.v3.prereg import PreregistroV3
from oncogo_crc.v3.rutas import exigir_salida_v3
from statsmodels.stats.multitest import multipletests

CHANGELOG = (
    "metrica",
    "valor_anterior",
    "valor_v3",
    "numerador",
    "denominador",
    "motivo",
    "correccion",
)


def _json(path: Path, payload: dict) -> None:
    exigir_salida_v3(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=_jsonable), encoding="utf-8")


def _jsonable(valor):
    if isinstance(valor, (np.floating, float)):
        numero = float(valor)
        return numero if np.isfinite(numero) else None
    if isinstance(valor, (np.integer,)):
        return int(valor)
    if isinstance(valor, (np.bool_,)):
        return bool(valor)
    raise TypeError(type(valor))


def _escribir(df: pd.DataFrame, path: Path) -> None:
    exigir_salida_v3(path)
    write_tsv(df, path)


def _clave_cptac(genes: pd.Series, gen_a_uni: dict[str, str]) -> pd.Series:
    return genes.astype(str).str.upper().map(lambda g: gen_a_uni.get(g) or np.nan)


def _clave_replica(identificadores: pd.Series) -> pd.Series:
    return identificadores.astype(str).map(lambda v: uniprot_base(v.split(";")[0]) or np.nan)


def _resumen_signo(fraccion: Fraccion, extra: dict) -> dict:
    payload = fraccion.to_dict()
    payload.update(extra)
    return payload


def _k1(settings: Settings, prereg: PreregistroV3, salida: Path, filas_log: list[dict]) -> dict:
    alpha = float(prereg.get("inferencia.alpha"))
    lfc = float(prereg.get("inferencia.treat_lfc"))
    min_pares = int(prereg.get("inferencia.min_pares_observados"))
    ds = settings.datasets["datasets"]["pride_crc_meta"]
    matriz_path = settings.resolve(ds["matrix_path"])
    if not matriz_path.is_file():
        return {"estado": "no_disponible", "motivo": f"no está {matriz_path}"}
    cruda = pd.read_csv(matriz_path, sep="\t", low_memory=False)
    meta = filas_desde_columnas([c for c in cruda.columns if c != "protein_id"])
    esperado = (
        prereg.get("k1_etiquetas.contraste_principal.grupo_a"),
        prereg.get("k1_etiquetas.contraste_principal.grupo_b"),
    )
    if CONTRASTE_PRINCIPAL != esperado:
        raise RuntimeError("el contraste principal del código no es el del preregistro")
    for destino in (
        settings.root / "metadata/sample_metadata/pride_pxd034575_v3.tsv",
        settings.resolve("metadata/sample_metadata/pride_pxd034575_v3.tsv"),
    ):
        destino.parent.mkdir(parents=True, exist_ok=True)
        meta.to_csv(destino, sep="\t", index=False)
    log2mat = matriz_log2_lfq(cruda)
    principal, n_pares = contraste_pareado(
        log2mat,
        meta,
        grupo_a="primary",
        grupo_b="mucosa",
        nombre="primario_frente_a_mucosa",
        lfc=lfc,
        alpha=alpha,
        min_pares=min_pares,
    )
    metastasis, n_met = contraste_pareado(
        log2mat,
        meta,
        grupo_a="liver_metastasis",
        grupo_b="mucosa",
        nombre="metastasis_hepatica_frente_a_mucosa",
        lfc=lfc,
        alpha=alpha,
        min_pares=min_pares,
    )
    frente_primario, n_mp = contraste_pareado(
        log2mat,
        meta,
        grupo_a="liver_metastasis",
        grupo_b="primary",
        nombre="metastasis_hepatica_frente_a_primario",
        lfc=lfc,
        alpha=alpha,
        min_pares=min_pares,
    )
    no_pareado = contraste_no_pareado(
        log2mat,
        meta,
        grupo_a="primary",
        grupo_b="mucosa",
        nombre="primario_frente_a_mucosa_no_pareado",
        lfc=lfc,
        alpha=alpha,
        min_n=min_pares,
    )
    cptac_path = settings.resolve("results/certified/cptac_unshared_treat.tsv")
    cptac = pd.read_csv(cptac_path, sep="\t")
    certificadas = cptac.loc[cptac["significant"].astype(bool)].copy()
    n_cert = int(certificadas["protein_id"].nunique())
    gen_a_uni = load_gene_uniprot_from_mapping(effective_protein_mapping(settings))
    resumenes = {}
    for nombre, tabla in (
        ("primario_frente_a_mucosa", principal),
        ("metastasis_hepatica_frente_a_mucosa", metastasis),
        ("primario_frente_a_mucosa_no_pareado", no_pareado),
    ):
        cruzado = alinear_por_clave(
            certificadas,
            tabla,
            clave_descubrimiento=_clave_cptac(certificadas["protein_id"], gen_a_uni),
            clave_replica=_clave_replica(tabla["protein_id"]),
        )
        fraccion = concordancia_de_signo(
            cruzado["log2fc_cptac"].to_numpy(),
            cruzado["log2fc_replica"].to_numpy(),
            nombre=f"R1_{nombre}",
            universo=n_cert,
            universo_nombre="proteinas_certificadas_cptac",
        )
        resumenes[nombre] = _resumen_signo(
            fraccion,
            {
                "spearman": spearman(cruzado["log2fc_cptac"].to_numpy(), cruzado["log2fc_replica"].to_numpy()),
                "lin": concordancia_de_lin(cruzado["log2fc_cptac"].to_numpy(), cruzado["log2fc_replica"].to_numpy()),
                "n_pares_diseno": int(n_pares if "no_pareado" not in nombre else meta["participant_id"].nunique()),
            },
        )
        _escribir(tabla, salida / f"k1_{nombre}.tsv")
    _escribir(frente_primario, salida / "k1_metastasis_hepatica_frente_a_primario.tsv")
    peptidos = _unique_peptide_counts(settings) or {}
    _frac, n_pep, n_den = _peptide_fraction(cptac, gen_a_uni, peptidos, alpha=alpha)
    r2 = Fraccion(
        nombre="R2_detectabilidad_peptidos",
        numerador=int(n_pep),
        denominador=int(n_den),
        unidad_denominador="proteinas_unicas",
        nota="Péptido único en proteinGroups. No depende de qué columna sea el tumor.",
    )
    cobertura = Fraccion(
        nombre="R1_cobertura_primario_frente_a_mucosa",
        numerador=int(resumenes["primario_frente_a_mucosa"]["denominador"]),
        denominador=n_cert,
        unidad_denominador="proteinas_unicas",
        nota="Certificadas con log2FC evaluable en el contraste principal.",
    )
    r1 = resumenes["primario_frente_a_mucosa"]
    umbral = float(prereg.get("hipotesis.h1.umbral_concordancia_proteina"))
    if r1["ci_high"] is None:
        clausula = "no_evaluable"
    elif float(r1["ci_high"]) < umbral:
        clausula = "sostiene"
    else:
        clausula = "rechaza"
    veredicto_h1 = "rechaza_h1" if clausula == "rechaza" else (
        "parcial_solo_proteina" if clausula == "sostiene" else "no_evaluable"
    )
    filas_log.append({
        "metrica": "R1_concordancia_de_signo",
        "valor_anterior": "104/119",
        "valor_v3": f"{r1['numerador']}/{r1['denominador']}",
        "numerador": r1["numerador"],
        "denominador": r1["denominador"],
        "motivo": "v1 y v2 usaron la metástasis hepática (M) como si fuera el tumor. El contraste principal de v3 es el primario (P) frente a la mucosa (N), pareado.",
        "correccion": "K1",
    })
    filas_log.append({
        "metrica": "R2_detectabilidad_peptidos",
        "valor_anterior": "226/688",
        "valor_v3": f"{r2.numerador}/{r2.denominador}",
        "numerador": r2.numerador,
        "denominador": r2.denominador,
        "motivo": "Cuenta péptidos únicos del proteinGroups. La etiqueta de la columna no entra en la cuenta.",
        "correccion": "K1",
    })
    filas_log.append({
        "metrica": "R3_homogeneidad_1_menos_I2",
        "valor_anterior": "0.094",
        "valor_v3": "no_evaluable",
        "numerador": "",
        "denominador": "",
        "motivo": "v2 metió a PXD019504 como tercera cohorte. K2 todavía no descarta que comparta personas con PXD034575 o PXD031705. Sin eso, 1 − I² de tres cohortes no se recalcula.",
        "correccion": "K1",
    })
    comparacion = pd.DataFrame([
        {"contraste": k, **{kk: vv for kk, vv in v.items() if not isinstance(vv, dict)},
         "spearman": (v.get("spearman") or {}).get("rho"),
         "lin_ccc": (v.get("lin") or {}).get("ccc")}
        for k, v in resumenes.items()
    ])
    _escribir(comparacion, salida / "k1_comparacion_contrastes.tsv")
    payload = {
        "estado": "calculado",
        "n_pares_primario_mucosa": n_pares,
        "n_pares_metastasis_mucosa": n_met,
        "n_pares_metastasis_primario": n_mp,
        "n_certificadas": n_cert,
        "contrastes": resumenes,
        "r2_peptidos": r2.to_dict(),
        "cobertura_principal": cobertura.to_dict(),
        "h1_clausula_proteina": clausula,
        "h1_clausula_via": "no_evaluable",
        "h1_veredicto": veredicto_h1,
        "h1_motivo_via": "La cláusula de vía no se ha vuelto a calcular con P frente a N.",
        "r3": {"estado": "no_evaluable", "motivo": filas_log[-1]["motivo"]},
    }
    _json(salida / "k1_etiquetas.json", payload)
    return payload


def _grupos_suero(columnas: list[str]) -> dict[str, list[str]]:
    meta = build_plasma_metadata_for_pxd055218_runs(columnas)
    grupos: dict[str, list[str]] = {}
    for condicion, bloque in meta.groupby("condition"):
        grupos[str(condicion)] = bloque["sample_id"].astype(str).tolist()
    return grupos


def _contraste_suero(mat: pd.DataFrame, cols_a: list[str], cols_b: list[str], *, lfc: float, alpha: float, nombre: str):
    a = mat[cols_a].to_numpy(dtype=float)
    b = mat[cols_b].to_numpy(dtype=float)
    coef, s2, df, stdev, n_obs = unpaired_observed_stats(a, b, min_n=3)
    tabla, _d0, _s02 = moderate_table(
        protein_ids=mat.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method=f"moderated_t_treat_unpaired_{nombre}",
    )
    return tabla


def _k4(mat: pd.DataFrame, grupos: dict[str, list[str]], prereg: PreregistroV3, salida: Path) -> dict:
    alpha = float(prereg.get("inferencia.alpha"))
    lfc = float(prereg.get("inferencia.treat_lfc"))
    rng_alfa = np.random.default_rng(prereg.semilla("k4_alfa_efectivo"))
    rng_mde = np.random.default_rng(prereg.semilla("k4_efecto_minimo"))
    pares = {
        "crc_frente_a_colonoscopia_negativa": ("crc", "colonoscopy_negative"),
        "crc_frente_a_sanos": ("crc", "healthy"),
        "colonoscopia_negativa_frente_a_sanos": ("colonoscopy_negative", "healthy"),
    }
    conteos = {}
    for nombre, (a, b) in pares.items():
        tabla = _contraste_suero(mat, grupos[a], grupos[b], lfc=lfc, alpha=alpha, nombre=nombre)
        p = p_tost(tabla["log2fc"].to_numpy(), tabla["se"].to_numpy(), tabla["df_total"].to_numpy(), lfc)
        tabla["p_tost"] = p
        padj = np.full(len(tabla), np.nan)
        ok = np.isfinite(p)
        if ok.any():
            _, q, _, _ = multipletests(p[ok], alpha=alpha, method="fdr_bh")
            padj[ok] = q
        tabla["padj_tost"] = padj
        se = tabla["se"].to_numpy(dtype=float)
        df = tabla["df_total"].to_numpy(dtype=float)
        alfa = alfa_efectivo(
            se, df, alpha=alpha, lfc=lfc,
            n_sim=int(prereg.get("k4_suero.n_sim_alfa")), rng=rng_alfa,
        )
        mde = efecto_minimo_detectable(
            se, df, alpha_eff=alfa, lfc=lfc,
            potencia=float(prereg.get("k4_suero.potencia_objetivo")),
            delta_max=float(prereg.get("k4_suero.delta_maximo_busqueda")),
            n_sim=int(prereg.get("k4_suero.n_sim_efecto")),
            n_delta=20,
            rng=rng_mde,
        )
        tabla["alfa_efectivo"] = alfa
        tabla["efecto_minimo_detectable"] = mde
        tabla["clase"] = clasificar(
            tabla["significant"].fillna(False).to_numpy(dtype=bool),
            (tabla["padj_tost"] < alpha).fillna(False).to_numpy(dtype=bool),
        )
        _escribir(tabla, salida / f"k4_{nombre}.tsv")
        conteos[nombre] = {
            "alfa_efectivo": alfa,
            "n": int(len(tabla)),
            **tabla["clase"].value_counts().to_dict(),
        }
    _json(salida / "k4_potencia.json", {"estado": "calculado", "contrastes": conteos})
    return conteos


def _valores_grupo(mat: pd.DataFrame, columnas: list[str], indice) -> np.ndarray:
    return mat.loc[indice, columnas].to_numpy(dtype=float)


def _k5(mat: pd.DataFrame, grupos: dict[str, list[str]], prereg: PreregistroV3, salida: Path) -> dict:
    min_n = int(prereg.get("k5_auc.minimo_observados_por_grupo"))
    n_boot = int(prereg.get("k5_auc.bootstrap.n"))
    rng_boot = np.random.default_rng(prereg.semilla("k5_bootstrap_auc"))
    rng_dif = np.random.default_rng(prereg.semilla("k5_bootstrap_diferencia"))
    rng_perm = np.random.default_rng(prereg.semilla("k5_permutacion_etiquetas"))
    pares = {
        "crc_frente_a_colonoscopia_negativa": ("crc", "colonoscopy_negative", 1),
        "colonoscopia_negativa_frente_a_sanos": ("colonoscopy_negative", "healthy", 1),
    }
    filas_auc = []
    for nombre, (ga, gb, _signo) in pares.items():
        for proteina in mat.index:
            pos = _valores_grupo(mat, grupos[ga], proteina)
            neg = _valores_grupo(mat, grupos[gb], proteina)
            if np.isfinite(pos).sum() < min_n or np.isfinite(neg).sum() < min_n:
                continue
            auc, lo, hi = delong_ic(pos, neg)
            blo, bhi = bootstrap_auc(pos, neg, n_boot=n_boot, rng=rng_boot)
            filas_auc.append({
                "contraste": nombre,
                "protein_id": proteina,
                "auc_con_signo": auc,
                "delong_bajo": lo,
                "delong_alto": hi,
                "bootstrap_bajo": blo,
                "bootstrap_alto": bhi,
                "n_positivo": int(np.isfinite(pos).sum()),
                "n_negativo": int(np.isfinite(neg).sum()),
            })
    aucs = pd.DataFrame(filas_auc)
    _escribir(aucs, salida / "k5_auc.tsv")
    cols = grupos["crc"] + grupos["colonoscopy_negative"]
    etiquetas = np.array([1] * len(grupos["crc"]) + [0] * len(grupos["colonoscopy_negative"]))
    bloque = mat[cols].to_numpy(dtype=float)
    print(f"[v3] nulo del máximo AUC: {int(prereg.get('k5_auc.nulo_del_maximo.n_permutaciones'))} permutaciones")
    nulo = nulo_del_maximo_plegado(
        bloque, etiquetas,
        n_perm=int(prereg.get("k5_auc.nulo_del_maximo.n_permutaciones")),
        min_n=min_n,
        rng=rng_perm,
    )
    diferencias = []
    for proteina in mat.index:
        diff = bootstrap_diferencia(
            _valores_grupo(mat, grupos["crc"], proteina),
            _valores_grupo(mat, grupos["colonoscopy_negative"], proteina),
            _valores_grupo(mat, grupos["healthy"], proteina),
            n_boot=n_boot,
            rng=rng_dif,
        )
        if not np.isfinite(diff["diferencia"]):
            continue
        diff["protein_id"] = proteina
        diferencias.append(diff)
    dif = pd.DataFrame(diferencias)
    _escribir(dif, salida / "k5_diferencia_auc.tsv")
    _json(salida / "k5_auc.json", {"estado": "calculado", "nulo_maximo_plegado": nulo, "n_auc": int(len(aucs))})
    return {"nulo": nulo, "n_auc": int(len(aucs))}


def _k6(mat: pd.DataFrame, grupos: dict[str, list[str]], prereg: PreregistroV3, salida: Path) -> dict:
    alpha = float(prereg.get("inferencia.alpha"))
    umbral = float(prereg.get("k6_faltantes.umbral_observacion"))
    columnas = grupos["crc"] + grupos["colonoscopy_negative"] + grupos["healthy"]
    etiquetas = (
        ["crc"] * len(grupos["crc"])
        + ["colonoscopy_negative"] * len(grupos["colonoscopy_negative"])
        + ["healthy"] * len(grupos["healthy"])
    )
    bloque = mat[columnas].to_numpy(dtype=float)
    fraccion = np.isfinite(bloque).mean(axis=1)
    mascara = fraccion < umbral
    sub = bloque[mascara]
    ids = mat.index.to_numpy()[mascara]
    if sub.size == 0:
        tabla = pd.DataFrame()
    else:
        presencia = fisher_presencia(
            sub, np.array(etiquetas), grupo_a="crc", grupo_b="colonoscopy_negative", alpha=alpha,
        )
        tabla = presencia.copy()
        tabla.insert(0, "protein_id", ids)
        tabla.insert(1, "fraccion_observada", fraccion[mascara])
        en_sano = np.array(etiquetas) == "healthy"
        tabla["presentes_healthy"] = np.isfinite(sub[:, en_sano]).sum(axis=1)
        tabla["ausentes_healthy"] = int(en_sano.sum()) - tabla["presentes_healthy"]
    _escribir(tabla, salida / "k6_presencia.tsv")
    return {"estado": "calculado", "n_por_debajo_del_umbral": int(mascara.sum()), "umbral": umbral}


def _puntuacion(mat: pd.DataFrame, genes: list[str], gen_de_fila: pd.Series) -> pd.Series:
    filas = gen_de_fila[gen_de_fila.isin(genes)].index
    if len(filas) == 0:
        return pd.Series(np.nan, index=mat.columns)
    return mat.loc[filas].mean(axis=0, skipna=True)


def _matriz_lineal_para_fraccion(mat: pd.DataFrame, escala: str) -> tuple[pd.DataFrame, str]:
    """La fracción de albúmina e inmunoglobulinas se suma en intensidad lineal.

    Sumar log2 comprime el rango y deja una fracción que no es de masa.
    """
    if escala == "log2_de_intensidad_lineal":
        valores = mat.to_numpy(dtype=float)
        lineal = np.where(np.isfinite(valores), np.exp2(valores), np.nan)
        return pd.DataFrame(lineal, index=mat.index, columns=mat.columns), "intensidad_lineal_recuperada_con_exp2"
    return mat, escala


def _k7(mat: pd.DataFrame, grupos: dict[str, list[str]], prereg: PreregistroV3, settings: Settings, salida: Path, escala: str) -> dict:
    gen_a_uni = load_gene_uniprot_from_mapping(effective_protein_mapping(settings))
    uni_a_gen = {}
    for gen, uni in gen_a_uni.items():
        uni_a_gen.setdefault(uni, gen)
    genes = []
    for pid in mat.index.astype(str):
        base = uniprot_base(pid.split(";")[0]) or ""
        genes.append(uni_a_gen.get(base, pid.upper()))
    gen_de_fila = pd.Series(genes, index=mat.index)
    por_muestra = pd.DataFrame({"sample_id": mat.columns})
    por_muestra["n_grupos_proteicos"] = np.isfinite(mat.to_numpy(dtype=float)).sum(axis=0)
    por_muestra["mediana_intensidad"] = np.nanmedian(mat.to_numpy(dtype=float), axis=0)
    alta = list(prereg.get("k7_control_tecnico.alta_abundancia"))
    presentes = gen_de_fila[gen_de_fila.isin(alta)].index
    lineal, escala_fraccion = _matriz_lineal_para_fraccion(mat, escala)
    if len(presentes):
        intensidad = lineal.loc[presentes].to_numpy(dtype=float)
        total = np.nansum(lineal.to_numpy(dtype=float), axis=0)
        alta_suma = np.nansum(intensidad, axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            por_muestra["fraccion_alta_abundancia"] = np.where(total > 0, alta_suma / total, np.nan)
    else:
        por_muestra["fraccion_alta_abundancia"] = np.nan
    indicadores = {
        "hemolisis": list(prereg.get("k7_control_tecnico.hemolisis")),
        "plaquetas": list(prereg.get("k7_control_tecnico.plaquetas")),
        "coagulacion": list(prereg.get("k7_control_tecnico.coagulacion_baja_en_suero")),
    }
    for nombre, lista in indicadores.items():
        por_muestra[nombre] = _puntuacion(mat, lista, gen_de_fila).to_numpy()
    etiqueta = {}
    for condicion, cols in grupos.items():
        for col in cols:
            etiqueta[col] = condicion
    por_muestra["grupo"] = por_muestra["sample_id"].map(etiqueta)
    _escribir(por_muestra, salida / "k7_por_muestra.tsv")
    alpha = float(prereg.get("k7_control_tecnico.alpha_exploratorio"))
    pruebas = []
    for columna in ("n_grupos_proteicos", "mediana_intensidad", "fraccion_alta_abundancia", "hemolisis", "plaquetas", "coagulacion"):
        muestras = [por_muestra.loc[por_muestra["grupo"] == g, columna].dropna().to_numpy() for g in ("crc", "colonoscopy_negative", "healthy")]
        muestras = [m for m in muestras if m.size >= 3]
        if len(muestras) < 2:
            p = float("nan")
        else:
            p = float(kruskal(*muestras).pvalue)
        fila = {
            "indicador": columna,
            "p_kruskal": p,
            "confusion_tecnica": bool(np.isfinite(p) and p < alpha),
        }
        for g in ("crc", "colonoscopy_negative", "healthy"):
            valores_g = por_muestra.loc[por_muestra["grupo"] == g, columna].dropna()
            fila[f"mediana_{g}"] = float(valores_g.median()) if len(valores_g) else float("nan")
        pruebas.append(fila)
    informe = []
    if settings.data_workspace is not None:
        plasma_dir = settings.resolve("data/raw/plasma")
        if plasma_dir.is_dir():
            informe = list(plasma_dir.glob("DIA-NN_report.tsv"))
    payload = {
        "estado": "calculado",
        "pruebas": pruebas,
        "precursores": "no_disponible" if not informe else "informe_encontrado",
        "motivo_precursores": "DIA-NN_report.tsv no está en el workspace: la matriz ancha ya sustituyó al informe de precursores." if not informe else str(informe[0]),
        "orden_de_corrida": "limitacion_declarada",
        "motivo_orden": "El orden de corrida no está en la matriz ancha. Recuperarlo exigiría los crudos, y no se reconstruye.",
        "escala_puntuaciones": "log2",
        "escala_fraccion_alta_abundancia": escala_fraccion,
        "n_filas_alta_abundancia": int(len(presentes)),
    }
    _json(salida / "k7_control_tecnico.json", payload)
    _escribir(pd.DataFrame(pruebas), salida / "k7_pruebas.tsv")
    return payload


def _fila_log(metrica: str, anterior: str, nuevo: str, numerador: str, denominador: str, motivo: str, correccion: str) -> dict:
    return {
        "metrica": metrica,
        "valor_anterior": anterior,
        "valor_v3": nuevo,
        "numerador": numerador,
        "denominador": denominador,
        "motivo": motivo,
        "correccion": correccion,
    }


def _changelog_suero(filas: list[dict], k4: dict, k5: dict, k6: dict, k7: dict) -> None:
    for nombre, bloque in k4.items():
        clases = {k: int(bloque.get(k, 0)) for k in ("diferente", "equivalente", "indeterminada", "conflicto")}
        texto = "; ".join(f"{k}={v}" for k, v in clases.items())
        filas.append(_fila_log(
            f"suero_{nombre}_clases",
            "0 proteinas en cancer frente a colonoscopia negativa, sin potencia ni TOST" if nombre == "crc_frente_a_colonoscopia_negativa" else "sin clasificacion TOST",
            texto,
            str(clases["diferente"]),
            str(bloque["n"]),
            "Clases excluyentes. La frase de ausencia solo vale para la clase equivalente. El numerador es el de diferente.",
            "K4",
        ))
    nulo = k5["nulo"]
    filas.append(_fila_log(
        "auc_nulo_del_maximo_plegado",
        "AUC plegado sin nulo de permutacion",
        f"maximo={nulo['maximo_observado']}; media_nula={nulo['media_nula']}; p95_nulo={nulo['percentil95_nulo']}; b={nulo['b']}",
        str(int(nulo["b"]) + 1),
        str(int(nulo["n_permutaciones"]) + 1),
        "p = (b+1)/(B+1) sobre el maximo del AUC plegado del catalogo cancer frente a colonoscopia negativa. El AUC de cada proteina se reporta con signo.",
        "K5",
    ))
    filas.append(_fila_log(
        "proteinas_con_observacion_bajo_0,70",
        "CHMP4A y otras con faltantes analizadas solo sobre los observados",
        f"n_bajo_umbral={k6['n_por_debajo_del_umbral']}",
        str(k6["n_por_debajo_del_umbral"]),
        "",
        "Fisher de presencia y ausencia, con BH, separado del contraste de abundancia. No se imputa.",
        "K6",
    ))
    confusas = [p["indicador"] for p in k7["pruebas"] if p["confusion_tecnica"]]
    filas.append(_fila_log(
        "control_tecnico_pxd055218",
        "orden de corrida no revisado; sin pantalla de hemolisis",
        f"confusion={','.join(confusas) if confusas else 'ninguna'}; precursores={k7['precursores']}; orden={k7['orden_de_corrida']}",
        "",
        "",
        "Kruskal-Wallis exploratorio a 0,05. La fraccion de alta abundancia se suma en intensidad lineal. Los precursores no estan en el disco.",
        "K7",
    ))
    filas.append(_fila_log(
        "unidad_de_las_fracciones",
        "fracciones de v2 sin unidad explicita de proteina unica",
        "proteinas_unicas",
        "",
        "",
        "Cada fraccion de indice_R_v3.tsv lleva numerador, denominador y unidad. No se cuentan filas de un cruce.",
        "K10",
    ))


def _hipotesis(salida: Path, k1: dict) -> None:
    filas = [
        {
            "hipotesis": "H1",
            "estado": k1.get("h1_veredicto", "no_evaluable"),
            "metrica": "concordancia de signo, primario frente a mucosa",
            "detalle": k1.get("h1_clausula_proteina", ""),
            "fuente": "results/v3/k1_etiquetas.json",
        },
        {"hipotesis": "H2", "estado": "no_evaluada", "metrica": "red", "detalle": "Sin nulo que conserve el grado.", "fuente": "results/v2/certification/certificado_v2.json"},
        {"hipotesis": "H3", "estado": "no_evaluable", "metrica": "", "detalle": "No hay H3 en el preregistro de v2 ni en el de v3.", "fuente": "config/preregistro_v3.yaml"},
        {"hipotesis": "H4", "estado": "no_evaluable", "metrica": "fraccion medible", "detalle": "v2 fijó el umbral 0,90 y no cerró un veredicto en el certificado.", "fuente": "results/v2/certification/certificado_v2.json"},
    ]
    for codigo, nombre in (
        ("H5", "especificidad de órgano"),
        ("H6", "origen celular"),
        ("H7", "visibilidad en sangre"),
        ("H8", "carga tumoral en suero"),
        ("H9", "dirección causal"),
        ("H10", "especificidad clínica en sangre"),
    ):
        filas.append({
            "hipotesis": codigo,
            "estado": "no_evaluable",
            "metrica": nombre,
            "detalle": "Análisis preregistrado. Todavía no está calculado. No se publica un número aproximado.",
            "fuente": "config/preregistro_v3.yaml",
        })
    _escribir(pd.DataFrame(filas), salida / "hipotesis_v3.tsv")


def run(settings: Settings) -> None:
    if CONTRASTE_PRINCIPAL != ("primary", "mucosa"):
        raise RuntimeError("el contraste principal tiene que ser primario frente a mucosa")
    prereg = PreregistroV3(settings.root)
    salida = settings.resolve("results") / "v3"
    exigir_salida_v3(salida / "marca.txt")
    salida.mkdir(parents=True, exist_ok=True)
    filas_log: list[dict] = []
    k1 = _k1(settings, prereg, salida, filas_log)
    plasma = settings.resolve(settings.datasets["datasets"]["serum_pxd055218"]["matrix_path"])
    if not plasma.is_file():
        _json(salida / "k4_potencia.json", {"estado": "no_disponible", "motivo": str(plasma)})
    else:
        cruda = pd.read_csv(plasma, sep="\t")
        cruda = cruda.drop_duplicates("protein_id", keep="first").set_index("protein_id")
        valores, escala = a_escala_log2(cruda.to_numpy(dtype=float))
        mat = pd.DataFrame(valores, index=cruda.index, columns=cruda.columns)
        grupos = _grupos_suero(list(mat.columns))
        _json(salida / "suero_escala.json", {"escala": escala, "grupos": {k: len(v) for k, v in grupos.items()}})
        print("[v3] K4 potencia y TOST")
        k4 = _k4(mat, grupos, prereg, salida)
        print("[v3] K5 AUC")
        k5 = _k5(mat, grupos, prereg, salida)
        print("[v3] K6 faltantes")
        k6 = _k6(mat, grupos, prereg, salida)
        print("[v3] K7 control técnico")
        k7 = _k7(mat, grupos, prereg, settings, salida, escala)
        _changelog_suero(filas_log, k4, k5, k6, k7)
    _escribir(pd.DataFrame(filas_log, columns=list(CHANGELOG)), salida / "changelog_numeros.tsv")
    _hipotesis(salida, k1)
    indice = []
    if k1.get("estado") == "calculado":
        for codigo, bloque in (
            ("R1_primario_frente_a_mucosa", k1["contrastes"]["primario_frente_a_mucosa"]),
            ("R1_metastasis_frente_a_mucosa", k1["contrastes"]["metastasis_hepatica_frente_a_mucosa"]),
        ):
            indice.append({
                "codigo": codigo,
                "valor": bloque["valor"],
                "numerador": bloque["numerador"],
                "denominador": bloque["denominador"],
                "unidad_denominador": "proteinas_unicas",
                "ci_low": bloque["ci_low"],
                "ci_high": bloque["ci_high"],
                "cobertura": bloque.get("cobertura"),
                "estado": "calculado",
            })
        indice.append({
            "codigo": "R1_linea_base_h5",
            "valor": None,
            "numerador": None,
            "denominador": None,
            "unidad_denominador": "proteinas_unicas",
            "estado": "no_evaluable",
            "motivo": "LinkedOmicsKB todavía no está leído.",
        })
        indice.append({
            "codigo": "R2_detectabilidad_peptidos",
            **{k: k1["r2_peptidos"][k] for k in ("valor", "numerador", "denominador", "unidad_denominador", "ci_low", "ci_high")},
            "estado": "calculado",
        })
        indice.append({
            "codigo": "R3_homogeneidad",
            "valor": None,
            "numerador": None,
            "denominador": None,
            "unidad_denominador": "proteinas_unicas",
            "estado": "no_evaluable",
            "motivo": k1["r3"]["motivo"],
        })
    _escribir(pd.DataFrame(indice), salida / "indice_R_v3.tsv")
    print(f"[v3] escrito en {salida}")
