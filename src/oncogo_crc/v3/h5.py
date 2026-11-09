"""H5. Especificidad de órgano. La regla queda fija antes de leer los efectos.

Las matrices de LinkedOmicsKB ya están en log2. No se aplica un segundo log2.
Una casilla vacía es ausencia. Un Ensembl repetido se elimina entero.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.moderated import moderate_table, paired_observed_stats
from oncogo_crc.v3.concordancia import concordancia_de_signo
from oncogo_crc.v3.genes_tanaka import gen_unico, mapa_genes
from oncogo_crc.v3.prereg import PreregistroV3
from oncogo_crc.v3.salida import volcar_json, volcar_tsv

BASE = "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
OTROS = ("CCRCC", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")
COLON = "COAD"
N_REPETICIONES = 200
MIN_CARCINOMAS = 6
UMBRAL_PAN = 4


def clase_de_proteina(treat_colon: bool, signo_colon: float, otros: list[tuple[float, bool]], *, min_otros: int = MIN_CARCINOMAS, umbral_pan: int = UMBRAL_PAN) -> str:
    """pan_carcinoma, especifica_colon o intermedia. otros tiene un elemento por carcinoma disponible."""
    if not np.isfinite(signo_colon) or signo_colon == 0 or len(otros) < min_otros:
        return "no_clasificable"
    mismo = 0
    for signo, treat in otros:
        if treat and np.isfinite(signo) and np.sign(signo) == np.sign(signo_colon):
            mismo += 1
    if mismo >= umbral_pan:
        return "pan_carcinoma"
    if treat_colon and mismo == 0:
        return "especifica_colon"
    return "intermedia"


def _bajar(url: str, destino: Path) -> None:
    destino.parent.mkdir(parents=True, exist_ok=True)
    if destino.is_file() and destino.stat().st_size > 1000:
        return
    subprocess.run(
        ["curl", "-fL", "--retry", "5", "--retry-delay", "2", "--max-time", "180", "-o", str(destino), url],
        check=True,
    )


def _mapa_ensembl(path: Path) -> dict[str, str]:
    tabla = pd.read_csv(path, sep="\t", usecols=["symbol", "status", "locus_type", "ensembl_gene_id"], dtype=str)
    tabla = tabla.loc[
        tabla["status"].eq("Approved")
        & tabla["locus_type"].eq("gene with protein product")
        & tabla["ensembl_gene_id"].fillna("").str.startswith("ENSG")
        & tabla["symbol"].fillna("").ne("")
    ]
    conteo_e = tabla.groupby("ensembl_gene_id")["symbol"].nunique()
    conteo_s = tabla.groupby("symbol")["ensembl_gene_id"].nunique()
    unicos = set(conteo_e[conteo_e == 1].index) & set(tabla.loc[tabla["symbol"].isin(conteo_s[conteo_s == 1].index), "ensembl_gene_id"])
    limpia = tabla.loc[tabla["ensembl_gene_id"].isin(unicos)]
    return dict(zip(limpia["ensembl_gene_id"], limpia["symbol"].str.upper(), strict=True))


def _leer_matriz(path: Path) -> pd.DataFrame:
    tabla = pd.read_csv(path, sep="\t", low_memory=False)
    ident = tabla.columns[0]
    tabla = tabla.rename(columns={ident: "ensembl"})
    tabla["ensembl"] = tabla["ensembl"].astype(str).str.split(".").str[0]
    repetidos = set(tabla.loc[tabla["ensembl"].duplicated(keep=False), "ensembl"])
    tabla = tabla.loc[~tabla["ensembl"].isin(repetidos)]
    return tabla.drop_duplicates("ensembl").set_index("ensembl")


def _diagnostico_pares(tumor: pd.DataFrame, normal: pd.DataFrame) -> dict:
    comunes = [c for c in tumor.columns if c in normal.columns]
    return {
        "n_tumor": int(tumor.shape[1]),
        "n_normal": int(normal.shape[1]),
        "n_comunes": len(comunes),
        "ejemplos_tumor": list(map(str, tumor.columns[:3])),
        "ejemplos_normal": list(map(str, normal.columns[:3])),
    }


def _efecto(tumor: pd.DataFrame, normal: pd.DataFrame, min_pares: int, lfc: float, alpha: float) -> pd.DataFrame | None:
    comunes = [c for c in tumor.columns if c in normal.columns]
    if len(comunes) < min_pares:
        return None
    genes = tumor.index.intersection(normal.index)
    dif = tumor.loc[genes, comunes].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    dif = dif - normal.loc[genes, comunes].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    coef, s2, df, stdev, n_obs = paired_observed_stats(dif, min_pairs=min_pares)
    tabla, _, _ = moderate_table(
        protein_ids=np.array(genes.astype(str)),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method="moderated_t_treat_paired_h5",
    )
    tabla["n_pares"] = len(comunes)
    return tabla


def _efecto_sub(tumor: pd.DataFrame, normal: pd.DataFrame, columnas: list[str], min_pares: int, lfc: float, alpha: float) -> pd.DataFrame:
    genes = tumor.index.intersection(normal.index)
    dif = tumor.loc[genes, columnas].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    dif = dif - normal.loc[genes, columnas].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    coef, s2, df, stdev, n_obs = paired_observed_stats(dif, min_pairs=min_pares)
    tabla, _, _ = moderate_table(
        protein_ids=np.array(genes.astype(str)),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method="moderated_t_treat_paired_h5_sub",
    )
    return tabla


def _pares_comunes(tumor: pd.DataFrame, normal: pd.DataFrame) -> list[str]:
    return [c for c in tumor.columns if c in normal.columns]


def concordancias_mismo_universo(
    tanaka: pd.Series, efectos: dict[str, pd.DataFrame], carcinomas: list[str]
) -> tuple[dict, list[str]]:
    """tanaka: log2FC P-N indexado por gen. El mismo universo para todos los carcinomas."""
    universo = pd.Index(tanaka.index[np.isfinite(tanaka.to_numpy())])
    for cancer in carcinomas:
        bloque = efectos[cancer]["log2fc"]
        universo = universo.intersection(bloque.index[np.isfinite(bloque.to_numpy())])
    salida: dict[str, dict] = {}
    for cancer in carcinomas:
        a = np.sign(tanaka.loc[universo].to_numpy(dtype=float))
        b = np.sign(efectos[cancer].loc[universo, "log2fc"].to_numpy(dtype=float))
        ok = (a != 0) & (b != 0)
        fraccion = concordancia_de_signo(
            a[ok], b[ok], nombre=f"H5_{cancer}", universo=int(ok.sum()), universo_nombre="tanaka_treat_comun"
        )
        bajo, alto = fraccion.intervalo()
        salida[cancer] = {
            "estimacion": fraccion.valor,
            "numerador": fraccion.numerador,
            "denominador": fraccion.denominador,
            "ci_low": bajo,
            "ci_high": alto,
        }
    return salida, list(universo)


def submuestreo(
    matrices, todos, genes, mapa, rng, *, n_rep, min_pares, lfc, alpha,
    umbral_pan, fraccion_pan, colon="COAD", umbral_clase=4,
) -> dict:
    """Submuestrea todos los carcinomas, colon incluido, al mismo número de pares."""
    n_sub = min(len(_pares_comunes(*matrices[k])) for k in todos)
    if n_sub < min_pares:
        return {"estado": "no_calculado", "motivo": f"n_sub={n_sub} < {min_pares}", "n_sub": int(n_sub)}
    otros = [k for k in todos if k != colon]
    fracs = []
    n_sig = {k: [] for k in todos}
    for _ in range(n_rep):
        sub = {}
        for k in todos:
            tumor, normal = matrices[k]
            elegidas = list(rng.choice(_pares_comunes(tumor, normal), size=n_sub, replace=False))
            tabla = _efecto_sub(tumor, normal, elegidas, min_pares, lfc, alpha)
            tabla["gen"] = tabla["protein_id"].map(mapa)
            sub[k] = tabla.dropna(subset=["gen"]).drop_duplicates("gen").set_index("gen")
            n_sig[k].append(int(sub[k]["significant"].sum()) if len(sub[k]) else 0)
        if len(sub) != len(todos):
            raise RuntimeError("submuestreo incompleto: faltan carcinomas")
        clases = []
        for gen in genes:
            if gen not in sub[colon].index or not np.isfinite(sub[colon].at[gen, "log2fc"]):
                continue
            observados = [
                (float(sub[k].at[gen, "log2fc"]), bool(sub[k].at[gen, "significant"]))
                if gen in sub[k].index else (float("nan"), False)
                for k in otros
            ]
            clases.append(clase_de_proteina(
                bool(sub[colon].at[gen, "significant"]), float(sub[colon].at[gen, "log2fc"]), observados,
                min_otros=len(otros), umbral_pan=umbral_clase,
            ))
        clases = [c for c in clases if c != "no_clasificable"]
        if clases:
            fracs.append(clases.count("pan_carcinoma") / len(clases))
    fracs_np = np.asarray(fracs, dtype=float)
    if len(fracs_np) < n_rep:
        raise RuntimeError(f"submuestreo incompleto: {len(fracs_np)}/{n_rep}")
    mismo = float(np.mean((fracs_np < umbral_pan) == (fraccion_pan < umbral_pan)))
    return {
        "estado": "calculado",
        "n_sub": int(n_sub),
        "n_repeticiones": int(n_rep),
        "fraccion_mismo_lado_0_20": mismo,
        "mediana_fraccion_pan": float(np.median(fracs_np)),
        "p2_5": float(np.quantile(fracs_np, 0.025)),
        "p97_5": float(np.quantile(fracs_np, 0.975)),
        "potencia_submuestreada": {
            k: {"n_pares": int(n_sub), "n_significativas_mediana": float(np.median(v))} for k, v in n_sig.items()
        },
    }


def ejecutar(settings: Settings, salida: Path, sello: dict, prereg: PreregistroV3) -> dict:
    destino = settings.resolve("external/linkedomics")
    hgnc = settings.resolve("external/hgnc/hgnc_complete_set.txt")
    try:
        _bajar(HGNC_URL, hgnc)
        for cancer in (COLON, *OTROS):
            for lado in ("Tumor", "Normal"):
                nombre = f"{cancer}_proteomics_gene_abundance_log2_reference_intensity_normalized_{lado}.txt"
                _bajar(f"{BASE}/{cancer}/{nombre}", destino / nombre)
    except (subprocess.CalledProcessError, OSError) as exc:
        payload = {"estado": "no_disponible", "motivo": f"la descarga no terminó: {exc}"}
        volcar_json(payload, salida / "h5_especificidad_organo.json", sello)
        return payload
    mapa = _mapa_ensembl(hgnc)
    alpha = float(prereg.get("inferencia.alpha"))
    lfc = float(prereg.get("inferencia.treat_lfc"))
    min_pares = int(prereg.get("enmienda_v3_2.tejido.min_pares_observados"))
    efectos: dict[str, pd.DataFrame] = {}
    matrices: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}
    excluidos: dict[str, str] = {}
    diagnostico: dict[str, dict] = {}
    disponibles: list[str] = []
    for cancer in (COLON, *OTROS):
        tumor = _leer_matriz(destino / f"{cancer}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt")
        normal = _leer_matriz(destino / f"{cancer}_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt")
        matrices[cancer] = (tumor, normal)
        diagnostico[cancer] = _diagnostico_pares(tumor, normal)
        tabla = _efecto(tumor, normal, min_pares, lfc, alpha)
        if tabla is None:
            excluidos[cancer] = f"no tiene {min_pares} pares con el mismo identificador"
            continue
        tabla["gen"] = tabla["protein_id"].map(mapa)
        efectos[cancer] = tabla.dropna(subset=["gen"]).drop_duplicates("gen").set_index("gen")
        if cancer != COLON:
            disponibles.append(cancer)
    if COLON in excluidos or len(disponibles) < MIN_CARCINOMAS:
        motivo = excluidos.get(COLON) or f"hay {len(disponibles)} carcinomas y el mínimo es {MIN_CARCINOMAS}"
        payload = {
            "estado": "no_disponible",
            "motivo": motivo,
            "carcinomas_disponibles": disponibles,
            "carcinomas_excluidos": excluidos,
            "diagnostico_pares": diagnostico,
        }
        volcar_json(payload, salida / "h5_especificidad_organo.json", sello)
        return payload
    certificadas = pd.read_csv(settings.resolve("results/certified/cptac_unshared_treat.tsv"), sep="\t")
    genes = sorted(certificadas.loc[certificadas["significant"].astype(bool), "protein_id"].astype(str).str.upper().unique())
    colon = efectos[COLON]
    filas = []
    for gen in genes:
        if gen not in colon.index or not np.isfinite(colon.at[gen, "log2fc"]) or colon.at[gen, "log2fc"] == 0:
            filas.append({"gen": gen, "clase": "no_clasificable", "log2fc_colon": colon.at[gen, "log2fc"] if gen in colon.index else np.nan})
            continue
        otros = []
        for cancer in disponibles:
            if gen not in efectos[cancer].index:
                otros.append((float("nan"), False))
            else:
                otros.append((float(efectos[cancer].at[gen, "log2fc"]), bool(efectos[cancer].at[gen, "significant"])))
        filas.append({
            "gen": gen,
            "clase": clase_de_proteina(bool(colon.at[gen, "significant"]), float(colon.at[gen, "log2fc"]), otros),
            "log2fc_colon": float(colon.at[gen, "log2fc"]),
            "treat_colon": bool(colon.at[gen, "significant"]),
        })
    tabla = pd.DataFrame(filas)
    clasificadas = tabla.loc[tabla["clase"] != "no_clasificable"]
    conteo = clasificadas["clase"].value_counts().to_dict()
    n = int(len(clasificadas))
    fraccion_pan = float(conteo.get("pan_carcinoma", 0) / n) if n else float("nan")
    import json
    k1 = salida / "k1_etiquetas.json"
    concordancia_k1 = None
    if k1.is_file():
        bloque_k1 = json.loads(k1.read_text(encoding="utf-8"))["contrastes"]["primario_frente_a_mucosa"]
        concordancia_k1 = float(bloque_k1["estimacion"])
    tanaka = pd.read_csv(salida / "k1_primario_frente_a_mucosa.tsv", sep="\t")
    mapa_tejido = mapa_genes(
        settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt"),
        settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv"),
    )
    tanaka = tanaka.assign(gen=tanaka["protein_id"].map(lambda pid: gen_unico(pid, mapa_tejido)))
    columna_n = next((c for c in ("n_obs", "n_pares_observados", "n_pares") if c in tanaka.columns), None)
    tanaka = tanaka.dropna(subset=["gen"])
    if columna_n is not None:
        tanaka = tanaka.sort_values(columna_n, ascending=False, kind="stable")
    tanaka = tanaka.drop_duplicates("gen", keep="first")
    tanaka["log2fc"] = pd.to_numeric(tanaka["log2fc"], errors="coerce")
    serie_amplia = tanaka.set_index("gen")["log2fc"]
    if "significant" in tanaka.columns:
        marca = tanaka["significant"].map(lambda v: str(v).lower() in {"true", "1"} or v is True)
        serie_principal = tanaka.loc[marca].set_index("gen")["log2fc"]
    else:
        serie_principal = serie_amplia
    carcinomas = [COLON, *disponibles]
    concordancias, universo = concordancias_mismo_universo(serie_principal, efectos, carcinomas)
    concordancias_amplias, universo_amplio = concordancias_mismo_universo(serie_amplia, efectos, carcinomas)
    concordancia_colon = concordancias.get(COLON, {}).get("estimacion")
    max_otro = max(
        (v["estimacion"] for k, v in concordancias.items() if k != COLON and v["estimacion"] is not None),
        default=float("nan"),
    )
    umbral_pan = float(prereg.get("hipotesis.h5.umbral_fraccion_pan_carcinoma"))
    diferencia = float(prereg.get("hipotesis.h5.diferencia_minima_frente_a_colon"))
    if concordancia_colon is None or not np.isfinite(max_otro):
        estado = "no_evaluable"
    elif fraccion_pan < umbral_pan and (concordancia_colon - max_otro) >= diferencia:
        estado = "refuta"
    else:
        estado = "sostiene"
    potencia = {}
    for cancer, bloque in efectos.items():
        n_sig = int(bloque["significant"].sum()) if len(bloque) else 0
        potencia[cancer] = {
            "n_pares": int(bloque["n_pares"].iloc[0]) if len(bloque) else 0,
            "n_significativas": n_sig,
            "fraccion_significativas": (n_sig / len(bloque)) if len(bloque) else None,
        }
    rng = np.random.default_rng(
        int(prereg.get("semillas.base")) + int(prereg.get("enmienda_v3_4.semillas.desplazamientos.h5_submuestreo_v3_4"))
    )
    sub = submuestreo(
        matrices, carcinomas, list(universo), mapa, rng,
        n_rep=N_REPETICIONES, min_pares=min_pares, lfc=lfc, alpha=alpha,
        umbral_pan=umbral_pan, fraccion_pan=fraccion_pan, colon=COLON, umbral_clase=UMBRAL_PAN,
    )
    sin_ov = [cancer for cancer in carcinomas if cancer != "OV"]
    sub_sin_ov = submuestreo(
        matrices, sin_ov, list(universo), mapa, rng,
        n_rep=N_REPETICIONES, min_pares=min_pares, lfc=lfc, alpha=alpha,
        umbral_pan=umbral_pan, fraccion_pan=fraccion_pan, colon=COLON, umbral_clase=4,
    )
    sub_sin_ov["decide_el_veredicto"] = False
    volcar_tsv(tabla, salida / "h5_especificidad_organo.tsv", sello)
    payload = {
        "estado": estado,
        "n_clasificadas": n,
        "conteo": {k: int(conteo.get(k, 0)) for k in ("pan_carcinoma", "especifica_colon", "intermedia")},
        "fraccion_pan_carcinoma": fraccion_pan,
        "concordancia_colon_tanaka": concordancia_colon,
        "concordancia_otros": {k: v for k, v in concordancias.items() if k != COLON},
        "concordancias_universo_principal": concordancias,
        "universo_principal_n": len(universo),
        "concordancias_sensibilidad": {"n": len(universo_amplio), "concordancias": concordancias_amplias},
        "concordancia_k1_referencia": {
            "estimacion": concordancia_k1,
            "nota": "seleccionada por CPTAC; no entra en la regla",
        },
        "maxima_concordancia_otro_organo": max_otro,
        "submuestreo": sub,
        "sensibilidad_sin_ov": sub_sin_ov,
        "potencia_por_carcinoma": potencia,
        "fraccion_de_submuestreos_con_pan_al_mismo_lado_de_0_20": sub.get("fraccion_mismo_lado_0_20"),
        "n_submuestreos": sub.get("n_repeticiones", 0),
        "n_pares_colon": potencia.get(COLON, {}).get("n_pares"),
        "sostiene_significa": "No se cumple la refutación preregistrada.",
        "carcinomas_disponibles": disponibles,
        "carcinomas_excluidos": excluidos,
        "diagnostico_pares": diagnostico,
    }
    volcar_json(payload, salida / "h5_especificidad_organo.json", sello)
    return payload
