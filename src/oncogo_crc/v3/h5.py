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
    k1 = salida / "k1_etiquetas.json"
    concordancia_colon = None
    if k1.is_file():
        import json
        bloque = json.loads(k1.read_text(encoding="utf-8"))["contrastes"]["primario_frente_a_mucosa"]
        concordancia_colon = float(bloque["estimacion"])
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
    concordancias = {}
    for cancer in disponibles:
        otro = efectos[cancer]["log2fc"].rename("log2fc_otro")
        cruzado = tanaka.merge(otro, left_on="gen", right_index=True, how="inner")
        fraccion = concordancia_de_signo(
            cruzado["log2fc"].to_numpy(),
            cruzado["log2fc_otro"].to_numpy(),
            nombre=f"H5_{cancer}",
            universo=n,
            universo_nombre="certificadas_clasificadas",
        )
        concordancias[cancer] = {"estimacion": fraccion.valor, "numerador": fraccion.numerador, "denominador": fraccion.denominador}
    max_otro = max((v["estimacion"] for v in concordancias.values() if v["estimacion"] is not None), default=float("nan"))
    umbral_pan = float(prereg.get("hipotesis.h5.umbral_fraccion_pan_carcinoma"))
    diferencia = float(prereg.get("hipotesis.h5.diferencia_minima_frente_a_colon"))
    if concordancia_colon is None or not np.isfinite(max_otro):
        estado = "no_evaluable"
    elif fraccion_pan < umbral_pan and (concordancia_colon - max_otro) >= diferencia:
        estado = "refuta"
    else:
        estado = "sostiene"
    rng = np.random.default_rng(int(prereg.get("semillas.base")) + int(prereg.get("enmienda_v3_2.semillas.desplazamientos.h5_submuestreo")))
    n_colon = int(colon["n_pares"].iloc[0]) if len(colon) else 0
    iguales = 0
    for _ in range(N_REPETICIONES):
        efectos_sub = {COLON: colon}
        for cancer in disponibles:
            tumor, normal = matrices[cancer]
            comunes = [c for c in tumor.columns if c in normal.columns]
            if len(comunes) < n_colon or n_colon < min_pares:
                break
            elegidas = list(rng.choice(comunes, size=n_colon, replace=False))
            sub = _efecto_sub(tumor, normal, elegidas, min_pares, lfc, alpha)
            sub["gen"] = sub["protein_id"].map(mapa)
            efectos_sub[cancer] = sub.dropna(subset=["gen"]).drop_duplicates("gen").set_index("gen")
        if len(efectos_sub) != 1 + len(disponibles):
            continue
        clases = []
        for gen in clasificadas["gen"]:
            if gen not in colon.index:
                continue
            otros = []
            for cancer in disponibles:
                bloque = efectos_sub[cancer]
                if gen not in bloque.index:
                    otros.append((float("nan"), False))
                else:
                    otros.append((float(bloque.at[gen, "log2fc"]), bool(bloque.at[gen, "significant"])))
            clases.append(clase_de_proteina(bool(colon.at[gen, "significant"]), float(colon.at[gen, "log2fc"]), otros))
        frac = clases.count("pan_carcinoma") / len(clases) if clases else float("nan")
        if np.isfinite(frac) and (frac < umbral_pan) == (fraccion_pan < umbral_pan):
            iguales += 1
    volcar_tsv(tabla, salida / "h5_especificidad_organo.tsv", sello)
    payload = {
        "estado": estado,
        "n_clasificadas": n,
        "conteo": {k: int(conteo.get(k, 0)) for k in ("pan_carcinoma", "especifica_colon", "intermedia")},
        "fraccion_pan_carcinoma": fraccion_pan,
        "concordancia_colon_tanaka": concordancia_colon,
        "concordancia_otros": concordancias,
        "maxima_concordancia_otro_organo": max_otro,
        "fraccion_de_submuestreos_con_pan_al_mismo_lado_de_0_20": iguales / N_REPETICIONES,
        "n_submuestreos": N_REPETICIONES,
        "n_pares_colon": n_colon,
        "sostiene_significa": "No se cumple la refutación preregistrada.",
        "carcinomas_disponibles": disponibles,
        "carcinomas_excluidos": excluidos,
        "diagnostico_pares": diagnostico,
    }
    volcar_json(payload, salida / "h5_especificidad_organo.json", sello)
    return payload
