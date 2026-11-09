"""Controles de la sección 8. Si un control positivo falla, la etapa se detiene."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.stages.v3_3_tejido import ids_descartados
from oncogo_crc.utils.moderated import moderate_table, paired_observed_stats
from oncogo_crc.v3.etiquetas import _pares, filas_desde_columnas, matriz_log2_lfq
from oncogo_crc.v3.genes_tanaka import filas_de_gen, mapa_genes
from oncogo_crc.v3.prereg import PreregistroV3
from oncogo_crc.utils.sello import exigir_arbol_limpio, hacer_sello
from oncogo_crc.v3.salida import carpeta, exigir_vigente, manifiesto_congelado, volcar_json, volcar_tsv

MUCOSA_BAJA = ("CA1", "AQP8", "SLC26A3", "GUCA2A")
N_PERMUTACIONES = 200
SEMILLA_PERMUTACION = 20260914 + 701


def _abrir(settings: Settings, etapa: str):
    prereg = PreregistroV3(settings.root)
    salida = carpeta(settings)
    exigir_vigente(salida / "marca.txt")
    salida.mkdir(parents=True, exist_ok=True)
    manifest = manifiesto_congelado(settings)
    sello = hacer_sello(settings.root, etapa=etapa, manifest=manifest)
    exigir_arbol_limpio(sello)
    return prereg, salida, sello


def _genes_bajos(settings: Settings, salida) -> list[dict]:
    contraste = pd.read_csv(salida / "k1_primario_frente_a_mucosa.tsv", sep="\t")
    contraste["log2fc"] = pd.to_numeric(contraste["log2fc"], errors="coerce")
    mapa = mapa_genes(
        settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt"),
        settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv"),
    )
    columna_n = next((c for c in ("n_obs", "n_pares_observados", "n_pares") if c in contraste.columns), None)
    filas = []
    for gen in MUCOSA_BAJA:
        ids = set(filas_de_gen(contraste["protein_id"].astype(str), mapa, gen))
        bloque = contraste.loc[contraste["protein_id"].astype(str).isin(ids) & np.isfinite(contraste["log2fc"])]
        base = {
            "analisis": "tumor_frente_a_mucosa",
            "control": gen,
            "tipo": "positivo",
            "esperado": "log2fc menor que 0",
            "n_filas_del_gen": len(ids),
        }
        if bloque.empty:
            filas.append({**base, "observado": "no_evaluable", "pasa": "falla"})
            continue
        if columna_n is not None:
            bloque = bloque.sort_values(columna_n, ascending=False, kind="stable")
        valor = float(bloque["log2fc"].iloc[0])
        filas.append({
            **base,
            "protein_id": str(bloque["protein_id"].iloc[0]),
            "observado": valor,
            "pasa": "pasa" if valor < 0 else "falla",
        })
    return filas


def _permutacion(settings: Settings, prereg: PreregistroV3) -> dict:
    ds = settings.datasets["datasets"]["pride_crc_meta"]
    cruda = pd.read_csv(settings.resolve(ds["matrix_path"]), sep="\t", low_memory=False)
    grupos = settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt")
    flags = None
    if grupos.is_file():
        flags = pd.read_csv(
            grupos,
            sep="\t",
            usecols=lambda c: c in {"Protein IDs", "Reverse", "Potential contaminant", "Only identified by site"},
            low_memory=False,
        )
    malos = ids_descartados(cruda["protein_id"], flags)
    cruda = cruda.loc[~cruda["protein_id"].astype(str).isin(malos)]
    meta = filas_desde_columnas([c for c in cruda.columns if c != "protein_id"])
    log2mat = matriz_log2_lfq(cruda)
    pares = _pares(meta, "primary", "mucosa")
    dif = log2mat[[a for a, _ in pares]].to_numpy(dtype=float) - log2mat[[b for _, b in pares]].to_numpy(dtype=float)
    lfc = float(prereg.get("inferencia.treat_lfc"))
    alpha = float(prereg.get("inferencia.alpha"))
    min_pares = int(prereg.get("enmienda_v3_2.tejido.min_pares_observados"))
    rng = np.random.default_rng(SEMILLA_PERMUTACION)
    recuentos = []
    for _ in range(N_PERMUTACIONES):
        signos = rng.choice(np.array([-1.0, 1.0]), size=dif.shape[1])
        coef, s2, df, stdev, n_obs = paired_observed_stats(dif * signos, min_pairs=min_pares)
        tabla, _, _ = moderate_table(
            protein_ids=log2mat.index.to_numpy(),
            coef=coef,
            s2=s2,
            df=df,
            stdev_unscaled=stdev,
            n_obs=n_obs,
            lfc=lfc,
            alpha=alpha,
            stat_method="permutacion_de_par",
        )
        recuentos.append(int(tabla["significant"].fillna(False).sum()))
    mediana = float(np.median(recuentos))
    return {
        "analisis": "tumor_frente_a_mucosa",
        "control": "permutacion_dentro_del_par",
        "tipo": "negativo",
        "esperado": "mediana de proteínas significativas igual a 0",
        "observado": mediana,
        "pasa": "pasa" if mediana == 0 else "falla",
        "n_permutaciones": N_PERMUTACIONES,
        "percentil_95": float(np.quantile(recuentos, 0.95)),
    }


def _suero(salida) -> list[dict]:
    traza_path = salida / "suero_trazabilidad_grupos.tsv"
    principal_path = salida / "suero_matriz_principal.tsv"
    k5 = salida / "k5_auc.json"
    if not traza_path.is_file() or not principal_path.is_file() or not k5.is_file():
        return [{
            "analisis": "suero",
            "control": "filtros_y_queratinas",
            "tipo": "positivo",
            "esperado": "el suero v3.2 ya corrido",
            "observado": "no_disponible",
            "pasa": "no_evaluable",
        }]
    traza = pd.read_csv(traza_path, sep="\t")
    principal = set(pd.read_csv(principal_path, sep="\t", usecols=["protein_id"])["protein_id"].astype(str))
    fuera = traza.loc[
        ~traza["humano"].astype(bool)
        | traza["contaminante"].astype(bool)
        | traza["contaminante_piel"].astype(bool)
        | traza["panel"].fillna("").astype(str).str.len().gt(0),
        "protein_group",
    ].astype(str)
    cruce = principal & set(fuera)
    piel = int(traza["contaminante_piel"].astype(bool).sum())
    nulos = json.loads(k5.read_text(encoding="utf-8")).get("nulos", {})
    antes = (nulos.get("filtrada_humana") or {}).get("maximo_observado")
    despues = (nulos.get("sin_contaminantes_ni_paneles") or {}).get("maximo_observado")
    baja = antes is not None and despues is not None and float(despues) < float(antes)
    return [
        {
            "analisis": "suero_filtro",
            "control": "sin_no_humano_ni_contaminante_en_la_matriz_principal",
            "tipo": "positivo",
            "esperado": "ningún identificador de la matriz principal está marcado",
            "observado": len(cruce),
            "pasa": "pasa" if len(cruce) == 0 else "falla",
        },
        {
            "analisis": "suero_contaminacion",
            "control": "queratinas_de_piel",
            "tipo": "positivo",
            "esperado": "hay contaminante_piel y el máximo AUC baja al quitarlas",
            "observado": {"n_piel": piel, "auc_antes": antes, "auc_despues": despues},
            "pasa": "pasa" if piel > 0 and baja else "falla",
        },
    ]


def run(settings: Settings) -> None:
    prereg, salida, sello = _abrir(settings, "v3_3_controles")
    if not (salida / "k1_primario_frente_a_mucosa.tsv").is_file():
        raise FileNotFoundError("falta el tejido v3.2; corre v3-2-tejido antes")
    filas = _genes_bajos(settings, salida)
    print(f"[v3.2] control negativo: {N_PERMUTACIONES} permutaciones del par")
    filas.append(_permutacion(settings, prereg))
    filas.extend(_suero(salida))
    tabla = pd.DataFrame(filas)
    volcar_tsv(tabla, salida / "controles.tsv", sello)
    fallos = tabla.loc[tabla["pasa"].eq("falla")]
    volcar_json({"estado": "calculado", "n_fallos": int(len(fallos))}, salida / "controles.json", sello)
    if len(fallos):
        raise RuntimeError(f"un control no salió como debía: {fallos['control'].tolist()}")
    print(f"[v3.2] controles escritos, fallos {len(fallos)}")
