"""Cierre ejecutable de la tesis. Escribe solo en results/tesis_cierre/.

No toca el certificado de v1 ni las tablas de v2. Si la matriz de suero no está,
o si sus etiquetas son el corte alfabético de PXD036594, deja el motivo por
escrito y no publica un contraste biológico.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from oncogo_crc.certification.contract import TREAT_LFC
from oncogo_crc.settings import Settings
from oncogo_crc.tesis.cierre import (
    a_escala_log2,
    cargar_catalogo,
    cruzar_catalogo,
    etiqueta_tanaka,
    grupos_de_muestras,
    intervalo_auc,
    lectura_con_ruido,
    localizar_candidata,
    puntuacion_enriquecimiento,
    tabla_neoantigenos,
    veredicto_especificidad,
)
from oncogo_crc.utils.go_utils import load_go_long
from oncogo_crc.utils.ids import effective_protein_mapping
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.moderated import moderate_table, paired_observed_stats, unpaired_observed_stats
from oncogo_crc.utils.reproducibility_overlap import load_gene_uniprot_from_mapping, uniprot_base
from oncogo_crc.stages.certified_inference import _numeric_block

MIN_POR_GRUPO = 8


def _contraste(mat: pd.DataFrame, cols_a: list[str], cols_b: list[str], *, lfc: float, alpha: float, metodo: str) -> pd.DataFrame:
    a_raw = _numeric_block(mat, cols_a)
    b_raw = _numeric_block(mat, cols_b)
    a, escala = a_escala_log2(a_raw)
    b, _ = a_escala_log2(b_raw)
    coef, s2, df, stdev, n_obs = unpaired_observed_stats(a, b, min_n=3)
    table, _d0, _s02 = moderate_table(
        protein_ids=mat.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method=metodo,
    )
    table["escala"] = escala
    return table


def _especificidad_de_candidatas(
    mat: pd.DataFrame,
    grupos: dict[str, list[str]],
    catalogo: pd.DataFrame,
    *,
    semilla: int,
    mapa_uniprot: dict[str, str] | None = None,
) -> pd.DataFrame:
    valores, _ = a_escala_log2(mat.to_numpy(dtype=float))
    bloque = pd.DataFrame(valores, index=mat.index, columns=mat.columns)
    mapa = {k.upper(): v.upper() for k, v in (mapa_uniprot or {}).items()}
    identificadores = pd.Series(bloque.index.astype(str), index=bloque.index)
    filas = []
    for _, cand in catalogo.iterrows():
        gen = str(cand["gene"]).upper()
        base = {
            "gene": cand["gene"],
            "origen": "literatura",
            "lista_ruido": cand.get("lista_ruido", ""),
        }
        idx = localizar_candidata(identificadores, gen, mapa.get(gen))
        if idx is None:
            base["veredicto"] = "no_evaluable"
            base["lectura"] = "no_evaluable"
            filas.append(base)
            continue
        fila = bloque.loc[idx]
        crc = fila[grupos["crc"]].to_numpy(dtype=float)
        cn = fila[grupos["colonoscopy_negative"]].to_numpy(dtype=float)
        sano = fila[grupos["healthy"]].to_numpy(dtype=float)
        auc_crc, bajo, alto = intervalo_auc(crc, cn, semilla=semilla)
        auc_enf = intervalo_auc(cn, sano, n_boot=0, semilla=semilla)[0]
        base.update(
            {
                "auc_crc_vs_enfermedad": auc_crc,
                "auc_crc_vs_enfermedad_ci_bajo": bajo,
                "auc_crc_vs_enfermedad_ci_alto": alto,
                "auc_enfermedad_vs_sano": auc_enf,
                "veredicto": veredicto_especificidad(auc_crc, bajo, auc_enf),
            }
        )
        base["lectura"] = lectura_con_ruido(base["veredicto"], cand.get("lista_ruido", ""))
        filas.append(base)
    return pd.DataFrame(filas)


def _leer_si_existe(path) -> pd.DataFrame | None:
    if path.is_file():
        return pd.read_csv(path, sep="\t")
    return None


def run(settings: Settings) -> None:
    cfg = settings.config.get("certification", {}) or {}
    alpha = float(cfg.get("alpha", 0.05))
    lfc = float(cfg.get("treat_lfc", TREAT_LFC))
    salida = settings.resolve("results") / "tesis_cierre"
    salida.mkdir(parents=True, exist_ok=True)

    catalogo = cargar_catalogo()
    mapa = load_gene_uniprot_from_mapping(effective_protein_mapping(settings))
    write_tsv(catalogo, salida / "candidatas_literatura.tsv")
    write_tsv(tabla_neoantigenos(), salida / "peptidos_neoantigenos.tsv")

    auditoria: dict = {
        "version": "tesis_cierre",
        "alpha": alpha,
        "treat_lfc": lfc,
        "min_por_grupo": MIN_POR_GRUPO,
        "origen_de_las_candidatas": "literatura",
        "regla_especificidad": (
            "especifica si AUC(CRC vs enfermedad) > 0,70, el límite inferior > 0,5 "
            "y ese AUC supera al de enfermedad vs sano"
        ),
        "neoantigenos": (
            "Los 19 péptidos de Rai 2025 y los neoepítopos de Xiang 2026 no se confirman "
            "con una matriz de grupos proteicos. La tabla peptidos_neoantigenos.tsv fija ese límite."
        ),
        "contrastes": {},
    }

    plasma_ds = settings.datasets["datasets"]["plasma_crc_2024"]
    matriz_path = settings.resolve(plasma_ds["matrix_path"])
    meta_path = settings.sample_metadata_path("plasma_crc_2024")
    metadata = _leer_si_existe(meta_path)

    if not matriz_path.is_file():
        auditoria["estado_suero"] = "matriz_ausente"
        auditoria["matriz"] = str(matriz_path)
        auditoria["motivo"] = (
            "No está la matriz de suero. Con el disco de datos montado, "
            "oncogo-crc run --stage tesis_cierre vuelve a intentar el contraste."
        )
    else:
        cabecera = pd.read_csv(matriz_path, sep="\t", nrows=1)
        if "Run" in cabecera.columns and "protein_id" not in cabecera.columns:
            auditoria["estado_suero"] = "informe_diann_largo"
            auditoria["motivo"] = (
                "El archivo sigue en formato largo de DIA-NN. Hay que pivotarlo a proteína × corrida "
                "antes de contrastar. No se sobrescribe el crudo desde esta etapa."
            )
        elif any(str(c).endswith(" Ratio To Standard") for c in cabecera.columns):
            auditoria["estado_suero"] = "panel_dirigido_sin_etiqueta_clinica"
            auditoria["motivo"] = (
                "PXD036594 es un panel dirigido y el nombre de la corrida no trae el diagnóstico. "
                "No se parte el archivo por orden alfabético."
            )
        else:
            ancha = pd.read_csv(matriz_path, sep="\t", low_memory=False)
            ancha = ancha.drop_duplicates("protein_id", keep="first").set_index("protein_id")
            columnas = [c for c in ancha.columns if c != "protein_id"]
            grupos = grupos_de_muestras(columnas, metadata)
            conteos = {k: len(v) for k, v in grupos.items()}
            auditoria["n_por_grupo"] = conteos
            if min(conteos["crc"], conteos["healthy"]) < MIN_POR_GRUPO:
                auditoria["estado_suero"] = "grupos_insuficientes_o_sin_etiqueta"
                auditoria["motivo"] = (
                    "No hay al menos 8 corridas de cáncer y 8 de sanos con etiqueta clínica. "
                    f"Conteos: {conteos}."
                )
            else:
                pares = {
                    "crc_vs_sano": ("crc", "healthy"),
                    "crc_vs_enfermedad": ("crc", "colonoscopy_negative"),
                    "enfermedad_vs_sano": ("colonoscopy_negative", "healthy"),
                }
                for nombre, (a, b) in pares.items():
                    if len(grupos[a]) < MIN_POR_GRUPO or len(grupos[b]) < MIN_POR_GRUPO:
                        auditoria["contrastes"][nombre] = {
                            "estado": "no_evaluable",
                            "n_a": len(grupos[a]),
                            "n_b": len(grupos[b]),
                        }
                        continue
                    tabla = _contraste(
                        ancha,
                        grupos[a],
                        grupos[b],
                        lfc=lfc,
                        alpha=alpha,
                        metodo=f"tesis_cierre_{nombre}",
                    )
                    write_tsv(tabla, salida / f"{nombre}.tsv")
                    cruce = cruzar_catalogo(catalogo, tabla, cohorte=nombre, mapa_uniprot=mapa)
                    write_tsv(cruce, salida / f"auditoria_{nombre}.tsv")
                    auditoria["contrastes"][nombre] = {
                        "estado": "calculado",
                        "n_a": len(grupos[a]),
                        "n_b": len(grupos[b]),
                        "n_significativas": int(tabla["significant"].sum()),
                        "n_probadas": int(tabla["p_value"].notna().sum()),
                        "alcance": "cribado" if min(len(grupos[a]), len(grupos[b])) < 30 else "cohorte",
                    }
                if all(len(grupos[g]) >= MIN_POR_GRUPO for g in grupos):
                    especificidad = _especificidad_de_candidatas(
                        ancha, grupos, catalogo, semilla=20261031, mapa_uniprot=mapa
                    )
                    write_tsv(especificidad, salida / "especificidad_candidatas.tsv")
                    auditoria["estado_suero"] = "tres_grupos"
                    auditoria["n_especificas"] = int((especificidad["veredicto"] == "especifica").sum())
                    auditoria["n_especificas_sin_bandera"] = int((especificidad["lectura"] == "especifica").sum())
                else:
                    auditoria["estado_suero"] = "sin_grupo_de_enfermedad"
                    auditoria["motivo"] = (
                        "Hay cáncer y sanos, pero no el grupo de colonoscopia negativa. "
                        "La especificidad frente a otra enfermedad intestinal queda no evaluable."
                    )

    tejido = _leer_si_existe(settings.resolve("results") / "certified" / "cptac_unshared_treat.tsv")
    if tejido is not None and {"protein_id", "log2fc", "padj"}.issubset(tejido.columns):
        cruce_tejido = cruzar_catalogo(catalogo, tejido, cohorte="cptac_tejido", mapa_uniprot=mapa)
        write_tsv(cruce_tejido, salida / "auditoria_cptac.tsv")
        auditoria["cptac"] = {
            "estado": "cruzado",
            "n_detectadas": int(cruce_tejido["detectada"].sum()),
            "n_replicadas": int((cruce_tejido["veredicto"] == "replicada").sum()),
        }
    else:
        auditoria["cptac"] = {"estado": "tabla_certificada_ausente"}

    auditoria["pxd034575"] = _tejido_tanaka(settings, mapa, tejido, salida, lfc=lfc, alpha=alpha)
    auditoria["vias_suero"] = _vias_suero(settings, mapa, salida)
    _fichas(settings, mapa, salida)

    (salida / "cierre_audit.json").write_text(
        json.dumps(auditoria, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"[tesis-cierre] {auditoria.get('estado_suero', 'sin_suero')} → {salida}")


def _acceso_a_gen(mapa: dict) -> dict[str, str]:
    salida: dict[str, str] = {}
    for gen, acc in mapa.items():
        if acc is None or (isinstance(acc, float) and np.isnan(acc)):
            continue
        base = str(acc).upper().split("-")[0].strip()
        if base and base not in salida:
            salida[base] = str(gen).upper()
    return salida


def _por_gen(tabla: pd.DataFrame, acceso: dict[str, str]) -> pd.DataFrame:
    filas = []
    for _, fila in tabla.iterrows():
        pid = str(fila["protein_id"])
        gen = acceso.get(uniprot_base(pid) or "")
        log2fc = fila.get("log2fc")
        if not gen or log2fc is None or not np.isfinite(log2fc):
            continue
        p = fila.get("p_value")
        n_obs = fila.get("n_obs")
        filas.append(
            {
                "gene": gen,
                "protein_id": pid,
                "log2fc": float(log2fc),
                "p_value": float(p) if p is not None and np.isfinite(p) else np.nan,
                "padj": float(fila["padj"]) if np.isfinite(fila.get("padj", np.nan)) else np.nan,
                "significant": bool(fila.get("significant", False)),
                "n_obs": float(n_obs) if n_obs is not None and np.isfinite(n_obs) else np.nan,
            }
        )
    if not filas:
        return pd.DataFrame(columns=["gene", "protein_id", "log2fc", "p_value", "padj", "significant", "n_obs"])
    salida = pd.DataFrame(filas).sort_values("n_obs", ascending=False)
    return salida.drop_duplicates("gene").set_index("gene")


def _concordancia(izquierda: pd.DataFrame, derecha: pd.DataFrame) -> dict:
    comun = izquierda.index.intersection(derecha.index)
    a = izquierda.loc[comun, "log2fc"].to_numpy(dtype=float)
    b = derecha.loc[comun, "log2fc"].to_numpy(dtype=float)
    mask = np.isfinite(a) & np.isfinite(b) & (a != 0) & (b != 0)
    if int(mask.sum()) < 3:
        return {"n": int(mask.sum()), "mismo_signo": 0, "fraccion": None, "pearson": None}
    mismo = np.sign(a[mask]) == np.sign(b[mask])
    return {
        "n": int(mask.sum()),
        "mismo_signo": int(mismo.sum()),
        "fraccion": float(mismo.mean()),
        "pearson": float(np.corrcoef(a[mask], b[mask])[0, 1]),
    }


def _pareado(bloque: pd.DataFrame, pacientes: dict, izquierda: str, derecha: str, *, lfc: float, alpha: float, metodo: str) -> pd.DataFrame:
    pares = [
        (tejidos[izquierda], tejidos[derecha])
        for tejidos in pacientes.values()
        if izquierda in tejidos and derecha in tejidos
    ]
    diffs = np.column_stack(
        [bloque[a].to_numpy(dtype=float) - bloque[b].to_numpy(dtype=float) for a, b in pares]
    )
    coef, s2, df, stdev, n_obs = paired_observed_stats(diffs, min_pairs=8)
    tabla, _, _ = moderate_table(
        protein_ids=bloque.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method=metodo,
    )
    return tabla


def _tejido_tanaka(settings: Settings, mapa: dict, cptac: pd.DataFrame | None, salida, *, lfc: float, alpha: float) -> dict:
    ruta = settings.resolve(settings.datasets["datasets"]["pride_crc_meta"]["matrix_path"])
    if not ruta.is_file():
        return {"estado": "matriz_ausente", "ruta": str(ruta)}
    crudo = pd.read_csv(ruta, sep="\t", low_memory=False)
    crudo = crudo.drop_duplicates("protein_id", keep="first").set_index("protein_id")
    pacientes: dict[str, dict[str, str]] = {}
    for columna in crudo.columns:
        etiqueta = etiqueta_tanaka(str(columna))
        if etiqueta is None:
            continue
        paciente, tejido = etiqueta
        pacientes.setdefault(paciente, {})[tejido] = columna
    if len(pacientes) < 8:
        return {"estado": "columnas_no_reconocidas", "n_pacientes": len(pacientes)}
    matriz = crudo.to_numpy(dtype=float)
    # En MaxQuant el cero es «no identificado», no abundancia cero.
    matriz = np.where(matriz > 0, matriz, np.nan)
    valores, escala = a_escala_log2(matriz)
    bloque = pd.DataFrame(valores, index=crudo.index, columns=crudo.columns)
    acceso = _acceso_a_gen(mapa)
    resumen = {
        "estado": "calculado",
        "articulo": "Tanaka et al. Cell Reports 2024. doi:10.1016/j.celrep.2024.113810",
        "leyenda": "P tumor primario, M metastasis hepatica, N mucosa normal. No es polipo.",
        "n_pacientes": len(pacientes),
        "escala": escala,
        "contrastes": {},
    }
    tablas = {}
    for nombre, (izq, der) in {
        "primario_vs_mucosa": ("tumor_primario", "mucosa_normal"),
        "metastasis_vs_mucosa": ("metastasis_hepatica", "mucosa_normal"),
        "metastasis_vs_primario": ("metastasis_hepatica", "tumor_primario"),
    }.items():
        tabla = _pareado(bloque, pacientes, izq, der, lfc=lfc, alpha=alpha, metodo=f"tesis_cierre_{nombre}")
        write_tsv(tabla, salida / f"tanaka_{nombre}.tsv")
        tablas[nombre] = _por_gen(tabla, acceso)
        resumen["contrastes"][nombre] = {
            "n_significativas": int(tabla["significant"].sum()),
            "n_probadas": int(tabla["p_value"].notna().sum()),
            "n_pares": len(pacientes),
        }
    if "primario_vs_mucosa" in tablas and "metastasis_vs_mucosa" in tablas:
        primario = tablas["primario_vs_mucosa"]
        metastasis = tablas["metastasis_vs_mucosa"]
        comun = primario.index.intersection(metastasis.index)
        fp = primario.loc[comun, "log2fc"].to_numpy(dtype=float)
        fm = metastasis.loc[comun, "log2fc"].to_numpy(dtype=float)
        mask = np.isfinite(fp) & np.isfinite(fm) & (fp != 0) & (fm != 0)
        mismo = np.sign(fp[mask]) == np.sign(fm[mask])
        menor = np.abs(fp[mask]) < np.abs(fm[mask])
        resumen["progresion_primario_hacia_metastasis"] = {
            "n": int(mask.sum()),
            "mismo_signo_que_la_metastasis": int(mismo.sum()),
            "fraccion_mismo_signo": float(mismo.mean()) if mask.sum() else None,
            "mismo_signo_y_magnitud_menor": int((mismo & menor).sum()),
            "fraccion_menor_entre_las_de_mismo_signo": float((menor[mismo]).mean()) if mismo.sum() else None,
            "lectura": (
                "El primario se compara con la metástasis del mismo paciente. "
                "No es un gradiente adenoma–carcinoma."
            ),
        }
    if cptac is not None and {"protein_id", "log2fc"}.issubset(cptac.columns):
        cptac_gen = cptac.copy()
        cptac_gen["protein_id"] = cptac_gen["protein_id"].astype(str).str.upper()
        cptac_gen = cptac_gen.drop_duplicates("protein_id").set_index("protein_id")
        resumen["concordancia_con_cptac"] = {
            "primario_vs_mucosa": _concordancia(tablas["primario_vs_mucosa"], cptac_gen),
            "metastasis_vs_mucosa": _concordancia(tablas["metastasis_vs_mucosa"], cptac_gen),
            "nota": (
                "CPTAC es tumor primario frente a mucosa. La comparación justa "
                "es primario_vs_mucosa. metastasis_vs_mucosa es la que v2 usó al llamar tumor a la columna M."
            ),
        }
    return resumen


def _ranking(tabla_gen: pd.DataFrame) -> pd.DataFrame:
    p = tabla_gen["p_value"].clip(lower=1e-300)
    score = np.sign(tabla_gen["log2fc"].to_numpy(dtype=float)) * (-np.log10(p.to_numpy(dtype=float)))
    salida = tabla_gen.copy()
    salida["score"] = score
    return salida.replace([np.inf, -np.inf], np.nan).dropna(subset=["score"]).sort_values("score", ascending=False)


def _correlacion_vias(ranking: pd.DataFrame, terminos: dict[str, np.ndarray], nes_cptac: pd.Series, *, n_perm: int = 200) -> dict:
    genes = ranking.index.to_numpy()
    scores = ranking["score"].to_numpy(dtype=float)
    comunes = [t for t in nes_cptac.index if t in terminos and len(set(terminos[t]) & set(genes)) >= 10]
    if len(comunes) < 8:
        return {"estado": "pocos_terminos", "n_terminos": len(comunes)}
    conjuntos = [np.array(sorted(set(terminos[t]) & set(genes))) for t in comunes]
    nes_ref = nes_cptac.loc[comunes].to_numpy(dtype=float)

    def _es_vector(g_ord: np.ndarray, s_ord: np.ndarray) -> np.ndarray:
        vector = np.empty(len(conjuntos), dtype=float)
        for i, conjunto in enumerate(conjuntos):
            vector[i] = puntuacion_enriquecimiento(s_ord, np.isin(g_ord, conjunto))
        return vector

    orden = np.argsort(-scores)
    es_obs = _es_vector(genes[orden], scores[orden])
    rng = np.random.default_rng(20261031)
    es_null = np.empty((n_perm, len(comunes)), dtype=float)
    for i in range(n_perm):
        shuf = rng.permutation(scores)
        o = np.argsort(-shuf)
        es_null[i] = _es_vector(genes[o], shuf[o])
    escala = np.nanmean(np.abs(es_null), axis=0)
    usable = np.isfinite(es_obs) & np.isfinite(nes_ref) & np.isfinite(escala) & (escala > 0)
    if int(usable.sum()) < 8:
        return {"estado": "sin_escala", "n_terminos": int(usable.sum())}
    nes = np.full(len(comunes), np.nan)
    nes[usable] = es_obs[usable] / escala[usable]
    r_obs = float(np.corrcoef(nes[usable], nes_ref[usable])[0, 1])
    nulos = []
    for i in range(n_perm):
        nes_i = es_null[i, usable] / escala[usable]
        if np.isfinite(nes_i).sum() < 8:
            continue
        nulos.append(float(np.corrcoef(nes_i, nes_ref[usable])[0, 1]))
    nulos_a = np.asarray(nulos, dtype=float)
    nulos_a = nulos_a[np.isfinite(nulos_a)]
    p = float((1 + np.sum(np.abs(nulos_a) >= abs(r_obs))) / (1 + nulos_a.size))
    return {
        "estado": "calculado" if int(usable.sum()) >= 10 else "exploratorio_pocos_terminos",
        "n_terminos": int(usable.sum()),
        "pearson_nes": r_obs,
        "p_nulo": p,
        "n_permutaciones": int(nulos_a.size),
        "terminos": "GO biologicos con FDR <= 0,25 en el GSEA de CPTAC y al menos 10 genes en el ranking",
    }


def _vias_suero(settings: Settings, mapa: dict, salida) -> dict:
    gsea_path = settings.resolve("results/gsea/cptac_tumor_vs_normal_gsea_go_signed_neglog10_p.tsv")
    go_path = settings.go_gene_long_path()
    if not gsea_path.is_file() or not go_path.is_file():
        return {"estado": "faltan_vias_o_go"}
    cptac = pd.read_csv(gsea_path, sep="\t")
    cptac = cptac.loc[cptac["FDR q-val"] <= 0.25, ["Term", "NES"]].drop_duplicates("Term")
    nes_cptac = cptac.set_index("Term")["NES"]
    go = load_go_long(go_path)
    go = go.loc[go["namespace"] == "BP", ["gene_symbol", "go_id"]]
    terminos = {term: genes.to_numpy() for term, genes in go.groupby("go_id")["gene_symbol"]}
    acceso = _acceso_a_gen(mapa)
    resultado = {"estado": "calculado", "contrastes": {}}
    for nombre in ("crc_vs_sano", "crc_vs_enfermedad", "enfermedad_vs_sano"):
        ruta = salida / f"{nombre}.tsv"
        if not ruta.is_file():
            resultado["contrastes"][nombre] = {"estado": "sin_tabla"}
            continue
        tabla = pd.read_csv(ruta, sep="\t")
        ranking = _ranking(_por_gen(tabla, acceso))
        resultado["contrastes"][nombre] = _correlacion_vias(ranking, terminos, nes_cptac)
    return resultado


def _fichas(settings: Settings, mapa: dict, salida) -> None:
    acceso = _acceso_a_gen(mapa)
    gen_a_acc = {gen: acc.upper().split("-")[0] for gen, acc in mapa.items() if isinstance(acc, str) and acc}
    piezas = []
    for nombre in ("crc_vs_sano", "crc_vs_enfermedad", "enfermedad_vs_sano"):
        ruta = salida / f"{nombre}.tsv"
        if not ruta.is_file():
            continue
        tabla = _por_gen(pd.read_csv(ruta, sep="\t"), acceso)
        tabla = tabla.reset_index()
        tabla["contraste"] = nombre
        piezas.append(tabla)
    if not piezas:
        return
    largo = pd.concat(piezas, ignore_index=True)
    especificidad = _leer_si_existe(salida / "especificidad_candidatas.tsv")
    interes = set(largo.loc[largo["significant"], "gene"])
    if especificidad is not None and "lectura" in especificidad.columns:
        interes |= set(especificidad.loc[especificidad["lectura"] != "no_evaluable", "gene"].astype(str).str.upper())
    if not interes:
        interes = set(largo["gene"].head(12))
    go_path = settings.go_gene_long_path()
    nombres_go: dict[str, list[str]] = {}
    if go_path.is_file():
        go = load_go_long(go_path)
        tamanos = go.groupby("go_id")["gene_symbol"].nunique()
        utiles = set(tamanos[(tamanos >= 15) & (tamanos <= 300)].index)
        sub = go[go["gene_symbol"].isin(interes) & go["go_id"].isin(utiles)]
        for gen, bloque in sub.groupby("gene_symbol"):
            elegidos = []
            for espacio, cupo in (("BP", 5), ("CC", 3)):
                trozo = bloque.loc[bloque["namespace"] == espacio].copy()
                trozo["n"] = trozo["go_id"].map(tamanos)
                trozo = trozo.sort_values("n").drop_duplicates("go_id").head(cupo)
                elegidos.extend(f"{fila.go_id}|{espacio}|{int(fila.n)}" for fila in trozo.itertuples(index=False))
            nombres_go[str(gen)] = elegidos
    hpa_path = settings.resolve("external/hpa/hpa_concentracion_sangre.tsv")
    hpa = _leer_si_existe(hpa_path)
    if hpa is not None and "gene" in hpa.columns:
        hpa = hpa.drop_duplicates("gene")
        hpa["gene"] = hpa["gene"].astype(str).str.upper()
        hpa = hpa.set_index("gene")
    filas = []
    for gen in sorted(interes):
        fila = {"gene": gen, "uniprot": gen_a_acc.get(gen, ""), "origen": "resultado_propio"}
        if especificidad is not None and gen in set(especificidad["gene"].astype(str).str.upper()):
            origen = especificidad.loc[especificidad["gene"].astype(str).str.upper() == gen].iloc[0]
            fila["origen"] = "literatura"
            fila["lista_ruido"] = origen.get("lista_ruido", "")
            fila["veredicto"] = origen.get("veredicto", "")
            fila["lectura"] = origen.get("lectura", "")
            fila["auc_crc_vs_enfermedad"] = origen.get("auc_crc_vs_enfermedad", "")
            fila["auc_enfermedad_vs_sano"] = origen.get("auc_enfermedad_vs_sano", "")
        else:
            fila["lista_ruido"] = ""
            fila["veredicto"] = ""
            fila["lectura"] = "pasa_fdr_en_algun_contraste"
        for nombre in ("crc_vs_sano", "crc_vs_enfermedad", "enfermedad_vs_sano"):
            trozo = largo[(largo["gene"] == gen) & (largo["contraste"] == nombre)]
            if trozo.empty:
                continue
            una = trozo.iloc[0]
            fila[f"log2fc_{nombre}"] = una["log2fc"]
            fila[f"padj_{nombre}"] = una["padj"]
            fila[f"n_obs_{nombre}"] = una["n_obs"]
            fila[f"significativa_{nombre}"] = una["significant"]
            if not fila["uniprot"]:
                fila["uniprot"] = uniprot_base(str(una["protein_id"])) or ""
        if hpa is not None and gen in hpa.index:
            fila["hpa_localizacion"] = hpa.at[gen, "localizacion_secretoma"] if "localizacion_secretoma" in hpa.columns else ""
            fila["hpa_ms_pg_l"] = hpa.at[gen, "conc_sangre_ms_pg_l"] if "conc_sangre_ms_pg_l" in hpa.columns else ""
            fila["hpa_inmuno_pg_l"] = hpa.at[gen, "conc_sangre_inmunoensayo_pg_l"] if "conc_sangre_inmunoensayo_pg_l" in hpa.columns else ""
        fila["go"] = ";".join(nombres_go.get(gen, []))
        acc = fila["uniprot"]
        fila["enlace_uniprot"] = f"https://www.uniprot.org/uniprotkb/{acc}" if acc else ""
        fila["enlace_alphafold"] = f"https://alphafold.ebi.ac.uk/entry/{acc}" if acc else ""
        filas.append(fila)
    if filas:
        write_tsv(pd.DataFrame(filas), salida / "fichas_proteinas.tsv")

