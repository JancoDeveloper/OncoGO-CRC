"""H8. Carga tumoral en suero de PXD064655, agregada con directLFQ."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, kruskal, spearmanr

from oncogo_crc.settings import Settings
from oncogo_crc.tesis.cierre import a_escala_log2
from oncogo_crc.utils.dia_nn_plasma import HUMANO_OX, acceso_base, leer_libreria
from oncogo_crc.v3.salida import volcar_json, volcar_tsv
from oncogo_crc.v3.suero_v32 import _contraste, _contraste_covariable, _hodges, _puntuacion_z

_ANOTACION = {
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
    "Modified.Sequence", "Precursor.Charge", "Precursor.Id",
}
_ORDINAL = {"1": 1, "2a": 2, "2b": 2, "3a": 3, "3b": 3, "4": 4}
_IV = {"4"}
_I_II = {"1", "2a", "2b"}
_HEMOLISIS = {"HBA1", "HBB", "CA1", "PRDX2"}
_CONTROLES = {"ALB", "APOA1", "C3", "TF", "SERPINA1", "A2M", "APOE", "FGA", "FGB", "FGG", "LRG1"}


def nombre_de_columna(columna: str) -> str:
    from pathlib import PurePosixPath, PureWindowsPath
    texto = str(columna)
    if "\\" in texto:
        return PureWindowsPath(texto).name
    return PurePosixPath(texto).name


def clasificar_corridas(columnas: list[str]) -> tuple[dict[str, str], list[str]]:
    """Une SH_NN.mzML y aparta los estándares ST_."""
    muestras: dict[str, str] = {}
    estandares: list[str] = []
    for col in columnas:
        if col in _ANOTACION:
            continue
        nombre = nombre_de_columna(col)
        hallado = re.fullmatch(r"(SH_\d+)\.mzML", nombre)
        if hallado:
            muestras[hallado.group(1)] = col
        elif re.fullmatch(r"ST_.+\.mzML", nombre):
            estandares.append(nombre)
    return muestras, estandares


def limpiar_estadio(meta: pd.DataFrame) -> pd.DataFrame:
    renombre = {c: str(c).replace("\r", "").strip() for c in meta.columns}
    tabla = meta.rename(columns=renombre).copy()
    for col in tabla.columns:
        tabla[col] = tabla[col].astype(str).str.replace("\r", "", regex=False).str.strip()
    id_col = next(c for c in tabla.columns if c.lower() in {"source name", "muestra", "sample"})
    stage = next(c for c in tabla.columns if re.search("stage", c, re.I))
    return pd.DataFrame({"muestra": tabla[id_col].tolist(), "estadio": tabla[stage].tolist()})


def unir_mapa(muestras: dict[str, str], meta: pd.DataFrame) -> dict:
    limpio = meta if {"muestra", "estadio"}.issubset(meta.columns) else limpiar_estadio(meta)
    ids = set(limpio["muestra"])
    claves = set(muestras)
    solo_matriz = sorted(claves - ids)
    solo_meta = sorted(ids - claves)
    return {
        "uno_a_uno": bool(claves) and not solo_matriz and not solo_meta,
        "n_unidas": len(claves & ids),
        "n_matriz": len(claves),
        "n_metadato": len(ids),
        "solo_matriz": solo_matriz,
        "solo_meta": solo_meta,
        "estadio": dict(zip(limpio["muestra"], limpio["estadio"], strict=True)),
    }


def fisher_carga(tabla) -> dict:
    arr = np.asarray(tabla, dtype=float)
    n_te = int(arr[0].sum())
    if n_te < 10 or (arr.sum(axis=1) == 0).any():
        return {
            "estado": "no_evaluable",
            "motivo": "menos de 10 tumorales epiteliales en el universo" if n_te < 10 else "una fila de la tabla está vacía",
            "odds_ratio": None,
            "p": None,
            "tabla": arr.astype(int).tolist(),
            "n_tumorales_epiteliales": n_te,
        }
    razon, p = fisher_exact(arr, alternative="greater")
    estado = "sostiene" if float(razon) > 1 and float(p) < 0.05 else "no_sostiene"
    return {
        "estado": estado,
        "motivo": "Fisher unilateral, alternative greater; sostiene si OR > 1 y p < 0,05",
        "odds_ratio": float(razon),
        "p": float(p),
        "tabla": arr.astype(int).tolist(),
        "n_tumorales_epiteliales": n_te,
    }


def _directlfq_acepta_arrays_de_solo_lectura() -> None:
    """Pandas 2.2+ entrega arrays de solo lectura. directLFQ 0.3.3 escribe encima."""
    import directlfq.config as config
    import directlfq.normalization as norm
    if getattr(norm.get_normfacts, "_oncogo", False):
        return
    original = norm.get_normfacts

    def con_copia(samples):
        return original(np.array(samples, dtype=float, copy=True))

    con_copia._oncogo = True
    norm.get_normfacts = con_copia
    config.check_wether_to_copy_numpy_arrays_derived_from_pandas()


def agregar_directlfq(entrada: Path) -> tuple[Path, str]:
    """Llama a directLFQ con los parámetros por defecto. Si rechaza NaN, reintenta con 0."""
    from directlfq.lfq_manager import run_lfq
    _directlfq_acepta_arrays_de_solo_lectura()
    usado = "nan"
    try:
        run_lfq(str(entrada))
    except Exception:
        tabla = pd.read_csv(entrada, sep="\t")
        tabla.to_csv(entrada, sep="\t", index=False, na_rep="0")
        usado = "cero"
        run_lfq(str(entrada))
    nuevos = sorted(
        (p for p in entrada.parent.glob("*protein_intensities.tsv") if not p.name.startswith("._")),
        key=lambda p: p.stat().st_mtime,
    )
    if not nuevos:
        raise RuntimeError("directLFQ no escribió la tabla de proteínas")
    return nuevos[-1], usado


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _simbolos(texto: object) -> set[str]:
    return {p.strip().upper() for p in str(texto).split(";") if p.strip() and p.strip().lower() != "nan"}


def _clasificar_grupo(fila: pd.Series, libreria) -> dict:
    nombres = str(fila.get("Protein.Names", ""))
    humano = "HUMAN" in nombres.upper() or any(a in libreria.accesos_humanos for a in _accesos(fila))
    accesos = _accesos(fila)
    simbolos = _simbolos(fila.get("Genes", ""))
    if humano:
        if accesos & set(libreria.accesos_humanos):
            motivo = "acceso_humano_en_libreria"
        elif not accesos and simbolos & set(libreria.genes_humanos):
            motivo = "simbolo_humano_sin_acceso"
        else:
            motivo = ""
    else:
        motivo = "acceso_no_humano_en_libreria" if accesos & set(libreria.accesos) else ""
    return {
        "humano": bool(humano),
        "contaminante": bool(motivo),
        "motivo": motivo,
        "genes": simbolos,
    }


def _accesos(fila: pd.Series) -> set[str]:
    texto = str(fila.get("Protein.Ids", fila.get("Protein.Group", "")))
    return {acceso_base(p) for p in texto.split(";") if p.strip() and p.strip().lower() != "nan"}


def _higado_secretado(settings: Settings) -> tuple[set[str], str]:
    import zipfile
    atlas = settings.resolve("external/hpa/proteinatlas.tsv.zip")
    sangre = settings.resolve("external/hpa/hpa_concentracion_sangre.tsv")
    if not atlas.is_file() or not sangre.is_file():
        return set(), f"falta {atlas.name if not atlas.is_file() else sangre.name}"
    with zipfile.ZipFile(atlas) as zf:
        tabla = pd.read_csv(zf.open(zf.namelist()[0]), sep="\t", usecols=["Gene", "RNA tissue specificity", "RNA tissue specific nTPM"], dtype=str)
    tabla["gen"] = tabla["Gene"].str.upper()
    especifico = tabla["RNA tissue specificity"].fillna("").str.lower().eq("tissue enriched")
    tejidos = tabla["RNA tissue specific nTPM"].fillna("").str.lower()
    en_higado = tejidos.map(lambda t: "liver" in {p.split(":")[0].strip() for p in t.split(";")})
    candidatos = set(tabla.loc[especifico & en_higado, "gen"])
    conc = pd.read_csv(sangre, sep="\t", dtype=str)
    secretadas = set(conc.loc[conc["localizacion_secretoma"].eq("Secreted to blood"), "gene"].str.upper())
    return candidatos & secretadas, "RNA tissue specificity; RNA tissue specific nTPM; localizacion_secretoma"


def _fase_aguda(settings: Settings) -> set[str]:
    go = settings.resolve("metadata/gene_ontology/go_gene_long.tsv")
    if not go.is_file():
        return set()
    anotacion = pd.read_csv(go, sep="\t", usecols=["gene_symbol", "go_id"], dtype=str)
    return set(anotacion.loc[anotacion["go_id"].eq("GO:0006953"), "gene_symbol"].str.upper())


def _gen_de_grupo(fila: pd.Series, uni: dict[str, str]) -> str | None:
    simbolos = _simbolos(fila.get("Genes", ""))
    if len(simbolos) == 1:
        return next(iter(simbolos))
    for acceso in _accesos(fila):
        if acceso in uni and uni[acceso]:
            return uni[acceso]
    return None


def _mapa_uniprot(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    tabla = pd.read_csv(path, sep="\t", usecols=["gene_symbol", "uniprot_accession"], dtype=str).dropna()
    salida: dict[str, set[str]] = {}
    for gen, acceso in zip(tabla["gene_symbol"].str.upper(), tabla["uniprot_accession"], strict=True):
        salida.setdefault(acceso_base(acceso), set()).add(gen)
    return {acceso: next(iter(genes)) for acceso, genes in salida.items() if len(genes) == 1}


def ejecutar(settings: Settings, salida: Path, sello: dict, prereg=None) -> dict:
    if prereg is None:
        from oncogo_crc.v3.prereg import PreregistroV3
        prereg = PreregistroV3(settings.root)
    import directlfq
    matriz_path = settings.resolve("external/pxd064655/report.pr_matrix.tsv")
    meta_path = settings.resolve("external/pxd064655/Sample_Metadata.txt")
    if not matriz_path.is_file() or not meta_path.is_file():
        payload = {"estado": "no_evaluable", "motivo": "no está la matriz o el metadato de PXD064655", "mapa_uno_a_uno": False}
        volcar_json(payload, salida / "h8_carga_tumoral.json", sello)
        return payload
    meta = pd.read_csv(meta_path, sep="\t", dtype=str)
    limpio = limpiar_estadio(meta)
    cabecera = pd.read_csv(matriz_path, sep="\t", nrows=0)
    muestras, estandares = clasificar_corridas(list(cabecera.columns))
    mapa = unir_mapa(muestras, limpio)
    tejido = next((c for c in meta.columns if re.search("organism part", str(c), re.I)), None)
    tejido_valor = str(meta[tejido].iloc[0]).strip() if tejido else ""
    base = {
        "mapa": mapa,
        "n_estandares_excluidos": len(estandares),
        "estandares": estandares,
        "discrepancia_tejido": {
            "metadato": tejido_valor,
            "protocolo_pride": "suero",
            "decision": "se anota y no se decide",
        },
        "pride_frente_a_matriz": "la matriz trae 100 corridas SH_; la revisión anota 98 crudos en PRIDE. No bloquea.",
        "directlfq_version": getattr(directlfq, "__version__", "desconocida"),
    }
    if not mapa["uno_a_uno"]:
        payload = base | {"estado": "no_evaluable", "motivo": "el mapa de muestras no es uno a uno", "mapa_uno_a_uno": False}
        volcar_json(payload, salida / "h8_carga_tumoral.json", sello)
        return payload
    libreria_path = settings.resolve(str(prereg.get("enmienda_v3_2.contaminantes.fichero")))
    if not libreria_path.is_file():
        payload = base | {"estado": "no_evaluable", "motivo": f"no está la librería de contaminantes {libreria_path}"}
        volcar_json(payload, salida / "h8_carga_tumoral.json", sello)
        return payload
    libreria = leer_libreria(libreria_path)
    columnas = ["Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "Precursor.Id", *muestras.values()]
    report = pd.read_csv(matriz_path, sep="\t", usecols=lambda c: c in columnas, low_memory=False)
    grupos = report.drop_duplicates("Protein.Group")
    clasificacion = {fila["Protein.Group"]: _clasificar_grupo(fila, libreria) for _, fila in grupos.iterrows()}
    marcados = {g for g, info in clasificacion.items() if info["genes"] & _CONTROLES and info["contaminante"]}
    if marcados:
        raise RuntimeError(f"un control positivo quedó como contaminante: {sorted(marcados)[:5]}")
    quedan = [g for g, info in clasificacion.items() if info["humano"] and not info["contaminante"]]
    report = report.loc[report["Protein.Group"].isin(quedan)].copy()
    renombre = {col: sid for sid, col in muestras.items()}
    bloque = report[["Protein.Group", "Precursor.Id", *muestras.values()]].rename(columns={"Protein.Group": "protein", "Precursor.Id": "ion", **renombre})
    derivados = settings.resolve("external/pxd064655/derivados")
    derivados.mkdir(parents=True, exist_ok=True)
    entrada = derivados / "pr_matrix_SH.aq_reformat.tsv"
    bloque.to_csv(entrada, sep="\t", index=False)
    proteina_path, convencion = agregar_directlfq(entrada)
    intensidades = pd.read_csv(proteina_path, sep="\t")
    id_col = "protein" if "protein" in intensidades.columns else intensidades.columns[0]
    intensidades = intensidades.set_index(id_col)
    columnas_sh = [c for c in intensidades.columns if re.fullmatch(r"SH_\d+", str(c))]
    valores = intensidades[columnas_sh].apply(pd.to_numeric, errors="coerce")
    valores = valores.replace(0, np.nan)
    log, escala = a_escala_log2(valores.to_numpy(dtype=float))
    mat = pd.DataFrame(log, index=valores.index.astype(str), columns=columnas_sh)
    n_proteinas_directlfq = int(len(mat))
    estadio = mapa["estadio"]
    iv = [c for c in mat.columns if estadio.get(c) in _IV]
    i_ii = [c for c in mat.columns if estadio.get(c) in _I_II]
    presencia = pd.DataFrame({
        "iv": mat[iv].notna().mean(axis=1),
        "i_ii": mat[i_ii].notna().mean(axis=1),
    })
    mat = mat.loc[(presencia["iv"] >= 0.70) & (presencia["i_ii"] >= 0.70)]
    n_presencia = int(len(mat))
    lfc = float(prereg.get("inferencia.treat_lfc"))
    alpha = float(prereg.get("inferencia.alpha"))
    contraste = _contraste(mat, iv, i_ii, lfc=lfc, alpha=alpha, nombre="pxd064655_iv_frente_a_i_ii")
    if "protein_id" in contraste.columns:
        contraste = contraste.set_index("protein_id")
    uni = _mapa_uniprot(settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv"))
    genes = {}
    for _, fila in grupos.iterrows():
        gen = _gen_de_grupo(fila, uni)
        if gen:
            genes[str(fila["Protein.Group"])] = gen
    from oncogo_crc.v3.genes_tanaka import gen_unico, mapa_genes
    tanaka_path = salida / "k1_primario_frente_a_mucosa.tsv"
    tanaka = pd.read_csv(tanaka_path, sep="\t")
    mapa_tejido = mapa_genes(
        settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt"),
        settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv"),
    )
    tanaka["gen"] = tanaka["protein_id"].map(lambda pid: gen_unico(pid, mapa_tejido))
    tanaka = tanaka.dropna(subset=["gen", "log2fc"])
    sube_tejido = set(tanaka.loc[tanaka["significant"].astype(bool) & (tanaka["log2fc"] > 0), "gen"].str.upper())
    con_log = set(tanaka["gen"].str.upper())
    h6_path = salida / "h6_origen_celular.tsv"
    h6 = pd.read_csv(h6_path, sep="\t")
    epitelial = set(h6.loc[h6["clase"].eq("epitelial"), "gen"].str.upper())
    aguda = _fase_aguda(settings)
    higado, columnas_hpa = _higado_secretado(settings)
    eritrocitos = {g.upper() for g in prereg.get("enmienda_v3_2.paneles_geyer_2019.eritrocitos")}
    ruido = {
        "fase_aguda": aguda,
        "higado_secretado": higado,
        "hemolisis": _HEMOLISIS,
        "eritrocitos_geyer": eritrocitos,
    }
    if not aguda or (not higado and columnas_hpa.startswith("falta")):
        payload = base | {
            "estado": "no_evaluable",
            "motivo": "falta una lista de ruido (GO:0006953 o el HPA de hígado y secreción)",
            "listas_de_ruido": {k: len(v) for k, v in ruido.items()},
        }
        volcar_json(payload, salida / "h8_carga_tumoral.json", sello)
        return payload
    filas = []
    for proteina, gen in genes.items():
        if proteina not in contraste.index:
            continue
        listas = [nombre for nombre, miembros in ruido.items() if gen in miembros]
        fila = contraste.loc[proteina]
        filas.append({
            "protein_id": proteina,
            "gen": gen,
            "log2fc": float(fila["log2fc"]) if np.isfinite(fila["log2fc"]) else None,
            "padj": float(fila["padj"]) if np.isfinite(fila["padj"]) else None,
            "clase": fila["clase"],
            "tumoral_epitelial": bool(gen in sube_tejido and gen in epitelial and gen in con_log and not listas),
            "en_universo": bool(gen in con_log),
            "lista_de_ruido": ";".join(listas),
        })
    tabla = pd.DataFrame(filas)
    if tabla.empty or "en_universo" not in tabla.columns:
        payload = base | {
            "estado": "no_evaluable",
            "motivo": "ninguna proteína cuantificada mapea a un gen único del contraste",
            "n_proteinas_directlfq": int(mat.shape[0]),
            "n_tras_presencia": int(len(contraste)),
        }
        volcar_json(payload, salida / "h8_carga_tumoral.json", sello)
        return payload
    universo = tabla.loc[tabla["en_universo"] & tabla["lista_de_ruido"].eq("")]
    te = universo["tumoral_epitelial"]
    sube = universo["clase"].eq("diferente") & (universo["log2fc"] > 0)
    cuenta = [
        [int((te & sube).sum()), int((te & ~sube).sum())],
        [int((~te & sube).sum()), int((~te & ~sube).sum())],
    ]
    prueba = fisher_carga(cuenta)
    ordinal = {m: _ORDINAL[e] for m, e in estadio.items() if e in _ORDINAL and m in mat.columns}
    rho = []
    if len(ordinal) >= 10 and len(tabla):
        orden = pd.Series(ordinal)
        for proteina in tabla["protein_id"]:
            if proteina not in mat.index:
                rho.append((None, None))
                continue
            comun = orden.index.intersection(mat.columns)
            x = orden.loc[comun].to_numpy(dtype=float)
            y = mat.loc[proteina, comun].to_numpy(dtype=float)
            ok = np.isfinite(y)
            if ok.sum() < 10:
                rho.append((None, None))
                continue
            r, p = spearmanr(x[ok], y[ok])
            rho.append((float(r), float(p)))
        tabla["rho_estadio"] = [a[0] for a in rho]
        tabla["p_estadio"] = [a[1] for a in rho]
    paneles = {
        "eritrocitos": list(prereg.get("enmienda_v3_2.paneles_geyer_2019.eritrocitos")),
        "plaquetas": list(prereg.get("enmienda_v3_2.paneles_geyer_2019.plaquetas")),
        "coagulacion": list(prereg.get("enmienda_v3_2.paneles_geyer_2019.coagulacion")),
    }
    piel = list(prereg.get("enmienda_v3_2.contaminante_piel.genes"))
    gen_por_proteina = pd.Series(genes)
    puntuaciones = {}
    for nombre, lista in {**paneles, "piel": piel}.items():
        miembros = gen_por_proteina[gen_por_proteina.isin({g.upper() for g in lista})].index.tolist()
        presentes = [m for m in miembros if m in mat.index]
        puntuaciones[nombre] = _puntuacion_z(mat, presentes)
    qc = pd.DataFrame({"muestra": list(mat.columns)})
    qc["estadio"] = qc["muestra"].map(estadio)
    qc["grupo"] = qc["estadio"].map(lambda e: "IV" if e in _IV else "I_II" if e in _I_II else "III" if e in {"3a", "3b"} else "otro")
    for nombre, serie in puntuaciones.items():
        qc[nombre] = serie.reindex(qc["muestra"]).to_numpy()
    volcar_tsv(qc, salida / "h8_qc_muestras.tsv", sello)
    comparacion = {}
    for nombre in puntuaciones:
        a = qc.loc[qc["grupo"].eq("IV"), nombre].dropna().to_numpy()
        b = qc.loc[qc["grupo"].eq("I_II"), nombre].dropna().to_numpy()
        p = float(kruskal(a, b).pvalue) if a.size >= 3 and b.size >= 3 else None
        comparacion[nombre] = {"p_kruskal": p, "hodges_iv_menos_i_ii": _hodges(a, b)["estimador"]}
    covariable = {"estado": "no_disponible", "motivo": "la puntuación de coagulación no tiene valores en los dos grupos"}
    coag = puntuaciones["coagulacion"]
    if coag.reindex(iv + i_ii).notna().sum() >= 8:
        ajustado = _contraste_covariable(mat, iv, i_ii, coag, lfc=lfc, alpha=alpha, nombre="pxd064655_coagulacion")
        covariable = {
            "estado": "calculado",
            "n": int(len(ajustado)),
            "n_diferente": int((ajustado["clase"] == "diferente").sum()) if "clase" in ajustado.columns else None,
            "decide": False,
        }
    volcar_tsv(tabla.drop(columns=["en_universo"]), salida / "h8_carga_tumoral.tsv", sello)
    payload = base | {
        "estado": prueba["estado"],
        "motivo": prueba["motivo"],
        "mapa_uno_a_uno": True,
        "n_por_grupo": {"IV": len(iv), "I_II": len(i_ii), "III": int(qc["grupo"].eq("III").sum())},
        "tabla": prueba["tabla"],
        "odds_ratio": prueba["odds_ratio"],
        "p": prueba["p"],
        "n_tumorales_epiteliales": prueba["n_tumorales_epiteliales"],
        "n_universo": int(len(universo)),
        "escala": escala,
        "convencion_ausentes": convencion,
        "ceros_de_salida_son_ausentes": True,
        "sha256_entrada": _sha(entrada),
        "sha256_salida": _sha(proteina_path),
        "fichero_proteinas": proteina_path.name,
        "columnas_hpa": columnas_hpa,
        "listas_de_ruido": {k: len(v) for k, v in ruido.items()},
        "qc": comparacion,
        "sensibilidad_coagulacion": covariable,
        "sensibilidad_ordinal": "rho_estadio en h8_carga_tumoral.tsv; no decide",
        "ox_humano": HUMANO_OX,
        "embudo_universo": {
            "n_proteinas_directlfq": n_proteinas_directlfq,
            "n_presencia_iv_y_i_ii": n_presencia,
            "n_gen_unico": int(tabla["gen"].nunique()),
            "n_con_log2fc_tanaka": int(tabla["en_universo"].sum()),
            "n_suben_en_tanaka": int(tabla["gen"].str.upper().isin(sube_tejido).sum()),
            "n_epiteliales": int(tabla["gen"].str.upper().isin(epitelial).sum()),
            "n_fuera_de_listas_de_ruido": int(tabla["lista_de_ruido"].eq("").sum()),
            "n_tumorales_epiteliales": int(tabla["tumoral_epitelial"].sum()),
        },
        "motivo_universo": (
            f"De {int(tabla['gen'].nunique())} proteínas con gen único, "
            f"{int(tabla['en_universo'].sum())} tienen log2FC en toda la tabla de Tanaka. "
            f"Tras las listas de ruido el universo de la prueba tiene {int(len(universo))}. "
            "El código no recorta a las proteínas concordantes ni a las 688 certificadas."
        ),
        "lectura": "Un sostiene diría que el estadio IV lleva más proteína epitelial del tumor a la sangre. Un no_sostiene con pocas proteínas tumorales es, sobre todo, un límite de profundidad de un DIA sin depleción.",
    }
    volcar_json(payload, salida / "h8_carga_tumoral.json", sello)
    return payload
