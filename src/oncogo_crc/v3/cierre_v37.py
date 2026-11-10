"""Piezas del cierre v3.7. Ninguna decide una hipótesis ni copia una cifra a mano."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import beta, spearmanr

from oncogo_crc.utils.sello import escribir_sello
from oncogo_crc.v3.salida import volcar_json, volcar_tsv

PERMITIDOS_AL_EMPEZAR = {"crudos_integridad.tsv", "crudos_integridad.tsv.sello.json"}
CAMPOS_SELLO = (
    "commit", "arbol_limpio", "sha256_preregistro", "sha256_enmiendas",
    "sha256_manifest_entradas", "etapa",
)
OBJETIVOS_RAI = (("BAIAP3", "O94812"), ("CEACAM5", "P06731"))
ENTRADAS = (
    "data/raw/plasma/DIA-NN_report.tsv",
    "results/certified/cptac_unshared_treat.tsv",
    "results/certified/inference_audit.json",
    "results/certification/certificado_resumen.json",
    "metadata/gene_ontology/go_gene_long.tsv",
    "metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv",
    "data/raw/pride/Cohort1_proteinGroups.txt",
    "external/hpa/hpa_concentracion_sangre.tsv",
    "config/preregistro_v3.yaml",
    "config/config.yaml",
    "config/entradas_v3_7.yaml",
    "config/esperado_v3_7.tsv",
)
PASOS = (
    "1. diferencial certificada en CPTAC (TREAT)",
    "2. medida en Tanaka P frente a N y con el mismo signo",
    "3. clase H5 (se informa; no elimina)",
    "4. HPA: secretada a sangre o con concentración plasmática conocida",
    "5. medida en la matriz principal del suero PXD055218",
    "6. K4 CP frente a CN: clase diferente con el mismo signo que el tejido",
    "7. H10: especifica_crc en el HPA Disease Blood Atlas",
)


def exigir_carpeta_inicial(salida: Path) -> None:
    if not salida.exists():
        return
    restos = [
        path.name for path in salida.iterdir()
        if not path.name.startswith(".") and path.name not in PERMITIDOS_AL_EMPEZAR
    ]
    if restos:
        raise RuntimeError(f"results/v3_7 tiene restos y la etapa no empieza: {restos[:8]}")


def bh(pares: list[tuple[str, float]]) -> dict[str, float]:
    ordenado = sorted(pares, key=lambda item: item[1])
    n = len(ordenado)
    q: dict[str, float] = {}
    previo = 1.0
    for rango in range(n, 0, -1):
        nombre, p = ordenado[rango - 1]
        previo = min(previo, p * n / rango)
        q[nombre] = float(min(1.0, previo))
    return q


def _sha256(path: Path) -> str:
    resumen = hashlib.sha256()
    with path.open("rb") as fichero:
        for bloque in iter(lambda: fichero.read(1 << 20), b""):
            resumen.update(bloque)
    return resumen.hexdigest()


def sellos_coinciden(salida: Path) -> dict:
    sellos = sorted(path for path in salida.glob("*.sello.json") if not path.name.startswith("._"))
    if not sellos:
        return {"estado": "falla", "motivo": "no hay sellos", "n": 0}
    claves = []
    for path in sellos:
        datos = json.loads(path.read_text(encoding="utf-8"))
        claves.append(tuple(datos.get(campo) for campo in CAMPOS_SELLO))
    if len(set(claves)) != 1 or claves[0][4] in (None, "") or claves[0][5] != "v3_7":
        return {"estado": "falla", "n": len(sellos), "motivo": "commit, manifiesto o etapa distintos"}
    return {"estado": "pasa", "n": len(sellos)}


def unificar_sellos(salida: Path, sello: dict) -> None:
    for path in salida.glob("*.sello.json"):
        if path.name.startswith("._"):
            continue
        escribir_sello(path, sello)
    escribir_sello(salida / "sello_corrida.json", sello)


def emitir(casillas: list[dict], prohibidas: list[str]) -> dict:
    salida = []
    alguna_falla = False
    for casilla in casillas:
        estado = casilla["estado"]
        ruta = casilla.get("ruta")
        sha = None
        if ruta is not None:
            path = Path(ruta)
            if not path.is_file() or path.stat().st_size == 0:
                estado = "falla"
            else:
                sha = _sha256(path)
        if estado == "falla":
            alguna_falla = True
        salida.append({
            "nombre": casilla["nombre"],
            "fuente": casilla.get("fuente"),
            "criterio": casilla.get("criterio"),
            "valor": casilla.get("valor"),
            "estado": estado,
            "sha256_fuente": sha,
        })
    estado_global = "no_publicable" if prohibidas or alguna_falla else "publicable"
    return {"estado": estado_global, "casillas": salida, "prohibidas": list(prohibidas)}


def _base(acceso: str) -> str:
    return acceso.split("-")[0].strip()


def _fila_rai(traza: pd.DataFrame, gen: str, acceso: str):
    if traza.empty or "protein_group" not in traza.columns:
        return None
    for fila in traza.itertuples(index=False):
        genes = {parte.strip().upper() for parte in str(getattr(fila, "genes", "")).split(";") if parte.strip()}
        accesos = {_base(parte) for parte in str(fila.protein_group).split(";") if parte.strip()}
        if gen in genes or acceso in accesos:
            return fila
    return None


def auditar_rai(traza: pd.DataFrame, matriz: pd.DataFrame, celda: pd.DataFrame) -> pd.DataFrame:
    """Recorre 897, el filtro global, la matriz y la celda. Los números salen de las tablas."""
    en_matriz = set(matriz["protein_id"].astype(str)) if "protein_id" in matriz.columns else set()
    en_celda = set(celda["protein_id"].astype(str).str.upper()) if "protein_id" in celda.columns else set()
    filas = []
    for gen, acceso in OBJETIVOS_RAI:
        fila = _fila_rai(traza, gen, acceso)
        if fila is None:
            filas.append({"gen": gen, "acceso": acceso, "paso": "no_esta_en_el_report",
                          "motivo": "no aparece entre los grupos del report",
                          "n_precursores": None, "fraccion_proteotipicos": None, "global_pg_q": None})
            continue
        global_q = float(fila.global_pg_q)
        pg = float(fila.pg_q_min)
        if global_q > 0.01:
            paso, motivo = "filtro_global_1pct", "global_pg_q por encima de 0.01"
        elif not bool(fila.humano) or pg > 0.01:
            paso, motivo = "humanos_pg_y_global", "no pasa humano con PG.Q.Value <= 0.01"
        elif acceso not in en_matriz and not any(acceso == _base(item) or item.startswith(acceso) for item in en_matriz):
            paso, motivo = "matriz_principal", "no está en la matriz principal"
        elif gen not in en_celda:
            paso, motivo = "celda", "no llega a la celda preregistrada"
        else:
            paso, motivo = "llega_a_la_celda", "está en la celda preregistrada"
        filas.append({
            "gen": gen, "acceso": acceso, "paso": paso, "motivo": motivo,
            "n_precursores": int(fila.n_precursores),
            "fraccion_proteotipicos": float(fila.fraccion_proteotipicos),
            "global_pg_q": global_q,
        })
    return pd.DataFrame(filas)


def filtrar_terminos(terminos: dict, universo: set[str], minimo: int, maximo: int) -> dict[str, set[str]]:
    universo = {str(gen).upper() for gen in universo}
    filtrados = {}
    for termino, genes in terminos.items():
        dentro = {str(gen).upper() for gen in genes} & universo
        if minimo <= len(dentro) <= maximo:
            filtrados[termino] = dentro
    return filtrados


def prerank_gsea(**kwargs):
    import gseapy
    return gseapy.prerank(**kwargs)


def gsea_prerank(ranking: pd.Series, terminos: dict, semilla: int, permutaciones: int, minimo: int, maximo: int):
    return prerank_gsea(
        rnk=ranking.sort_values(ascending=False),
        gene_sets={termino: sorted(genes) for termino, genes in terminos.items()},
        permutation_num=permutaciones,
        min_size=minimo,
        max_size=maximo,
        seed=semilla,
        outdir=None,
        no_plot=True,
        threads=1,
        verbose=False,
    )


def shap_de_un_ajuste(x: np.ndarray, y: np.ndarray, semilla: int = 0):
    """SHAP lineal exacto: coeficiente por (x escalada − media). Suma el cambio del logit."""
    from oncogo_crc.v3.ml_suero import _tubo

    tubo = _tubo(semilla, "proteinas").fit(x, y)
    escalado = tubo[:-1].transform(x)
    media = escalado.mean(axis=0)
    coef = tubo.named_steps["modelo"].coef_.ravel()
    valores = coef * (escalado - media)
    modelo = tubo.named_steps["modelo"]
    diferencia = modelo.decision_function(escalado) - float(modelo.decision_function(media.reshape(1, -1))[0])
    return valores, diferencia


def comprobar_con_shap(x: np.ndarray, y: np.ndarray, semilla: int = 0) -> float:
    import shap
    from oncogo_crc.v3.ml_suero import _tubo

    tubo = _tubo(semilla, "proteinas").fit(x, y)
    escalado = np.asarray(tubo[:-1].transform(x), dtype=float)
    modelo = tubo.named_steps["modelo"]
    media = escalado.mean(axis=0)
    exacto = modelo.coef_.ravel() * (escalado - media)
    explicado = shap.LinearExplainer(modelo, escalado).shap_values(escalado)
    if isinstance(explicado, list):
        explicado = explicado[-1]
    return float(np.max(np.abs(np.asarray(explicado, dtype=float) - exacto)))


def _parametros_gsea(root: Path) -> dict:
    cfg = yaml.safe_load((root / "config" / "config.yaml").read_text(encoding="utf-8"))
    gsea = cfg["gsea"]
    return {
        "fdr": float(gsea["fdr_threshold"]),
        "minimo": int(gsea["min_size"]),
        "maximo": int(gsea["max_size"]),
        "permutaciones": int(gsea["permutations"]),
        "semilla": int(gsea["seed"]),
        "ranking": "signed_neglog10_p",
    }


def _clopper(aciertos: int, n: int) -> tuple[float | None, float | None]:
    if n <= 0:
        return None, None
    bajo = 0.0 if aciertos <= 0 else float(beta.ppf(0.025, aciertos, n - aciertos + 1))
    alto = 1.0 if aciertos >= n else float(beta.ppf(0.975, aciertos + 1, n - aciertos))
    return bajo, alto


def _ranking(tabla: pd.DataFrame) -> pd.DataFrame:
    trabajo = tabla.copy()
    trabajo["gen"] = trabajo["gen"].astype(str).str.upper()
    trabajo = trabajo.dropna(subset=["gen", "log2fc", "p_value"]).drop_duplicates("gen")
    p = pd.to_numeric(trabajo["p_value"], errors="coerce").clip(lower=1e-300, upper=1)
    signo = np.sign(pd.to_numeric(trabajo["log2fc"], errors="coerce")).replace(0, 1)
    trabajo["ranking"] = signo * -np.log10(p)
    return trabajo


def _terminos_bp(path: Path) -> dict[str, set[str]]:
    tabla = pd.read_csv(path, sep="\t", usecols=["gene_symbol", "go_id", "namespace"])
    bp = tabla.loc[tabla["namespace"].eq("BP")]
    return {
        str(go_id): set(genes.astype(str).str.upper())
        for go_id, genes in bp.groupby("go_id")["gene_symbol"]
    }


def ejecutar_go(settings, salida: Path, sello: dict) -> dict:
    from oncogo_crc.utils.go_utils import ora_fisher
    from oncogo_crc.v3.genes_tanaka import gen_unico, mapa_genes

    params = _parametros_gsea(settings.root)
    try:
        import gseapy
        version = str(gseapy.__version__)
    except ImportError:
        version = "no_instalado"
    terminos = _terminos_bp(settings.resolve("metadata/gene_ontology/go_gene_long.tsv"))
    cohortes = {}
    nes: dict[str, dict[str, float]] = {}
    filas_ora = []
    filas_gsea = []
    fuentes = {
        "cptac_tumor_frente_a_normal": settings.resolve("results/certified/cptac_unshared_treat.tsv"),
        "tanaka_primario_frente_a_mucosa": salida / "k1_primario_frente_a_mucosa.tsv",
    }
    mapa = None
    for nombre, ruta in fuentes.items():
        if not ruta.is_file():
            cohortes[nombre] = {"estado": "no_evaluable", "motivo": f"falta {ruta.name}"}
            continue
        tabla = pd.read_csv(ruta, sep="\t")
        if nombre.startswith("tanaka"):
            if mapa is None:
                mapa = mapa_genes(
                    settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt"),
                    settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv"),
                )
            tabla["gen"] = tabla["protein_id"].map(lambda pid: gen_unico(pid, mapa))
        else:
            tabla["gen"] = tabla["protein_id"].astype(str)
        trabajo = _ranking(tabla)
        universo = set(trabajo["gen"])
        terminos_cohorte = filtrar_terminos(terminos, universo, params["minimo"], params["maximo"])
        listas = {termino: sorted(genes) for termino, genes in terminos_cohorte.items()}
        significativa = trabajo["significant"].map(lambda valor: str(valor).lower() in {"true", "1"})
        suben = set(trabajo.loc[significativa & (pd.to_numeric(trabajo["log2fc"]) > 0), "gen"])
        bajan = set(trabajo.loc[significativa & (pd.to_numeric(trabajo["log2fc"]) < 0), "gen"])
        for direccion, genes in (("suben", suben), ("bajan", bajan)):
            ora = ora_fisher(genes, universo, terminos_cohorte, fdr_threshold=params["fdr"], min_genes_in_term=params["minimo"])
            ora.insert(0, "cohorte", nombre)
            ora.insert(1, "direccion", direccion)
            filas_ora.append(ora)
        gsea_estado = "calculado"
        gsea_error = None
        n_gsea = 0
        if version == "no_instalado":
            gsea_estado, gsea_error = "no_evaluable", "gseapy no está instalado"
        else:
            try:
                resultado = gsea_prerank(
                    trabajo.set_index("gen")["ranking"], listas,
                    params["semilla"], params["permutaciones"], params["minimo"], params["maximo"],
                )
                res = resultado.res2d.copy()
                res.insert(0, "cohorte", nombre)
                filas_gsea.append(res)
                col_fdr = "FDR q-val" if "FDR q-val" in res.columns else "fdr"
                col_nes = "NES" if "NES" in res.columns else "nes"
                col_term = "Term" if "Term" in res.columns else "term"
                if col_fdr in res.columns:
                    n_gsea = int((pd.to_numeric(res[col_fdr], errors="coerce") <= params["fdr"]).sum())
                    nes[nombre] = {
                        str(fila[col_term]): float(fila[col_nes])
                        for _, fila in res.iterrows()
                        if pd.notna(fila.get(col_nes))
                    }
            except Exception as exc:
                gsea_estado, gsea_error = "no_evaluable", str(exc)
        ora_cat = pd.concat(filas_ora, ignore_index=True) if filas_ora else pd.DataFrame()
        de_esta = ora_cat.loc[ora_cat["cohorte"].eq(nombre)] if len(ora_cat) and "cohorte" in ora_cat.columns else ora_cat
        def _n_sig(direccion: str) -> int:
            if de_esta.empty or "significant" not in de_esta.columns:
                return 0
            return int(((de_esta["direccion"] == direccion) & de_esta["significant"]).sum())
        cohortes[nombre] = {
            "estado": "calculado",
            "n_genes": int(len(universo)),
            "n_terminos_probados": int(len(listas)),
            "n_ora_significativos_suben": _n_sig("suben"),
            "n_ora_significativos_bajan": _n_sig("bajan"),
            "n_gsea_fdr": n_gsea,
            "gsea_estado": gsea_estado,
            "gsea_error": gsea_error,
        }
    concordancia = {"numerador": None, "denominador": None, "ci_low": None, "ci_high": None,
                    "nota": "descriptiva, sin veredicto"}
    if "cptac_tumor_frente_a_normal" in nes and "tanaka_primario_frente_a_mucosa" in nes:
        gsea_cptac = pd.concat(filas_gsea, ignore_index=True)
        col_fdr = "FDR q-val" if "FDR q-val" in gsea_cptac.columns else "fdr"
        col_term = "Term" if "Term" in gsea_cptac.columns else "term"
        bajos = gsea_cptac.loc[
            gsea_cptac["cohorte"].eq("cptac_tumor_frente_a_normal")
            & (pd.to_numeric(gsea_cptac[col_fdr], errors="coerce") <= params["fdr"]),
            col_term,
        ].astype(str)
        mismos = 0
        n = 0
        for termino in bajos:
            a = nes["cptac_tumor_frente_a_normal"].get(termino)
            b = nes["tanaka_primario_frente_a_mucosa"].get(termino)
            if a is None or b is None or a == 0 or b == 0:
                continue
            n += 1
            mismos += int(np.sign(a) == np.sign(b))
        bajo, alto = _clopper(mismos, n)
        concordancia = {"numerador": mismos, "denominador": n, "ci_low": bajo, "ci_high": alto,
                        "nota": "descriptiva, sin veredicto"}
    payload = {
        "veredicto_h1": "no_emite",
        "gseapy": version,
        "parametros": params,
        "cohortes": cohortes,
        "concordancia_nes": concordancia,
        "limitaciones": (
            f"GO sin propagar, sin IEA, {params['permutaciones']} permutaciones, "
            f"p mínimo 1/{params['permutaciones'] + 1}"
        ),
    }
    if filas_ora:
        volcar_tsv(pd.concat(filas_ora, ignore_index=True), salida / "go_ora_tejido.tsv", sello)
    else:
        volcar_tsv(pd.DataFrame(columns=["cohorte", "direccion", "go_id"]), salida / "go_ora_tejido.tsv", sello)
    if filas_gsea:
        volcar_tsv(pd.concat(filas_gsea, ignore_index=True), salida / "gsea_tejido.tsv", sello)
    else:
        volcar_tsv(pd.DataFrame({"estado": ["no_evaluable"]}), salida / "gsea_tejido.tsv", sello)
    volcar_json(payload, salida / "go_tejido.json", sello)
    return payload


def _shap_de_pliegues(x: np.ndarray, y: np.ndarray, semilla: int, tipo: str) -> tuple[np.ndarray, float | None]:
    from oncogo_crc.v3.ml_suero import N_REPETICIONES, N_SPLITS, _tubo
    from sklearn.model_selection import RepeatedStratifiedKFold

    suma = np.zeros(x.shape[1], dtype=float)
    cuenta = np.zeros(x.shape[1], dtype=float)
    partes = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPETICIONES, random_state=semilla)
    comprobacion = None
    for entreno, prueba in partes.split(x, y):
        tubo = _tubo(semilla, tipo).fit(x[entreno], y[entreno])
        escalado_te = np.asarray(tubo[:-1].transform(x[prueba]), dtype=float)
        escalado_tr = np.asarray(tubo[:-1].transform(x[entreno]), dtype=float)
        media = escalado_tr.mean(axis=0)
        coef = tubo.named_steps["modelo"].coef_.ravel()
        valores = np.abs(coef * (escalado_te - media))
        suma += valores.sum(axis=0)
        cuenta += len(prueba)
        if comprobacion is None:
            try:
                import shap
                exacto = coef * (escalado_tr - media)
                explicado = shap.LinearExplainer(tubo.named_steps["modelo"], escalado_tr).shap_values(escalado_tr)
                if isinstance(explicado, list):
                    explicado = explicado[-1]
                comprobacion = float(np.max(np.abs(np.asarray(explicado, dtype=float) - exacto)))
            except ImportError:
                comprobacion = None
            except Exception as exc:
                comprobacion = str(exc)
    return suma / np.maximum(cuenta, 1), comprobacion


def ejecutar_explicacion(settings, salida: Path, sello: dict, prereg) -> dict:
    from oncogo_crc.v3.ml_suero import matrices_de_suero
    from oncogo_crc.v3.suero_v32 import semilla_v33

    principal = pd.read_csv(salida / "suero_matriz_principal.tsv", sep="\t")
    bloques = matrices_de_suero(
        principal,
        pd.read_csv(salida / "suero_matriz_humana_fdr.tsv", sep="\t"),
        pd.read_csv(salida / "suero_calidad_muestra.tsv", sep="\t"),
        prereg,
    )
    nombres = {
        "M2_proteinas_principal": principal["protein_id"].astype(str).tolist(),
        "M3_solo_qc": ["eritrocitos", "plaquetas", "coagulacion", "piel", "log2_n_precursores"],
    }
    traza = pd.read_csv(salida / "suero_trazabilidad_grupos.tsv", sep="\t")
    orden = pd.read_csv(salida / "preanalitica_orden.tsv", sep="\t")
    orden["sample_id"] = orden["crudo"].astype(str).str.removesuffix(".raw")
    calidad = pd.read_csv(salida / "suero_calidad_muestra.tsv", sep="\t")
    filas = []
    comprobaciones = {}
    semilla = semilla_v33(prereg, "ml_cv")
    for modelo in ("M2_proteinas_principal", "M3_solo_qc"):
        x, y, tipo = bloques[modelo]
        media, comprobacion = _shap_de_pliegues(x, y, semilla, tipo)
        comprobaciones[modelo] = comprobacion
        orden_idx = np.argsort(-media)[:20]
        for indice in orden_idx:
            filas.append({"modelo": modelo, "variable": nombres[modelo][indice], "media_abs_shap": float(media[indice])})
    tabla = pd.DataFrame(filas)
    tabla = _anotar_shap(tabla, principal, traza, orden, calidad)
    volcar_tsv(tabla, salida / "ml_explicacion.tsv", sello)
    de_m2 = tabla.loc[tabla["modelo"].eq("M2_proteinas_principal")]
    q_columnas = [col for col in de_m2.columns if col.startswith("q_")]
    correlacionan = 0
    if q_columnas and len(de_m2):
        correlacionan = int(de_m2[q_columnas].apply(pd.to_numeric, errors="coerce").le(0.05).any(axis=1).sum())
    payload = {
        "comprobacion_shap": comprobaciones,
        "n_piel_en_las_20": int(de_m2["contaminante_piel"].map(lambda v: str(v).lower() == "true").sum()) if len(de_m2) else 0,
        "n_coagulacion_en_las_20": int(de_m2["panel"].astype(str).str.contains("coagulacion").sum()) if len(de_m2) else 0,
        "n_correlacionan_con_calidad_bh": correlacionan,
        "lectura": "M3 usa solo la calidad. Si separa los grupos, la separación de M2 no puede atribuirse a la biología del tumor.",
    }
    volcar_json(payload, salida / "ml_explicacion.json", sello)
    return payload


def _anotar_shap(tabla: pd.DataFrame, principal: pd.DataFrame, traza: pd.DataFrame, orden: pd.DataFrame, calidad: pd.DataFrame) -> pd.DataFrame:
    anotacion = {}
    for fila in traza.itertuples(index=False):
        for acceso in str(fila.protein_group).split(";"):
            base = _base(acceso)
            if base:
                anotacion[base] = fila
    muestras = [c for c in principal.columns if str(c).startswith(("CP_", "CN_"))]
    matriz = principal.assign(protein_id=principal["protein_id"].astype(str)).set_index("protein_id")
    calidad = calidad.assign(sample_id=calidad["sample_id"].astype(str)).set_index("sample_id")
    orden = orden.drop_duplicates("sample_id").set_index("sample_id")
    comun = [m for m in muestras if m in calidad.index and m in orden.index]
    pares_p = []
    extra = []
    for fila in tabla.itertuples(index=False):
        info = anotacion.get(_base(str(fila.variable)))
        registro = {
            "gen": str(getattr(info, "genes", fila.variable)) if info is not None else str(fila.variable),
            "contaminante": bool(getattr(info, "contaminante", False)) if info is not None else False,
            "contaminante_piel": bool(getattr(info, "contaminante_piel", False)) if info is not None else False,
            "panel": str(getattr(info, "panel", "")) if info is not None else "",
        }
        if fila.variable in matriz.index and comun:
            valores = pd.to_numeric(matriz.loc[fila.variable, comun], errors="coerce")
            if isinstance(valores, pd.DataFrame):
                valores = valores.iloc[0]
            valores.index = comun
            for metrica in ("eritrocitos", "plaquetas", "coagulacion", "piel", "n_precursores", "orden"):
                if metrica == "orden":
                    referencia = pd.to_numeric(orden.loc[comun, "orden"], errors="coerce")
                elif metrica == "n_precursores":
                    referencia = np.log2(pd.to_numeric(calidad.loc[comun, "n_precursores"], errors="coerce"))
                else:
                    referencia = pd.to_numeric(calidad.loc[comun, metrica], errors="coerce")
                referencia.index = comun
                ok = valores.notna().to_numpy() & referencia.notna().to_numpy()
                clave = f"{fila.modelo}:{fila.variable}:{metrica}"
                if int(ok.sum()) < 5:
                    registro[f"rho_{metrica}"] = None
                    registro[f"p_{metrica}"] = None
                else:
                    rho, p = spearmanr(valores.to_numpy()[ok], referencia.to_numpy()[ok])
                    registro[f"rho_{metrica}"] = float(rho)
                    registro[f"p_{metrica}"] = float(p)
                    pares_p.append((clave, float(p)))
                registro[f"_clave_{metrica}"] = clave
        extra.append(registro)
    q = bh(pares_p) if pares_p else {}
    for registro in extra:
        for metrica in ("eritrocitos", "plaquetas", "coagulacion", "piel", "n_precursores", "orden"):
            clave = registro.pop(f"_clave_{metrica}", None)
            registro[f"q_{metrica}"] = q.get(clave) if clave else None
    return pd.concat([tabla.reset_index(drop=True), pd.DataFrame(extra)], axis=1)


def _signo(valor) -> int:
    numero = pd.to_numeric(pd.Series([valor]), errors="coerce").iloc[0]
    if pd.isna(numero) or numero == 0:
        return 0
    return int(np.sign(numero))


def ejecutar_priorizacion(settings, salida: Path, sello: dict) -> dict:
    from oncogo_crc.v3.genes_tanaka import gen_unico, mapa_genes

    cptac = pd.read_csv(settings.resolve("results/certified/cptac_unshared_treat.tsv"), sep="\t")
    cptac["gen"] = cptac["protein_id"].astype(str).str.upper()
    tanaka = pd.read_csv(salida / "k1_primario_frente_a_mucosa.tsv", sep="\t")
    mapa = mapa_genes(
        settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt"),
        settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv"),
    )
    tanaka["gen"] = tanaka["protein_id"].map(lambda pid: gen_unico(pid, mapa))
    tanaka = tanaka.dropna(subset=["gen"]).drop_duplicates("gen")
    tanaka_fc = dict(zip(tanaka["gen"].str.upper(), tanaka["log2fc"], strict=False))
    h5 = pd.read_csv(salida / "h5_especificidad_organo.tsv", sep="\t")
    clase_h5 = dict(zip(h5["gen"].astype(str).str.upper(), h5["clase"], strict=False))
    hpa = pd.read_csv(settings.resolve("external/hpa/hpa_concentracion_sangre.tsv"), sep="\t")
    hpa["gen"] = hpa["gene"].astype(str).str.upper()
    en_sangre = set()
    for fila in hpa.itertuples(index=False):
        secretada = str(getattr(fila, "localizacion_secretoma", "")) == "Secreted to blood"
        concentracion = pd.notna(getattr(fila, "conc_sangre_ms_pg_l", None)) or pd.notna(getattr(fila, "conc_sangre_inmunoensayo_pg_l", None))
        if secretada or concentracion:
            en_sangre.add(fila.gen)
    mapa_acc = pd.read_csv(settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv"), sep="\t", usecols=["gene_symbol", "uniprot_base"])
    principal = set(pd.read_csv(salida / "suero_matriz_principal.tsv", sep="\t", usecols=["protein_id"])["protein_id"].astype(str))
    genes_suero = set(mapa_acc.loc[mapa_acc["uniprot_base"].astype(str).isin(principal), "gene_symbol"].astype(str).str.upper())
    k4 = pd.read_csv(salida / "k4_genes_maxlfq_crc_frente_a_colonoscopia_negativa.tsv", sep="\t")
    k4["gen"] = k4["protein_id"].astype(str).str.upper()
    k4_clase = dict(zip(k4["gen"], k4["clase"], strict=False))
    k4_fc = dict(zip(k4["gen"], k4["log2fc"], strict=False))
    h10 = pd.read_csv(salida / "h10_especificidad.tsv", sep="\t")
    h10_clase = dict(zip(h10["gen"].astype(str).str.upper(), h10["clase"], strict=False))
    filas = []
    for fila in cptac.itertuples(index=False):
        gen = str(fila.gen)
        ultimo = 0
        motivo = "no es diferencial certificada en CPTAC"
        if str(fila.significant).lower() in {"true", "1"}:
            ultimo = 1
            fc_t = tanaka_fc.get(gen)
            if fc_t is None:
                motivo = "no está medida en Tanaka primario frente a mucosa"
            elif _signo(fc_t) != _signo(fila.log2fc):
                motivo = "el signo en Tanaka no coincide con el de CPTAC"
            else:
                ultimo = 3
                if gen not in en_sangre:
                    motivo = "HPA no la anota como secretada a sangre ni con concentración plasmática"
                elif gen not in genes_suero:
                    ultimo = 4
                    motivo = "no está en la matriz principal del suero"
                else:
                    ultimo = 5
                    if k4_clase.get(gen) != "diferente" or _signo(k4_fc.get(gen)) != _signo(fila.log2fc):
                        motivo = "K4 CP frente a CN no es diferente con el mismo signo que el tejido"
                    elif h10_clase.get(gen) != "especifica_crc":
                        ultimo = 6
                        motivo = "H10 no la clasifica como específica de CRC"
                    else:
                        ultimo = 7
                        motivo = ""
        filas.append({
            "gen": gen,
            "padj_cptac": float(fila.padj),
            "log2fc_cptac": float(fila.log2fc),
            "ultimo_paso": ultimo,
            "motivo_no_pasa_el_siguiente": motivo,
            "clase_h5": clase_h5.get(gen, "no_informada"),
        })
    tabla = pd.DataFrame(filas).sort_values(["ultimo_paso", "padj_cptac"], ascending=[False, True])
    volcar_tsv(tabla, salida / "priorizacion.tsv", sello)
    pasos = [{"paso": i, "nombre": PASOS[i - 1], "n": int((tabla["ultimo_paso"] >= i).sum())} for i in range(1, 8)]
    llegan_5 = tabla.loc[tabla["ultimo_paso"] >= 5, "gen"].tolist()
    payload = {
        "pasos": pasos,
        "nadie_pasa_el_paso_6": bool((tabla["ultimo_paso"] >= 6).sum() == 0),
        "hipotesis_para_cohortes_futuras": llegan_5,
        "lectura": "Lista de hipótesis para cohortes futuras, no biomarcadores.",
    }
    volcar_json(payload, salida / "priorizacion_embudo.json", sello)
    return payload


def construir_manifiesto(settings, salida: Path) -> Path:
    lineas = []
    for relativa in ENTRADAS:
        ruta = settings.root / relativa if relativa.startswith("config/") else settings.resolve(relativa)
        if ruta.is_file() and ruta.stat().st_size > 0:
            lineas.append(f"{_sha256(ruta)}  {relativa}")
        else:
            lineas.append(f"AUSENTE  {relativa}")
    crudos = settings.root / "results" / "v3_7" / "crudos_integridad.tsv"
    if crudos.is_file():
        lineas.append(f"{_sha256(crudos)}  results/v3_7/crudos_integridad.tsv")
    else:
        lineas.append("AUSENTE  results/v3_7/crudos_integridad.tsv")
    destino = salida / "manifest_entradas.sha256"
    destino.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    return destino


def marcar_declaradas(eventos: list[dict], manifiesto: Path, salida: Path) -> list[dict]:
    declaradas = set()
    for linea in manifiesto.read_text(encoding="utf-8").splitlines():
        partes = linea.split("  ", 1)
        if len(partes) == 2:
            declaradas.add(partes[1])
    salida_res = str(salida.resolve())
    for evento in eventos:
        texto = str(evento["ruta"]).replace("\\", "/")
        if str(Path(evento["ruta"]).resolve()).startswith(salida_res):
            evento["declarada"] = "salida"
        elif any(texto.endswith(item) or f"/{item}" in texto for item in declaradas):
            evento["declarada"] = "entrada"
        else:
            evento["declarada"] = "no_declarada"
    return eventos


def version_trfp(settings) -> str:
    binario = settings.resolve("external/bin/trfp/ThermoRawFileParser")
    if not binario.is_file():
        return "no_registrada"
    try:
        proc = subprocess.run([str(binario), "--version"], capture_output=True, text=True, timeout=60, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return "no_registrada"
    texto = (proc.stdout or proc.stderr or "").strip()
    return texto.splitlines()[0] if texto else "no_registrada"


def _leer_campo(path: Path, campo: str):
    if not path.is_file():
        return None
    if path.suffix == ".json":
        datos = json.loads(path.read_text(encoding="utf-8"))
        for parte in campo.split("."):
            if isinstance(datos, dict) and parte in datos:
                datos = datos[parte]
            else:
                return None
        return datos
    filtro, columna = campo.split("|", 1)
    clave, valor = filtro.split("=", 1)
    tabla = pd.read_csv(path, sep="\t")
    filas = tabla.loc[tabla[clave].astype(str).eq(valor)]
    if filas.empty or columna not in filas.columns:
        return None
    return filas.iloc[0][columna]


def comparar_esperado(salida: Path, esperado: Path) -> tuple[pd.DataFrame, bool]:
    tabla = pd.read_csv(esperado, sep="\t", dtype=str)
    filas = []
    todo = True
    for fila in tabla.itertuples(index=False):
        valor = _leer_campo(salida / fila.fichero, fila.campo)
        try:
            dentro = abs(float(valor) - float(fila.valor)) <= float(fila.tolerancia)
        except (TypeError, ValueError):
            dentro = str(valor) == str(fila.valor)
        todo = todo and dentro
        filas.append({
            "metrica": fila.metrica, "valor_esperado": fila.valor, "valor_esta_corrida": valor,
            "fichero": fila.fichero, "campo": fila.campo, "tolerancia": fila.tolerancia,
            "dentro_de_tolerancia": dentro,
        })
    return pd.DataFrame(filas), todo


def embudo_es_encajado(tabla: pd.DataFrame, n_pasos: int = 7) -> bool:
    alcanzados = [set(tabla.loc[tabla["ultimo_paso"] >= paso, "gen"]) for paso in range(0, n_pasos + 1)]
    return all(alcanzados[k] <= alcanzados[k - 1] for k in range(1, len(alcanzados)))


def escribir_corrida(salida: Path, sello: dict, inicio: str, codigo: int | None) -> None:
    destino = salida / "corrida.json"
    previo = json.loads(destino.read_text(encoding="utf-8")) if destino.is_file() else {}
    payload = {
        "pid": os.getpid(),
        "inicio": previo.get("inicio", inicio),
        "fin": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) if codigo is not None else None,
        "codigo": codigo,
    }
    volcar_json(payload, destino, sello)


def correr_tests(root: Path, salida: Path) -> int:
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "--tb=line", "-p", "no:cacheprovider"],
        cwd=root, capture_output=True, text=True, check=False,
    )
    (salida / "tests.txt").write_text(proc.stdout + proc.stderr, encoding="utf-8")
    return int(proc.returncode)


def casillas_de_esta_corrida(salida: Path, settings, prereg, prohibidas: list[str], lecturas: list[dict], tests_codigo: int, sello: dict) -> list[dict]:
    def caja(nombre, fuente, criterio, valor, estado):
        if fuente.startswith("results/v3_7/"):
            ruta = salida / fuente.split("results/v3_7/", 1)[1]
        else:
            ruta = settings.root / fuente
        return {"nombre": nombre, "fuente": fuente, "ruta": ruta, "criterio": criterio, "valor": valor, "estado": estado}

    casillas = []
    manifiesto = salida / "manifest_entradas.sha256"
    ausentes = [linea for linea in manifiesto.read_text(encoding="utf-8").splitlines() if linea.startswith("AUSENTE")] if manifiesto.is_file() else ["falta el manifiesto"]
    no_declaradas = [ev["ruta"] for ev in lecturas if ev.get("declarada") == "no_declarada" and ev.get("clase") != "prohibida"]
    casillas.append(caja(
        "entradas", "results/v3_7/manifest_entradas.sha256",
        "cada entrada declarada tiene sha256 y ninguna lectura del vigilante queda fuera del manifiesto",
        {"ausentes": ausentes, "no_declaradas": len(no_declaradas)},
        "pasa" if not ausentes and not no_declaradas else "falla",
    ))
    casillas.append({
        "nombre": "sellos", "fuente": "results/v3_7/sello_corrida.json", "ruta": salida / "sello_corrida.json",
        "criterio": "todos los sellos comparten commit, árbol, preregistro y manifiesto",
        "valor": sellos_coinciden(salida), "estado": sellos_coinciden(salida)["estado"],
    })
    corrida = json.loads((salida / "corrida.json").read_text(encoding="utf-8")) if (salida / "corrida.json").is_file() else {}
    casillas.append(caja(
        "corrida_unica", "results/v3_7/corrida.json",
        "el mismo proceso escribe el pid, el inicio, el fin y el código 0",
        corrida, "pasa" if corrida.get("pid") == os.getpid() and corrida.get("inicio") and corrida.get("fin") and corrida.get("codigo") == 0 else "falla",
    ))
    controles = pd.read_csv(salida / "controles.tsv", sep="\t") if (salida / "controles.tsv").is_file() else pd.DataFrame()
    fallos = controles.loc[controles["pasa"].astype(str).eq("falla")] if "pasa" in controles.columns else pd.DataFrame()
    tejido = json.loads((salida / "controles_tejido.json").read_text(encoding="utf-8")) if (salida / "controles_tejido.json").is_file() else {}
    no_medibles = bool(len(fallos)) and set(fallos["observado"].astype(str)) <= {"no_evaluable"} and tejido.get("veredicto_v3_4") == "pasa"
    casillas.append(caja(
        "controles_suero", "results/v3_7/controles.tsv",
        "los fallos son controles no medibles en la réplica; controles_tejido.json veredicto_v3_4 es pasa",
        fallos[["control", "observado", "n_filas_del_gen"]].to_dict("records") if len(fallos) else [],
        "pasa" if no_medibles else "falla",
    ))
    k1 = json.loads((salida / "k1_etiquetas.json").read_text(encoding="utf-8")) if (salida / "k1_etiquetas.json").is_file() else {}
    cobertura = k1.get("cobertura") or {}
    casillas.append(caja(
        "h1", "results/v3_7/k1_etiquetas.json",
        "la cláusula de proteína da el veredicto; C15 condiciona la concordancia a las proteínas que la réplica detecta",
        {"h1_veredicto": k1.get("h1_veredicto"), "cobertura": cobertura},
        "pasa" if k1.get("h1_veredicto") and cobertura.get("numerador") is not None else "falla",
    ))
    casillas.append(caja(
        "auditoria_rai2025", "results/v3_7/auditoria_rai2025.tsv",
        "BAIAP3 y CEACAM5 tienen un paso leído de suero_trazabilidad_grupos.tsv",
        None, "pasa" if (salida / "auditoria_rai2025.tsv").is_file() else "falla",
    ))
    enmienda = prereg.get("enmienda_v3_7_cierre.no_se_ejecutan_en_esta_tesis") or {}
    casillas.append({
        "nombre": "alcance", "fuente": "config/preregistro_v3.yaml", "ruta": settings.root / "config" / "preregistro_v3.yaml",
        "criterio": "lo que la enmienda declara no ejecutado es desviación, no un resultado",
        "valor": sorted(enmienda), "estado": "desviacion_declarada" if enmienda else "falla",
    })
    changelog, reproduce = comparar_esperado(salida, settings.root / "config" / "esperado_v3_7.tsv")
    volcar_tsv(changelog, salida / "changelog_numeros.tsv", sello)
    casillas.append(caja(
        "reproduce", "config/esperado_v3_7.tsv",
        "cada métrica de config/esperado_v3_7.tsv coincide dentro de su tolerancia",
        {"n": int(len(changelog)), "fuera": int((~changelog["dentro_de_tolerancia"]).sum()) if len(changelog) else None},
        "pasa" if reproduce else "falla",
    ))
    texto_tests = (salida / "tests.txt").read_text(encoding="utf-8") if (salida / "tests.txt").is_file() else ""
    saltados = [
        linea for linea in texto_tests.splitlines()
        if "SKIPPED" in linea and "ONCOGO_REGRESION_V1" not in linea
    ]
    casillas.append(caja(
        "tests", "results/v3_7/tests.txt",
        "pytest termina en 0. La única omisión admitida es reejecutar v1, que pide ONCOGO_REGRESION_V1=1",
        {"codigo": tests_codigo}, "pasa" if tests_codigo == 0 and not saltados else "falla",
    ))
    for nombre, fichero, campo in (
        ("suero", "suero_embudo.json", "celda_preregistrada"),
        ("ml", "ml_celda_preregistrada.json", "auc_oof_media"),
        ("multiverso", "multiverso_resumen.json", "n_con_diferente"),
        ("h5", "h5_especificidad_organo.json", "estado"),
        ("h6", "h6_origen_celular.json", "estado"),
        ("h7", "h7_visibilidad_en_sangre.json", "estado"),
        ("h8", "h8_carga_tumoral.json", "estado"),
        ("h9", "h9_wald.json", "estado"),
        ("h10", "h10_especificidad.json", "estado"),
        ("preanalitica", "preanalitica.json", "estado"),
        ("go", "go_tejido.json", "veredicto_h1"),
        ("ml_explicacion", "ml_explicacion.json", "comprobacion_shap"),
        ("priorizacion", "priorizacion_embudo.json", "nadie_pasa_el_paso_6"),
    ):
        ruta = salida / fichero
        valor = None
        estado = "pasa"
        if ruta.is_file():
            datos = json.loads(ruta.read_text(encoding="utf-8"))
            valor = datos.get(campo)
            if nombre == "ml_explicacion" and isinstance(valor, dict):
                for item in valor.values():
                    if isinstance(item, str) or (isinstance(item, (int, float)) and not isinstance(item, bool) and item > 1e-8):
                        estado = "falla"
            if nombre == "go" and valor != "no_emite":
                estado = "falla"
        casillas.append(caja(nombre, f"results/v3_7/{fichero}", "el fichero de esta corrida existe", valor, estado))
    if prohibidas:
        casillas.append({
            "nombre": "sin_mezcla", "fuente": "results/v3_7/lecturas_v3_7.json", "ruta": salida / "lecturas_v3_7.json",
            "criterio": "el vigilante no vio lecturas prohibidas", "valor": prohibidas[:5], "estado": "falla",
        })
    return casillas
