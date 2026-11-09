"""H7. Las variables conocidas antes de mirar la sangre predicen qué proteínas tumorales cambian en sangre."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from oncogo_crc.settings import Settings
from oncogo_crc.v3.salida import volcar_json, volcar_tsv
from oncogo_crc.v3.suero_v32 import semilla_v33

N_GRUPOS_PROTEINA = 5
N_BOOTSTRAP = 2000
MINIMO_POR_CLASE = 10
GO_COMPONENTE = {
    "go_extracelular": "GO:0005576",
    "go_membrana": "GO:0005886",
    "go_exosoma": "GO:0070062",
}


def tubo(variables: list[str], columnas_esm: list[str], semilla: int) -> Pipeline:
    bloques = [("var", Pipeline([
        ("imputar", SimpleImputer(strategy="median")),
        ("escalar", StandardScaler()),
    ]), variables)]
    if columnas_esm:
        bloques.append(("esm", Pipeline([
            ("imputar", SimpleImputer(strategy="median")),
            ("escalar", StandardScaler()),
            ("pca", PCA(n_components=5, random_state=semilla)),
        ]), columnas_esm))
    return Pipeline([
        ("columnas", ColumnTransformer(bloques)),
        ("modelo", LogisticRegressionCV(
            Cs=np.logspace(-1, 4, 10),
            cv=StratifiedKFold(5, shuffle=True, random_state=semilla),
            penalty="l2",
            scoring="roc_auc",
            max_iter=5000,
        )),
    ])


def fuera_de_recurso(tabla: pd.DataFrame, recursos: list[str], recurso_test: str,
                     variables: list[str], columnas_esm: list[str], semilla: int) -> pd.Series:
    """tabla tiene una fila por proteína y una columna y_<recurso> con 1, 0 o NaN."""
    rng = np.random.default_rng(semilla)
    grupo = pd.Series(rng.integers(0, N_GRUPOS_PROTEINA, size=len(tabla)), index=tabla.index)
    pred = pd.Series(np.nan, index=tabla.index)
    otros = [r for r in recursos if r != recurso_test]
    for g in range(N_GRUPOS_PROTEINA):
        entrenar = grupo.index[grupo != g]
        largo = pd.concat([
            pd.DataFrame({"proteina": entrenar, "y": tabla.loc[entrenar, f"y_{r}"].to_numpy()}) for r in otros
        ]).dropna()
        if largo["y"].nunique() < 2:
            continue
        X = tabla.loc[largo["proteina"], variables + columnas_esm]
        modelo = tubo(variables, columnas_esm, semilla).fit(X, largo["y"].astype(int))
        probar = [p for p in grupo.index[grupo == g] if pd.notna(tabla.at[p, f"y_{recurso_test}"])]
        if probar:
            pred.loc[probar] = modelo.predict_proba(tabla.loc[probar, variables + columnas_esm])[:, 1]
    return pred


def auc_bootstrap(y, p, semilla: int) -> dict:
    y = np.asarray(y, dtype=int)
    p = np.asarray(p, dtype=float)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    if len(pos) < MINIMO_POR_CLASE or len(neg) < MINIMO_POR_CLASE:
        return {"estado": "no_evaluable", "n_pos": int(len(pos)), "n_neg": int(len(neg))}
    rng = np.random.default_rng(semilla)
    valores = []
    for _ in range(N_BOOTSTRAP):
        i = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        valores.append(roc_auc_score(y[i], p[i]))
    logit = np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))
    cal = LogisticRegression(C=1e6).fit(logit.reshape(-1, 1), y)
    return {
        "estado": "calculado",
        "auc": float(roc_auc_score(y, p)),
        "ic95_bajo": float(np.quantile(valores, 0.025)),
        "ic95_alto": float(np.quantile(valores, 0.975)),
        "pendiente_calibracion": float(cal.coef_[0, 0]),
        "intercepto_calibracion": float(cal.intercept_[0]),
        "n_pos": int(len(pos)),
        "n_neg": int(len(neg)),
    }


def veredicto(por_recurso: dict[str, dict], umbral: float = 0.55, minimo: int = 2) -> str:
    calculados = [r for r in por_recurso.values() if r.get("estado") == "calculado"]
    if len(calculados) < minimo:
        return "no_evaluable"
    if sum(r["ic95_bajo"] > umbral for r in calculados) >= minimo:
        return "sostiene"
    if all(r["ic95_bajo"] <= 0.5 <= r["ic95_alto"] for r in calculados):
        return "refuta"
    return "no_concluyente"


def _etiqueta(padj, signo_recurso, signo_cptac) -> float:
    if not np.isfinite(padj) or not np.isfinite(signo_recurso) or not np.isfinite(signo_cptac) or signo_cptac == 0:
        return float("nan")
    if padj < 0.05 and np.sign(signo_recurso) == np.sign(signo_cptac):
        return 1.0
    return 0.0


def construir_tabla_h7(settings: Settings, salida: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Devuelve la tabla, la ficha de variables y la lista de las que faltan."""
    cptac = pd.read_csv(settings.resolve("results/certified/cptac_unshared_treat.tsv"), sep="\t")
    cptac = cptac.loc[cptac["significant"].astype(bool)].copy()
    cptac["gen"] = cptac["protein_id"].astype(str).str.upper()
    cptac = cptac.drop_duplicates("gen")
    tabla = cptac.set_index("gen")[["log2fc"]].rename(columns={"log2fc": "log2fc_cptac"})
    tabla["magnitud_en_tejido"] = tabla["log2fc_cptac"].abs()
    variables = [{
        "variable": "magnitud_en_tejido",
        "fuente": "results/certified/cptac_unshared_treat.tsv",
        "definicion": "valor absoluto del log2FC de CPTAC en las proteínas certificadas",
    }]
    faltantes = []

    go = settings.resolve("metadata/gene_ontology/go_gene_long.tsv")
    if go.is_file():
        anotacion = pd.read_csv(go, sep="\t", usecols=["gene_symbol", "go_id"], dtype=str)
        for nombre, termino in GO_COMPONENTE.items():
            genes = set(anotacion.loc[anotacion["go_id"].eq(termino), "gene_symbol"].str.upper())
            tabla[nombre] = tabla.index.isin(genes).astype(float)
            variables.append({"variable": nombre, "fuente": str(go), "definicion": termino})
    else:
        faltantes.append("componente_celular_go")

    hpa = settings.resolve("external/hpa/hpa_concentracion_sangre.tsv")
    if hpa.is_file():
        sangre = pd.read_csv(hpa, sep="\t", dtype=str).drop_duplicates("gene")
        sangre["gen"] = sangre["gene"].str.upper()
        sangre = sangre.set_index("gen")
        tabla["hpa_secreted_to_blood"] = sangre["localizacion_secretoma"].eq("Secreted to blood").reindex(tabla.index).fillna(False).astype(float)
        conc = pd.to_numeric(sangre["conc_sangre_ms_pg_l"], errors="coerce").reindex(tabla.index)
        tabla["concentracion_plasmatica_hpa_log10"] = np.where(conc > 0, np.log10(conc), np.nan)
        variables.append({"variable": "hpa_secreted_to_blood", "fuente": str(hpa), "definicion": "localizacion_secretoma == Secreted to blood"})
        variables.append({"variable": "concentracion_plasmatica_hpa_log10", "fuente": str(hpa), "definicion": "log10 de conc_sangre_ms_pg_l"})
    else:
        faltantes.extend(["hpa_secreted_to_blood", "concentracion_plasmatica_hpa_log10"])
    faltantes.append("especificidad_tisular")

    h5 = salida / "h5_especificidad_organo.tsv"
    if h5.is_file():
        clases = pd.read_csv(h5, sep="\t", dtype=str)
        if "gen" in clases.columns and "clase" in clases.columns:
            clases["gen"] = clases["gen"].str.upper()
            for clase in ("pan_carcinoma", "especifica_colon", "intermedia"):
                miembros = set(clases.loc[clases["clase"].eq(clase), "gen"])
                tabla[f"clase_h5_{clase}"] = tabla.index.isin(miembros).astype(float)
                variables.append({"variable": f"clase_h5_{clase}", "fuente": str(h5), "definicion": "one-hot de la clase H5"})
        else:
            faltantes.append("clase_h5")
    else:
        faltantes.append("clase_h5")

    k1 = salida / "k1_primario_frente_a_mucosa.tsv"
    grupos = settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt")
    uni = settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv")
    matriz = settings.resolve("data/raw/pride/pride_crc_matrix.tsv")
    if k1.is_file() and grupos.is_file() and uni.is_file() and matriz.is_file():
        from oncogo_crc.v3.etiquetas import matriz_log2_lfq
        from oncogo_crc.v3.genes_tanaka import gen_unico, mapa_genes
        cruda = pd.read_csv(matriz, sep="\t", low_memory=False)
        log2mat = matriz_log2_lfq(cruda)
        mapa = mapa_genes(grupos, uni)
        abundancia = {}
        for pid in log2mat.index.astype(str):
            gen = gen_unico(pid, mapa)
            if gen and gen not in abundancia:
                abundancia[gen] = float(np.nanmedian(log2mat.loc[pid].to_numpy(dtype=float)))
        tabla["abundancia_en_tejido"] = pd.Series(abundancia).reindex(tabla.index)
        variables.append({"variable": "abundancia_en_tejido", "fuente": str(matriz), "definicion": "mediana log2 en Tanaka, una fila por gen único"})
    else:
        faltantes.append("abundancia_en_tejido")

    faltantes.extend(["origen_celular_h6", "presencia_en_vesiculas_pxd033236", "esm2_t12_35M_pca5", "tmbed_helices_tm", "tmbed_peptido_senal"])
    return tabla, pd.DataFrame(variables), faltantes


def _recursos(settings: Settings, tabla: pd.DataFrame) -> list[str]:
    alvez = settings.resolve("external/alvez_2023/alvez_2023_colorrectal_por_proteina.tsv")
    if alvez.is_file():
        bloque = pd.read_csv(alvez, sep="\t")
        bloque["gen"] = bloque["proteina"].astype(str).str.upper()
        bloque = bloque.drop_duplicates("gen").set_index("gen")
        signo = bloque["direccion"].map({"sube": 1.0, "baja": -1.0})
        padj = pd.to_numeric(bloque["p_ajustado"], errors="coerce")
        tabla["y_alvez_2023"] = [
            _etiqueta(padj.get(g, np.nan), signo.get(g, np.nan), tabla.at[g, "log2fc_cptac"])
            if g in bloque.index else np.nan
            for g in tabla.index
        ]
    de = settings.resolve("external/hpa/blood_pea_disease_de.tsv.zip")
    if de.is_file():
        import zipfile
        from oncogo_crc.v3.h10 import tipo_de_control
        with zipfile.ZipFile(de) as zf:
            hpa = pd.read_csv(zf.open(zf.namelist()[0]), sep="\t")
        hpa = hpa.rename(columns={"Gene": "gen", "p-value adjusted": "padj"})
        hpa = hpa.loc[hpa["Disease"].eq("Colorectal cancer")].copy()
        hpa["tipo_control"] = hpa["Control"].map(tipo_de_control)
        hpa = hpa.loc[hpa["tipo_control"].eq("sanos")].drop_duplicates("gen")
        hpa["gen"] = hpa["gen"].astype(str).str.upper()
        hpa = hpa.set_index("gen")
        tabla["y_hpa_disease_blood_atlas"] = [
            _etiqueta(
                float(hpa.at[g, "padj"]) if g in hpa.index else np.nan,
                float(hpa.at[g, "logFC"]) if g in hpa.index else np.nan,
                tabla.at[g, "log2fc_cptac"],
            )
            for g in tabla.index
        ]
    return [r for r in ("alvez_2023", "hpa_disease_blood_atlas", "papier_2024") if f"y_{r}" in tabla.columns]


def ejecutar(settings: Settings, salida: Path, sello: dict, prereg) -> dict:
    tabla, ficha, faltantes = construir_tabla_h7(settings, salida)
    recursos = _recursos(settings, tabla)
    if "papier_2024" not in recursos:
        faltantes.append("papier_2024")
    variables = [c for c in ficha["variable"] if c in tabla.columns]
    volcar_tsv(ficha, salida / "h7_variables.tsv", sello)
    volcar_tsv(tabla.reset_index(names="gen"), salida / "h7_tabla_proteinas.tsv", sello)
    if len(recursos) < 2 or not variables:
        payload = {
            "estado": "no_evaluable",
            "motivo": "hacen falta dos recursos de sangre o una variable del preregistro",
            "recursos": recursos,
            "variables_faltantes": faltantes,
            "modelo_secundario": "no_disponible",
            "nota_tmbed": "TMbed no se ejecutó. El sustituto de UniProt no entra en el modelo principal.",
            "nota_tripod": "TRIPOD+AI: esto no es una prueba diagnóstica. La validación es fuera de recurso y fuera de proteína.",
        }
        volcar_json(payload, salida / "h7_visibilidad_en_sangre.json", sello)
        return payload
    semilla = semilla_v33(prereg, "h7_modelo")
    semilla_b = semilla_v33(prereg, "h7_bootstrap")
    predicciones = []
    por_recurso = {}
    for recurso in recursos:
        pred = fuera_de_recurso(tabla, recursos, recurso, variables, [], semilla)
        y = tabla[f"y_{recurso}"]
        ok = y.notna() & pred.notna()
        por_recurso[recurso] = auc_bootstrap(y.loc[ok].to_numpy(), pred.loc[ok].to_numpy(), semilla_b)
        for gen in tabla.index[ok]:
            predicciones.append({"gen": gen, "recurso": recurso, "y": float(y.at[gen]), "p": float(pred.at[gen])})
    base_vars = [c for c in ("concentracion_plasmatica_hpa_log10", "magnitud_en_tejido") if c in variables]
    incremental = {"estado": "no_disponible", "motivo": "el modelo secundario de secuencia no se ajustó en esta corrida"}
    if len(base_vars) == 2 and recursos:
        pred_base = fuera_de_recurso(tabla, recursos, recursos[0], base_vars, [], semilla)
        y0 = tabla[f"y_{recursos[0]}"]
        ok = y0.notna() & pred_base.notna()
        base = auc_bootstrap(y0.loc[ok].to_numpy(), pred_base.loc[ok].to_numpy(), semilla_b)
        completo = por_recurso[recursos[0]]
        if base.get("estado") == "calculado" and completo.get("estado") == "calculado":
            incremental = {
                "estado": "calculado",
                "recurso": recursos[0],
                "auc_base": base["auc"],
                "auc_completo": completo["auc"],
                "incremento": completo["auc"] - base["auc"],
                "nota": "El completo usa las variables que existen. ESM-2 y TMbed no entraron: quedan fuera del principal.",
            }
    minimo = 2 if "papier_2024" not in recursos else 2
    estado = veredicto(por_recurso, minimo=minimo)
    payload = {
        "estado": estado,
        "por_recurso": por_recurso,
        "variables_usadas": variables,
        "variables_faltantes": faltantes,
        "recursos": recursos,
        "modelo_secundario": incremental,
        "nota_alvez_hpa": "El solapamiento de casos entre Álvez 2023 y el atlas HPA no está comprobado. Si lo hay, no son dos recursos independientes.",
        "nota_tmbed": "TMbed y ESM-2 no se ejecutaron en el modelo principal. ProtT5 queda en el disco externo si se descarga; no cambia la regla.",
        "nota_tripod": "TRIPOD+AI: no es una prueba diagnóstica. Validación fuera de recurso y fuera de proteína, 5 grupos, semilla h7_modelo.",
        "fase_pepe": 1,
    }
    volcar_tsv(pd.DataFrame(predicciones), salida / "h7_predicciones.tsv", sello)
    volcar_json(payload, salida / "h7_visibilidad_en_sangre.json", sello)
    return payload
