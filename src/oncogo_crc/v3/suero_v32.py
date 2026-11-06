"""S2 a S9. El análisis oficial es treat y TOST sobre la matriz filtrada.

Los contaminantes y los paneles salen del contraste principal. Vuelven solo
en las sensibilidades declaradas en la enmienda. No se imputa.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

from oncogo_crc.settings import Settings
from oncogo_crc.tesis.cierre import a_escala_log2
from oncogo_crc.utils.dia_nn_plasma import build_plasma_metadata_for_pxd055218_runs, matriz_filtrada
from oncogo_crc.utils.moderated import moderate_table, unpaired_observed_stats
from oncogo_crc.v3.auc import auc_firmado, bootstrap_auc, delong_ic, nulo_del_maximo_plegado
from oncogo_crc.v3.faltantes import fisher_presencia
from oncogo_crc.v3.potencia import alfa_efectivo, clasificar, efecto_minimo_detectable, p_tost
from oncogo_crc.v3.prereg import PreregistroV3
from oncogo_crc.v3.salida_v32 import volcar_json, volcar_tsv

REACTIVOS = frozenset({"TRY1", "TRY2", "PRSS1", "PRSS2", "CSN1S1", "CSN1S2", "CSN1S2A", "CSN2", "CSN3"})
PARES = {
    "crc_frente_a_colonoscopia_negativa": ("crc", "colonoscopy_negative"),
    "crc_frente_a_sanos": ("crc", "healthy"),
    "colonoscopia_negativa_frente_a_sanos": ("colonoscopy_negative", "healthy"),
}
Z95 = 1.959963984540054


def semilla_v32(prereg: PreregistroV3, nombre: str) -> int:
    base = int(prereg.get("enmienda_v3_2.semillas.base"))
    return base + int(prereg.get(f"enmienda_v3_2.semillas.desplazamientos.{nombre}"))


def _log2(ancha: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    mat = ancha.drop_duplicates("protein_id", keep="first").set_index("protein_id")
    mat.index = mat.index.astype(str)
    valores, escala = a_escala_log2(mat.to_numpy(dtype=float))
    return pd.DataFrame(valores, index=mat.index, columns=mat.columns.astype(str)), escala


def _grupos(columnas: list[str]) -> dict[str, list[str]]:
    meta = build_plasma_metadata_for_pxd055218_runs(columnas)
    grupos: dict[str, list[str]] = {}
    for condicion, bloque in meta.groupby("condition"):
        grupos[str(condicion)] = bloque["sample_id"].astype(str).tolist()
    return grupos


def _simbolos(texto: object) -> set[str]:
    return {parte.strip() for parte in str(texto).split(";") if parte.strip() and parte.strip() != "nan"}


def ids_fuera_de_principal(traza: pd.DataFrame) -> set[str]:
    """Contaminante de la librería, piel o panel de Geyer. El humano ya está filtrado."""
    panel = traza["panel"].fillna("").astype(str).str.len() > 0
    fuera = traza["contaminante"].astype(bool) | traza["contaminante_piel"].astype(bool) | panel
    return set(traza.loc[fuera, "protein_group"].astype(str))


def _recortar(mat: pd.DataFrame, fuera: set[str]) -> pd.DataFrame:
    keep = [i for i in mat.index.astype(str) if i not in fuera]
    return mat.loc[keep]


def _clases(tabla: pd.DataFrame, p_tost_vals: np.ndarray, alpha: float) -> pd.DataFrame:
    tabla = tabla.copy()
    tabla["p_tost"] = p_tost_vals
    padj = np.full(len(tabla), np.nan)
    ok = np.isfinite(p_tost_vals)
    if ok.any():
        _, q, _, _ = multipletests(p_tost_vals[ok], alpha=alpha, method="fdr_bh")
        padj[ok] = q
    tabla["padj_tost"] = padj
    tabla["clase"] = clasificar(
        tabla["significant"].fillna(False).to_numpy(dtype=bool),
        (tabla["padj_tost"] < alpha).fillna(False).to_numpy(dtype=bool),
    )
    return tabla


def _contraste(mat: pd.DataFrame, cols_a: list[str], cols_b: list[str], *, lfc: float, alpha: float, nombre: str, min_n: int = 3):
    coef, s2, df, stdev, n_obs = unpaired_observed_stats(
        mat[cols_a].to_numpy(dtype=float),
        mat[cols_b].to_numpy(dtype=float),
        min_n=min_n,
    )
    tabla, _, _ = moderate_table(
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
    return _clases(tabla, p_tost(tabla["log2fc"].to_numpy(), tabla["se"].to_numpy(), tabla["df_total"].to_numpy(), lfc), alpha)


def _contraste_covariable(mat, cols_a, cols_b, covariable: pd.Series, *, lfc, alpha, nombre, min_n=3):
    cols = list(cols_a) + list(cols_b)
    grupo = np.array([1.0] * len(cols_a) + [0.0] * len(cols_b))
    x = covariable.reindex(cols).to_numpy(dtype=float)
    y = mat[cols].to_numpy(dtype=float)
    n_prot = y.shape[0]
    coef = np.full(n_prot, np.nan)
    s2 = np.full(n_prot, np.nan)
    df = np.full(n_prot, np.nan)
    stdev = np.full(n_prot, np.nan)
    n_obs = np.zeros(n_prot)
    for i in range(n_prot):
        ok = np.isfinite(y[i]) & np.isfinite(x)
        n = int(ok.sum())
        if n < max(min_n * 2, 4):
            continue
        gg = grupo[ok]
        if int((gg == 1).sum()) < min_n or int((gg == 0).sum()) < min_n:
            continue
        xx = x[ok]
        X = np.column_stack([np.ones(n), gg, xx - float(xx.mean())])
        beta, *_ = np.linalg.lstsq(X, y[i, ok], rcond=None)
        resid = y[i, ok] - X @ beta
        dof = n - 3
        if dof < 1:
            continue
        var = float(resid @ resid / dof)
        try:
            inv = np.linalg.inv(X.T @ X)
        except np.linalg.LinAlgError:
            continue
        coef[i] = float(beta[1])
        s2[i] = var
        df[i] = dof
        stdev[i] = float(np.sqrt(max(inv[1, 1], 0.0)))
        n_obs[i] = n
    tabla, _, _ = moderate_table(
        protein_ids=mat.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method=f"moderated_t_treat_unpaired_coagulacion_{nombre}",
    )
    return _clases(tabla, p_tost(tabla["log2fc"].to_numpy(), tabla["se"].to_numpy(), tabla["df_total"].to_numpy(), lfc), alpha)


def _potencia(tabla: pd.DataFrame, prereg: PreregistroV3, rng_alfa, rng_mde) -> pd.DataFrame:
    alpha = float(prereg.get("inferencia.alpha"))
    lfc = float(prereg.get("inferencia.treat_lfc"))
    se = tabla["se"].to_numpy(dtype=float)
    df = tabla["df_total"].to_numpy(dtype=float)
    alfa = alfa_efectivo(se, df, alpha=alpha, lfc=lfc, n_sim=int(prereg.get("k4_suero.n_sim_alfa")), rng=rng_alfa)
    mde = efecto_minimo_detectable(
        se, df, alpha_eff=alfa, lfc=lfc,
        potencia=float(prereg.get("k4_suero.potencia_objetivo")),
        delta_max=float(prereg.get("k4_suero.delta_maximo_busqueda")),
        n_sim=int(prereg.get("k4_suero.n_sim_efecto")),
        n_delta=20,
        rng=rng_mde,
    )
    tabla = tabla.copy()
    tabla["alfa_efectivo"] = alfa
    tabla["efecto_minimo_detectable"] = mde
    return tabla


def _conteo(tabla: pd.DataFrame) -> dict:
    clases = tabla["clase"].value_counts().to_dict() if len(tabla) else {}
    return {"n": int(len(tabla)), "alfa_efectivo": float(tabla["alfa_efectivo"].iloc[0]) if len(tabla) else float("nan"), **{k: int(clases.get(k, 0)) for k in ("diferente", "equivalente", "indeterminada", "conflicto")}}


def _k4_bloque(mat, grupos, prereg, *, ajustar, coagulacion: pd.Series | None):
    alpha = float(prereg.get("inferencia.alpha"))
    lfc = float(prereg.get("inferencia.treat_lfc"))
    rng_alfa = np.random.default_rng(semilla_v32(prereg, "k4_alfa_efectivo"))
    rng_mde = np.random.default_rng(semilla_v32(prereg, "k4_efecto_minimo"))
    salida = {}
    tablas = {}
    for nombre, (a, b) in PARES.items():
        if a not in grupos or b not in grupos or not grupos[a] or not grupos[b]:
            salida[nombre] = {"estado": "no_disponible", "motivo": f"falta el grupo {a} o {b}"}
            continue
        if ajustar and coagulacion is not None:
            tabla = _contraste_covariable(mat, grupos[a], grupos[b], coagulacion, lfc=lfc, alpha=alpha, nombre=nombre)
        else:
            tabla = _contraste(mat, grupos[a], grupos[b], lfc=lfc, alpha=alpha, nombre=nombre)
        tabla = _potencia(tabla, prereg, rng_alfa, rng_mde)
        tablas[nombre] = tabla
        salida[nombre] = _conteo(tabla)
    return salida, tablas


def _puntuacion_z(mat: pd.DataFrame, filas: list[str]) -> pd.Series:
    presentes = [f for f in filas if f in mat.index]
    if not presentes:
        return pd.Series(np.nan, index=mat.columns)
    bloque = mat.loc[presentes]
    sd = bloque.std(axis=1, ddof=1).replace(0, np.nan)
    z = bloque.sub(bloque.mean(axis=1), axis=0).div(sd, axis=0)
    return z.median(axis=0, skipna=True)


def _hodges(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return {"estimador": float("nan"), "ic_bajo": float("nan"), "ic_alto": float("nan")}
    diffs = np.sort((a[:, None] - b[None, :]).ravel())
    k = int(diffs.size)
    se = float(np.sqrt(a.size * b.size * (a.size + b.size + 1) / 12.0))
    m = int(np.floor(k / 2.0 - Z95 * se))
    m = min(max(m, 0), k - 1)
    return {"estimador": float(np.median(diffs)), "ic_bajo": float(diffs[m]), "ic_alto": float(diffs[k - 1 - m])}


def _peptidos_humanos(path: Path) -> set[str] | None:
    if not path.is_file() or path.stat().st_size < 1000:
        return None
    abrir = gzip.open if path.suffix == ".gz" else open
    peptidos: set[str] = set()
    seq: list[str] = []

    def volcar(texto: str) -> None:
        amino = "".join(ch for ch in texto if ch.isalpha())
        if not amino:
            return
        frag: list[str] = []
        buf: list[str] = []
        for i, aa in enumerate(amino):
            buf.append(aa)
            if aa in "KR" and (i + 1 == len(amino) or amino[i + 1] != "P"):
                frag.append("".join(buf))
                buf = []
        if buf:
            frag.append("".join(buf))
        for i, trozo in enumerate(frag):
            if 7 <= len(trozo) <= 30:
                peptidos.add(trozo)
            if i + 1 < len(frag):
                duo = trozo + frag[i + 1]
                if 7 <= len(duo) <= 30:
                    peptidos.add(duo)

    with abrir(path, "rt", encoding="utf-8", errors="replace") as entrada:
        for linea in entrada:
            if linea.startswith(">"):
                if seq:
                    volcar("".join(seq))
                seq = []
            else:
                seq.append(linea.strip())
        if seq:
            volcar("".join(seq))
    return peptidos


def _especies(tabla: pd.DataFrame, traza: pd.DataFrame, peptidos: set[str] | None) -> pd.DataFrame:
    secuencias: dict[str, set[str]] = {}
    if "Stripped.Sequence" in tabla.columns:
        for grupo, bloque in tabla.groupby("Protein.Group", sort=False):
            secs = {str(s) for s in bloque["Stripped.Sequence"].dropna().astype(str) if 7 <= len(str(s)) <= 30}
            secuencias[str(grupo)] = secs
    filas = []
    for _, fila in traza.iterrows():
        grupo = str(fila["protein_group"])
        especies = _simbolos(str(fila["especies"]).replace(";", ",")) if False else set(filter(None, str(fila["especies"]).split(";")))
        genes = _simbolos(fila["genes"])
        secs = secuencias.get(grupo, set())
        if peptidos is None or not secs:
            ortologo = "no_disponible"
        else:
            ortologo = "si" if secs and secs <= peptidos else "no"
        reactivo = (not bool(fila["humano"])) and bool(fila["contaminante"]) and bool(especies & {"BOVIN", "PIG"} or genes & REACTIVOS)
        if bool(fila["humano"]):
            tipo = "humano"
        elif reactivo:
            tipo = "reactivo"
        elif ortologo == "si":
            tipo = "ortologo"
        elif ortologo == "no":
            tipo = "microorganismo"
        else:
            tipo = "no_disponible"
        filas.append({
            "protein_group": grupo,
            "especies": fila["especies"],
            "humano": bool(fila["humano"]),
            "ortologo_compartido": ortologo,
            "tipo": tipo,
            "n_peptidos_7_30": len(secs),
        })
    return pd.DataFrame(filas)


def _fraccion_no_humana(tabla: pd.DataFrame, especies: pd.DataFrame) -> pd.DataFrame:
    if "Precursor.Quantity" not in tabla.columns:
        return pd.DataFrame()
    no_humanos = set(especies.loc[~especies["humano"].astype(bool), "protein_group"].astype(str))
    trabajo = tabla.dropna(subset=["Precursor.Quantity"]).copy()
    trabajo["Run"] = trabajo["Run"].astype(str)
    trabajo["no_humano"] = trabajo["Protein.Group"].astype(str).isin(no_humanos)
    filas = []
    for corrida, bloque in trabajo.groupby("Run", sort=False):
        total = float(bloque["Precursor.Quantity"].sum())
        ajena = float(bloque.loc[bloque["no_humano"], "Precursor.Quantity"].sum())
        nombre = corrida.replace("\\", "/").split("/")[-1]
        if nombre.startswith("CP_"):
            grupo = "crc"
        elif nombre.startswith("CN_"):
            grupo = "colonoscopy_negative"
        elif nombre.startswith("Normal_"):
            grupo = "healthy"
        else:
            grupo = "otro"
        filas.append({"Run": corrida, "grupo": grupo, "fraccion_intensidad_no_humana": ajena / total if total else float("nan")})
    return pd.DataFrame(filas)


def _entrapment(especies: pd.DataFrame, traza: pd.DataFrame) -> dict:
    """Cota de Wen 2025. r no está en el depósito: no se inventa un cociente de bases."""
    pasa = traza["global_pg_q"] <= 0.01
    unidos = traza.loc[pasa, ["protein_group"]].merge(especies, on="protein_group", how="left")
    n_t = int(unidos["humano"].astype(bool).sum())
    entran = unidos.loc[~unidos["humano"].astype(bool) & (unidos["ortologo_compartido"] == "no") & (unidos["tipo"] != "reactivo")]
    n_e = int(len(entran))
    cota = float(n_e / (n_t + n_e)) if (n_t + n_e) else float("nan")
    return {
        "estado": "cota",
        "estimador": "FDP ≈ N_e * (1 + 1/r) / (N_t + N_e)",
        "nivel": "grupo_de_proteina_con_Global.PG.Q.Value_<=_0.01",
        "N_t": n_t,
        "N_e": n_e,
        "r": "no_disponible",
        "motivo_r": "PRIDE lista 60 crudos Thermo y el report. No hay FASTA ni report.log.txt, así que el cociente de péptidos únicos de las dos bases no se conoce.",
        "cota_inferior_si_r_es_grande": cota,
        "advertencia": "La base no se diseñó como entrapment. La cota inferior es el estimador cuando r tiende a infinito. Si r es menor, el estimador es mayor. No es un FDR oficial.",
    }


def _k5(mat: pd.DataFrame, grupos: dict[str, list[str]], prereg: PreregistroV3, etiqueta: str):
    min_n = int(prereg.get("k5_auc.minimo_observados_por_grupo"))
    n_boot = int(prereg.get("k5_auc.bootstrap.n"))
    rng_boot = np.random.default_rng(semilla_v32(prereg, "k5_bootstrap_auc"))
    rng_dif = np.random.default_rng(semilla_v32(prereg, "k5_bootstrap_diferencia"))
    rng_perm = np.random.default_rng(semilla_v32(prereg, "k5_permutacion_etiquetas"))
    filas = []
    for nombre, (ga, gb) in (
        ("crc_frente_a_colonoscopia_negativa", ("crc", "colonoscopy_negative")),
        ("colonoscopia_negativa_frente_a_sanos", ("colonoscopy_negative", "healthy")),
    ):
        if ga not in grupos or gb not in grupos:
            continue
        for proteina in mat.index:
            pos = mat.loc[proteina, grupos[ga]].to_numpy(dtype=float)
            neg = mat.loc[proteina, grupos[gb]].to_numpy(dtype=float)
            if np.isfinite(pos).sum() < min_n or np.isfinite(neg).sum() < min_n:
                continue
            auc, lo, hi = delong_ic(pos, neg)
            blo, bhi = bootstrap_auc(pos, neg, n_boot=n_boot, rng=rng_boot)
            filas.append({
                "version": etiqueta,
                "contraste": nombre,
                "protein_id": proteina,
                "auc_con_signo": auc,
                "delong_bajo": lo,
                "delong_alto": hi,
                "bootstrap_bajo": blo,
                "bootstrap_alto": bhi,
            })
    aucs = pd.DataFrame(filas)
    cols = grupos.get("crc", []) + grupos.get("colonoscopy_negative", [])
    etiquetas = np.array([1] * len(grupos.get("crc", [])) + [0] * len(grupos.get("colonoscopy_negative", [])))
    nulo = nulo_del_maximo_plegado(
        mat[cols].to_numpy(dtype=float),
        etiquetas,
        n_perm=int(prereg.get("k5_auc.nulo_del_maximo.n_permutaciones")),
        min_n=min_n,
        rng=rng_perm,
    )
    nulo["version"] = etiqueta
    diferencias = _diferencias(mat, grupos, n_boot, rng_dif, min_n)
    return aucs, nulo, diferencias


def _diferencias(mat, grupos, n_boot, rng, min_n) -> pd.DataFrame:
    if not all(k in grupos and grupos[k] for k in ("crc", "colonoscopy_negative", "healthy")):
        return pd.DataFrame()
    filas = []
    for proteina in mat.index:
        a = mat.loc[proteina, grupos["crc"]].to_numpy(dtype=float)
        b = mat.loc[proteina, grupos["colonoscopy_negative"]].to_numpy(dtype=float)
        c = mat.loc[proteina, grupos["healthy"]].to_numpy(dtype=float)
        a, b, c = a[np.isfinite(a)], b[np.isfinite(b)], c[np.isfinite(c)]
        if min(a.size, b.size, c.size) < min_n:
            continue
        obs = float(auc_firmado(a, b) - auc_firmado(b, c))
        diffs = np.empty(n_boot, dtype=float)
        for i in range(n_boot):
            aa = a[rng.integers(0, a.size, a.size)]
            bb = b[rng.integers(0, b.size, b.size)]
            cc = c[rng.integers(0, c.size, c.size)]
            diffs[i] = auc_firmado(aa, bb) - auc_firmado(bb, cc)
        bajo, alto = np.nanquantile(diffs, [0.025, 0.975])
        cola = float(np.mean(diffs * np.sign(obs) <= 0)) if obs != 0 else float(np.mean(np.abs(diffs) <= 0))
        filas.append({
            "protein_id": proteina,
            "diferencia": obs,
            "ci_bajo": float(bajo),
            "ci_alto": float(alto),
            "p_bootstrap": min(1.0, 2.0 * cola),
        })
    dif = pd.DataFrame(filas)
    if len(dif):
        ok = np.isfinite(dif["p_bootstrap"].to_numpy())
        padj = np.full(len(dif), np.nan)
        if ok.any():
            _, q, _, _ = multipletests(dif.loc[ok, "p_bootstrap"], method="fdr_bh")
            padj[np.flatnonzero(ok)] = q
        dif["padj_bh"] = padj
    return dif


def _frase(principal: dict, con_paneles: dict, traza: pd.DataFrame, tabla_con_paneles: pd.DataFrame | None) -> str:
    crc = principal.get("crc_frente_a_colonoscopia_negativa", {})
    if crc.get("diferente", 0) == 0:
        return f"No se detecta diferencia; {int(crc.get('equivalente', 0))} proteínas son equivalentes"
    coag = set(traza.loc[traza["panel"].fillna("").astype(str).str.contains("coagulacion"), "protein_group"].astype(str))
    y = int(con_paneles.get("crc_frente_a_colonoscopia_negativa", {}).get("diferente", 0))
    z = 0
    if tabla_con_paneles is not None and len(tabla_con_paneles):
        z = int(tabla_con_paneles.loc[tabla_con_paneles["clase"] == "diferente", "protein_id"].astype(str).isin(coag).sum())
    return f"Tras quitar contaminantes, {y} proteínas difieren; {z} de ellas pertenecen al panel de coagulación"


def analizar_suero(tabla, traza, ancha, prereg: PreregistroV3, salida: Path, sello: dict, settings: Settings) -> None:
    mat, escala = _log2(ancha)
    grupos = _grupos(list(mat.columns))
    fuera = ids_fuera_de_principal(traza)
    solo_contaminante = set(traza.loc[traza["contaminante"].astype(bool) | traza["contaminante_piel"].astype(bool), "protein_group"].astype(str))
    principal = _recortar(mat, fuera)
    con_paneles = _recortar(mat, solo_contaminante)
    gen_de = traza.assign(protein_group=traza["protein_group"].astype(str)).drop_duplicates("protein_group").set_index("protein_group")["genes"]
    genes_de_fila = pd.Series({idx: next(iter(_simbolos(gen_de.get(idx, ""))), "") for idx in mat.index})
    paneles = {
        "eritrocitos": list(prereg.get("enmienda_v3_2.paneles_geyer_2019.eritrocitos")),
        "plaquetas": list(prereg.get("enmienda_v3_2.paneles_geyer_2019.plaquetas")),
        "coagulacion": list(prereg.get("enmienda_v3_2.paneles_geyer_2019.coagulacion")),
    }
    piel = list(prereg.get("enmienda_v3_2.contaminante_piel.genes"))
    puntuaciones = {nombre: _puntuacion_z(mat, genes_de_fila[genes_de_fila.isin(lista)].index.tolist()) for nombre, lista in {**paneles, "piel": piel}.items()}
    coag = puntuaciones["coagulacion"]

    fasta = settings.resolve("external/uniprot/human_reviewed_isoforms.fasta.gz")
    peptidos = _peptidos_humanos(fasta)
    especies = _especies(tabla, traza, peptidos)
    fraccion = _fraccion_no_humana(tabla, especies)
    volcar_tsv(especies, salida / "suero_especies.tsv", sello)
    if len(fraccion) >= 3 and fraccion["grupo"].nunique() >= 2:
        muestras = [fraccion.loc[fraccion["grupo"] == g, "fraccion_intensidad_no_humana"].dropna().to_numpy() for g in ("crc", "colonoscopy_negative", "healthy")]
        muestras = [m for m in muestras if m.size >= 3]
        p_kw = float(kruskal(*muestras).pvalue) if len(muestras) >= 2 else float("nan")
    else:
        p_kw = float("nan")
    volcar_json(_entrapment(especies, traza) | {"p_kruskal_fraccion_no_humana": p_kw, "ortologos": "calculado" if peptidos is not None else "no_disponible", "motivo_ortologos": None if peptidos is not None else f"no está {fasta}"}, salida / "suero_entrapment.json", sello)

    por_muestra = pd.DataFrame({"sample_id": mat.columns})
    for nombre, serie in puntuaciones.items():
        por_muestra[nombre] = serie.reindex(por_muestra["sample_id"]).to_numpy()
    if len(fraccion):
        por_muestra = por_muestra.merge(fraccion.rename(columns={"Run": "sample_id"}), on="sample_id", how="left")
    por_muestra["n_precursores"] = tabla.groupby(tabla["Run"].astype(str)).size().reindex(por_muestra["sample_id"]).to_numpy()
    por_muestra["n_grupos"] = np.isfinite(mat.to_numpy(dtype=float)).sum(axis=0)
    etiqueta = {col: grupo for grupo, cols in grupos.items() for col in cols}
    por_muestra["grupo"] = por_muestra["sample_id"].map(etiqueta)
    volcar_tsv(por_muestra, salida / "suero_calidad_muestra.tsv", sello)
    pruebas = []
    for columna in ("eritrocitos", "plaquetas", "coagulacion", "piel", "fraccion_intensidad_no_humana", "n_precursores", "n_grupos"):
        if columna not in por_muestra.columns:
            continue
        bloques = {}
        for g in ("crc", "colonoscopy_negative", "healthy"):
            bloques[g] = por_muestra.loc[por_muestra["grupo"] == g, columna].dropna().to_numpy()
        validos = [v for v in bloques.values() if v.size >= 3]
        p = float(kruskal(*validos).pvalue) if len(validos) >= 2 else float("nan")
        fila = {"indicador": columna, "p_kruskal": p}
        for g, valores in bloques.items():
            fila[f"mediana_{g}"] = float(np.median(valores)) if valores.size else float("nan")
        pares = (("crc", "colonoscopy_negative"), ("crc", "healthy"), ("colonoscopy_negative", "healthy"))
        for ga, gb in pares:
            hl = _hodges(bloques[ga], bloques[gb])
            fila[f"hl_{ga}_menos_{gb}"] = hl["estimador"]
            fila[f"hl_ic_bajo_{ga}_{gb}"] = hl["ic_bajo"]
            fila[f"hl_ic_alto_{ga}_{gb}"] = hl["ic_alto"]
        pruebas.append(fila)
    volcar_tsv(pd.DataFrame(pruebas), salida / "suero_calidad_pruebas.tsv", sello)

    conteos = {}
    tablas_principal = {}
    for etiqueta_version, matriz, ajustar in (
        ("principal", principal, False),
        ("con_contaminantes", mat, False),
        ("sin_contaminantes_con_paneles", con_paneles, False),
        ("ajustado_coagulacion", principal, True),
    ):
        resumen, tablas = _k4_bloque(matriz, grupos, prereg, ajustar=ajustar, coagulacion=coag if ajustar else None)
        conteos[etiqueta_version] = resumen
        for nombre, tabla in tablas.items():
            nombre_fichero = f"k4_{nombre}.tsv" if etiqueta_version == "principal" else f"k4_{etiqueta_version}_{nombre}.tsv"
            volcar_tsv(tabla, salida / nombre_fichero, sello)
            if etiqueta_version == "principal":
                tablas_principal[nombre] = tabla
    genes = matriz_filtrada(tabla, valor="Genes.MaxLFQ", indice="Genes", q_de_grupo="GG.Q.Value") if "Genes.MaxLFQ" in tabla.columns and "GG.Q.Value" in tabla.columns else pd.DataFrame()
    if len(genes):
        mat_g, _ = _log2(genes)
        genes_fuera: set[str] = set()
        for _, fila in traza.loc[traza["contaminante"].astype(bool) | traza["contaminante_piel"].astype(bool) | (traza["panel"].fillna("").astype(str).str.len() > 0)].iterrows():
            genes_fuera |= _simbolos(fila["genes"])
        resumen_g, tablas_g = _k4_bloque(_recortar(mat_g, genes_fuera), _grupos(list(mat_g.columns)), prereg, ajustar=False, coagulacion=None)
        conteos["genes_maxlfq"] = resumen_g
        for nombre, tabla in tablas_g.items():
            volcar_tsv(tabla, salida / f"k4_genes_maxlfq_{nombre}.tsv", sello)
    else:
        conteos["genes_maxlfq"] = {"estado": "no_disponible", "motivo": "el report no trae Genes.MaxLFQ o GG.Q.Value"}
    volcar_json({"estado": "calculado", "escala": escala, "versiones": conteos}, salida / "k4_potencia.json", sello)

    nulos = {}
    aucs = []
    diferencias = None
    for etiqueta_version, matriz in (
        ("filtrada_humana", mat),
        ("sin_contaminantes", con_paneles),
        ("sin_contaminantes_ni_paneles", principal),
    ):
        print(f"[v3.2] K5 {etiqueta_version}: {int(prereg.get('k5_auc.nulo_del_maximo.n_permutaciones'))} permutaciones")
        auc, nulo, dif = _k5(matriz, grupos, prereg, etiqueta_version)
        aucs.append(auc)
        nulos[etiqueta_version] = nulo
        if etiqueta_version == "sin_contaminantes_ni_paneles":
            diferencias = dif
    volcar_tsv(pd.concat(aucs, ignore_index=True) if aucs else pd.DataFrame(), salida / "k5_auc.tsv", sello)
    if diferencias is not None:
        volcar_tsv(diferencias, salida / "k5_diferencia_auc.tsv", sello)
    volcar_json({"estado": "calculado", "nulos": nulos}, salida / "k5_auc.json", sello)

    alpha = float(prereg.get("inferencia.alpha"))
    umbral = float(prereg.get("k6_faltantes.umbral_observacion"))
    columnas = grupos.get("crc", []) + grupos.get("colonoscopy_negative", []) + grupos.get("healthy", [])
    etiquetas_k6 = (
        ["crc"] * len(grupos.get("crc", []))
        + ["colonoscopy_negative"] * len(grupos.get("colonoscopy_negative", []))
        + ["healthy"] * len(grupos.get("healthy", []))
    )
    bloque = principal[columnas].to_numpy(dtype=float) if all(c in principal.columns for c in columnas) else np.empty((0, 0))
    if bloque.size:
        fr = np.isfinite(bloque).mean(axis=1)
        mascara = fr < umbral
        sub = bloque[mascara]
        ids = principal.index.to_numpy()[mascara]
        if sub.size:
            presencia = fisher_presencia(sub, np.array(etiquetas_k6), grupo_a="crc", grupo_b="colonoscopy_negative", alpha=alpha)
            presencia.insert(0, "protein_id", ids)
            presencia.insert(1, "fraccion_observada", fr[mascara])
        else:
            presencia = pd.DataFrame()
        k6 = {"estado": "calculado", "n_por_debajo_del_umbral": int(mascara.sum()), "umbral": umbral}
    else:
        presencia = pd.DataFrame()
        k6 = {"estado": "no_disponible", "motivo": "la matriz principal no tiene las tres columnas de grupo"}
    volcar_tsv(presencia, salida / "k6_presencia.tsv", sello)
    chmp = traza.loc[traza["genes"].astype(str).str.contains("CHMP4A", regex=False)]
    k6["chmp4a"] = chmp.to_dict(orient="records")
    volcar_json(k6, salida / "k6_presencia.json", sello)

    precursores = tabla.groupby(tabla["Run"].astype(str)).size()
    med = {}
    for g, cols in grupos.items():
        vals = precursores.reindex(cols).dropna().to_numpy()
        med[g] = float(np.median(vals)) if vals.size else float("nan")
    confusiones = []
    for indicador in ("coagulacion", "piel"):
        fila = next((p for p in pruebas if p["indicador"] == indicador), None)
        if fila and np.isfinite(fila["p_kruskal"]) and fila["p_kruskal"] < float(prereg.get("k7_control_tecnico.alpha_exploratorio")):
            confusiones.append({
                "indicador": indicador,
                "lectura": "confusión preanalítica",
                "p_kruskal": fila["p_kruskal"],
                "hl_crc_menos_colonoscopy_negative": fila.get("hl_crc_menos_colonoscopy_negative"),
            })
    k7 = {
        "estado": "calculado",
        "precursores_min": int(precursores.min()) if len(precursores) else None,
        "precursores_max": int(precursores.max()) if len(precursores) else None,
        "precursores_mediana": med,
        "orden_de_corrida": "limitacion_declarada",
        "motivo_orden": "El nombre de la corrida no trae la fecha de adquisición. Los crudos no se descargaron. Queda la limitación de Leek 2010.",
        "confusion_preanalitica": confusiones,
        "escala": escala,
    }
    volcar_json(k7, salida / "k7_control_tecnico.json", sello)
    volcar_tsv(pd.DataFrame(pruebas), salida / "k7_pruebas.tsv", sello)

    tabla_paneles = None
    ruta_paneles = salida / "k4_sin_contaminantes_con_paneles_crc_frente_a_colonoscopia_negativa.tsv"
    if ruta_paneles.is_file():
        tabla_paneles = pd.read_csv(ruta_paneles, sep="\t")
    frase = _frase(conteos.get("principal", {}), conteos.get("sin_contaminantes_con_paneles", {}), traza, tabla_paneles)
    resumen = {
        "escala": escala,
        "n_grupos_report": int(traza["protein_group"].nunique()),
        "n_matriz_humana_fdr": int(len(mat)),
        "n_principal": int(len(principal)),
        "clases": conteos,
        "k5_nulos": {k: {"maximo_observado": v.get("maximo_observado"), "p": v.get("p")} for k, v in nulos.items()},
        "frase": frase,
    }
    volcar_json(resumen, salida / "suero_resumen.json", sello)
