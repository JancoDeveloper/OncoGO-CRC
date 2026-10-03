"""Contraste certificado v2: C1, C3, C7, C9 y C15.

No toca v1. Todo se escribe bajo ``results/v2/certified/``.

Qué cambia respecto a la corrida del 1 de octubre de 2026:

- **C1**: 10.000 permutaciones en CPTAC en vez de 100, y la enumeración completa
  de las 2^16 inversiones de signo en la réplica. El p deja de estar en su piso.
- **C3**: cada fracción sale con numerador, denominador y unidad. En v1 convivían
  226/688 y 119/702 como la misma métrica, y el 702 eran filas del cruce.
- **C7**: al ajustar una proteína de matriz extracelular, el índice estromal se
  recalcula sin ella. Además se ajusta por firmas de músculo liso, inmune y
  epitelial, todas preregistradas.
- **C9**: R1 estratificado por abundancia, por péptidos únicos, por dirección del
  efecto y por si sobrevive al ajuste de composición.
- **C15**: R1 arrastra la cobertura sobre los hallazgos. En v1 era 104/119, y 119
  es el 17,3 % de las 688 certificadas.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.stages.certified_inference import (
    MIN_COMPLETE,
    _cptac_contrast,
    _pride_contrast,
    _unique_peptide_counts,
)
from oncogo_crc.utils.go_utils import load_go_long
from oncogo_crc.utils.ids import effective_protein_mapping
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.moderated import moderate_table
from oncogo_crc.utils.quant_columns import columns_for_family, is_summary_protein
from oncogo_crc.utils.reproducibility_overlap import (
    cptac_pride_overlap_tables,
    load_gene_uniprot_from_mapping,
)
from oncogo_crc.v2 import prereg as prereg_mod
from oncogo_crc.v2.composicion import (
    construir_indice,
    covariables_leave_one_out,
    delta_por_par,
    regresion_multiple_por_proteina,
    regresion_por_proteina,
)
from oncogo_crc.v2.metricas import Fraccion, p_empirico
from oncogo_crc.v2.nulos import nulo_aleatorio, nulo_exacto

ECM_GO_TERM = "GO:0031012"


def _matriz_limpia(raw: pd.DataFrame) -> pd.DataFrame:
    """Quita las filas de resumen (Mean, Median, StdDev) y los duplicados."""
    mat = raw.loc[~raw["protein_id"].map(is_summary_protein)].copy()
    return mat.drop_duplicates("protein_id", keep="first").set_index("protein_id")


def _genes_de_termino(go_long: pd.DataFrame, termino: str) -> set[str]:
    return {
        str(g).strip().upper()
        for g in go_long.loc[go_long["go_id"].astype(str) == termino, "gene_symbol"]
    }


def _ajuste_por_indice(
    *,
    nombre: str,
    indice,
    pares: list[tuple[str, str]],
    diffs: np.ndarray,
    protein_ids: np.ndarray,
    certificadas: set[str],
    leave_one_out: bool,
    lfc: float,
    alpha: float,
) -> tuple[pd.DataFrame, dict]:
    """Ajusta el contraste por un índice de composición y mide qué sobrevive."""
    if leave_one_out:
        covariables, n_excluidas = covariables_leave_one_out(indice, pares, protein_ids)
    else:
        base = delta_por_par(indice.score, pares)
        covariables = np.tile(base, (len(protein_ids), 1))
        n_excluidas = 0
    coef, s2, df, stdev, n_used, pendiente, pendiente_t = regresion_por_proteina(
        diffs, covariables, min_pares=MIN_COMPLETE, evaluar_en=0.0
    )
    tabla, prior_df, _s02 = moderate_table(
        protein_ids=protein_ids,
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_used,
        lfc=lfc,
        alpha=alpha,
        stat_method=f"moderated_t_treat_paired_ajustado_{nombre}",
    )
    tabla["pendiente_composicion"] = pendiente
    tabla["pendiente_composicion_t"] = pendiente_t
    tabla["certificada_sin_ajustar"] = tabla["protein_id"].astype(str).isin(certificadas)
    cert = tabla.loc[tabla["certificada_sin_ajustar"]]
    sobreviven = int(cert["significant"].sum())
    base_delta = delta_por_par(indice.score, pares)
    frac = Fraccion(
        nombre=f"R4_composicion_{nombre}",
        numerador=sobreviven,
        denominador=len(certificadas),
        unidad_denominador="proteinas_unicas",
        nota=(
            "Proteínas certificadas que siguen pasando treat y BH al evaluar el efecto "
            "en diferencia de composición cero. Es una cota, no una pureza medida."
        ),
    )
    resumen = {
        "indice": nombre,
        "leave_one_out": leave_one_out,
        "n_proteinas_en_firma": indice.m,
        "n_proteinas_con_indice_propio": n_excluidas,
        "delta_media": float(np.nanmean(base_delta)),
        "delta_sd": float(np.nanstd(base_delta, ddof=1)),
        "evaluado_en_delta": 0.0,
        "prior_df": float(prior_df),
        "pendiente_dependiente_n": int((cert["pendiente_composicion_t"].abs() > 2.0).sum()),
        "mediana_abs_log2fc_ajustado": float(cert["log2fc"].abs().median()),
        "n_significativas_total_ajustado": int(tabla["significant"].sum()),
        "fraccion": frac.to_dict(),
    }
    return tabla, resumen


def _ajuste_multiindice(
    *,
    indices: list,
    nombres: list[str],
    pares: list[tuple[str, str]],
    protein_ids: np.ndarray,
    diffs: np.ndarray,
    certificadas: set[str],
    lfc: float,
    alpha: float,
) -> tuple[pd.DataFrame, dict]:
    """Ajusta por varios compartimentos a la vez y evalúa el efecto con todos en cero."""
    covs = []
    for idx in indices:
        base = delta_por_par(idx.score, pares)
        covs.append(np.tile(base, (len(protein_ids), 1)))
    coef, s2, df, stdev, n_used, pendientes = regresion_multiple_por_proteina(
        diffs, covs, min_pares=MIN_COMPLETE
    )
    tabla, prior_df, _s02 = moderate_table(
        protein_ids=protein_ids,
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_used,
        lfc=lfc,
        alpha=alpha,
        stat_method="moderated_t_treat_paired_ajustado_conjunto_multiindice",
    )
    for j, nombre in enumerate(nombres):
        tabla[f"pendiente_{nombre}"] = pendientes[:, j]
    tabla["certificada_sin_ajustar"] = tabla["protein_id"].astype(str).isin(certificadas)
    cert = tabla.loc[tabla["certificada_sin_ajustar"]]
    frac = Fraccion(
        nombre="R4_composicion_conjunto_multiindice",
        numerador=int(cert["significant"].sum()),
        denominador=len(certificadas),
        unidad_denominador="proteinas_unicas",
        nota=(
            "Efecto evaluado con todos los compartimentos medidos en diferencia cero. "
            "Es la cota más estricta: quita la varianza compartida por los índices."
        ),
    )
    resumen = {
        "indice": "conjunto_multiindice",
        "componentes": list(nombres),
        "leave_one_out": False,
        "n_proteinas_en_firma": sum(i.m for i in indices),
        "n_proteinas_con_indice_propio": 0,
        "delta_media": None,
        "delta_sd": None,
        "evaluado_en_delta": 0.0,
        "prior_df": float(prior_df),
        "pendiente_dependiente_n": None,
        "mediana_abs_log2fc_ajustado": float(cert["log2fc"].abs().median()),
        "n_significativas_total_ajustado": int(tabla["significant"].sum()),
        "pendiente_mediana_por_componente": {
            nombre: float(np.nanmedian(pendientes[:, j])) for j, nombre in enumerate(nombres)
        },
        "fraccion": frac.to_dict(),
    }
    return tabla, resumen


def _control_permutado(
    *,
    indices: list,
    pares: list[tuple[str, str]],
    protein_ids: np.ndarray,
    diffs: np.ndarray,
    certificadas: set[str],
    lfc: float,
    alpha: float,
    n_permutaciones: int,
    semilla: int,
) -> dict:
    """Control negativo de C7 (enmienda E6): ¿cuánto se pierde solo por ajustar?

    Permuta el orden de los pares en cada covariable. La distribución de la
    covariable no cambia, pero deja de corresponder al par al que pertenece. Lo
    que el modelo pierde aquí es pérdida de potencia, no confusión.
    """
    rng = np.random.default_rng(semilla)
    bases = [delta_por_par(idx.score, pares) for idx in indices]
    n_prot = len(protein_ids)
    cert_mask = pd.Series(protein_ids).astype(str).isin(certificadas).to_numpy()
    fracciones: list[float] = []
    for _ in range(n_permutaciones):
        covs = [np.tile(rng.permutation(b), (n_prot, 1)) for b in bases]
        coef, s2, df, stdev, n_used, _p = regresion_multiple_por_proteina(
            diffs, covs, min_pares=MIN_COMPLETE
        )
        tabla, _d0, _s02 = moderate_table(
            protein_ids=protein_ids,
            coef=coef,
            s2=s2,
            df=df,
            stdev_unscaled=stdev,
            n_obs=n_used,
            lfc=lfc,
            alpha=alpha,
            stat_method="control_permutado",
        )
        sobreviven = int((tabla["significant"].to_numpy() & cert_mask).sum())
        fracciones.append(sobreviven / len(certificadas))
    arr = np.asarray(fracciones)
    return {
        "n_permutaciones": int(n_permutaciones),
        "semilla": int(semilla),
        "fraccion_mediana": float(np.median(arr)),
        "fraccion_p2_5": float(np.percentile(arr, 2.5)),
        "fraccion_p97_5": float(np.percentile(arr, 97.5)),
        "nota": (
            "Covariables con el orden de los pares permutado. Conservan su "
            "distribución y pierden su relación con el par. Lo que se cae aquí es "
            "pérdida de potencia por gastar grados de libertad, no confusión."
        ),
    }


def _r1_estratificado(
    overlap: pd.DataFrame,
    *,
    alpha: float,
    n_certificadas: int,
    peptidos: dict[str, int] | None,
    sobreviven_composicion: set[str],
    intensidad_pride: pd.Series | None,
) -> dict:
    """C9: la concordancia de dirección, partida por los estratos preregistrados."""
    padj = pd.to_numeric(overlap["padj_cptac"], errors="coerce")
    sig = overlap.loc[padj < alpha].copy()
    a = pd.to_numeric(sig["log2fc_cptac"], errors="coerce")
    b = pd.to_numeric(sig["log2fc_pride"], errors="coerce")
    sig["_a"], sig["_b"] = a, b
    evaluable = np.isfinite(a) & np.isfinite(b) & (a != 0) & (b != 0)
    ev = sig.loc[evaluable].drop_duplicates("protein_id_cptac").copy()
    ev["mismo_signo"] = np.sign(ev["_a"]) == np.sign(ev["_b"])

    global_frac = Fraccion(
        nombre="R1_direccion",
        numerador=int(ev["mismo_signo"].sum()),
        denominador=int(len(ev)),
        unidad_denominador="proteinas_unicas",
        universo=n_certificadas,
        universo_nombre="n_certificadas",
        nota=(
            "El denominador son los hallazgos que la réplica alcanza a medir, no los "
            "hallazgos. La cobertura dice qué parte de los certificados entra."
        ),
    )

    def _estrato(frame: pd.DataFrame, nombre: str) -> dict:
        return Fraccion(
            nombre=nombre,
            numerador=int(frame["mismo_signo"].sum()),
            denominador=int(len(frame)),
            unidad_denominador="proteinas_unicas",
        ).to_dict()

    estratos: dict[str, list[dict]] = {}

    # Dirección del efecto en tumor.
    estratos["direccion_del_efecto"] = [
        _estrato(ev.loc[ev["_a"] > 0], "sube_en_tumor"),
        _estrato(ev.loc[ev["_a"] < 0], "baja_en_tumor"),
    ]

    # Sobrevive el ajuste por composición.
    pid = ev["protein_id_cptac"].astype(str)
    estratos["sobrevive_ajuste_composicion"] = [
        _estrato(ev.loc[pid.isin(sobreviven_composicion)], "sobrevive"),
        _estrato(ev.loc[~pid.isin(sobreviven_composicion)], "se_cae"),
    ]

    # Péptidos únicos en la réplica.
    if peptidos:
        uni = ev["_join_uni"].astype(str).map(lambda v: peptidos.get(v))
        n_pep = pd.to_numeric(pd.Series(uni.to_numpy(), index=ev.index), errors="coerce")
        cortes = [1, 2, 5, 10, np.inf]
        filas = []
        for lo, hi in zip(cortes[:-1], cortes[1:]):
            sel = ev.loc[(n_pep >= lo) & (n_pep < hi)]
            etiqueta = f"peptidos_{int(lo)}_a_{'inf' if hi == np.inf else int(hi) - 1}"
            filas.append(_estrato(sel, etiqueta))
        filas.append(_estrato(ev.loc[n_pep.isna()], "peptidos_sin_dato"))
        estratos["n_peptidos_unicos"] = filas

    # Cuartiles de intensidad en la réplica.
    if intensidad_pride is not None and len(intensidade := intensidad_pride.dropna()) > 0:
        inten = ev["protein_id_pride"].astype(str).map(intensidade.to_dict())
        inten = pd.to_numeric(pd.Series(inten.to_numpy(), index=ev.index), errors="coerce")
        validos = inten.dropna()
        filas = []
        if len(validos) >= 8:
            bordes = np.nanpercentile(validos, [0, 25, 50, 75, 100])
            for q in range(4):
                lo, hi = bordes[q], bordes[q + 1]
                sel = ev.loc[(inten >= lo) & (inten <= hi)] if q == 3 else ev.loc[(inten >= lo) & (inten < hi)]
                filas.append(_estrato(sel, f"intensidad_q{q + 1}"))
        filas.append(_estrato(ev.loc[inten.isna()], "intensidad_sin_dato"))
        estratos["cuartiles_de_intensidad_pride"] = filas

    return {"global": global_frac.to_dict(), "estratos": estratos}


def _intensidad_mediana_pride(path) -> pd.Series | None:
    """Mediana de log2 de las intensidades LFQ positivas, por proteína."""
    try:
        raw = pd.read_csv(path, sep="\t", low_memory=False)
    except OSError:
        return None
    raw = raw.loc[~raw["protein_id"].map(is_summary_protein)].drop_duplicates("protein_id")
    mat = raw.set_index("protein_id")
    vals = mat.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        logged = np.where(vals > 0, np.log2(vals), np.nan)
    return pd.Series(np.nanmedian(logged, axis=1), index=mat.index.astype(str))


def _pride_pares(settings: Settings, *, lfc: float, alpha: float):
    """Brazo pareado tumor frente a normal del mismo paciente, para el nulo exacto.

    En PXD034575 las columnas se llaman ``LFQ intensity CRC<n><M|N|P>``: M es el
    tumor (maligno), N el tejido normal y P el pólipo. No es ``T``.
    """
    from oncogo_crc.stages.translational_axes import _paired_arm, _pride_log2

    ds = settings.datasets["datasets"]["pride_crc_meta"]
    log2mat, pacientes = _pride_log2(settings.resolve(ds["matrix_path"]))
    tabla, n_pares = _paired_arm(
        log2mat, pacientes, target="M", label="tumor", lfc=lfc, alpha=alpha
    )
    pares = [(v["M"], v["N"]) for v in pacientes.values() if "M" in v and "N" in v]
    if not pares:
        raise RuntimeError(
            "PXD034575: no se encontró ningún par M/N. Revise el patrón de columnas LFQ."
        )
    diffs = (
        log2mat[[a for a, _ in pares]].to_numpy() - log2mat[[b for _, b in pares]].to_numpy()
    )
    return tabla, diffs, n_pares


def run(settings: Settings) -> None:
    prereg = prereg_mod.cargar(settings.root)
    alpha = float(prereg.get("inferencia.alpha"))
    lfc = float(prereg.get("inferencia.treat_lfc"))
    out = settings.resolve("results") / "v2" / "certified"
    out.mkdir(parents=True, exist_ok=True)

    cptac_ds = settings.datasets["datasets"]["cptac_colorectal"]
    pride_ds = settings.datasets["datasets"]["pride_crc_meta"]
    raw = pd.read_csv(settings.resolve(cptac_ds["matrix_path"]), sep="\t", low_memory=False)
    meta = pd.read_csv(settings.sample_metadata_path("cptac_colorectal"), sep="\t")

    cptac, diffs, pares, prior_df = _cptac_contrast(
        raw,
        meta,
        family=str(prereg.get("inferencia.familia_columnas_cptac")),
        stat_method="moderated_t_treat_paired_unshared_logratio",
        lfc=lfc,
        alpha=alpha,
    )
    certificadas = set(cptac.loc[cptac["significant"], "protein_id"].astype(str))
    n_cert = len(certificadas)
    print(f"[v2-certified] CPTAC: pares={len(pares)} d0={prior_df:.2f} certificadas={n_cert}")
    write_tsv(cptac, out / "cptac_unshared_treat.tsv")

    pride_meta = pd.read_csv(settings.sample_metadata_path("pride_crc_meta"), sep="\t")
    pride = _pride_contrast(
        settings.resolve(pride_ds["matrix_path"]), pride_meta, lfc=lfc, alpha=alpha
    )
    write_tsv(pride, out / "pride_lfq_treat.tsv")

    mapping = effective_protein_mapping(settings)
    gene_to_uni = load_gene_uniprot_from_mapping(mapping, exclude_ambiguous=True)
    overlap, resumen_cruce = cptac_pride_overlap_tables(
        cptac, pride, gene_to_uni=gene_to_uni, alpha=alpha
    )
    write_tsv(overlap, out / "cptac_pride_overlap.tsv")

    # ------------------------------------------------------------------ C1
    n_perm = int(prereg.get("c1_permutaciones.cptac_n_permutaciones"))
    semilla = prereg.semilla("r5_permutacion_cptac")
    print(f"[v2-certified] C1: {n_perm} permutaciones de signo en CPTAC (semilla {semilla})")
    nulo_cptac = nulo_aleatorio(
        diffs,
        lfc=lfc,
        alpha=alpha,
        min_pares=MIN_COMPLETE,
        n_permutaciones=n_perm,
        semilla=semilla,
    )
    p_cptac = p_empirico(nulo_cptac["b_al_menos_observado"], nulo_cptac["n_permutaciones"])
    write_tsv(
        pd.DataFrame({"n_significativas": nulo_cptac["nulo"]}),
        out / "r5_nulo_cptac.tsv",
    )
    print(
        f"[v2-certified] C1 CPTAC: observado={nulo_cptac['observado']} "
        f"mediana nula={nulo_cptac['nulo_mediana']:.0f} p={p_cptac['p']:.6f} "
        f"(piso {p_cptac['p_minimo_alcanzable']:.6f}, tope={p_cptac['es_tope']})"
    )

    pride_pareado, diffs_pride, n_pares_pride = _pride_pares(settings, lfc=lfc, alpha=alpha)
    write_tsv(pride_pareado, out / "pride_paired_treat.tsv")
    nulo_pride = None
    p_pride = None
    if bool(prereg.get("c1_permutaciones.pride_exacto")) and n_pares_pride <= 20:
        print(
            f"[v2-certified] C1: enumerando 2^{n_pares_pride} inversiones de signo en la réplica"
        )
        nulo_pride = nulo_exacto(diffs_pride, lfc=lfc, alpha=alpha, min_pares=MIN_COMPLETE)
        p_pride = p_empirico(nulo_pride["b_al_menos_observado"], nulo_pride["n_permutaciones"])
        write_tsv(
            pd.DataFrame({"n_significativas": nulo_pride["nulo"]}),
            out / "r5_nulo_pride_exacto.tsv",
        )
        print(
            f"[v2-certified] C1 réplica exacta: observado={nulo_pride['observado']} "
            f"de {nulo_pride['n_permutaciones']} inversiones, p={p_pride['p']:.6f}"
        )

    # ------------------------------------------------------------------ C7
    go_path = settings.go_gene_long_path()
    ajustes: list[dict] = []
    sobreviven_principal: set[str] = set()
    if go_path.exists():
        go_long = load_go_long(go_path)
        cols = columns_for_family(
            [c for c in raw.columns if c != "protein_id"],
            str(prereg.get("inferencia.familia_columnas_cptac")),
        )
        mat = _matriz_limpia(raw)
        protein_ids = cptac["protein_id"].to_numpy()
        firmas = {"ecm": _genes_de_termino(go_long, ECM_GO_TERM)}
        for nombre, genes in (prereg.get("c7_composicion.firmas") or {}).items():
            firmas[nombre] = {str(g).strip().upper() for g in genes}

        indices = {}
        for nombre, genes in firmas.items():
            idx = construir_indice(mat, cols, genes, nombre=nombre)
            if idx is None:
                print(f"[v2-certified] C7: firma {nombre} sin proteínas suficientes, se omite")
                continue
            indices[nombre] = idx

        # Método principal preregistrado: matriz extracelular con exclusión propia.
        principal = str(prereg.get("c7_composicion.metodo_principal"))
        plan = [("ecm_leave_one_out", "ecm", True), ("ecm_sin_exclusion", "ecm", False)]
        plan += [(n, n, False) for n in indices if n != "ecm"]
        # Ajuste conjunto: matriz extracelular y músculo liso comparten proteínas de
        # pared intestinal, así que por separado no dicen cuál compartimento manda.
        conjunto = [n for n in ("ecm", "musculo_liso", "inmune", "epitelial") if n in indices]
        if len(conjunto) >= 2:
            tabla_m, resumen_m = _ajuste_multiindice(
                indices=[indices[n] for n in conjunto],
                nombres=conjunto,
                pares=pares,
                protein_ids=protein_ids,
                diffs=diffs,
                certificadas=certificadas,
                lfc=lfc,
                alpha=alpha,
            )
            resumen_m["es_principal"] = False
            ajustes.append(resumen_m)
            write_tsv(tabla_m, out / "cptac_treat_ajustado_conjunto_multiindice.tsv")
            fm = resumen_m["fraccion"]
            print(
                f"[v2-certified] C7 conjunto_multiindice ({'+'.join(conjunto)}): "
                f"sobreviven {fm['numerador']}/{fm['denominador']} = {fm['valor']:.4f}"
            )
            control = _control_permutado(
                indices=[indices[n] for n in conjunto],
                pares=pares,
                protein_ids=protein_ids,
                diffs=diffs,
                certificadas=certificadas,
                lfc=lfc,
                alpha=alpha,
                n_permutaciones=int(prereg.get("c7_composicion.control_permutado.n_permutaciones")),
                semilla=prereg.semilla("c7_control_permutado"),
            )
            resumen_m["control_permutado"] = control
            print(
                f"[v2-certified] C7 control permutado: mediana {control['fraccion_mediana']:.4f} "
                f"(IC {control['fraccion_p2_5']:.4f}–{control['fraccion_p97_5']:.4f}). "
                f"La caída atribuible a composición es "
                f"{control['fraccion_mediana'] - fm['valor']:+.4f}"
            )

        for etiqueta, clave, loo in plan:
            if clave not in indices:
                continue
            tabla, resumen = _ajuste_por_indice(
                nombre=etiqueta,
                indice=indices[clave],
                pares=pares,
                diffs=diffs,
                protein_ids=protein_ids,
                certificadas=certificadas,
                leave_one_out=loo,
                lfc=lfc,
                alpha=alpha,
            )
            resumen["es_principal"] = etiqueta == principal
            ajustes.append(resumen)
            write_tsv(tabla, out / f"cptac_treat_ajustado_{etiqueta}.tsv")
            if etiqueta == principal:
                sobreviven_principal = set(
                    tabla.loc[tabla["certificada_sin_ajustar"] & tabla["significant"], "protein_id"]
                    .astype(str)
                    .tolist()
                )
            f = resumen["fraccion"]
            print(
                f"[v2-certified] C7 {etiqueta}: sobreviven {f['numerador']}/{f['denominador']} "
                f"= {f['valor']:.4f} (firma={resumen['n_proteinas_en_firma']}, "
                f"con índice propio={resumen['n_proteinas_con_indice_propio']})"
            )
    else:
        print(f"[v2-certified] sin anotación GO en {go_path}: C7 no se calcula")

    # -------------------------------------------------------------- C3, C9, C15
    peptidos = _unique_peptide_counts(settings)
    r2 = None
    if peptidos:
        minimo = 1
        ok = 0
        for gen in cptac.loc[cptac["significant"], "protein_id"].astype(str):
            uni = gene_to_uni.get(gen.strip().upper())
            if uni is not None and (peptidos.get(uni) or 0) >= minimo:
                ok += 1
        r2 = Fraccion(
            nombre="R2_detectabilidad",
            numerador=ok,
            denominador=n_cert,
            unidad_denominador="proteinas_unicas",
            nota=(
                f"UniProt con al menos {minimo} péptido único en el proteinGroups de la "
                "réplica. No demuestra que la misma secuencia se haya visto en las dos "
                "plataformas."
            ),
        ).to_dict()
        print(f"[v2-certified] C3 R2: {ok}/{n_cert} = {r2['valor']:.4f}")

    intensidad = _intensidad_mediana_pride(settings.resolve(pride_ds["matrix_path"]))
    r1 = _r1_estratificado(
        overlap,
        alpha=alpha,
        n_certificadas=n_cert,
        peptidos=peptidos,
        sobreviven_composicion=sobreviven_principal,
        intensidad_pride=intensidad,
    )
    g = r1["global"]
    print(
        f"[v2-certified] C15 R1: {g['numerador']}/{g['denominador']} = {g['valor']:.4f}; "
        f"cobertura {g['cobertura']:.4f} de {g['universo']} certificadas"
    )
    filas_estratos = []
    for estrato, filas in r1["estratos"].items():
        for fila in filas:
            filas_estratos.append({"estrato": estrato, **fila})
    write_tsv(pd.DataFrame(filas_estratos), out / "r1_estratificado.tsv")

    # ----------------------------------------------------------------- auditoría
    auditoria = {
        "version": 2,
        "metodo": "moderated_t_treat_paired_unshared_logratio",
        "alpha": alpha,
        "treat_lfc": lfc,
        "semilla_base": int(prereg.get("semillas.base")),
        "min_pares_observados": MIN_COMPLETE,
        "cptac": {
            "pares": len(pares),
            "prior_df": float(prior_df),
            "n_certificadas": n_cert,
            "n_evaluadas": int(cptac["p_value"].notna().sum()),
            "n_filas": int(len(cptac)),
        },
        "pride": {
            "n_significativas_no_pareado": int(pride["significant"].sum()),
            "n_evaluadas_no_pareado": int(pride["p_value"].notna().sum()),
            "n_significativas_pareado": int(pride_pareado["significant"].sum()),
            "n_pares": int(n_pares_pride),
        },
        "c1_nulo_cptac": {
            **{k: v for k, v in nulo_cptac.items() if k != "nulo"},
            "p_empirico": p_cptac,
            "tabla_nulo": "r5_nulo_cptac.tsv",
        },
        "c1_nulo_pride_exacto": (
            {
                **{k: v for k, v in nulo_pride.items() if k != "nulo"},
                "p_empirico": p_pride,
                "tabla_nulo": "r5_nulo_pride_exacto.tsv",
            }
            if nulo_pride is not None
            else {"estado": "no_calculado", "motivo": "la réplica no tiene brazo pareado enumerable"}
        ),
        "c3_r2_detectabilidad": r2,
        "c15_r1_direccion": r1,
        "c7_composicion": {
            "metodo_principal": str(prereg.get("c7_composicion.metodo_principal")),
            "evaluado_en_delta": float(prereg.get("c7_composicion.evaluar_efecto_en")),
            "lenguaje_obligatorio": str(prereg.get("c7_composicion.lenguaje_obligatorio")),
            "ajustes": ajustes,
        },
        "resumen_cruce": resumen_cruce.to_dict(orient="records"),
    }
    (out / "inference_audit_v2.json").write_text(
        json.dumps(auditoria, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[v2-certified] escrito en {out}")
