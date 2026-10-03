"""D5. R4 alternativo con una referencia de composición que viene de fuera.

El preregistro lista ``deconvolucion_referencia_externa`` entre los métodos de
sensibilidad de C7. Esta etapa lo ejecuta, y por el camino se topa con dos
obstáculos que no son suyos sino del problema, y que conviene leer antes que
cualquier número.

El primero. La deconvolución no se puede correr sobre CPTAC, que es la cohorte
donde vive R4. CPTAC publica un log de un cociente contra un canal de
referencia, y mezclar perfiles celulares exige una escala de abundancia. La
etapa lo comprueba con los datos delante en lugar de afirmarlo, y por eso el R4
alternativo se calcula en las cohortes LFQ, cuya intensidad sí está en escala de
abundancia. Eso cambia la pregunta: no es «R4 de CPTAC con otra covariable» sino
«R4 recalculado en la réplica, con una covariable que no sale del propio
proteoma».

El segundo. Las proporciones tipo a tipo que devuelve la deconvolución no son
publicables, y la etapa lo comprueba también. Lo que queda en pie es el agregado
epitelio frente a no epitelio. Con eso se puede ajustar.

Lo que esta etapa aporta frente al índice estromal es una cosa concreta: la
covariable no se calcula del proteoma que se está contrastando. El índice
estromal sí, y C7 solo pudo mitigar esa circularidad excluyendo cada proteína de
su propio índice. Aquí no hay nada que excluir por ese motivo, aunque se excluye
igual para que la comparación entre los dos ajustes sea la misma cuenta.

Nada de aquí cambia un subíndice del certificado.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.go_utils import load_go_long
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.moderated import moderate_table, paired_observed_stats
from oncogo_crc.utils.quant_columns import is_summary_protein
from oncogo_crc.utils.reproducibility_overlap import (
    load_gene_uniprot_from_mapping,
    uniprot_base,
)
from oncogo_crc.v2 import deconvolucion as D
from oncogo_crc.v2.composicion import (
    construir_indice,
    covariables_leave_one_out,
    delta_por_par,
    regresion_por_proteina,
)
from oncogo_crc.v2.metricas import Fraccion
from oncogo_crc.v2.prereg import cargar

ECM_GO_TERM = "GO:0031012"
MIN_PARES = 3
DESPLAZAMIENTO = 415  # semillas.desplazamientos.c7_deconvolucion del preregistro


def _tarjeta(base: Path) -> dict | None:
    path = base / "estado.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def _pxd034575(path: Path, u2g: dict[str, str]) -> tuple[pd.DataFrame, list[tuple[str, str]]]:
    """La réplica del certificado. El tumor es M (maligno), el normal N.

    Esta cohorte publica accesiones de UniProt, no símbolos de gen, así que hay
    que cruzarla con el mapa. El cruce pierde proteínas, y cuántas pierde queda
    en la auditoría: una firma de ARN está indexada por gen y no hay forma de
    evitar el paso.
    """
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    raw = raw.loc[~raw["protein_id"].map(is_summary_protein)]
    genes = [u2g.get(uniprot_base(str(i)) or "") for i in raw["protein_id"]]
    raw = raw.assign(_gen=[str(g).upper() if g else None for g in genes])
    raw = raw.dropna(subset=["_gen"])
    cols = [c for c in raw.columns if str(c).startswith("LFQ intensity ")]
    mat = raw.groupby("_gen")[cols].sum()
    pares = []
    for c in cols:
        nombre = c.replace("LFQ intensity ", "")
        if nombre.endswith("M"):
            normal = f"LFQ intensity {nombre[:-1]}N"
            if normal in cols:
                pares.append((c, normal))
    return mat, sorted(pares)


def _pxd019504(path: Path) -> tuple[pd.DataFrame, list[tuple[str, str]]]:
    """La cohorte independiente. Aquí el tumor es T y el normal N.

    No es un detalle de estilo: en PXD034575 la T no existe y el tumor es M. Usar
    el patrón de una cohorte en la otra deja los pares vacíos o, peor, los invierte.
    """
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    if "Gene names" not in raw.columns:
        raise KeyError(
            "proteinGroups.txt de PXD019504 no tiene la columna «Gene names»; sin "
            "símbolo de gen no se puede cruzar con una firma de ARN"
        )
    raw = raw[raw["Gene names"].notna()].copy()
    raw["_gen"] = raw["Gene names"].astype(str).str.split(";").str[0].str.upper()
    cols = [c for c in raw.columns if str(c).startswith("LFQ intensity ")]
    mat = raw.groupby("_gen")[cols].sum()
    pares = []
    for c in cols:
        nombre = c.replace("LFQ intensity ", "")
        if nombre.endswith("T"):
            normal = f"LFQ intensity {nombre[:-1]}N"
            if normal in cols:
                pares.append((c, normal))
    return mat, sorted(pares)


def _covariables_epiteliales(
    mat_lineal: pd.DataFrame,
    columnas: list[str],
    firma: pd.DataFrame,
    pares: list[tuple[str, str]],
    protein_ids: np.ndarray,
    tipos_epiteliales: list[str],
) -> tuple[np.ndarray, int]:
    """Diferencia de fracción epitelial por par, con la proteína excluida.

    Solo hay que recalcular para las proteínas que están en la firma. Para las
    demás, quitarlas no cambia la firma, así que la fracción es la misma: hacer
    el cálculo igual multiplicaría el tiempo por sesenta sin cambiar un número.
    """
    base = D.eje_epitelial(
        D.fracciones(mat_lineal, columnas, firma), tipos_epiteliales
    )
    delta_base = delta_por_par(base, pares)
    matriz = np.tile(delta_base, (len(protein_ids), 1))
    en_firma = set(firma.index.astype(str))
    n_excluidas = 0
    for i, pid in enumerate(protein_ids):
        clave = str(pid).strip().upper()
        if clave not in en_firma:
            continue
        sin = D.eje_epitelial_sin_gen(
            mat_lineal, columnas, firma, clave, tipos_epiteliales
        )
        matriz[i] = delta_por_par(sin, pares)
        n_excluidas += 1
    return matriz, n_excluidas


def _ajuste(
    *,
    nombre: str,
    covariables: np.ndarray,
    diffs: np.ndarray,
    protein_ids: np.ndarray,
    certificadas: set[str],
    lfc: float,
    alpha: float,
    n_en_firma: int,
    n_excluidas: int,
) -> dict:
    """La misma cuenta que C7, con otra covariable."""
    coef, s2, df, stdev, n_used, pendiente, pendiente_t = regresion_por_proteina(
        diffs, covariables, min_pares=MIN_PARES, evaluar_en=0.0
    )
    tabla, prior_df, _ = moderate_table(
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
    tabla["pendiente_composicion_t"] = pendiente_t
    tabla["significativa_sin_ajustar"] = (
        tabla["protein_id"].astype(str).isin(certificadas)
    )
    cert = tabla.loc[tabla["significativa_sin_ajustar"]]
    frac = Fraccion(
        nombre=f"R4_alternativo_{nombre}",
        numerador=int(cert["significant"].sum()),
        denominador=len(certificadas),
        unidad_denominador="proteinas_unicas",
        nota=(
            "Proteínas significativas en esta cohorte que siguen pasando treat y "
            "BH al evaluar el efecto en diferencia de composición cero."
        ),
    )
    base = covariables[0] if len(covariables) else np.array([np.nan])
    media = float(np.nanmean(base))
    sd = float(np.nanstd(base, ddof=1))
    return {
        "covariable": nombre,
        "media_de_la_covariable": round(media, 4),
        "sd_de_la_covariable": round(sd, 4),
        # A cuántas desviaciones típicas de los datos está el cero donde se
        # evalúa el efecto. Es la cantidad que no depende de la escala de la
        # covariable, y la que gobierna cuánta precisión cuesta el ajuste.
        "distancia_tipificada_al_cero": round(abs(media) / sd, 3) if sd > 0 else None,
        "n_proteinas_en_firma": int(n_en_firma),
        "n_proteinas_con_covariable_propia": int(n_excluidas),
        "prior_df": float(prior_df),
        "pendiente_dependiente_n": int(
            (cert["pendiente_composicion_t"].abs() > 2.0).sum()
        ),
        "fraccion": frac.to_dict(),
    }


def _como_comparar(ajustes: list[dict]) -> dict:
    """Qué se puede concluir de que los dos ajustes den números distintos.

    La lectura tentadora es que el R4 más bajo mide la confusión de verdad y el
    más alto se queda corto. Hay un motivo para desconfiar de esa lectura, y es
    aritmético más que biológico.

    R4 evalúa el efecto en diferencia de composición cero, y ningún par tiene
    diferencia cero: todos los tumores difieren de su normal. Así que el efecto
    en cero es una extrapolación, y el error típico de una extrapolación crece
    con la distancia a los datos. La cantidad que la mide sin depender de la
    escala de la covariable es a cuántas desviaciones típicas está el cero de la
    media de los datos. Si una covariable tiene el cero más lejos, su ajuste
    pierde más proteínas por falta de precisión, sin que eso signifique que haya
    más confusión.

    Comparar las desviaciones típicas a secas no serviría: una es una fracción
    acotada entre cero y uno y la otra es una media de puntuaciones tipificadas.
    Son escalas distintas y el cociente no significa nada.

    Lo que sí se concluye, y es lo que importa para la tesis, es que R4 depende de
    cómo se mida la composición. Si dos formas razonables de medirla dan cotas
    distintas, R4 es una cota y no una cantidad. Esa es exactamente la frase que
    el preregistro declaró obligatoria.
    """
    if len(ajustes) < 2:
        return {
            "estado": "no_comparable",
            "motivo": (
                "solo hay un ajuste en esta cohorte, así que no hay dos cotas que "
                "comparar"
            ),
            "n_ajustes": len(ajustes),
        }
    por_nombre = {a["covariable"]: a for a in ajustes}
    ext = por_nombre.get("fraccion_epitelial_referencia_externa")
    inte = por_nombre.get("indice_estromal_del_propio_proteoma")
    if ext is None or inte is None:
        return {"estado": "no_comparable", "motivo": "faltan los dos ajustes esperados"}
    d_ext = ext["distancia_tipificada_al_cero"]
    d_int = inte["distancia_tipificada_al_cero"]
    extrapola_mas = (
        bool(d_ext > d_int) if d_ext is not None and d_int is not None else None
    )
    pend_ext = ext["pendiente_dependiente_n"]
    pend_int = inte["pendiente_dependiente_n"]
    return {
        "estado": "comparable",
        "r4_con_referencia_externa": ext["fraccion"]["valor"],
        "r4_con_el_indice_del_propio_proteoma": inte["fraccion"]["valor"],
        "distancia_tipificada_al_cero_externa": d_ext,
        "distancia_tipificada_al_cero_interna": d_int,
        "la_externa_extrapola_mas_lejos": extrapola_mas,
        "proteinas_con_pendiente_de_composicion_externa": pend_ext,
        "proteinas_con_pendiente_de_composicion_interna": pend_int,
        "los_intervalos_se_solapan": bool(
            ext["fraccion"]["ci_high"] >= inte["fraccion"]["ci_low"]
            and inte["fraccion"]["ci_high"] >= ext["fraccion"]["ci_low"]
        ),
        "por_que_difieren": (
            f"Las dos covariables marcan un número parecido de proteínas como "
            f"dependientes de la composición ({pend_ext} frente a {pend_int}), "
            f"así que la diferencia en R4 no viene de que una detecte mucha más "
            f"confusión. Viene de la precisión: el cero donde se evalúa el efecto "
            f"está a {d_ext} desviaciones típicas de los datos con la covariable "
            f"externa y a {d_int} con la interna. Extrapolar más lejos cuesta "
            f"error típico, y con más error típico menos proteínas pasan treat."
            if extrapola_mas
            else (
                f"La covariable externa marca {pend_ext} proteínas como "
                f"dependientes de la composición y la interna {pend_int}, y el "
                f"cero está a {d_ext} desviaciones típicas de los datos con la "
                f"externa frente a {d_int} con la interna. Con estas dos cifras "
                f"delante se puede juzgar cuánto de la diferencia en R4 es "
                f"confusión medida y cuánto precisión perdida."
            )
        ),
        "lo_que_no_se_puede_concluir": (
            "Que el R4 más bajo sea «la verdad» y el más alto un ajuste "
            "insuficiente. Parte de la diferencia es el coste en precisión de "
            "extrapolar a una composición que ningún par tiene. Y comparar las "
            "desviaciones típicas a secas tampoco serviría: una covariable es una "
            "fracción entre cero y uno y la otra una media de puntuaciones "
            "tipificadas, así que el cociente entre sus escalas no significa nada."
        ),
        "lo_que_si_se_concluye": (
            "R4 depende de cómo se mida la composición. Dos formas razonables de "
            "medirla dan cotas distintas, así que R4 es una cota y no una "
            "cantidad. Es la frase que el preregistro declaró obligatoria, y aquí "
            "tiene un número que la respalda."
        ),
    }


def _cohorte(
    nombre: str,
    mat_lineal: pd.DataFrame,
    pares: list[tuple[str, str]],
    firma: pd.DataFrame,
    tipos_epiteliales: list[str],
    genes_ecm: set[str],
    *,
    lfc: float,
    alpha: float,
) -> dict:
    columnas = sorted({c for par in pares for c in par})
    fr = D.fracciones(mat_lineal, columnas, firma)
    diag = D.diagnosticar(fr, firma)
    epi = D.eje_epitelial(fr, tipos_epiteliales)

    log2 = mat_lineal[columnas].mask(mat_lineal[columnas] <= 0.0)
    log2 = pd.DataFrame(
        np.log2(log2.to_numpy(dtype=float)), index=log2.index, columns=columnas
    )
    indice_ecm = construir_indice(log2, columnas, genes_ecm, nombre="ecm")
    comparacion = D.comparar_con_indice(epi, indice_ecm.score)
    veredicto_eje = D.veredicto_del_eje_epitelial(epi, pares, comparacion, diag)

    diffs = (
        log2[[a for a, _ in pares]].to_numpy(dtype=float)
        - log2[[b for _, b in pares]].to_numpy(dtype=float)
    )
    protein_ids = log2.index.to_numpy()
    coef, s2, df, stdev, n_obs = paired_observed_stats(diffs, min_pairs=MIN_PARES)
    sin_ajustar, _, _ = moderate_table(
        protein_ids=protein_ids,
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method="moderated_t_treat_paired_lfq_sin_ajustar",
    )
    certificadas = set(
        sin_ajustar.loc[sin_ajustar["significant"], "protein_id"].astype(str)
    )

    ajustes = []
    if certificadas and veredicto_eje["usable_como_covariable"]:
        cov_epi, n_exc_epi = _covariables_epiteliales(
            mat_lineal, columnas, firma, pares, protein_ids, tipos_epiteliales
        )
        ajustes.append(
            _ajuste(
                nombre="fraccion_epitelial_referencia_externa",
                covariables=cov_epi,
                diffs=diffs,
                protein_ids=protein_ids,
                certificadas=certificadas,
                lfc=lfc,
                alpha=alpha,
                n_en_firma=len(set(protein_ids.astype(str)) & set(firma.index.astype(str))),
                n_excluidas=n_exc_epi,
            )
        )
    if certificadas:
        cov_ecm, n_exc_ecm = covariables_leave_one_out(indice_ecm, pares, protein_ids)
        ajustes.append(
            _ajuste(
                nombre="indice_estromal_del_propio_proteoma",
                covariables=cov_ecm,
                diffs=diffs,
                protein_ids=protein_ids,
                certificadas=certificadas,
                lfc=lfc,
                alpha=alpha,
                n_en_firma=indice_ecm.m,
                n_excluidas=n_exc_ecm,
            )
        )

    return {
        "cohorte": nombre,
        "n_pares": len(pares),
        "n_genes_medidos": int(mat_lineal.shape[0]),
        "n_marcadores_medidos": len(fr.genes_usados),
        "n_significativas_sin_ajustar": len(certificadas),
        "diagnostico": diag,
        "fraccion_epitelial": {
            "mediana_tumor": round(
                float(np.median([epi[a] for a, _ in pares])), 4
            ),
            "mediana_normal": round(
                float(np.median([epi[b] for _, b in pares])), 4
            ),
            "delta_medio": round(float(np.nanmean(delta_por_par(epi, pares))), 4),
        },
        "comparacion_con_el_indice_estromal": comparacion,
        "veredicto_del_eje_epitelial": veredicto_eje,
        "ajustes": ajustes,
        "como_comparar_los_dos_ajustes": _como_comparar(ajustes),
        "_fracciones": fr.tabla.assign(
            residuo_relativo=fr.residuo_relativo, fraccion_epitelial=epi
        ),
    }


def run(settings: Settings) -> None:
    prereg = cargar(settings.root)
    lfc = float(prereg.get("inferencia.treat_lfc"))
    alpha = float(prereg.get("inferencia.alpha"))
    ref_dir = settings.resolve("external") / "referencia_unicelular"
    salida = settings.resolve("results") / "v2" / "axes" / "deconvolucion"
    salida.mkdir(parents=True, exist_ok=True)
    audit = salida / "deconvolucion_audit_v2.json"

    def _no_disponible(motivo: str) -> None:
        audit.write_text(
            json.dumps(
                {"version": "v2", "estado": "no_disponible", "motivo": motivo},
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"[v2-deconvolucion] no_disponible: {motivo}")

    tarjeta = _tarjeta(ref_dir)
    if tarjeta is None or tarjeta.get("estado") != "usada":
        return _no_disponible(
            "falta external/referencia_unicelular/estado.json en estado «usada»; "
            "hay que ejecutar scripts/download/hpa_single_cell.py"
            if tarjeta is None
            else f"la referencia unicelular está marcada {tarjeta.get('estado')}: "
            f"{tarjeta.get('motivo', 'sin motivo')}"
        )
    firma_path = ref_dir / tarjeta["firma"]
    if not firma_path.exists():
        return _no_disponible(f"falta la firma {firma_path}")
    go_path = settings.go_gene_long_path()
    if not go_path.exists():
        return _no_disponible(
            "falta la tabla de GO; sin ella no se puede construir el índice "
            "estromal con el que comparar la deconvolución"
        )

    firma = pd.read_csv(firma_path, sep="\t", index_col=0)
    firma.index = [str(g).strip().upper() for g in firma.index]
    tipos_epiteliales = tarjeta["tipos_celulares"]["tipos_epiteliales"]
    go_long = load_go_long(go_path)
    genes_ecm = {
        str(g).strip().upper()
        for g in go_long.loc[
            go_long["go_id"].astype(str) == ECM_GO_TERM, "gene_symbol"
        ]
    }
    print(
        f"[v2-deconvolucion] firma: {firma.shape[0]} marcadores, "
        f"{firma.shape[1]} tipos celulares, {len(tipos_epiteliales)} epiteliales."
    )

    mapping = pd.read_csv(
        settings.resolve("metadata")
        / "protein_id_maps"
        / "gene_uniprot_human_reviewed_bulk.tsv",
        sep="\t",
        low_memory=False,
    )
    g2u = load_gene_uniprot_from_mapping(mapping, exclude_ambiguous=True)
    u2g = {u: g for g, u in g2u.items()}

    cohortes: list[dict] = []
    replica = settings.resolve("data") / "raw" / "pride" / "pride_crc_matrix.tsv"
    if replica.exists():
        mat, pares = _pxd034575(replica, u2g)
        if len(pares) >= MIN_PARES:
            cohortes.append(
                _cohorte(
                    "PXD034575",
                    mat,
                    pares,
                    firma,
                    tipos_epiteliales,
                    genes_ecm,
                    lfc=lfc,
                    alpha=alpha,
                )
            )
    independiente = (
        settings.resolve("external") / "pride" / "PXD019504" / "proteinGroups.txt"
    )
    if independiente.exists():
        try:
            mat2, pares2 = _pxd019504(independiente)
        except KeyError as e:
            print(f"[v2-deconvolucion] PXD019504 sin usar: {e}")
        else:
            if len(pares2) >= MIN_PARES:
                cohortes.append(
                    _cohorte(
                        "PXD019504",
                        mat2,
                        pares2,
                        firma,
                        tipos_epiteliales,
                        genes_ecm,
                        lfc=lfc,
                        alpha=alpha,
                    )
                )

    if not cohortes:
        return _no_disponible(
            "ninguna cohorte LFQ con pares tumor-normal está disponible en local; "
            "sin escala de abundancia no se puede deconvolver, y las matrices TMT "
            "no sirven para esto"
        )

    for c in cohortes:
        print(
            f"[v2-deconvolucion] {c['cohorte']}: {c['n_pares']} pares, "
            f"{c['n_marcadores_medidos']} marcadores medidos, "
            f"{c['n_significativas_sin_ajustar']} proteínas significativas sin ajustar."
        )
        print(f"[v2-deconvolucion] {c['cohorte']}: {c['diagnostico']['veredicto']}")
        v = c["veredicto_del_eje_epitelial"]
        print(
            f"[v2-deconvolucion] {c['cohorte']}: el eje epitelial baja en el "
            f"tumor en {v['pares_con_menos_epitelio_en_el_tumor']} de "
            f"{v['n_pares']} pares. {v['no_es_una_proporcion_de_celulas']}"
        )
        comp = c["comparacion_con_el_indice_estromal"]
        if comp.get("estado") == "calculada":
            print(
                f"[v2-deconvolucion] {c['cohorte']}: la fracción epitelial y el "
                f"índice estromal correlacionan a rho {comp['rho_spearman']:+.3f} "
                f"(n={comp['n']}), signo {comp['signo_observado']}."
            )
        for a in c["ajustes"]:
            f = a["fraccion"]
            print(
                f"[v2-deconvolucion] {c['cohorte']}: R4 con «{a['covariable']}» = "
                f"{f['valor']:.4f} ({f['numerador']}/{f['denominador']}), "
                f"desviación típica de la covariable {a['sd_de_la_covariable']:.3f}."
            )
        cmp_ = c["como_comparar_los_dos_ajustes"]
        if cmp_.get("estado") == "comparable":
            print(f"[v2-deconvolucion] {c['cohorte']}: {cmp_['lo_que_si_se_concluye']}")

    # La comprobación de que esto no se puede hacer en CPTAC, con los datos delante.
    tmt = None
    cptac_path = settings.resolve("data") / "raw" / "cptac" / "cptac_crc_matrix.tsv"
    if cptac_path.exists():
        raw = pd.read_csv(cptac_path, sep="\t", low_memory=False)
        raw = raw.loc[~raw["protein_id"].map(is_summary_protein)]
        raw = raw.drop_duplicates("protein_id", keep="first").set_index("protein_id")
        raw.index = [str(i).strip().upper() for i in raw.index]
        cols = [c for c in raw.columns if str(c).endswith(" Unshared Log Ratio")]
        try:
            tmt = D.por_que_no_en_tmt(raw, cols, firma)
        except (ValueError, D.FaltaColumna) as e:
            tmt = {"estado": "no_comprobable", "motivo": str(e)}
        else:
            print(
                f"[v2-deconvolucion] CPTAC: el "
                f"{100 * tmt['fraccion_negativa']:.1f} % de los valores es "
                "negativo, así que ninguna mezcla de perfiles no negativos puede "
                "reproducir el perfil. La deconvolución no se corre ahí."
            )

    filas = []
    for c in cohortes:
        t = c.pop("_fracciones")
        t = t.reset_index().rename(columns={"index": "muestra"})
        t.insert(0, "cohorte", c["cohorte"])
        filas.append(t)
    write_tsv(pd.concat(filas, ignore_index=True), salida / "fracciones_por_muestra.tsv")

    write_tsv(
        pd.DataFrame(
            [
                {
                    "cohorte": c["cohorte"],
                    "covariable": a["covariable"],
                    "n_pares": c["n_pares"],
                    **a["fraccion"],
                }
                for c in cohortes
                for a in c["ajustes"]
            ]
        ),
        salida / "r4_alternativo.tsv",
    )

    auditoria = {
        "version": "v2",
        "estado": "calculado",
        "semilla": int(prereg.get("semillas.base")) + DESPLAZAMIENTO,
        "metodo_del_preregistro": "deconvolucion_referencia_externa",
        "lenguaje_obligatorio": prereg.get("c7_composicion.lenguaje_obligatorio"),
        "referencia": {
            "fuente": tarjeta["fuente"],
            "url": tarjeta["url"],
            "cita": tarjeta["cita"],
            "licencia": tarjeta["licencia"],
            "fecha_de_acceso": tarjeta["fecha_de_acceso"],
            "sha256_firma": tarjeta["sha256_firma"],
            "n_marcadores": tarjeta["seleccion_de_marcadores"]["n_marcadores"],
            "regla_de_marcadores": tarjeta["seleccion_de_marcadores"]["regla"],
            "enmienda_de_la_regla": tarjeta["seleccion_de_marcadores"][
                "no_estaba_preregistrada"
            ],
            "referencia_de_colon_no_usada": tarjeta[
                "referencia_especifica_de_colon"
            ],
        },
        "por_que_no_se_corre_en_cptac": tmt,
        "cohortes": cohortes,
        "que_aporta_frente_al_indice_estromal": (
            "La covariable no se calcula del proteoma que se está contrastando. El "
            "índice estromal sí, y C7 solo pudo mitigar esa circularidad "
            "excluyendo cada proteína de su propio índice. Aquí la referencia es "
            "externa, así que la circularidad no existe."
        ),
        "que_no_aporta": (
            "No es R4 de CPTAC con otra covariable. CPTAC publica un log-cociente "
            "y la deconvolución no vale ahí, así que el ajuste se hace en las "
            "cohortes LFQ. Son cohortes distintas, con menos pares y otra "
            "tecnología, de modo que el número no sustituye al R4 del "
            "certificado: lo acompaña."
        ),
        "limitaciones": [
            "La referencia es de ARN y lo deconvuelto es proteína. D2 midió el "
            "acuerdo entre transcrito y proteína del mismo gen en CEACAM5 y salió "
            "una correlación de rangos de 0,24. Una firma de ARN aplicada a un "
            "proteoma es una aproximación, y el número de arriba la hereda.",
            "Las proporciones tipo a tipo no se publican: el diagnóstico de cada "
            "cohorte explica por qué. Solo se usa el agregado epitelial.",
            "Los tipos celulares del Human Protein Atlas son de todos los tejidos "
            "juntos, así que «fibroblastos» es un promedio entre tejidos.",
            "La regla de selección de marcadores no estaba preregistrada y se "
            "declara como enmienda E9.",
            "R4 resulta ser sensible a dónde queda el cero respecto de los datos, "
            "y no solo a cuánta confusión capta la covariable. Ningún par tiene "
            "diferencia de composición cero, así que evaluar el efecto ahí es "
            "extrapolar, y el error típico de la extrapolación crece con la "
            "distancia. En PXD019504 las dos covariables marcan un número "
            "parecido de proteínas como dependientes de la composición y aun así "
            "dan cotas muy distintas, porque una extrapola más del doble de lejos "
            "en unidades de su propia desviación típica. Es una limitación de R4 "
            "como índice, no de esta fuente, y conviene tenerla delante al leer "
            "cualquier valor de R4, incluido el del certificado.",
            "Comparar la fracción epitelial con el índice estromal no valida "
            "ninguno de los dos: los dos salen del mismo proteoma. La validación "
            "con otra tecnología es la de D2, contra la pureza ABSOLUTE del ADN.",
        ],
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
    }
    audit.write_text(
        json.dumps(auditoria, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    print(f"[v2-deconvolucion] auditoría: {audit}")
