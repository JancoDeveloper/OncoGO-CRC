"""Cláusula de vías de H1 recalculada (C4) y su nulo de etiquetas (C13).

Qué cambia respecto de v1. La cláusula de vía de v1 se declaró satisfecha con
cuatro términos y un intervalo de 0,398 a 1,0. Aquí se exige un mínimo de diez
términos evaluables, la prueba de conjunto es competitiva en lugar de una
sobrerrepresentación de Fisher, se repite en tres colecciones independientes y
la cláusula solo se sostiene si el límite inferior del intervalo supera el
umbral. Después se aplica la tabla de verdad preregistrada de H1.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.go_utils import load_go_long
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.reproducibility_overlap import (
    load_gene_uniprot_from_mapping,
    uniprot_base,
)
from oncogo_crc.utils.ids import effective_protein_mapping
from oncogo_crc.stages.certified_inference import _cptac_contrast
from oncogo_crc.stages.v2_certified import _pride_pares
from oncogo_crc.v2 import vias
from oncogo_crc.v2.conjuntos import correlacion_intergenica
from oncogo_crc.v2.prereg import cargar


def _estadisticos_por_gen(tabla: pd.DataFrame, *, a_gen=None) -> pd.Series:
    """t moderado por símbolo de gen. Si hay duplicados, gana el |t| mayor."""
    t = tabla[["protein_id", "t_stat"]].dropna(subset=["t_stat"]).copy()
    if a_gen is None:
        t["gen"] = t["protein_id"].astype(str).str.upper()
    else:
        t["gen"] = t["protein_id"].map(a_gen)
    t = t.loc[t["gen"].astype(str).str.len() > 0]
    t["_abs"] = t["t_stat"].abs()
    t = t.sort_values("_abs", ascending=False, kind="stable").drop_duplicates("gen")
    return t.set_index("gen")["t_stat"].astype(float)


def _rho_y_df(tabla: pd.DataFrame, diffs: np.ndarray) -> tuple[float, float]:
    """Correlación intergénica de los residuos y grados de libertad del contraste.

    El preregistro exige estimar rho de los residuos. Dejarlo en cero haría
    VIF = 1 y la prueba saldría anticonservadora: justo el problema que la
    corrección por correlación existe para evitar.
    """
    residuos = diffs - np.nanmean(diffs, axis=1, keepdims=True)
    rho = correlacion_intergenica(np.nan_to_num(residuos, nan=0.0))
    df = float(np.nanmedian(tabla["df_total"])) if "df_total" in tabla.columns else 10.0
    if not np.isfinite(df) or df <= 0:
        df = 10.0
    return rho, df


def run(settings: Settings) -> None:
    prereg = cargar(settings.root)
    alpha = float(prereg.get("inferencia.alpha"))
    n_min = int(prereg.get("c4_vias.n_min_terminos"))
    min_prot = int(prereg.get("c4_vias.min_proteinas_por_termino"))
    umbral = float(prereg.get("c4_vias.umbral_h1_via_min"))
    salida = settings.resolve("results") / "v2" / "vias"
    salida.mkdir(parents=True, exist_ok=True)

    cert = settings.resolve("results") / "v2" / "certified"
    cptac = pd.read_csv(cert / "cptac_unshared_treat.tsv", sep="\t")
    pride = pd.read_csv(cert / "pride_paired_treat.tsv", sep="\t")

    # Las matrices de diferencias por par hacen falta para estimar rho.
    cptac_ds = settings.datasets["datasets"]["cptac_colorectal"]
    raw = pd.read_csv(settings.resolve(cptac_ds["matrix_path"]), sep="\t", low_memory=False)
    meta = pd.read_csv(settings.sample_metadata_path("cptac_colorectal"), sep="\t")
    _tabla_c, diffs_c, _pares_c, _prior = _cptac_contrast(
        raw,
        meta,
        family=str(prereg.get("inferencia.familia_columnas_cptac")),
        stat_method="moderated_t_treat_paired_unshared_logratio",
        lfc=float(prereg.get("inferencia.treat_lfc")),
        alpha=alpha,
    )
    _tabla_p, diffs_p, _n_pares = _pride_pares(
        settings, lfc=float(prereg.get("inferencia.treat_lfc")), alpha=alpha
    )

    mapping = effective_protein_mapping(settings)
    gene_to_uni = load_gene_uniprot_from_mapping(mapping, exclude_ambiguous=True)
    uni_to_gene = {u: g for g, u in gene_to_uni.items()}

    go_path = settings.go_gene_long_path()
    go_long = load_go_long(go_path) if go_path.exists() else None
    colecciones = vias.cargar_colecciones(
        go_long=go_long, base_externa=settings.resolve("external") / "vias"
    )
    for nombre, c in colecciones.items():
        print(f"[v2-vias] colección {nombre}: {c.estado}, {c.n_terminos} términos. {c.motivo}")

    est_c = _estadisticos_por_gen(cptac)
    est_p = _estadisticos_por_gen(
        pride, a_gen=lambda v: uni_to_gene.get(uniprot_base(str(v)) or "", "")
    )
    print(f"[v2-vias] universo: CPTAC {est_c.size} genes, réplica {est_p.size} genes.")

    rho_c, df_c = _rho_y_df(cptac, diffs_c)
    rho_p, df_p = _rho_y_df(pride, diffs_p)
    rho_cfg = prereg.get("c4_vias.correlacion_intergenica")
    print(f"[v2-vias] rho intergénica: CPTAC {rho_c:.4f}, réplica {rho_p:.4f}")

    por_coleccion: dict[str, dict] = {}
    clausulas: dict[str, str] = {}
    for nombre, c in colecciones.items():
        if c.estado != "usada":
            por_coleccion[nombre] = {
                "estado": "no_disponible",
                "motivo": c.motivo,
                "n_terminos_en_la_coleccion": c.n_terminos,
            }
            clausulas[nombre] = "no_evaluable"
            continue
        tabla_c = vias.prueba_por_cohorte(
            est_c, c.terminos, rho=rho_c, df=df_c, min_proteinas=min_prot, alpha=alpha
        )
        tabla_p = vias.prueba_por_cohorte(
            est_p, c.terminos, rho=rho_p, df=df_p, min_proteinas=min_prot, alpha=alpha
        )
        write_tsv(tabla_c, salida / f"competitiva_cptac_{nombre}.tsv")
        write_tsv(tabla_p, salida / f"competitiva_replica_{nombre}.tsv")
        conc = vias.concordancia_entre_cohortes(
            tabla_c, tabla_p, coleccion=nombre, n_min_terminos=n_min, umbral=umbral
        )
        detalle = conc.pop("detalle", None)
        if detalle is not None and len(detalle):
            write_tsv(detalle, salida / f"concordancia_{nombre}.tsv")
        nulo = vias.nulo_de_etiquetas(
            est_c,
            c.terminos,
            rho=rho_c,
            df=df_c,
            min_proteinas=min_prot,
            n_permutaciones=int(prereg.get("c13_nulo_vias.n_permutaciones")),
            semilla=prereg.semilla("r5_permutacion_vias"),
        )
        # El detalle por término va a su tabla, no al JSON.
        por_termino = nulo.pop("por_termino", {})
        if por_termino:
            write_tsv(
                pd.DataFrame(
                    {"termino": list(por_termino), "p_nulo_etiquetas": list(por_termino.values())}
                ).sort_values("p_nulo_etiquetas", kind="stable"),
                salida / f"nulo_etiquetas_cptac_{nombre}.tsv",
            )
        por_coleccion[nombre] = {
            "fuente": c.fuente,
            "n_terminos_en_la_coleccion": c.n_terminos,
            "n_terminos_probados_cptac": int(len(tabla_c)),
            "n_terminos_probados_replica": int(len(tabla_p)),
            "rho_cptac": rho_c,
            "rho_replica": rho_p,
            "concordancia_entre_cohortes": conc,
            "nulo_etiquetas_de_gen_c13": nulo,
        }
        clausulas[nombre] = conc["estado"]
        print(
            f"[v2-vias] {nombre}: {conc.get('n_terminos_significativos_descubrimiento')} "
            f"significativos en CPTAC, {conc['n_terminos_evaluados']} evaluables en la "
            f"réplica → {conc['estado']}"
        )

    # La cláusula de vía global. Preregistrado: tres colecciones, y una conclusión
    # que solo aparece en una no es robusta. Se exige que al menos una sea
    # evaluable; si las evaluables discrepan, la cláusula se rechaza.
    evaluables = {k: v for k, v in clausulas.items() if v in ("sostiene", "rechaza")}
    if not evaluables:
        clausula_via = "no_evaluable"
        motivo_via = (
            "Ninguna de las tres colecciones alcanzó el mínimo de "
            f"{n_min} términos evaluables en la réplica."
        )
    elif all(v == "sostiene" for v in evaluables.values()):
        clausula_via = "sostiene"
        motivo_via = f"Las {len(evaluables)} colecciones evaluables coinciden."
    else:
        clausula_via = "rechaza"
        motivo_via = (
            "Las colecciones evaluables no coinciden o alguna rechaza: "
            + ", ".join(f"{k}={v}" for k, v in evaluables.items())
        )
    margenes = {
        k: por_coleccion[k]["concordancia_entre_cohortes"].get(
            "margen_del_ic_inferior_sobre_el_umbral"
        )
        for k in evaluables
    }
    robustez = {
        "n_colecciones_pedidas": len(colecciones),
        "n_colecciones_evaluables": len(evaluables),
        "una_sola_coleccion": len(evaluables) == 1,
        "margen_del_ic_inferior_sobre_el_umbral": margenes,
        "holgado": {
            k: por_coleccion[k]["concordancia_entre_cohortes"].get("holgado")
            for k in evaluables
        },
        "advertencia": (
            "La cláusula se apoya en una sola de las tres colecciones y el límite "
            "inferior de su intervalo supera el umbral por poco. No es un resultado "
            "robusto: es el mínimo que la regla preregistrada acepta."
            if len(evaluables) == 1
            and any(m is not None and m <= 0.05 for m in margenes.values())
            else ""
        ),
        "cuello_de_botella": (
            "Lo que deja no evaluables a GO BP y a Hallmark no es el resultado, es "
            "la cobertura: la réplica solo mapea una fracción de sus proteínas a "
            "símbolo de gen, así que casi ningún término llega al mínimo de "
            "proteínas medidas allí."
        ),
    }
    if robustez["advertencia"]:
        print(f"[v2-vias] ADVERTENCIA: {robustez['advertencia']}")

    # La cláusula de proteína viene del stage certificado de v2, no se recalcula.
    auditoria_v2 = json.loads(
        (cert / "inference_audit_v2.json").read_text(encoding="utf-8")
    )
    r1 = (auditoria_v2.get("c15_r1_direccion") or {}).get("global") or {}
    umbral_prot = float(prereg.get("hipotesis.h1.clausula_proteina.umbral"))
    if r1.get("ci_high") is None:
        clausula_proteina = "no_evaluable"
        motivo_prot = "R1 no es calculable."
    elif float(r1["ci_high"]) < umbral_prot:
        clausula_proteina = "sostiene"
        motivo_prot = (
            f"El intervalo de R1 queda entero por debajo de {umbral_prot}: "
            f"{r1['ci_low']:.3f}–{r1['ci_high']:.3f}."
        )
    else:
        clausula_proteina = "rechaza"
        motivo_prot = (
            f"El intervalo de R1 llega hasta {float(r1['ci_high']):.3f}, que no está "
            f"por debajo de {umbral_prot}."
        )

    h1 = vias.tabla_verdad_h1(clausula_proteina, clausula_via)
    h1["motivo_clausula_proteina"] = motivo_prot
    h1["motivo_clausula_via"] = motivo_via
    h1["r1_usado"] = r1
    print(f"[v2-vias] H1: proteína={clausula_proteina}, vía={clausula_via} → {h1['veredicto']}")

    auditoria = {
        "version": "v2",
        "alpha": alpha,
        "n_min_terminos": n_min,
        "min_proteinas_por_termino": min_prot,
        "umbral_concordancia_via": umbral,
        "correlacion_intergenica": rho_cfg,
        "prueba": prereg.get("c4_vias.prueba"),
        "colecciones_pedidas": prereg.get("c4_vias.colecciones"),
        "por_coleccion": por_coleccion,
        "universo": {
            "n_genes_cptac": int(est_c.size),
            "n_genes_replica": int(est_p.size),
            "n_proteinas_en_la_tabla_de_la_replica": int(len(pride)),
            "fraccion_de_la_replica_con_simbolo_de_gen": (
                float(est_p.size / len(pride)) if len(pride) else None
            ),
        },
        "clausula_via": {
            "estado": clausula_via,
            "motivo": motivo_via,
            "por_coleccion": clausulas,
            "robustez": robustez,
        },
        "h1": h1,
        "h2": {
            "estado": prereg.get("hipotesis.h2.estado"),
            "motivo": prereg.get("hipotesis.h2.motivo"),
            "requisito_para_evaluar": prereg.get("hipotesis.h2.requisito_para_evaluar"),
        },
        "comparacion_con_v1": {
            "v1_concordancia_via": 1.0,
            "v1_n_terminos": 4,
            "v1_ci_low": 0.3976353643835142,
            "v1_estado": "sostiene",
            "lectura": (
                "v1 dio la cláusula por satisfecha con cuatro términos y un intervalo "
                "que bajaba a 0,398. Con cuatro términos ese intervalo no distingue "
                "concordancia perfecta de una moneda, así que la cláusula no estaba "
                "satisfecha: estaba sin evaluar."
            ),
        },
        "limitaciones": [
            "La prueba competitiva compara cada conjunto con el resto de proteínas "
            "medidas. El resultado depende del universo medido, que no es el proteoma.",
            "Las dos cohortes usan tecnologías distintas y universos distintos. Un "
            "término puede ser no evaluable en la réplica solo por cobertura.",
            "El nulo de etiquetas de gen rompe la pertenencia pero conserva la "
            "correlación entre proteínas, así que no es un nulo de la correlación.",
            "Ningún resultado de este fichero sirve para decidir sobre un paciente.",
        ],
    }
    (salida / "vias_audit_v2.json").write_text(
        json.dumps(auditoria, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    print(f"[v2-vias] escrito en {salida}")
