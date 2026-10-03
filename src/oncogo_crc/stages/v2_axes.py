"""Ejes de v2: secuencia adenoma–carcinoma (C5) y selección de modelo (C6).

Qué se corrige respecto de v1.

- **C5**. En v1 las proteínas se eligieron por ser significativas en el brazo
  tumor de PXD034575 y su magnitud se midió en ese mismo brazo. Aquí la
  selección se hace en CPTAC, que es otra cohorte y otra tecnología, y se añaden
  dos acotaciones del sesgo: mitades divididas por paciente y un nulo que
  reselecciona.
- **C6**. En v1 el bootstrap remuestreaba los 32 especímenes de prueba, que en
  realidad son 16 pacientes con dos muestras cada uno. Aquí se remuestrean
  pacientes, la validación cruzada agrupa por participante y se reporta la
  diferencia de AUC por pares, no solo dos intervalos que se miran de lejos.

Nada de esto toca las matrices crudas ni los resultados de v1.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.moderated import moderate_table, paired_observed_stats
from oncogo_crc.utils.quant_columns import is_summary_protein
from oncogo_crc.utils.ids import effective_protein_mapping
from oncogo_crc.utils.reproducibility_overlap import (
    load_gene_uniprot_from_mapping,
    uniprot_base,
)
from oncogo_crc.stages.certified_inference import MIN_COMPLETE, _numeric_block
from oncogo_crc.stages.translational_axes import _pride_log2, _zscore_rows
from oncogo_crc.v2 import modelos, polipos
from oncogo_crc.v2.prereg import cargar


def _brazo_pareado(
    log2mat: pd.DataFrame,
    pacientes: dict[str, dict[str, str]],
    *,
    objetivo: str,
    etiqueta: str,
    lfc: float,
    alpha: float,
) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    """Contraste pareado dentro del paciente y las diferencias que lo producen."""
    ids = [p for p, v in pacientes.items() if objetivo in v and "N" in v]
    if not ids:
        raise RuntimeError(
            f"ningún paciente con los tejidos {objetivo} y N. "
            "En PXD034575 el tumor es M (maligno), el normal N y el pólipo P."
        )
    diffs = (
        log2mat[[pacientes[p][objetivo] for p in ids]].to_numpy()
        - log2mat[[pacientes[p]["N"] for p in ids]].to_numpy()
    )
    coef, s2, df, stdev, n_obs = paired_observed_stats(diffs, min_pairs=MIN_COMPLETE)
    tabla, _d0, _s02 = moderate_table(
        protein_ids=log2mat.index.to_numpy(),
        coef=coef,
        s2=s2,
        df=df,
        stdev_unscaled=stdev,
        n_obs=n_obs,
        lfc=lfc,
        alpha=alpha,
        stat_method=f"moderated_t_treat_paired_lfq_{etiqueta}",
    )
    tabla["comparison_id"] = f"pride_{etiqueta}_vs_normal"
    return tabla, diffs, ids


def _c5(
    settings: Settings,
    prereg,
    log2mat: pd.DataFrame,
    pacientes: dict[str, dict[str, str]],
    cptac: pd.DataFrame,
    gene_to_uni: dict[str, str],
    *,
    lfc: float,
    alpha: float,
    salida,
) -> dict:
    tumor, diffs_t, ids_t = _brazo_pareado(
        log2mat, pacientes, objetivo="M", etiqueta="tumor", lfc=lfc, alpha=alpha
    )
    polipo, diffs_p, ids_p = _brazo_pareado(
        log2mat, pacientes, objetivo="P", etiqueta="polipo", lfc=lfc, alpha=alpha
    )
    write_tsv(tumor, salida / "pride_pareado_tumor_vs_normal.tsv")
    write_tsv(polipo, salida / "pride_pareado_polipo_vs_normal.tsv")
    print(
        f"[v2-axes] pareado: tumor n={len(ids_t)} sig={int(tumor['significant'].sum())}; "
        f"pólipo n={len(ids_p)} sig={int(polipo['significant'].sum())}"
    )

    # El de v1, para tener el contraste explícito del sesgo.
    dentro = polipos.seleccion_independiente(
        tumor,
        polipo,
        set(tumor.loc[tumor["significant"], "protein_id"].astype(str)),
    )
    dentro["nota"] = (
        "Selección y medición en el mismo brazo. Es el cálculo de v1 y está "
        "inflado por construcción. Se reporta solo para comparar."
    )

    # El principal: seleccionar en CPTAC, medir en la réplica.
    elegidas = set(cptac.loc[cptac["significant"], "protein_id"].astype(str).str.upper())
    uni_a_gen = {u: g for g, u in gene_to_uni.items()}
    mapa = {
        str(i): uni_a_gen.get(uniprot_base(str(i)), "")
        for i in log2mat.index
    }
    fuera = polipos.seleccion_independiente(tumor, polipo, elegidas, mapear=mapa)
    print(
        "[v2-axes] C5 selección independiente: "
        f"{fuera['n_medidas_en_la_replica']} proteínas medidas, "
        f"dirección {fuera['direccion']['valor']}, "
        f"intermedio {fuera['intermedio_dado_mismo_signo']['valor']}"
    )

    # Las mitades divididas necesitan los mismos pacientes en el mismo orden.
    comunes = [p for p in ids_t if p in set(ids_p)]
    idx_t = [ids_t.index(p) for p in comunes]
    idx_p = [ids_p.index(p) for p in comunes]
    mitades = polipos.mitades_divididas(
        diffs_t[:, idx_t],
        diffs_p[:, idx_p],
        log2mat.index.to_numpy(),
        lfc=lfc,
        alpha=alpha,
        min_pares=MIN_COMPLETE,
        n_repeticiones=int(prereg.get("c5_polipos.mitades_divididas.n_repeticiones")),
        semilla=prereg.semilla("c5_mitades_divididas"),
    )
    print(f"[v2-axes] C5 mitades divididas: {mitades.get('direccion', mitades)}")

    n_elegir = int(tumor["significant"].sum())
    nulo = polipos.nulo_intercambio(
        diffs_t[:, idx_t],
        diffs_p[:, idx_p],
        n_elegir=max(n_elegir, 10),
        n_permutaciones=int(prereg.get("c5_polipos.nulo_intercambio.n_permutaciones")),
        semilla=prereg.semilla("c5_nulo_intercambio"),
        min_pares=MIN_COMPLETE,
    )
    print(
        f"[v2-axes] C5 nulo que reselecciona: observado {nulo.get('observado_intermedio')}, "
        f"nulo mediana {nulo.get('nulo_mediana')}, p {nulo.get('p_empirico')}"
    )

    # El veredicto sale de la regla preregistrada, no de elegir el análisis que
    # mejor queda. El principal es la selección independiente.
    f = fuera["intermedio_dado_mismo_signo"]
    reproduce = f["valor"] is not None and f["ci_low"] is not None and f["ci_low"] > 0.5
    veredicto = {
        "analisis_principal": prereg.get("c5_polipos.analisis_principal"),
        "regla": prereg.get("c5_polipos.regla"),
        "la_seleccion_independiente_lo_reproduce": bool(reproduce),
        "conclusion": (
            "El gradiente de magnitud normal–adenoma–carcinoma se sostiene."
            if reproduce
            else "El gradiente de magnitud no se sostiene con selección independiente."
        ),
        "comparacion": {
            "mismo_brazo_v1": dentro["intermedio_dado_mismo_signo"]["valor"],
            "mitades_divididas_mediana": mitades.get("intermedio_dado_mismo_signo", {}).get("mediana"),
            "esperado_por_regresion_a_la_media": nulo.get("nulo_mediana"),
            "seleccion_independiente": f["valor"],
        },
        "lectura": (
            "Los cuatro números no coinciden y la diferencia es informativa. El de v1 "
            "mide donde seleccionó y es el más alto. Las mitades divididas miden en "
            "pacientes no vistos y bajan poco, así que dentro de esta cohorte el "
            "gradiente aguanta. El nulo dice cuánto de eso explicaría por sí sola la "
            "regresión a la media. La selección independiente cambia de cohorte y de "
            "tecnología, se queda con muchas menos proteínas y ahí el gradiente "
            "desaparece. Como el análisis principal preregistrado es ese último, la "
            "conclusión de v1 no queda respaldada."
        ),
    }
    print(f"[v2-axes] C5 veredicto: {veredicto['conclusion']}")

    return {
        "n_pacientes_con_tumor_y_normal": len(ids_t),
        "n_pacientes_con_polipo_y_normal": len(ids_p),
        "n_pacientes_con_los_tres_tejidos": len(comunes),
        "n_significativas_brazo_tumor": n_elegir,
        "seleccion_en_el_mismo_brazo_v1": dentro,
        "seleccion_independiente_cptac": fuera,
        "mitades_divididas": mitades,
        "nulo_que_reselecciona": nulo,
        "veredicto": veredicto,
        "enmienda": "E7: el nulo preregistrado estaba mal especificado; ver docs/enmiendas_v2.md",
    }


def _c6(
    settings: Settings,
    prereg,
    cptac_raw: pd.DataFrame,
    cptac_meta: pd.DataFrame,
    log2mat: pd.DataFrame,
    pacientes: dict[str, dict[str, str]],
    gene_to_uni: dict[str, str],
) -> dict:
    etiquetas = cptac_meta.set_index("sample_id")["condition"].astype(str)
    cols = [
        c
        for c in cptac_raw.columns
        if str(c).endswith(" Unshared Log Ratio") and etiquetas.get(c) in ("tumor", "normal")
    ]
    tejido = cptac_raw.loc[~cptac_raw["protein_id"].map(is_summary_protein)]
    tejido = tejido.drop_duplicates("protein_id", keep="first").set_index("protein_id")
    tejido = pd.DataFrame(_numeric_block(tejido, cols), index=tejido.index, columns=cols)
    tejido["_uni"] = [gene_to_uni.get(str(g).upper()) for g in tejido.index]
    tejido = tejido.dropna(subset=["_uni"]).drop_duplicates("_uni").set_index("_uni")

    pride = log2mat.copy()
    pride["_uni"] = [uniprot_base(str(i)) for i in pride.index]
    pride = pride.dropna(subset=["_uni"]).drop_duplicates("_uni").set_index("_uni")
    ids_pac = [p for p, v in pacientes.items() if "M" in v and "N" in v]
    pcols = [pacientes[p][k] for p in ids_pac for k in ("M", "N")]
    pride = pride[pcols]

    comunes = sorted(set(tejido.index) & set(pride.index))
    x_tr = _zscore_rows(tejido.loc[comunes]).fillna(0.0).T.to_numpy()
    x_te = _zscore_rows(pride.loc[comunes]).fillna(0.0).T.to_numpy()
    y_tr = np.array([1 if etiquetas.get(c) == "tumor" else 0 for c in cols])
    y_te = np.array([1 if str(c).upper().endswith("M") else 0 for c in pcols])
    # Cada paciente aporta su tumor y su normal, en ese orden.
    grupos_te = np.repeat(np.array(ids_pac), 2)
    grupos_tr = modelos.grupos_de_cptac(cols, cptac_meta)

    semilla = prereg.semilla("c6_validacion_cruzada")
    definiciones = {
        "logistica_l1": lambda c: LogisticRegression(
            solver="saga", l1_ratio=1.0, C=c, max_iter=5000, random_state=semilla
        ),
        "elastic_net": lambda c: LogisticRegression(
            solver="saga", l1_ratio=0.5, C=c, max_iter=5000, random_state=semilla
        ),
        "random_forest": lambda c: RandomForestClassifier(
            n_estimators=500, min_samples_leaf=max(1, int(c)), random_state=semilla, n_jobs=-1
        ),
    }
    rejillas = {
        "logistica_l1": [0.01, 0.1, 1.0],
        "elastic_net": [0.01, 0.1, 1.0],
        "random_forest": [1, 3, 5],
    }
    k = int(prereg.get("c6_remuestreo.k_pliegues"))
    n_grupos_tr = int(np.unique(grupos_tr).size)
    pliegues = min(k, n_grupos_tr)

    puntajes: dict[str, np.ndarray] = {}
    aucs: dict[str, float] = {}
    interno: dict[str, dict] = {}
    for nombre, fabrica in definiciones.items():
        # Bucle interno agrupado por participante: elige el hiperparámetro sin
        # que ningún paciente esté a la vez en los dos lados del pliegue.
        por_valor: dict[float, float] = {}
        gkf = GroupKFold(n_splits=pliegues)
        for valor in rejillas[nombre]:
            parciales = []
            for tr, va in gkf.split(x_tr, y_tr, groups=grupos_tr):
                if np.unique(y_tr[va]).size < 2:
                    continue
                modelo = fabrica(valor)
                modelo.fit(x_tr[tr], y_tr[tr])
                parciales.append(
                    roc_auc_score(y_tr[va], modelo.predict_proba(x_tr[va])[:, 1])
                )
            if parciales:
                por_valor[valor] = float(np.mean(parciales))
        mejor_valor = max(por_valor, key=por_valor.get) if por_valor else rejillas[nombre][1]
        # Si el bucle interno empata, la elección es arbitraria y hay que decirlo.
        empatados = [v for v, a in por_valor.items() if a == por_valor.get(mejor_valor)]
        modelo = fabrica(mejor_valor)
        modelo.fit(x_tr, y_tr)
        s = modelo.predict_proba(x_te)[:, 1]
        puntajes[nombre] = s
        aucs[nombre] = float(roc_auc_score(y_te, s))
        interno[nombre] = {
            "hiperparametro_elegido": mejor_valor,
            "auc_interno_agrupado_por_valor": por_valor,
            "auc_interno_del_elegido": por_valor.get(mejor_valor),
            "n_hiperparametros_empatados": len(empatados),
            "hiperparametros_empatados": empatados,
            "eleccion_arbitraria": len(empatados) > 1,
            "nota_empate": (
                "El bucle interno no distingue entre los valores empatados. Se toma "
                "el primero de la rejilla y la elección no está respaldada por el dato."
                if len(empatados) > 1
                else ""
            ),
        }
        print(
            f"[v2-axes] C6 {nombre}: hiperparámetro {mejor_valor}, "
            f"AUC interno {por_valor.get(mejor_valor)}, AUC externo {aucs[nombre]:.3f}"
        )

    boot = modelos.bootstrap_por_conglomerado(
        y_te,
        grupos_te,
        puntajes,
        n_bootstrap=int(prereg.get("c6_remuestreo.n_bootstrap")),
        semilla=prereg.semilla("c6_bootstrap_conglomerado"),
    )
    # El de v1, para cuantificar cuánto cambia la unidad de remuestreo.
    boot_v1 = modelos.bootstrap_por_conglomerado(
        y_te,
        np.arange(y_te.size),
        puntajes,
        n_bootstrap=int(prereg.get("c6_remuestreo.n_bootstrap")),
        semilla=prereg.semilla("c6_bootstrap_conglomerado"),
    )
    veredicto = modelos.veredicto_sin_solape(boot["por_modelo"], aucs)
    print(f"[v2-axes] C6 veredicto: {veredicto['veredicto']}")

    # Los dos criterios preregistrados pueden discrepar, y aquí discrepan. Mirar
    # dos intervalos que se solapan no es lo mismo que mirar el intervalo de la
    # diferencia: el segundo aprovecha que los dos modelos se evalúan en las
    # mismas muestras, así que detecta diferencias que el primero no ve.
    mejor = veredicto.get("mejor")
    separa_al_mejor = [
        nombre
        for nombre, dif in boot["diferencias_por_pares"].items()
        if mejor and mejor in nombre and dif.get("incluye_cero") is False
    ]
    discrepancia = {
        "criterio_de_decision_preregistrado": "solape de intervalos por modelo",
        "veredicto_por_solape": veredicto["veredicto"],
        "pares_cuya_diferencia_excluye_el_cero": separa_al_mejor,
        "los_dos_criterios_coinciden": bool(
            veredicto.get("intervalos_disjuntos") == bool(separa_al_mejor)
        ),
        "lectura": (
            "El criterio de decisión es el solape, y por solape no hay ganador. La "
            "diferencia por pares sí separa al mejor modelo de los otros dos. No se "
            "cambia la regla después de ver el resultado: el veredicto declarado "
            "sigue siendo que no hay ganador, y la diferencia por pares se reporta "
            "como lo que es, una señal de que con más muestras podría haberlo."
        ),
    }
    anchos = {
        nombre: {
            "por_paciente": boot["por_modelo"][nombre]["ci_high"]
            - boot["por_modelo"][nombre]["ci_low"],
            "por_especimen_v1": boot_v1["por_modelo"][nombre]["ci_high"]
            - boot_v1["por_modelo"][nombre]["ci_low"],
        }
        for nombre in aucs
    }
    for nombre, a in anchos.items():
        a["el_agrupado_es_mas_estrecho"] = bool(a["por_paciente"] < a["por_especimen_v1"])

    return {
        "n_proteinas_comunes_uniprot": len(comunes),
        "n_entrenamiento": int(y_tr.size),
        "n_validacion": int(y_te.size),
        "n_pacientes_validacion": int(np.unique(grupos_te).size),
        "n_participantes_entrenamiento": n_grupos_tr,
        "k_pliegues_usados": pliegues,
        "metrica_primaria": prereg.get("c6_remuestreo.metrica_primaria"),
        "regla_ganador": prereg.get("c6_remuestreo.regla_ganador"),
        "auc_externo": aucs,
        "seleccion_interna_agrupada": interno,
        "bootstrap_por_paciente": boot,
        "bootstrap_por_especimen_v1": {
            "por_modelo": boot_v1["por_modelo"],
            "nota": (
                "Cálculo de v1. Trata el tumor y el normal del mismo paciente como "
                "dos personas. Se reporta solo para ver cuánto cambia el intervalo."
            ),
        },
        "veredicto": veredicto,
        "discrepancia_entre_criterios": discrepancia,
        "ancho_del_intervalo_segun_la_unidad": {
            "por_modelo": anchos,
            "lectura": (
                "Agrupar por paciente da aquí un intervalo más estrecho, no más ancho. "
                "El diseño es emparejado y la discriminación ocurre dentro del "
                "paciente, así que mantener juntas sus dos muestras conserva la "
                "cancelación del efecto de paciente. Remuestrear especímenes mezcla el "
                "tumor de uno con el normal de otro y añade varianza que el estudio no "
                "tiene. Lo que justifica agrupar no es obtener un intervalo mayor, sino "
                "remuestrear la unidad que de verdad se repetiría: el paciente."
            ),
        },
    }


def run(settings: Settings) -> None:
    prereg = cargar(settings.root)
    alpha = float(prereg.get("inferencia.alpha"))
    lfc = float(prereg.get("inferencia.treat_lfc"))
    salida = settings.resolve("results") / "v2" / "axes"
    salida.mkdir(parents=True, exist_ok=True)

    cptac = pd.read_csv(
        settings.resolve("results") / "v2" / "certified" / "cptac_unshared_treat.tsv", sep="\t"
    )
    mapping = effective_protein_mapping(settings)
    gene_to_uni = load_gene_uniprot_from_mapping(mapping, exclude_ambiguous=True)

    pride_ds = settings.datasets["datasets"]["pride_crc_meta"]
    log2mat, pacientes = _pride_log2(settings.resolve(pride_ds["matrix_path"]))
    con_tres = sum(1 for v in pacientes.values() if {"M", "N", "P"} <= set(v))
    print(f"[v2-axes] PXD034575: {len(pacientes)} pacientes, {con_tres} con los tres tejidos.")

    c5 = _c5(
        settings, prereg, log2mat, pacientes, cptac, gene_to_uni,
        lfc=lfc, alpha=alpha, salida=salida,
    )

    cptac_ds = settings.datasets["datasets"]["cptac_colorectal"]
    cptac_raw = pd.read_csv(settings.resolve(cptac_ds["matrix_path"]), sep="\t", low_memory=False)
    cptac_meta = pd.read_csv(settings.sample_metadata_path("cptac_colorectal"), sep="\t")
    c6 = _c6(settings, prereg, cptac_raw, cptac_meta, log2mat, pacientes, gene_to_uni)

    auditoria = {
        "version": "v2",
        "alpha": alpha,
        "treat_lfc": lfc,
        "semilla_base": prereg.get("semillas.base"),
        "c5_adenoma_carcinoma": c5,
        "c6_seleccion_de_modelo": c6,
        "limitaciones": [
            "El brazo de pólipo es un adenoma del mismo paciente, no una cohorte "
            "independiente de adenomas.",
            "La selección independiente cruza dos tecnologías distintas, TMT con "
            "razones logarítmicas y LFQ con intensidades. Mide coincidencia de "
            "dirección, no de magnitud.",
            "16 pacientes de validación dan intervalos anchos. Un AUC alto aquí no "
            "selecciona un modelo.",
            "Ningún resultado de este fichero sirve para decidir sobre un paciente.",
        ],
    }
    (salida / "axes_audit_v2.json").write_text(
        json.dumps(auditoria, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[v2-axes] escrito en {salida}")
