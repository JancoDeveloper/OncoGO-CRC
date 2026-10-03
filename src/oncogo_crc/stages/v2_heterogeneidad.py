"""D3. Variación dentro de un mismo tumor, y la búsqueda de cohortes que la trajo.

Lee la matriz que descargó ``scripts/download/pride_heterogeneidad.py``. Esa
cohorte no cumple el mínimo de pacientes del preregistro, así que **no entra en
el meta-análisis**: lo que aporta es la única cosa que ninguna cohorte del
certificado tiene, que son ocho regiones del mismo tumor.

Si la descarga no está, la etapa lo dice y no calcula nada.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.v2 import heterogeneidad as H
from oncogo_crc.v2.prereg import cargar

# La semilla se deriva de la base del preregistro con un desplazamiento propio,
# igual que el resto de los análisis estocásticos.
DESPLAZAMIENTO = 303


def _tarjeta(base: Path) -> dict | None:
    path = base / "estado.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def _certificadas(settings: Settings) -> set[str]:
    path = (
        settings.resolve("results")
        / "v2"
        / "certified"
        / "cptac_unshared_treat.tsv"
    )
    if not path.exists():
        return set()
    d = pd.read_csv(path, sep="\t", low_memory=False)
    if "significant" not in d.columns or "protein_id" not in d.columns:
        return set()
    sig = d["significant"]
    if sig.dtype == object:
        sig = sig.astype(str).str.lower().isin(("true", "1"))
    return {str(p).strip() for p in d.loc[sig.astype(bool), "protein_id"]}


def run(settings: Settings) -> None:
    prereg = cargar(settings.root)
    fuente = settings.resolve("external") / "heterogeneidad_intratumoral"
    salida = settings.resolve("results") / "v2" / "axes" / "heterogeneidad"
    salida.mkdir(parents=True, exist_ok=True)
    audit = salida / "heterogeneidad_audit_v2.json"

    tarjeta = _tarjeta(fuente)
    if tarjeta is None or tarjeta.get("estado") != "usada":
        motivo = (
            "falta external/heterogeneidad_intratumoral/estado.json; hay que "
            "ejecutar scripts/download/pride_heterogeneidad.py"
            if tarjeta is None
            else f"la fuente está marcada {tarjeta.get('estado')}: "
                 f"{tarjeta.get('motivo', 'sin motivo')}"
        )
        audit.write_text(
            json.dumps(
                {"version": "v2", "estado": "no_disponible", "motivo": motivo},
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"[v2-heterogeneidad] no_disponible: {motivo}")
        return

    matriz = pd.read_csv(
        fuente / tarjeta["matriz_limpia"], sep="\t", index_col=0, low_memory=False
    )
    diseno = tarjeta["parseo"]["diseno"]
    regiones = list(diseno["regiones_del_mismo_tumor"])
    otros = list(diseno["tumores_de_un_paciente"])
    mucosas = list(diseno["mucosas_normales"])
    metastasis = list(diseno["metastasis_hepaticas"])
    t4 = str(diseno["tumor_elegido_para_el_paciente_4"])
    semilla = int(prereg.get("semillas.base")) + DESPLAZAMIENTO

    print(
        f"[v2-heterogeneidad] {matriz.shape[0]} genes; {len(regiones)} regiones del "
        f"mismo tumor y {len(otros) + 1} pacientes."
    )

    comparacion = H.variacion_dentro_y_entre(
        matriz,
        regiones_del_mismo_tumor=regiones,
        tumores_de_otros_pacientes=otros,
        semilla=semilla,
    )
    if comparacion["estado"] == "calculado":
        print(
            f"[v2-heterogeneidad] con el número de observaciones igualado: "
            f"desviación mediana dentro del tumor "
            f"{comparacion['sd_mediana_dentro_del_tumor']:.3f}, entre pacientes "
            f"{comparacion['sd_mediana_entre_pacientes']:.3f}; cociente "
            f"{comparacion['cociente_mediano']:.3f} "
            f"({comparacion['cociente_ci_low']:.3f}–"
            f"{comparacion['cociente_ci_high']:.3f})."
        )
        print(f"[v2-heterogeneidad] {comparacion['veredicto']}")
    else:
        print(f"[v2-heterogeneidad] {comparacion['motivo']}")

    cert = _certificadas(settings)
    por_proteina = H.proteinas_mas_heterogeneas(
        matriz,
        regiones_del_mismo_tumor=regiones,
        tumores_de_otros_pacientes=otros,
        certificadas=cert,
    )
    tabla = por_proteina.pop("tabla", None)
    if tabla is not None:
        write_tsv(tabla, salida / "variacion_por_proteina.tsv")
    if por_proteina["estado"] == "calculado":
        c = por_proteina.get("certificadas", {})
        if "cociente_mediano" in c:
            print(
                f"[v2-heterogeneidad] de las {len(cert)} proteínas certificadas, "
                f"{c['n_presentes']} están aquí: cociente mediano "
                f"{c['cociente_mediano']:.3f} frente a "
                f"{c['cociente_mediano_del_resto']:.3f} en el resto."
            )

    pares_tumor = [(t, n) for t, n in zip(otros + [t4], mucosas, strict=False)]
    pares_meta = [(lm, t) for lm, t in zip(metastasis, otros + [t4], strict=False)]
    tumor_vs_mucosa = H.contraste_pareado(
        matriz, pares=pares_tumor, nombre="tumor_frente_a_mucosa_normal"
    )
    meta_vs_tumor = H.contraste_pareado(
        matriz, pares=pares_meta, nombre="metastasis_hepatica_frente_a_tumor_primario"
    )
    for c in (tumor_vs_mucosa, meta_vs_tumor):
        if c["estado"] == "calculado":
            print(
                f"[v2-heterogeneidad] {c['contraste']}: "
                f"{c['n_con_todos_los_pares_en_la_misma_direccion']} de "
                f"{c['n_proteinas']} proteínas con los {c['n_pares']} pares en la "
                f"misma dirección ({100 * c['fraccion_con_todos_los_pares_en_la_misma_direccion']:.0f} %), "
                f"frente al {100 * c['fraccion_esperada_por_azar']:.0f} % que daría "
                f"el azar: {c['veces_por_encima_del_azar']:.1f} veces más."
            )

    write_tsv(
        pd.DataFrame(tarjeta["busqueda_de_cohortes_d3"]),
        salida / "busqueda_de_cohortes_d3.tsv",
    )

    auditoria = {
        "version": "v2",
        "estado": "calculado",
        "fuente": {
            "nombre": tarjeta["nombre"],
            "accesion": tarjeta["accesion"],
            "fecha_de_acceso": tarjeta["fecha_de_acceso"],
            "cita": tarjeta["cita"],
            "licencia": tarjeta["licencia"],
            "sha256_matriz_limpia": tarjeta["sha256_matriz_limpia"],
            "n_pacientes": tarjeta["n_pacientes"],
            "parseo": tarjeta["parseo"],
        },
        "busqueda_de_cohortes_d3": tarjeta["busqueda_de_cohortes_d3"],
        "criterios_del_preregistro": {
            "min_pacientes": int(prereg.get("d3_inclusion_cohortes.min_pacientes")),
            "requiere_tumor_y_normal": bool(
                prereg.get("d3_inclusion_cohortes.requiere_tumor_y_normal")
            ),
            "requiere_matriz_descargable": bool(
                prereg.get("d3_inclusion_cohortes.requiere_matriz_descargable")
            ),
            "resultado_de_la_busqueda": (
                "Ningún depósito nuevo cumple los cuatro criterios. Los que "
                "publican matriz no llegan a 10 pacientes o no traen tejido "
                "normal; los que tienen el diseño adecuado no publican matriz de "
                "proteínas. El meta-análisis de C10 se queda con las cohortes que "
                "ya tenía."
            ),
        },
        "no_entra_en_el_meta_analisis": tarjeta["no_entra_en_el_meta_analisis"],
        "que_se_pregunta": (
            "Si un trozo de tumor representa al tumor. Las tres cohortes del "
            "certificado miden una biopsia por paciente y suponen que sí, sin "
            "decirlo, porque ninguna tiene con qué comprobarlo."
        ),
        "variacion_dentro_frente_a_entre": comparacion,
        "variacion_por_proteina": por_proteina,
        "tumor_frente_a_mucosa_normal": tumor_vs_mucosa,
        "metastasis_frente_a_tumor_primario": meta_vs_tumor,
        "consecuencia_para_la_tesis": (
            "Si dos trozos del mismo tumor se parecen menos entre sí que dos "
            "tumores de pacientes distintos, el ruido de la biopsia no se reduce "
            "añadiendo pacientes: entra en cada medida. No invalida los efectos "
            "del certificado, porque el contraste es dentro de cada paciente y el "
            "ruido no tiene dirección. Lo que cambia es la interpretación: lo "
            "medido es un trozo de tumor, no el tumor."
        ),
        "no_entra_en_el_certificado": (
            "Ningún número de esta etapa es un subíndice. Son cuatro pacientes y "
            "un solo tumor con regiones múltiples."
        ),
        "limitaciones": list(tarjeta.get("limitaciones", [])),
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
        "semilla": semilla,
        "derivacion_de_la_semilla": (
            f"semillas.base ({prereg.get('semillas.base')}) + {DESPLAZAMIENTO}"
        ),
    }
    audit.write_text(
        json.dumps(auditoria, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[v2-heterogeneidad] escrito en {salida}")
