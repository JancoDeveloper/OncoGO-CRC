"""D2. Validación del índice estromal de R4 contra una pureza medida.

Lee el proteoma retrospectivo de TCGA que descargó
``scripts/download/tcga_cptac_retro.py`` y la tabla clínica de TCGA que descargó
``scripts/download/tcga_coadread.py``, y cruza las dos por paciente. Son los
mismos pacientes, que es justo lo que no ocurre con PDC000116.

Si falta cualquiera de las dos fuentes, la etapa lo dice y no calcula nada.

Nada de aquí cambia un subíndice del certificado. Lo que decide es si el lenguaje
obligatorio del preregistro sobre R4 —«una cota del efecto de composición, no una
pureza medida»— es una cautela o una necesidad.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.go_utils import load_go_long
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.v2 import pureza_medida as P
from oncogo_crc.v2.prereg import cargar

# Por encima de esto el índice se puede presentar como un sustituto de la pureza.
# Por debajo, no: sigue siendo un análisis de sensibilidad legítimo, pero no un
# control de pureza.
#
# Este umbral **no** está en el preregistro. El preregistro listaba
# `pureza_absolute_tcga` como método de sensibilidad pero no fijaba cuánto tenía
# que valer la correlación para que el índice pasara de cota a sustituto. Se fija
# aquí y queda declarado como enmienda E8, con la fecha y el motivo.
RHO_PARA_SUSTITUIR_LA_PUREZA = 0.7
ENMIENDA_DEL_UMBRAL = "E8"


def _tarjeta(base: Path) -> dict | None:
    path = base / "estado.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def _veredicto(validacion: dict) -> dict:
    """Qué se puede decir de R4 a la vista de la correlación con la pureza."""
    ecm = validacion.get("ecm", {})
    if ecm.get("estado") != "calculado":
        return {
            "estado": "no_calculable",
            "motivo": ecm.get("motivo", "el índice de matriz extracelular no salió"),
        }
    rho = float(ecm["rho"])
    sustituye = abs(rho) >= RHO_PARA_SUSTITUIR_LA_PUREZA and rho < 0
    return {
        "estado": "calculado",
        "rho_con_la_pureza_medida": rho,
        "varianza_compartida": float(ecm["varianza_compartida"]),
        "umbral": RHO_PARA_SUSTITUIR_LA_PUREZA,
        "el_umbral_no_estaba_preregistrado": (
            f"El preregistro contemplaba usar la pureza ABSOLUTE como método de "
            f"sensibilidad, pero no fijaba cuánto tenía que valer la correlación. "
            f"Este umbral se fija después y queda declarado como enmienda "
            f"{ENMIENDA_DEL_UMBRAL} en docs/enmiendas_v2.md. No tiene el mismo "
            f"valor que un umbral preregistrado."
        ),
        "el_indice_sustituye_a_la_pureza": bool(sustituye),
        "lectura": (
            "El índice estromal correlaciona con la pureza medida por ADN en el "
            f"sentido esperado (rho {rho:+.3f}), así que no está midiendo ruido. "
            f"Pero comparte solo el {100 * float(ecm['varianza_compartida']):.1f} % "
            "de la varianza: de cada cien partes de variación en pureza, el "
            "índice captura unas pocas. No es un sustituto de la pureza."
            if not sustituye
            else (
                f"El índice correlaciona con la pureza medida a rho {rho:+.3f}, "
                "por encima del umbral de la enmienda E8. Se puede presentar como "
                "un sustituto de la pureza."
            )
        ),
        "consecuencia_para_r4": (
            "R4 sigue siendo un análisis de sensibilidad válido: ajustar por el "
            "índice responde a «cuánto del efecto sobrevive si los pares no "
            "difieren en composición medida así». Lo que no se puede decir es que "
            "R4 controle la pureza. El lenguaje obligatorio del preregistro no es "
            "una cautela de más: es la descripción correcta."
        )
        if not sustituye
        else (
            "R4 puede describirse como un ajuste por pureza, porque el índice y "
            "la pureza medida coinciden por encima del umbral de E8."
        ),
    }


def run(settings: Settings) -> None:
    prereg = cargar(settings.root)
    proteoma = settings.resolve("external") / "tcga_cptac_retro"
    clinica_dir = settings.resolve("external") / "tcga_coadread"
    salida = settings.resolve("results") / "v2" / "axes" / "pureza_medida"
    salida.mkdir(parents=True, exist_ok=True)
    audit = salida / "pureza_medida_audit_v2.json"

    def _no_disponible(motivo: str) -> None:
        audit.write_text(
            json.dumps(
                {"version": "v2", "estado": "no_disponible", "motivo": motivo},
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"[v2-pureza] no_disponible: {motivo}")

    t_prot = _tarjeta(proteoma)
    t_clin = _tarjeta(clinica_dir)
    if t_prot is None or t_prot.get("estado") != "usada":
        return _no_disponible(
            "falta external/tcga_cptac_retro/estado.json en estado «usada»; hay "
            "que ejecutar scripts/download/tcga_cptac_retro.py"
            if t_prot is None
            else f"el proteoma retrospectivo está marcado {t_prot.get('estado')}: "
                 f"{t_prot.get('motivo', 'sin motivo')}"
        )
    if t_clin is None or t_clin.get("estado") != "usada":
        return _no_disponible(
            "falta external/tcga_coadread/estado.json en estado «usada»; sin la "
            "clínica de TCGA no hay pureza medida contra la que validar"
        )

    go_path = settings.go_gene_long_path()
    if not go_path.exists():
        return _no_disponible(
            "falta la tabla de GO; sin ella no se puede construir la misma firma "
            "de matriz extracelular que usa R4, y construir otra no validaría R4"
        )

    matriz = pd.read_csv(
        proteoma / t_prot["matriz_limpia"], sep="\t", index_col=0, low_memory=False
    )
    alicuotas = pd.read_csv(
        proteoma / "tcga_cptac_retro_alicuotas.tsv", sep="\t"
    )
    clin = pd.read_csv(
        clinica_dir / t_clin["tabla_limpia"], sep="\t", low_memory=False
    ).set_index("paciente")
    clin = clin[~clin.index.duplicated(keep="first")]
    print(
        f"[v2-pureza] proteoma: {matriz.shape[0]} genes x {matriz.shape[1]} "
        f"alícuotas de {alicuotas['paciente'].nunique()} pacientes."
    )

    go_long = load_go_long(go_path)
    firmas = {"ecm": P.genes_del_termino(go_long, P.ECM_GO_TERM)}
    for nombre, genes in (prereg.get("c7_composicion.firmas") or {}).items():
        firmas[nombre] = {str(g).strip().upper() for g in genes}

    try:
        validacion = P.validar_el_indice(matriz, alicuotas, clin, firmas)
    except P.FaltaColumna as e:
        return _no_disponible(str(e))

    for nombre, v in validacion.items():
        if v.get("estado") == "calculado":
            print(
                f"[v2-pureza] firma {nombre}: {v['n_proteinas_de_la_firma']} "
                f"proteínas, n={v['n']}, rho con la pureza {v['rho']:+.3f} "
                f"(p={v['p_rho']:.2e})"
            )
        else:
            print(f"[v2-pureza] firma {nombre}: {v.get('motivo')}")

    veredicto = _veredicto(validacion)
    if veredicto["estado"] == "calculado":
        print(f"[v2-pureza] {veredicto['lectura']}")

    repro = P.reproducibilidad(matriz, alicuotas)
    if repro["estado"] == "calculado":
        print(
            f"[v2-pureza] reproducibilidad del ensayo: r mediana "
            f"{repro['r_mediana']:.3f} entre {repro['n_pares_de_replicas']} pares "
            "de alícuotas del mismo tumor."
        )

    try:
        cea = P.cea_de_proteina(matriz, alicuotas, clin)
    except P.FaltaColumna as e:
        cea = {"estado": "no_calculable", "motivo": str(e)}

    repro_genes = {
        g: P.reproducibilidad_de_un_gen(matriz, alicuotas, g) for g in P.GENES_CEA
    }
    r5 = repro_genes.get("CEACAM5", {})
    if r5.get("estado") == "calculado":
        print(
            f"[v2-pureza] CEACAM5 medido dos veces en el mismo tumor difiere "
            f"{r5['veces_de_diferencia_mediana']:.1f} veces (mediana), y la "
            f"variación entre pacientes distintos es "
            f"{r5['desviacion_entre_pacientes_log2']:.2f} en log2: el ruido es el "
            f"{100 * r5['ruido_frente_a_senal']:.0f} % de la señal."
        )

    filas = [
        {"firma": k, **{kk: vv for kk, vv in v.items() if not isinstance(vv, str)}}
        for k, v in validacion.items()
    ]
    write_tsv(pd.DataFrame(filas), salida / "indice_frente_a_pureza.tsv")

    prot_cea = cea.get("proteina", {})
    write_tsv(
        pd.DataFrame(
            [
                {"gen": g, **{kk: vv for kk, vv in v.items() if not isinstance(vv, str)}}
                for g, v in prot_cea.items()
            ]
        ),
        salida / "ceacam_proteina_frente_a_cea.tsv",
    )

    auditoria = {
        "version": "v2",
        "estado": "calculado",
        "fuentes": {
            "proteoma": {
                "nombre": t_prot["nombre"],
                "pdc_study_id": t_prot["pdc_study_id"],
                "study_submitter_id": t_prot["study_submitter_id"],
                "metodo": t_prot["metodo"],
                "fecha_de_acceso": t_prot["fecha_de_acceso"],
                "cita": t_prot["cita"],
                "licencia": t_prot["licencia"],
                "sha256_matriz_limpia": t_prot["sha256_matriz_limpia"],
                "muestras": t_prot["muestras"],
                "discrepancias_con_el_enunciado": t_prot[
                    "discrepancias_con_el_enunciado"
                ],
                "independencia_del_estudio_hermano": t_prot[
                    "independencia_del_estudio_hermano"
                ],
            },
            "clinica": {
                "nombre": t_clin["nombre"],
                "fecha_de_acceso": t_clin["fecha_de_acceso"],
                "cita": t_clin["cita"],
            },
        },
        "que_se_pregunta": (
            "Si el índice estromal que R4 usa como covariable se parece a la "
            "pureza del tumor medida por ADN en los mismos pacientes."
        ),
        "por_que_aqui_se_puede": (
            "Son pacientes de TCGA, así que tienen proteoma y pureza ABSOLUTE en "
            "la misma persona. Los de PDC000116 no son de TCGA."
        ),
        "validacion_del_indice": validacion,
        "veredicto_sobre_r4": veredicto,
        "lenguaje_obligatorio_del_preregistro": str(
            prereg.get("c7_composicion.lenguaje_obligatorio")
        ),
        "reproducibilidad_del_ensayo": repro,
        "reproducibilidad_por_gen": repro_genes,
        "cea_de_proteina_frente_a_cea_serico": cea,
        "relacion_con_d1": (
            "D1 midió el CEACAM5 del tejido con transcrito y encontró que explica "
            "el 2,4 % de la varianza del CEA en sangre. Aquí está la proteína en "
            "los mismos pacientes, y antes de leer su correlación hay que leer su "
            "reproducibilidad."
        ),
        "no_entra_en_el_certificado": (
            "Ningún número de esta etapa es un subíndice. R4 no se recalcula aquí: "
            "esta cohorte no tiene normales, así que no hay contraste tumor-normal. "
            "Lo que se decide es cómo hay que describir R4."
        ),
        "limitaciones": list(t_prot.get("limitaciones", []))
        + [
            "La pureza ABSOLUTE es a su vez una estimación, hecha con datos de "
            "copias y de variantes. Es una medida de ADN y no un recuento de "
            "células al microscopio, así que la comparación es entre dos medidas "
            "indirectas y no contra un patrón de oro.",
            "La pieza del proteoma y la pieza del ADN no son el mismo trozo de "
            "tumor. Un tumor no es homogéneo, así que parte de la discrepancia es "
            "muestreo y no error de ninguna de las dos medidas.",
            "Las correlaciones son de rangos y el intervalo es la aproximación de "
            "Fisher. Con 83 pacientes va bien, pero es una aproximación.",
            "Las cuatro firmas se correlacionan con la misma pureza, así que son "
            "cuatro pruebas sobre la misma variable. Ninguna está corregida por "
            "multiplicidad: la que decide es la de matriz extracelular, que es la "
            "que R4 usa, y las otras tres son descriptivas.",
        ],
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
        "semilla_base": int(prereg.get("semillas.base")),
        "nota_de_semilla": (
            "Nada de esta etapa es estocástico: son correlaciones y desviaciones "
            "típicas. La semilla se registra por coherencia, no porque se use."
        ),
    }
    audit.write_text(
        json.dumps(auditoria, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[v2-pureza] escrito en {salida}")
