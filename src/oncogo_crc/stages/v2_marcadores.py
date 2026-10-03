"""D8. Meta-análisis bivariado del CEA frente a la metilación y al ctDNA.

Qué pregunta responde esta etapa. Todo este trabajo busca proteínas que separen
tumor de tejido normal, y C8 ya mostró que lo que se ve en tejido no llega al
plasma. Pero falta la pregunta de arriba: el marcador que sí se usa en clínica,
el CEA, ¿qué rinde de verdad? Y las alternativas de ADN, ¿rinden más?

La respuesta es un punto con su intervalo, no un número. Y sale de tablas 2x2
copiadas de la literatura, con la referencia, la tabla y la página de cada fila,
y con una columna que una persona tiene que firmar.

Lo que esta etapa no hace: no combina el CEA con la metilación en un mismo
modelo. El CEA de estas tablas detecta recaída después de una cirugía; la
metilación detecta cáncer en cribado. Son dos preguntas distintas, y restar sus
sensibilidades no significaría nada aunque la aritmética funcionara.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.v2 import reitsma as R

# La tabla de extracción vive en el repositorio, no en el workspace de datos: es
# un artefacto curado a mano que necesita revisión y historia de cambios.
TABLA = Path(__file__).resolve().parents[3] / "data" / "literatura" / "extraccion_2x2.tsv"
FICHA = TABLA.parent / "extraccion_2x2_ficha.json"


def _especificidad_desde_fpr(fpr: float, fpr_ci: list | tuple) -> tuple[float, tuple]:
    """Pasa una tasa de falsos positivos a especificidad.

    Al restar de uno los extremos del intervalo se invierten de orden: el límite
    superior de la tasa de falsos positivos da el límite inferior de la
    especificidad. Es un detalle de dos líneas que, si se olvida, produce un
    intervalo del revés que sigue teniendo buen aspecto.
    """
    bajo, alto = float(fpr_ci[0]), float(fpr_ci[1])
    return 1.0 - float(fpr), (1.0 - alto, 1.0 - bajo)


def run(settings: Settings) -> None:
    salida = settings.resolve("results") / "v2" / "marcadores"
    salida.mkdir(parents=True, exist_ok=True)
    destino_json = salida / "marcadores_audit_v2.json"

    if not TABLA.exists():
        auditoria = {
            "version": "v2",
            "estado": "no_disponible",
            "motivo": (
                f"no existe {TABLA.name}; hay que ejecutar "
                "scripts/literatura/extraer_2x2.py"
            ),
        }
        destino_json.write_text(
            json.dumps(auditoria, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        print(f"[v2-marcadores] no_disponible: {auditoria['motivo']}")
        return

    d = pd.read_csv(TABLA, sep="\t")
    ficha = json.loads(FICHA.read_text(encoding="utf-8")) if FICHA.exists() else {}
    pendientes = int(d["verificado_por"].isna().sum() + (d["verificado_por"] == "").sum())

    ajustes: dict[str, dict] = {}
    filas = []
    for (prueba, indicacion), g in d.groupby(["prueba", "indicacion"]):
        clave = f"{prueba}|{indicacion}"
        if len(g) < 4:
            ajustes[clave] = {
                "estado": "no_calculable",
                "k": int(len(g)),
                "motivo": (
                    f"el modelo bivariado tiene cinco parámetros y hay "
                    f"{len(g)} estudios"
                ),
            }
            print(f"[v2-marcadores] {clave}: no_calculable ({len(g)} estudios)")
            continue
        res = R.reitsma(g["tp"], g["fp"], g["fn"], g["tn"])
        ajustes[clave] = {"estado": "calculado", **res.to_dict()}
        filas.append({"prueba": prueba, "indicacion": indicacion, **res.to_dict()})
        print(
            f"[v2-marcadores] {clave}: k={res.k}, "
            f"sensibilidad {res.sensibilidad:.3f} "
            f"({res.sensibilidad_ci_low:.3f}-{res.sensibilidad_ci_high:.3f}), "
            f"especificidad {res.especificidad:.3f} "
            f"({res.especificidad_ci_low:.3f}-{res.especificidad_ci_high:.3f}), "
            f"correlacion {res.correlacion:+.3f}"
        )
    if filas:
        write_tsv(pd.DataFrame(filas), salida / "reitsma_por_prueba.tsv")

    cea = ajustes.get("CEA|deteccion_de_recaida", {})
    sept9 = ajustes.get("metilacion_SEPT9|cribado_o_diagnostico", {})

    # La comparación entre las dos pruebas con tablas propias queda marcada como
    # no comparable, porque las indicaciones son distintas. El argumento se pasa
    # explícito y no por defecto.
    cruce = None
    if cea.get("estado") == "calculado" and sept9.get("estado") == "calculado":
        cruce = R.diferencia(
            R.ResultadoReitsma(**{k: v for k, v in cea.items() if k != "estado"}),
            R.ResultadoReitsma(**{k: v for k, v in sept9.items() if k != "estado"}),
            comparable=False,
        )

    # El ctDNA sí es la misma indicación que el CEA de estas tablas: detectar
    # recaída tras la cirugía. Pero de él no hay tablas 2x2 por estudio, así que
    # se compara con su estimación combinada publicada y se dice que no es un
    # contraste.
    ctdna = {}
    publicado = (ficha.get("ctdna_solo_combinado") or {})
    if cea.get("estado") == "calculado":
        objeto_cea = R.ResultadoReitsma(
            **{k: v for k, v in cea.items() if k != "estado"}
        )
        for nombre, v in publicado.items():
            espec, espec_ci = _especificidad_desde_fpr(v["fpr"], v["fpr_ci"])
            ctdna[nombre] = {
                "que_es": v["que_es"],
                **R.comparar_con_publicado(
                    objeto_cea,
                    sensibilidad=v["sensibilidad"],
                    sensibilidad_ci=tuple(v["sensibilidad_ci"]),
                    especificidad=espec,
                    especificidad_ci=espec_ci,
                    fuente=f"{(ficha.get('fuentes') or {}).get('zhu_2024', {}).get('pmc', 'zhu_2024')}, tabla 2, fila CRC",
                    k_publicado=v["k"],
                ),
                "misma_indicacion_que_el_cea": True,
            }
            c = ctdna[nombre]
            print(
                f"[v2-marcadores] ctDNA {nombre}: sensibilidad "
                f"{v['sensibilidad']:.2f} frente a "
                f"{objeto_cea.sensibilidad:.2f} del CEA; los intervalos "
                f"{'se solapan' if c['intervalos_de_sensibilidad_se_solapan'] else 'no se solapan'}"
            )

    auditoria = {
        "version": "v2",
        "estado": "calculado",
        "modelo": (
            "Reitsma bivariado de efectos aleatorios sobre el logit de la "
            "sensibilidad y el logit de la especificidad, con la varianza dentro "
            "del estudio tomada del binomio. Implementado en Python porque no "
            "hay R (enmienda E2)."
        ),
        "tabla_de_extraccion": str(TABLA.relative_to(TABLA.parents[2])),
        "n_filas": int(len(d)),
        "pendiente_de_verificacion_humana": pendientes,
        "ajustes": ajustes,
        "comparacion_entre_cea_y_metilacion": cruce,
        "comparacion_entre_cea_y_ctdna": ctdna,
        "fuentes": ficha.get("fuentes"),
        "regla_de_algoritmo_en_la_metilacion": ficha.get(
            "regla_de_algoritmo_en_la_metilacion"
        ),
        "filas_excluidas_a_proposito": ficha.get("filas_excluidas_a_proposito"),
        "que_aporta_a_la_tesis": (
            "Pone un techo a la pregunta de partida. Antes de buscar una "
            "proteína nueva conviene saber qué rinde la que ya se usa, y qué "
            "rinde la alternativa de ADN en la misma indicación."
        ),
        "no_entra_en_el_certificado": (
            "Ninguno de estos números es un subíndice. Son estimaciones de la "
            "literatura, no resultados de las cohortes de este trabajo."
        ),
        "limitaciones": [
            "La varianza dentro de cada estudio se trata como conocida, que es "
            "la aproximación normal de Reitsma 2005 y no el binomio exacto. Con "
            "celdas pequeñas empeora, y se cuenta cuántos estudios tienen alguna "
            "celda por debajo de cinco.",
            "El intervalo del punto combinado sale de la covarianza de mínimos "
            "cuadrados generalizados, que trata la varianza entre estudios como "
            "conocida. Eso lo estrecha un poco.",
            "De los 25 estudios de metilación, la mayoría son de casos y "
            "controles, donde la proporción de enfermos la fija el diseño. La "
            "revisión de origen califica de alto el riesgo de sesgo de selección "
            "en 21 de ellos. La especificidad de un diseño así no es la que se "
            "vería en una población de cribado.",
            "Del ctDNA no hay tablas 2x2 por estudio en la fuente, así que se "
            "usa su estimación combinada publicada. Comparar un ajuste propio "
            "con el de otro grupo arrastra diferencias de criterio de inclusión "
            "que no entran en ningún intervalo.",
            "No se calcula el área bajo la curva resumen. Obtenerla exige "
            "reparametrizar a HSROC y extrapolar a umbrales que ningún estudio "
            "midió, y aquí no hace falta: la pregunta es qué rinde cada prueba "
            "en su punto de corte habitual.",
            "Las filas están pendientes de verificación humana mientras la "
            "columna verificado_por siga vacía.",
        ],
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
    }
    destino_json.write_text(
        json.dumps(auditoria, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[v2-marcadores] escrito en {salida}")
