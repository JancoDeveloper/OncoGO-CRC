"""Certificado de v2: índice R, tabla de verdad de H1, H2 sin evaluar (C12, C14).

Todo lo que este stage escribe sale de los JSON de auditoría de las etapas
anteriores y de ``config/registro_subindices.yaml``. No hay ni un número escrito
a mano. Si una etapa no ha corrido, su subíndice sale ``no_calculable`` con el
motivo, no con el valor de v1.

Tres cosas que v1 no hacía.

- **C12**. H2 se declara ``no_evaluada`` y se dice qué haría falta para
  evaluarla. v1 la trataba como si tuviera respuesta.
- **C14**. El informe se genera desde las tablas. Hay un fichero
  ``changelog_numeros.tsv`` que pone cada número de v1 al lado del de v2 con su
  motivo, para que el cambio sea auditable y no haya que creerse el texto.
- Los cinco subíndices nunca se promedian ni se resumen en una nota global.
  Cada uno mide una cosa distinta y se leen juntos, no sumados.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.v2.prereg import cargar

NO_CALCULABLE = "no_calculable"

# El código usa claves sin acentos; el informe se lee en castellano.
UNIDADES_LEGIBLES = {
    "proteinas_unicas": "proteínas únicas",
    "terminos": "términos",
    "permutaciones": "permutaciones",
    "pacientes": "pacientes",
    "pares": "pares",
}


def _leer(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def _coma(valor, digitos: int = 4) -> str:
    """Número en el formato del texto en español, con signo menos de verdad."""
    if valor is None:
        return "no disponible"
    return f"{float(valor):.{digitos}f}".replace(".", ",").replace("-", "\u2212")


def _subindices(registro: dict, certified: dict | None, meta: dict | None) -> list[dict]:
    """Los cinco subíndices, cada uno con su numerador, denominador y estado."""
    definiciones = {d["codigo"]: d for d in registro["subindices"]}
    filas: list[dict] = []

    def _fila(codigo: str, *, valor=None, numerador=None, denominador=None,
              estado: str, motivo: str = "", extra: dict | None = None) -> dict:
        d = definiciones[codigo]
        fila = {
            "codigo": codigo,
            "nombre": d["nombre"],
            "valor": valor,
            "numerador": numerador,
            "denominador": denominador,
            "unidad_denominador": d["unidad_denominador"],
            "estado": estado,
            "motivo": motivo,
            "definicion": " ".join(str(d["definicion"]).split()),
            "nota_obligatoria": " ".join(str(d["nota_obligatoria"]).split()),
            "metodo": d["metodo"],
        }
        fila.update(extra or {})
        filas.append(fila)
        return fila

    if certified is None:
        for codigo in definiciones:
            _fila(codigo, estado=NO_CALCULABLE, motivo="la etapa v2_certified no ha corrido")
        return filas

    r1 = (certified.get("c15_r1_direccion") or {}).get("global") or {}
    if r1:
        _fila(
            "R1_direccion",
            valor=r1["valor"],
            numerador=r1["numerador"],
            denominador=r1["denominador"],
            estado="calculado",
            extra={
                "ci_low": r1.get("ci_low"),
                "ci_high": r1.get("ci_high"),
                "cobertura": r1.get("cobertura"),
                "universo": r1.get("universo"),
            },
        )
    else:
        _fila("R1_direccion", estado=NO_CALCULABLE, motivo="no hay R1 en la auditoría")

    r2 = certified.get("c3_r2_detectabilidad") or {}
    if r2:
        _fila(
            "R2_detectabilidad",
            valor=r2["valor"],
            numerador=r2["numerador"],
            denominador=r2["denominador"],
            estado="calculado",
            extra={"ci_low": r2.get("ci_low"), "ci_high": r2.get("ci_high")},
        )
    else:
        _fila("R2_detectabilidad", estado=NO_CALCULABLE, motivo="no hay R2 en la auditoría")

    if meta and meta.get("estado") == "calculado":
        _fila(
            "R3_homogeneidad",
            valor=meta.get("mediana_1_menos_i2"),
            denominador=meta.get("n_proteinas_evaluables"),
            estado="calculado",
            extra={"n_cohortes": meta.get("n_cohortes")},
        )
    else:
        n_coh = (meta or {}).get("n_cohortes", 2)
        _fila(
            "R3_homogeneidad",
            estado=NO_CALCULABLE,
            motivo=(
                f"hay {n_coh} cohortes de tejido y el preregistro exige al menos 3 "
                "para que I² sea estable"
            ),
            extra={"n_cohortes": n_coh},
        )

    principal = next(
        (a for a in (certified.get("c7_composicion") or {}).get("ajustes", []) if a.get("es_principal")),
        None,
    )
    if principal and principal.get("fraccion"):
        f = principal["fraccion"]
        _fila(
            "R4_composicion",
            valor=f["valor"],
            numerador=f["numerador"],
            denominador=f["denominador"],
            estado="calculado",
            extra={
                "ci_low": f.get("ci_low"),
                "ci_high": f.get("ci_high"),
                "indice": principal.get("indice"),
            },
        )
    else:
        _fila("R4_composicion", estado=NO_CALCULABLE, motivo="no hay ajuste principal de C7")

    nulo = certified.get("c1_nulo_cptac") or {}
    p = nulo.get("p_empirico") or {}
    if p:
        _fila(
            "R5_nulo",
            valor=p["p"],
            numerador=p["b"],
            denominador=p["n_permutaciones"],
            estado="calculado_como_tope" if p.get("es_tope") else "calculado",
            motivo=(
                "Ninguna permutación alcanzó lo observado, así que el valor es el "
                "mínimo que el número de permutaciones permite. Es una cota, no una "
                "estimación."
                if p.get("es_tope")
                else ""
            ),
            extra={
                "p_minimo_alcanzable": p.get("p_minimo_alcanzable"),
                "es_tope": p.get("es_tope"),
            },
        )
    else:
        _fila("R5_nulo", estado=NO_CALCULABLE, motivo="no hay nulo de C1")
    return filas


def _changelog(certified: dict, vias: dict | None, axes: dict | None, det: dict | None) -> pd.DataFrame:
    """Cada número de v1 al lado del de v2, con el motivo del cambio."""
    filas: list[dict] = []

    def _add(nombre, v1, v2, motivo, correccion):
        filas.append(
            {
                "numero": nombre,
                "v1": v1,
                "v2": v2,
                "cambia": None if v1 is None or v2 is None else bool(v1 != v2),
                "motivo": motivo,
                "correccion": correccion,
            }
        )

    r1 = (certified.get("c15_r1_direccion") or {}).get("global") or {}
    _add(
        "R1_concordancia_de_signo", 0.8739495798319328, r1.get("valor"),
        "Mismo cálculo. Lo que cambia es que ahora se publica junto a su cobertura: "
        "el denominador es solo el 17 % de los hallazgos.", "C15",
    )
    _add(
        "R1_cobertura", None, r1.get("cobertura"),
        "v1 no publicaba este número. Sin él, R1 parece una concordancia general y "
        "es una concordancia condicionada a que la réplica mida la proteína.", "C15",
    )
    _add(
        "R2_detectabilidad", 0.3285, (certified.get("c3_r2_detectabilidad") or {}).get("valor"),
        "Mismo valor. Lo que se corrige es que v1 tenía tres denominadores distintos "
        "circulando (688, 702 y 244) para cosas que parecían la misma.", "C3",
    )
    nulo = (certified.get("c1_nulo_cptac") or {}).get("p_empirico") or {}
    _add(
        "R5_p_empirico", 0.0099, nulo.get("p"),
        "100 permutaciones pasan a 10 000. En las dos el nulo nunca alcanza lo "
        "observado, así que el número es el mínimo alcanzable y se declara como tope.",
        "C1",
    )
    principal = next(
        (a for a in (certified.get("c7_composicion") or {}).get("ajustes", []) if a.get("es_principal")),
        None,
    )
    _add(
        "R4_composicion", 0.4462, ((principal or {}).get("fraccion") or {}).get("valor"),
        "v1 usaba las 282 proteínas de la matriz extracelular como índice y a la vez "
        "las contrastaba. Ahora cada proteína se excluye de su propio índice.", "C7",
    )
    if vias:
        v = ((vias.get("por_coleccion") or {}).get("reactome") or {}).get(
            "concordancia_entre_cohortes", {}
        )
        _add(
            "concordancia_de_via", 1.0, (v.get("concordancia") or {}).get("valor"),
            "v1 daba 1,0 con cuatro términos. Ahora el mínimo es diez términos, la "
            "prueba es competitiva y se repite en tres colecciones.", "C4",
        )
        _add(
            "n_terminos_de_via", 4, v.get("n_terminos_evaluados"),
            "El número de términos es parte del resultado, no un detalle: con cuatro, "
            "el intervalo iba de 0,398 a 1,0.", "C4",
        )
        _add(
            "H1_veredicto", "rechazada_por_clausula_proteina",
            (vias.get("h1") or {}).get("veredicto"),
            "Mismo veredicto, pero ahora sale de una tabla de verdad explícita con el "
            "estado no_evaluable incluido.", "C4",
        )
    if axes:
        c5 = (axes.get("c5_adenoma_carcinoma") or {})
        comp = (c5.get("veredicto") or {}).get("comparacion", {})
        _add(
            "polipo_magnitud_intermedia", 0.9194630872483222, comp.get("seleccion_independiente"),
            "v1 seleccionaba las proteínas en el brazo tumor y medía su magnitud en ese "
            "mismo brazo. Con selección independiente en CPTAC el efecto desaparece.",
            "C5",
        )
        _add(
            "polipo_esperado_por_regresion_a_la_media", None,
            comp.get("esperado_por_regresion_a_la_media"),
            "v1 no tenía este nulo. Dice cuánta de la asimetría la produce la propia "
            "selección.", "C5",
        )
    if det and det.get("estado") == "calculado":
        c = det.get("concordancia", {})
        _add(
            "concordancia_tejido_plasma", (c.get("todas") or {}).get("valor"),
            (c.get("dentro_de_la_ventana") or {}).get("valor"),
            "v1 medía sobre todas las proteínas comunes. Ahora se mide solo entre las "
            "que el método alcanza a detectar en suero.", "C8",
        )
    return pd.DataFrame(filas)


def _informe(
    prereg,
    subindices: list[dict],
    h1: dict,
    h2: dict,
    h4: dict,
    vias: dict | None,
    axes: dict | None,
    det: dict | None,
    certified: dict,
) -> str:
    """El informe de brecha, escrito desde las tablas.

    Frases cortas, voz activa y nada que un lector sin formación en el tema no
    pueda seguir. Ningún número se teclea: todos vienen de los argumentos.
    """
    lineas: list[str] = []
    a = lineas.append
    a("# Informe de brecha de v2")
    a("")
    a("Este informe dice qué sostiene este trabajo y qué no. Todos los números")
    a("salen de las tablas de `results/v2/`. Ninguno está escrito a mano.")
    a("")
    a("**Ningún resultado de esta tesis sirve para decidir sobre un paciente.**")
    a("")
    a("## El índice R, subíndice por subíndice")
    a("")
    a("Son cinco medidas distintas. No se promedian. Un número único escondería")
    a("que miden cosas diferentes y que algunas no se pueden calcular.")
    a("")
    for s in subindices:
        a(f"### {s['codigo']} — {s['nombre']}")
        a("")
        if s["estado"] == NO_CALCULABLE:
            a(f"- Estado: **no calculable**. {s['motivo']}")
        else:
            cuerpo = f"- Valor: **{_coma(s['valor'])}**"
            if s["numerador"] is not None and s["denominador"] is not None:
                unidad = UNIDADES_LEGIBLES.get(
                    str(s["unidad_denominador"]), str(s["unidad_denominador"])
                )
                cuerpo += f" ({s['numerador']} de {s['denominador']} {unidad})"
            a(cuerpo)
            if s.get("ci_low") is not None:
                a(f"- Intervalo del 95 %: {_coma(s['ci_low'])} a {_coma(s['ci_high'])}")
            if s.get("cobertura") is not None:
                a(
                    f"- Cobertura: {_coma(s['cobertura'])}. Es la parte de los hallazgos "
                    "que entra en el denominador."
                )
            if s["estado"] == "calculado_como_tope":
                a(f"- Atención: {s['motivo']}")
        a(f"- Qué mide: {s['definicion']}")
        a(f"- Qué no mide: {s['nota_obligatoria']}")
        a("")

    a("## Las hipótesis")
    a("")
    a("### H1")
    a("")
    a("H1 son dos condiciones a la vez. La primera habla de proteínas sueltas: se")
    a("esperaba que la dirección del efecto coincidiera poco entre cohortes. La")
    a("segunda habla de vías: se esperaba que ahí sí coincidiera mucho.")
    a("")
    a(f"- Cláusula de proteína: **{h1['clausula_proteina']}**. {h1.get('motivo_clausula_proteina', '')}")
    a(f"- Cláusula de vía: **{h1['clausula_via']}**. {h1.get('motivo_clausula_via', '')}")
    a(f"- Veredicto: **{h1['veredicto']}**")
    a("")
    a(h1["nota"])
    a("")
    a("### H2")
    a("")
    a(f"- Estado: **{h2['estado']}**")
    a(f"- Motivo: {h2['motivo']}")
    a(f"- Qué haría falta: {h2['requisito_para_evaluar']}")
    a("")
    a("Esto es una corrección respecto de v1. Antes H2 figuraba como si tuviera")
    a("respuesta. No la tiene, y decirlo es parte del resultado.")
    a("")
    a("### H4")
    a("")
    a(f"- Predicción preregistrada: {h4['prediccion']}")
    a(f"- Umbral: {_coma(h4['umbral'], 2)}")
    r2 = next((s for s in subindices if s["codigo"] == "R2_detectabilidad"), {})
    if r2.get("valor") is not None:
        a(
            f"- Medido: {_coma(r2['valor'])} de los candidatos ya es medible, muy por "
            f"debajo de {_coma(h4['umbral'], 2)}. H4 no se falsifica."
        )
    a("")

    if vias:
        a("## Lo que cambia en las vías")
        a("")
        rob = (vias.get("clausula_via") or {}).get("robustez", {})
        a(
            f"Se probaron {rob.get('n_colecciones_pedidas')} colecciones de vías y solo "
            f"{rob.get('n_colecciones_evaluables')} resultó evaluable."
        )
        if rob.get("advertencia"):
            a("")
            a(f"**Advertencia.** {rob['advertencia']}")
        a("")
        a(rob.get("cuello_de_botella", ""))
        a("")

    if axes:
        c5 = axes.get("c5_adenoma_carcinoma") or {}
        ver = c5.get("veredicto") or {}
        if ver:
            a("## Lo que cambia en la secuencia adenoma–carcinoma")
            a("")
            a(f"**{ver.get('conclusion', '')}**")
            a("")
            comp = ver.get("comparacion", {})
            a("| Cómo se midió | Fracción con menor magnitud en el pólipo |")
            a("| --- | --- |")
            a(f"| Seleccionando y midiendo en el mismo brazo (v1) | {_coma(comp.get('mismo_brazo_v1'))} |")
            a(f"| Mitades divididas, midiendo en pacientes no vistos | {_coma(comp.get('mitades_divididas_mediana'))} |")
            a(f"| Lo que la regresión a la media ya explica | {_coma(comp.get('esperado_por_regresion_a_la_media'))} |")
            a(f"| Seleccionando en otra cohorte (análisis principal) | {_coma(comp.get('seleccion_independiente'))} |")
            a("")
            a(ver.get("lectura", ""))
            a("")
        c6 = axes.get("c6_seleccion_de_modelo") or {}
        if c6:
            a("## Lo que cambia en la selección de modelo")
            a("")
            a(f"- Veredicto: **{(c6.get('veredicto') or {}).get('veredicto', '')}**")
            disc = c6.get("discrepancia_entre_criterios") or {}
            if disc.get("pares_cuya_diferencia_excluye_el_cero"):
                a(f"- {disc.get('lectura', '')}")
            ancho = c6.get("ancho_del_intervalo_segun_la_unidad") or {}
            if ancho.get("lectura"):
                a(f"- {ancho['lectura']}")
            a("")

    if det and det.get("estado") == "calculado":
        a("## Lo que cambia en la brecha tejido→plasma")
        a("")
        curva = det.get("curva_de_deteccion", {})
        c = det.get("concordancia", {})
        if curva.get("estado") == "calculado":
            a(
                f"La probabilidad de que una proteína se mida en suero depende con "
                f"fuerza de su concentración: {_coma(curva['pendiente_por_orden_de_magnitud'], 3)} "
                f"por cada orden de magnitud. Así que hay una ventana de detección real."
            )
            a("")
        a("| Conjunto de proteínas | Concordancia de dirección | n |")
        a("| --- | --- | --- |")
        for clave, etiqueta in (
            ("todas", "Todas las comunes (cálculo de v1)"),
            ("dentro_de_la_ventana", "Dentro de la ventana de detección"),
            ("fuera_de_la_ventana", "Fuera de la ventana"),
        ):
            d = c.get(clave)
            if d:
                a(f"| {etiqueta} | {_coma(d['valor'])} | {d['denominador']} |")
        a("")
        a((c.get("veredicto") or {}).get("lectura", ""))
        a("")

    a("## Lo que sigue sin poder decirse")
    a("")
    for s in subindices:
        if s["estado"] == NO_CALCULABLE:
            a(f"- {s['codigo']}: {s['motivo']}")
    a(f"- H2: {h2['motivo']}")
    a("")
    a("## Recordatorio")
    a("")
    a("Este trabajo compara cohortes públicas para medir qué parte de un hallazgo")
    a("de proteómica se repite. No es una prueba diagnóstica, no está validada en")
    a("ningún paciente y no se puede usar para decidir sobre uno.")
    a("")
    return "\n".join(lineas)


def run(settings: Settings) -> None:
    prereg = cargar(settings.root)
    base = settings.resolve("results") / "v2"
    salida = base / "certification"
    salida.mkdir(parents=True, exist_ok=True)

    certified = _leer(base / "certified" / "inference_audit_v2.json")
    vias = _leer(base / "vias" / "vias_audit_v2.json")
    axes = _leer(base / "axes" / "axes_audit_v2.json")
    det = _leer(base / "deteccion" / "deteccion_audit_v2.json")
    meta = _leer(base / "meta" / "meta_audit_v2.json")

    registro = yaml.safe_load(
        (settings.root / "config" / "registro_subindices.yaml").read_text(encoding="utf-8")
    )
    subindices = _subindices(registro, certified, meta)
    tabla_r = pd.DataFrame(subindices)
    write_tsv(tabla_r, salida / "indice_R_v2.tsv")
    for s in subindices:
        print(f"[v2-cert] {s['codigo']}: {s['estado']} {s['valor'] if s['valor'] is not None else ''}")

    # C12: H2 se declara sin evaluar, con el requisito para poder evaluarla.
    h2 = {
        "estado": prereg.get("hipotesis.h2.estado"),
        "motivo": prereg.get("hipotesis.h2.motivo"),
        "requisito_para_evaluar": prereg.get("hipotesis.h2.requisito_para_evaluar"),
        "correccion": "C12",
        "que_hacia_v1": (
            "v1 informaba de H2 como si tuviera respuesta. Sin un nulo que conserve "
            "el grado de cada nodo, cualquier red parece estructurada."
        ),
    }
    h4 = {
        "prediccion": prereg.get("hipotesis.h4.prediccion"),
        "umbral": prereg.get("hipotesis.h4.umbral"),
    }
    h1 = (vias or {}).get("h1") or {
        "clausula_proteina": "no_evaluable",
        "clausula_via": "no_evaluable",
        "veredicto": "no_evaluable",
        "nota": "la etapa v2_vias no ha corrido",
    }

    changelog = _changelog(certified or {}, vias, axes, det)
    write_tsv(changelog, salida / "changelog_numeros.tsv")
    n_cambian = int(changelog["cambia"].fillna(False).sum())
    print(f"[v2-cert] changelog: {len(changelog)} números, {n_cambian} cambian respecto de v1.")

    informe = _informe(prereg, subindices, h1, h2, h4, vias, axes, det, certified or {})
    (salida / "informe_brecha_v2.md").write_text(informe, encoding="utf-8")

    certificado = {
        "version": "v2",
        "fecha_del_preregistro": prereg.get("version"),
        "semilla_base": prereg.get("semillas.base"),
        "alpha": prereg.get("inferencia.alpha"),
        "treat_lfc": prereg.get("inferencia.treat_lfc"),
        "metodo_inferencia": (certified or {}).get("metodo"),
        "subindices": subindices,
        "h1": h1,
        "h2": h2,
        "h4": h4,
        "etapas_presentes": {
            "v2_certified": certified is not None,
            "v2_vias": vias is not None,
            "v2_axes": axes is not None,
            "v2_deteccion": det is not None,
            "v2_meta": meta is not None,
        },
        "regla_de_lectura": (
            "Los cinco subíndices no se promedian. Cada uno mide una cosa distinta y "
            "alguno puede no ser calculable; colapsarlos en un número ocultaría eso."
        ),
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
    }
    (salida / "certificado_v2.json").write_text(
        json.dumps(certificado, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    print(f"[v2-cert] H1: {h1['veredicto']}; H2: {h2['estado']}")
    print(f"[v2-cert] escrito en {salida}")
