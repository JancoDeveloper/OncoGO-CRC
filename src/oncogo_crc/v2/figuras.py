"""Figuras de v2. Cada una sale de una tabla y escribe la tabla que la genera.

Tres reglas para todas.

1. Nunca se dibuja un punto sin su intervalo. Un punto solo invita a comparar
   magnitudes que el dato no distingue, y eso es lo que hizo v1 al declarar
   ganador un modelo con intervalos solapados.
2. Lo que no se puede evaluar se dibuja como no evaluable, no se omite. Si GO BP
   desaparece de la figura, el lector entiende que no se probó.
3. Cada figura devuelve la ruta de la tabla con la que se dibujó. El verificador
   de números y el informe leen esa misma tabla, así que el texto y la figura no
   pueden desincronizarse.

Formato: PNG a 300 dpi para el PDF y SVG para la web. Los decimales van con coma,
que es la convención del idioma del informe.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

DPI = 300
GRIS = "#5a5a5a"
AZUL = "#2a5d8f"
NARANJA = "#c1651a"
ROJO = "#a4243b"
VERDE = "#3f6b4a"


@dataclass
class Figura:
    """Una figura con su tabla de origen y su pie."""

    nombre: str
    titulo: str
    pie: str
    png: Path
    svg: Path
    tsv: Path


def _coma(x, _pos=None) -> str:
    return f"{x:.2f}".replace(".", ",").replace("-", "\u2212")


def _coma3(x, _pos=None) -> str:
    return f"{x:.3f}".replace(".", ",").replace("-", "\u2212")


def _coma_fino(x) -> str:
    """Como ``_coma3``, salvo que redondear a tres decimales borraría el número.

    Un p de 9,999·10⁻⁵ impreso con tres decimales sale «0,000», que es falso y
    además hace imposible distinguirlo del 0,0099 de v1, que es cien veces mayor.
    """
    v = float(x)
    if v != 0 and abs(v) < 0.005:
        mantisa, exp = f"{abs(v):.1e}".split("e")
        signo = "\u2212" if v < 0 else ""
        return f"{signo}{mantisa.replace('.', ',')}·10{_superindice(int(exp))}"
    return _coma3(v)


def _superindice(n: int) -> str:
    tabla = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
    return str(n).translate(tabla)


def _estilo(ax, *, formato_x=None, formato_y=None) -> None:
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    ax.tick_params(labelsize=8, colors=GRIS)
    ax.grid(axis="both", color="#e4e4e4", lw=0.6)
    ax.set_axisbelow(True)
    if formato_x:
        ax.xaxis.set_major_formatter(FuncFormatter(formato_x))
    if formato_y:
        ax.yaxis.set_major_formatter(FuncFormatter(formato_y))


def _guardar(fig, salida: Path, nombre: str) -> tuple[Path, Path]:
    salida.mkdir(parents=True, exist_ok=True)
    png = salida / f"{nombre}.png"
    svg = salida / f"{nombre}.svg"
    fig.savefig(png, dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png, svg


def _tabla(df: pd.DataFrame, salida: Path, nombre: str) -> Path:
    salida.mkdir(parents=True, exist_ok=True)
    path = salida / f"{nombre}.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def _linea_vertical(ax, x, texto, *, color, estilo="--", arriba=0.97) -> None:
    """Línea de referencia con su etiqueta dentro de los ejes.

    Etiquetarla con coordenadas de datos la deja fuera del lienzo en cuanto el
    número de filas cambia. En fracción de ejes siempre cae dentro.
    """
    ax.axvline(x, color=color, ls=estilo, lw=1.1, zorder=1)
    ax.annotate(
        texto,
        xy=(x, arriba),
        xycoords=ax.get_xaxis_transform(),
        xytext=(3, 0),
        textcoords="offset points",
        fontsize=7.5,
        color=color,
        ha="left",
        va="top",
    )


def _filas_a_la_derecha(ax) -> None:
    """Etiquetas del eje y a la derecha, para que no invadan el panel vecino."""
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.tick_params(axis="y", length=0, pad=4)


def _holgura(ax, n_filas: int) -> None:
    """Margen vertical para que la fila de arriba no choque con el título."""
    ax.set_ylim(-0.6, n_filas - 0.1)


def _barras_de_error(ax, y, valor, low, high, *, color, etiqueta=None) -> None:
    """Punto con intervalo. Los faltantes no se dibujan: se marcan aparte."""
    valor = np.atleast_1d(np.asarray(valor, dtype=float))
    low = np.atleast_1d(np.asarray(low, dtype=float))
    high = np.atleast_1d(np.asarray(high, dtype=float))
    # xerr tiene que ser (2, n) aunque n sea 1, o matplotlib lo lee como dos puntos.
    izq = np.where(np.isfinite(low), valor - low, 0.0)
    der = np.where(np.isfinite(high), high - valor, 0.0)
    ax.errorbar(
        valor,
        np.atleast_1d(y),
        xerr=np.vstack([izq, der]),
        fmt="o",
        color=color,
        ecolor=color,
        elinewidth=1.4,
        capsize=3,
        markersize=6,
        label=etiqueta,
    )


# --------------------------------------------------------------------- C5
def fig_c5_gradiente(axes: dict, salida: Path) -> Figura:
    """Los cuatro números del gradiente adenoma-carcinoma, en la misma escala."""
    c5 = axes["c5_adenoma_carcinoma"]
    v1 = c5["seleccion_en_el_mismo_brazo_v1"]["intermedio_dado_mismo_signo"]
    mitades = c5["mitades_divididas"]["intermedio_dado_mismo_signo"]
    indep = c5["seleccion_independiente_cptac"]["intermedio_dado_mismo_signo"]
    nulo = c5["nulo_que_reselecciona"]

    filas = [
        {
            "calculo": "Mismo brazo (v1)",
            "valor": v1["valor"],
            "ci_low": v1["ci_low"],
            "ci_high": v1["ci_high"],
            "n": v1["denominador"],
            "que_mide": "selección y medición en el mismo brazo",
        },
        {
            "calculo": "Mitades divididas",
            "valor": mitades["mediana"],
            "ci_low": mitades["p2_5"],
            "ci_high": mitades["p97_5"],
            "n": int(c5["mitades_divididas"]["n_proteinas_evaluables_mediana"]),
            "que_mide": "selección en unos pacientes, medición en otros",
        },
        {
            "calculo": "Esperado por azar",
            "valor": nulo["nulo_mediana"],
            "ci_low": nulo["nulo_p2_5"],
            "ci_high": nulo["nulo_p97_5"],
            "n": nulo["n_elegibles_en_los_dos_brazos"],
            "que_mide": "nulo que vuelve a seleccionar",
        },
        {
            "calculo": "Selección independiente",
            "valor": indep["valor"],
            "ci_low": indep["ci_low"],
            "ci_high": indep["ci_high"],
            "n": indep["denominador"],
            "que_mide": "selección en CPTAC, medición en la réplica",
        },
    ]
    df = pd.DataFrame(filas)
    tsv = _tabla(df, salida, "c5_gradiente_adenoma_carcinoma")

    fig, ax = plt.subplots(figsize=(7.4, 3.3))
    colores = [GRIS, AZUL, NARANJA, ROJO]
    ys = np.arange(len(df))[::-1]
    for y, (_, r), c in zip(ys, df.iterrows(), colores):
        _barras_de_error(
            ax, y, r["valor"], r["ci_low"], r["ci_high"], color=c
        )
        ax.text(
            r["valor"],
            y + 0.26,
            f"{_coma3(r['valor'])}  (n={int(r['n'])})",
            fontsize=8,
            ha="center",
            color=c,
        )
    _linea_vertical(ax, 0.5, "azar", color="#999999", estilo=":", arriba=0.06)
    ax.set_yticks(ys)
    ax.set_yticklabels(df["calculo"], fontsize=8.5)
    ax.set_xlim(0.3, 1.04)
    _holgura(ax, len(df))
    ax.set_xlabel(
        "Fracción de proteínas cuyo pólipo queda entre el normal y el tumor",
        fontsize=8.5,
    )
    _estilo(ax, formato_x=_coma)
    ax.set_title(
        "C5. El gradiente se desploma cuando la selección es independiente",
        fontsize=10,
        loc="left",
        color="#222222",
    )
    png, svg = _guardar(fig, salida, "c5_gradiente_adenoma_carcinoma")
    return Figura(
        nombre="c5_gradiente_adenoma_carcinoma",
        titulo="El gradiente normal-adenoma-carcinoma con selección independiente",
        pie=(
            "Cada fila mide lo mismo con un diseño distinto. La primera elige y mide "
            "las proteínas en el mismo brazo, que es lo que hizo v1. La última las "
            "elige en CPTAC y las mide en la réplica, que es el análisis principal "
            "preregistrado. La tercera es el valor que saldría solo por regresión a "
            "la media. La selección independiente cae hasta el azar."
        ),
        png=png,
        svg=svg,
        tsv=tsv,
    )


# --------------------------------------------------------------------- C6
def fig_c6_modelos(axes: dict, salida: Path) -> Figura:
    """AUC externo por modelo con las dos unidades de remuestreo, y las diferencias."""
    c6 = axes["c6_seleccion_de_modelo"]
    pac = c6["bootstrap_por_paciente"]["por_modelo"]
    esp = c6["bootstrap_por_especimen_v1"]["por_modelo"]
    difs = c6["bootstrap_por_paciente"]["diferencias_por_pares"]
    nombres = {
        "logistica_l1": "Logística L1",
        "elastic_net": "Elastic net",
        "random_forest": "Bosque aleatorio",
    }

    filas = []
    for clave, bonito in nombres.items():
        for unidad, fuente in (("paciente", pac), ("espécimen (v1)", esp)):
            d = fuente[clave]
            filas.append(
                {
                    "modelo": bonito,
                    "unidad_de_remuestreo": unidad,
                    "auc_mediana": d["mediana"],
                    "ci_low": d["ci_low"],
                    "ci_high": d["ci_high"],
                    "ancho_del_intervalo": d["ci_high"] - d["ci_low"],
                    "auc_en_la_cohorte_externa": c6["auc_externo"][clave],
                }
            )
    df = pd.DataFrame(filas)

    filas_dif = []
    for clave, d in difs.items():
        a, b = clave.split("_menos_")
        filas_dif.append(
            {
                "comparacion": f"{nombres[a]} − {nombres[b]}",
                "diferencia_mediana": d["mediana"],
                "ci_low": d["ci_low"],
                "ci_high": d["ci_high"],
                "incluye_cero": d["incluye_cero"],
            }
        )
    df_dif = pd.DataFrame(filas_dif)
    tsv = _tabla(
        pd.concat(
            [
                df.assign(panel="auc_por_modelo"),
                df_dif.assign(panel="diferencia_por_pares"),
            ],
            ignore_index=True,
        ),
        salida,
        "c6_modelos_auc_externo",
    )

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(10.4, 3.4), gridspec_kw={"width_ratios": [1.1, 1], "wspace": 0.08}
    )

    ys = np.arange(len(nombres))[::-1]
    for i, (clave, bonito) in enumerate(nombres.items()):
        y = ys[i]
        for desp, unidad, color, marca in (
            (0.16, "paciente", AZUL, "o"),
            (-0.16, "espécimen (v1)", GRIS, "s"),
        ):
            d = (pac if unidad == "paciente" else esp)[clave]
            ax1.errorbar(
                d["mediana"],
                y + desp,
                xerr=[[d["mediana"] - d["ci_low"]], [d["ci_high"] - d["mediana"]]],
                fmt=marca,
                color=color,
                ecolor=color,
                elinewidth=1.3,
                capsize=3,
                markersize=5,
                label=unidad if i == 0 else None,
            )
    _linea_vertical(ax1, 0.5, "azar", color="#999999", estilo=":", arriba=0.07)
    ax1.set_yticks(ys)
    ax1.set_yticklabels(list(nombres.values()), fontsize=8.5)
    ax1.set_xlim(0.45, 1.02)
    _holgura(ax1, len(nombres))
    ax1.set_xlabel("AUC en la cohorte externa", fontsize=8.5)
    ax1.legend(fontsize=7.5, frameon=False, loc="upper left")
    _estilo(ax1, formato_x=_coma)
    ax1.set_title(
        "Agrupar por paciente no ensancha el intervalo", fontsize=9.5, loc="left"
    )

    ys2 = np.arange(len(df_dif))[::-1]
    for y, (_, r) in zip(ys2, df_dif.iterrows()):
        color = GRIS if r["incluye_cero"] else VERDE
        _barras_de_error(ax2, y, r["diferencia_mediana"], r["ci_low"], r["ci_high"], color=color)
        ax2.text(r["diferencia_mediana"], y + 0.22, _coma3(r["diferencia_mediana"]),
                 fontsize=7.5, ha="center", color=color)
    _linea_vertical(ax2, 0.0, "cero", color=ROJO, arriba=0.07)
    ax2.set_yticks(ys2)
    ax2.set_yticklabels(df_dif["comparacion"], fontsize=7.5)
    _filas_a_la_derecha(ax2)
    _holgura(ax2, len(df_dif))
    ax2.set_xlabel("Diferencia de AUC en el mismo remuestreo", fontsize=8.5)
    _estilo(ax2, formato_x=_coma)
    ax2.set_title("Dos de las tres diferencias excluyen el cero", fontsize=9.5, loc="left")

    fig.suptitle(
        "C6. Comparar modelos: remuestreo por paciente y diferencia por pares",
        fontsize=10,
        x=0.02,
        y=1.04,
        ha="left",
    )
    png, svg = _guardar(fig, salida, "c6_modelos_auc_externo")
    return Figura(
        nombre="c6_modelos_auc_externo",
        titulo="Selección de modelo con remuestreo por paciente",
        pie=(
            "Izquierda: el AUC de cada modelo con su intervalo, remuestreando "
            "pacientes y remuestreando especímenes. En un diseño emparejado el "
            "intervalo por paciente sale más estrecho, porque el efecto del paciente "
            "se cancela dentro del par. Derecha: la diferencia de AUC calculada en el "
            "mismo remuestreo. La regla de decisión preregistrada es el solape de los "
            "intervalos de la izquierda, y por solape no hay ganador."
        ),
        png=png,
        svg=svg,
        tsv=tsv,
    )


# --------------------------------------------------------------------- C4
def fig_c4_vias(vias: dict, salida: Path) -> Figura:
    """La cláusula de vías por colección: cuántos términos quedan y qué concuerda."""
    bonito = {"go_bp": "GO procesos", "reactome": "Reactome", "msigdb_hallmark": "Hallmark"}
    umbral = float(vias["umbral_concordancia_via"])
    n_min = int(vias["n_min_terminos"])

    filas = []
    for clave, etiqueta in bonito.items():
        col = vias["por_coleccion"][clave]
        conc = col["concordancia_entre_cohortes"]
        c = conc.get("concordancia")
        filas.append(
            {
                "coleccion": etiqueta,
                "n_terminos_en_la_coleccion": col["n_terminos_en_la_coleccion"],
                "n_significativos_en_cptac": conc["n_terminos_significativos_descubrimiento"],
                "n_evaluables_en_la_replica": conc["n_terminos_evaluados"],
                "estado": conc["estado"],
                "concordancia": c["valor"] if c else np.nan,
                "ci_low": c["ci_low"] if c else np.nan,
                "ci_high": c["ci_high"] if c else np.nan,
                "margen_del_ic_sobre_el_umbral": conc.get(
                    "margen_del_ic_inferior_sobre_el_umbral", np.nan
                ),
            }
        )
    df = pd.DataFrame(filas)
    tsv = _tabla(df, salida, "c4_clausula_de_vias")

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(10.0, 3.2), gridspec_kw={"width_ratios": [1, 1.05], "wspace": 0.1}
    )

    ys = np.arange(len(df))[::-1]
    alto = 0.3
    ax1.barh(ys + alto / 2, df["n_significativos_en_cptac"], height=alto, color=AZUL,
             label="significativos en CPTAC")
    ax1.barh(ys - alto / 2, df["n_evaluables_en_la_replica"], height=alto, color=NARANJA,
             label="evaluables en la réplica")
    ax1.set_xscale("symlog")
    ax1.set_xlim(0, 1400)
    _linea_vertical(ax1, n_min, f"mínimo = {n_min}", color=ROJO, arriba=0.99)
    for y, (_, r) in zip(ys, df.iterrows()):
        ax1.annotate(str(int(r["n_significativos_en_cptac"])),
                     xy=(r["n_significativos_en_cptac"], y + alto / 2),
                     xytext=(4, 0), textcoords="offset points",
                     fontsize=7.5, va="center", color=AZUL)
        ax1.annotate(str(int(r["n_evaluables_en_la_replica"])),
                     xy=(r["n_evaluables_en_la_replica"], y - alto / 2),
                     xytext=(4, 0), textcoords="offset points",
                     fontsize=7.5, va="center", color=NARANJA)
    ax1.set_yticks(ys)
    ax1.set_yticklabels(df["coleccion"], fontsize=8.5)
    _holgura(ax1, len(df))
    ax1.set_xlabel("Número de términos (escala logarítmica)", fontsize=8.5)
    ax1.legend(fontsize=7.5, frameon=False, loc="lower right")
    _estilo(ax1)
    ax1.set_title("La réplica deja muy pocos términos evaluables", fontsize=9.5, loc="left")

    for y, (_, r) in zip(ys, df.iterrows()):
        if np.isfinite(r["concordancia"]):
            _barras_de_error(ax2, y, r["concordancia"], r["ci_low"], r["ci_high"], color=VERDE)
            ax2.text(r["concordancia"], y + 0.22, _coma3(r["concordancia"]),
                     fontsize=8, ha="center", color=VERDE)
        else:
            ax2.text(0.56, y, "no evaluable", fontsize=8.5, ha="center", va="center",
                     color=GRIS, style="italic")
    _linea_vertical(ax2, umbral, f"umbral = {_coma(umbral)}", color=ROJO, arriba=0.99)
    ax2.set_yticks(ys)
    ax2.set_yticklabels(df["coleccion"], fontsize=8.5)
    _filas_a_la_derecha(ax2)
    ax2.set_xlim(0.3, 1.04)
    _holgura(ax2, len(df))
    ax2.set_xlabel("Concordancia de dirección entre cohortes", fontsize=8.5)
    _estilo(ax2, formato_x=_coma)
    ax2.set_title("Solo Reactome se puede evaluar, y por poco", fontsize=9.5, loc="left")

    fig.suptitle(
        "C4. La cláusula de vías con un mínimo de términos y prueba competitiva",
        fontsize=10, x=0.02, y=1.04, ha="left",
    )
    png, svg = _guardar(fig, salida, "c4_clausula_de_vias")
    return Figura(
        nombre="c4_clausula_de_vias",
        titulo="La cláusula de vías exige un mínimo de términos",
        pie=(
            "v1 declaró la cláusula satisfecha con cuatro términos y un intervalo que "
            "bajaba hasta 0,398. v2 exige diez términos como mínimo y que el límite "
            "inferior del intervalo supere el umbral. De las tres colecciones solo "
            "Reactome llega al mínimo, y lo supera por 0,011. El cuello de botella no "
            "es el resultado: es que la réplica mapea a símbolo de gen una parte "
            "pequeña de sus proteínas."
        ),
        png=png,
        svg=svg,
        tsv=tsv,
    )


# --------------------------------------------------------------------- C8
def fig_c8_deteccion(deteccion: dict, tabla_por_proteina: pd.DataFrame, salida: Path) -> Figura:
    """La curva de detección y la concordancia dentro y fuera de la ventana."""
    curva = deteccion["curva_de_deteccion"]
    conc = deteccion["concordancia"]
    umbral = float(curva["umbral_log10_pg_l"])
    b0 = float(curva["intercepto"])
    b1 = float(curva["pendiente_por_orden_de_magnitud"])

    datos = tabla_por_proteina.copy()
    col_conc = next(
        (c for c in ("log10_concentracion", "log10_concentracion_pg_l") if c in datos.columns),
        None,
    )
    col_det = next(
        (c for c in ("detectada_en_plasma", "detectada_en_suero") if c in datos.columns),
        None,
    )
    tiene_curva = col_conc is not None and col_det is not None

    filas = [
        {
            "grupo": "Todas",
            "valor": conc["todas"]["valor"],
            "ci_low": conc["todas"]["ci_low"],
            "ci_high": conc["todas"]["ci_high"],
            "numerador": conc["todas"]["numerador"],
            "denominador": conc["todas"]["denominador"],
        },
        {
            "grupo": "Dentro de la ventana",
            "valor": conc["dentro_de_la_ventana"]["valor"],
            "ci_low": conc["dentro_de_la_ventana"]["ci_low"],
            "ci_high": conc["dentro_de_la_ventana"]["ci_high"],
            "numerador": conc["dentro_de_la_ventana"]["numerador"],
            "denominador": conc["dentro_de_la_ventana"]["denominador"],
        },
        {
            "grupo": "Fuera de la ventana",
            "valor": conc["fuera_de_la_ventana"]["valor"],
            "ci_low": conc["fuera_de_la_ventana"]["ci_low"],
            "ci_high": conc["fuera_de_la_ventana"]["ci_high"],
            "numerador": conc["fuera_de_la_ventana"]["numerador"],
            "denominador": conc["fuera_de_la_ventana"]["denominador"],
        },
    ]
    df = pd.DataFrame(filas)
    tsv = _tabla(
        df.assign(
            umbral_log10_pg_l=umbral,
            pendiente_por_orden_de_magnitud=b1,
            z_de_la_pendiente=curva["z_pendiente"],
        ),
        salida,
        "c8_ventana_de_deteccion",
    )

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(10.0, 3.2), gridspec_kw={"wspace": 0.12}
    )

    if tiene_curva:
        d = datos.dropna(subset=[col_conc, col_det])
        x = pd.to_numeric(d[col_conc], errors="coerce")
        # El TSV trae los booleanos como texto; mapear es más seguro que castear,
        # porque un valor inesperado tiene que salir NaN y no colarse como 0.
        y = (
            d[col_det]
            .astype(str)
            .str.strip()
            .str.lower()
            .map({"true": 1.0, "false": 0.0, "1": 1.0, "0": 0.0, "1.0": 1.0, "0.0": 0.0})
        )
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        bordes = np.linspace(float(x.min()), float(x.max()), 13)
        centros, fracs, tams = [], [], []
        for lo, hi in zip(bordes[:-1], bordes[1:]):
            sel = (x >= lo) & (x < hi)
            if sel.sum() >= 5:
                centros.append((lo + hi) / 2)
                fracs.append(float(y[sel].mean()))
                tams.append(int(sel.sum()))
        if centros:
            ax1.scatter(centros, fracs, s=np.clip(np.array(tams) / 4, 10, 120),
                        color=AZUL, alpha=0.75, zorder=3,
                        label="observado (por tramo)")
        rejilla = np.linspace(float(x.min()), float(x.max()), 200)
        ax1.plot(rejilla, 1 / (1 + np.exp(-(b0 + b1 * rejilla))), color=ROJO, lw=1.6,
                 label="logística ajustada")
    _linea_vertical(ax1, umbral, "umbral", color=VERDE, arriba=0.52)
    ax1.set_xlabel("Concentración en sangre, log₁₀ de pg/L", fontsize=8.5)
    ax1.set_ylabel("Probabilidad de cuantificarse en suero", fontsize=8.5)
    ax1.set_ylim(-0.04, 1.08)
    ax1.legend(fontsize=7.5, frameon=False, loc="upper left")
    _estilo(ax1, formato_y=_coma)
    ax1.set_title(
        f"La detección depende de la concentración (z = {_coma(curva['z_pendiente'])})",
        fontsize=9.5, loc="left",
    )

    ys = np.arange(len(df))[::-1]
    for y_, (_, r) in zip(ys, df.iterrows()):
        _barras_de_error(ax2, y_, r["valor"], r["ci_low"], r["ci_high"], color=NARANJA)
        ax2.text(r["valor"], y_ + 0.24,
                 f"{_coma3(r['valor'])}  ({int(r['numerador'])}/{int(r['denominador'])})",
                 fontsize=8, ha="center", color=NARANJA)
    _linea_vertical(ax2, 0.5, "azar", color="#999999", estilo=":", arriba=0.07)
    ax2.set_yticks(ys)
    ax2.set_yticklabels(df["grupo"], fontsize=8.5)
    _filas_a_la_derecha(ax2)
    ax2.set_xlim(0.05, 0.85)
    _holgura(ax2, len(df))
    ax2.set_xlabel("Concordancia de dirección tejido-plasma", fontsize=8.5)
    _estilo(ax2, formato_x=_coma)
    ax2.set_title("La ventana no rescata la concordancia", fontsize=9.5, loc="left")

    fig.suptitle(
        "C8. La brecha tejido-plasma no se explica por el límite de detección",
        fontsize=10, x=0.02, y=1.04, ha="left",
    )
    png, svg = _guardar(fig, salida, "c8_ventana_de_deteccion")
    return Figura(
        nombre="c8_ventana_de_deteccion",
        titulo="La ventana de detección no explica la brecha tejido-plasma",
        pie=(
            "Izquierda: la probabilidad de que una proteína se cuantifique en suero "
            "crece con su concentración en sangre, así que el experimento tiene una "
            "ventana de detección real. Derecha: restringir la concordancia a las "
            "proteínas que caen dentro de esa ventana no la separa del azar. Eso "
            "convierte la brecha en un resultado y no en un artefacto, que es una "
            "afirmación más fuerte que la de v1 y mejor acotada."
        ),
        png=png,
        svg=svg,
        tsv=tsv,
    )


# --------------------------------------------------------------------- C10
def fig_c10_heterogeneidad(
    meta: dict, por_proteina: pd.DataFrame, solo_lfq: pd.DataFrame, salida: Path
) -> Figura:
    """Distribución de I² con las tres cohortes y sin la de TMT."""
    diag = meta["diagnostico_de_tecnologia"]
    col_lfq = "i2_solo_lfq" if "i2_solo_lfq" in solo_lfq.columns else "i2"
    # Los dos histogramas tienen que ser de las mismas proteínas. Dibujar las 185
    # con al menos dos cohortes contra las 101 con tres daría una mediana distinta
    # de la del certificado, y el pie de figura dejaría de cuadrar con el texto.
    juntas = por_proteina.merge(
        solo_lfq[["uniprot", col_lfq]], on="uniprot", how="inner"
    )
    tres = pd.to_numeric(juntas.get("i2"), errors="coerce").dropna()
    lfq = pd.to_numeric(juntas.get(col_lfq), errors="coerce").dropna()

    df = pd.DataFrame(
        [
            {
                "conjunto": "Las tres cohortes",
                "n_proteinas": int(len(tres)),
                "i2_mediano": float(tres.median()) if len(tres) else np.nan,
                "i2_p25": float(tres.quantile(0.25)) if len(tres) else np.nan,
                "i2_p75": float(tres.quantile(0.75)) if len(tres) else np.nan,
            },
            {
                "conjunto": "Solo las dos de LFQ",
                "n_proteinas": int(len(lfq)),
                "i2_mediano": float(lfq.median()) if len(lfq) else np.nan,
                "i2_p25": float(lfq.quantile(0.25)) if len(lfq) else np.nan,
                "i2_p75": float(lfq.quantile(0.75)) if len(lfq) else np.nan,
            },
        ]
    ).assign(
        caida_al_quitar_cptac=diag["caida_al_quitar_cptac"],
        r3_del_certificado=meta["mediana_1_menos_i2"],
    )
    tsv = _tabla(df, salida, "c10_heterogeneidad_i2")

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(10.0, 3.2),
        gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.12},
    )

    bins = np.linspace(0, 1, 26)
    ax1.hist(tres, bins=bins, color=AZUL, alpha=0.72,
             label=f"las tres cohortes (n={len(tres)})")
    ax1.hist(lfq, bins=bins, color=NARANJA, alpha=0.62,
             label=f"las mismas, sin CPTAC (n={len(lfq)})")
    for serie, color in ((tres, AZUL), (lfq, NARANJA)):
        if len(serie):
            ax1.axvline(float(serie.median()), color=color, ls="--", lw=1.2)
    ax1.set_xlabel("I² por proteína", fontsize=8.5)
    ax1.set_ylabel("Número de proteínas", fontsize=8.5)
    ax1.legend(fontsize=7.5, frameon=False, loc="upper left")
    _estilo(ax1, formato_x=_coma)
    ax1.set_title("La heterogeneidad se concentra cerca de 1", fontsize=9.5, loc="left")

    ys = np.array([1, 0])
    for y_, (_, r) in zip(ys, df.iterrows()):
        _barras_de_error(ax2, y_, r["i2_mediano"], r["i2_p25"], r["i2_p75"],
                         color=AZUL if y_ == 1 else NARANJA)
        ax2.text(r["i2_mediano"], y_ + 0.2, _coma3(r["i2_mediano"]), fontsize=8,
                 ha="center", color=AZUL if y_ == 1 else NARANJA)
    ax2.set_yticks(ys)
    ax2.set_yticklabels(df["conjunto"], fontsize=8.5)
    _filas_a_la_derecha(ax2)
    ax2.set_xlim(0.4, 1.02)
    ax2.set_ylim(-0.6, 1.5)
    ax2.set_xlabel("I² mediano, con el rango intercuartílico", fontsize=8.5)
    _estilo(ax2, formato_x=_coma)
    ax2.set_title(
        f"Quitar la cohorte de TMT solo baja {_coma3(diag['caida_al_quitar_cptac'])}",
        fontsize=9.5, loc="left",
    )

    fig.suptitle(
        "C10. R3 ya es calculable, y dice que las cohortes no concuerdan en magnitud",
        fontsize=10, x=0.02, y=1.04, ha="left",
    )
    png, svg = _guardar(fig, salida, "c10_heterogeneidad_i2")
    return Figura(
        nombre="c10_heterogeneidad_i2",
        titulo="Heterogeneidad entre cohortes y de dónde viene",
        pie=(
            "Con tres cohortes de tejido R3 pasa de no calculable a tener un valor, y "
            "ese valor es bajo: el I² mediano está cerca de 1, es decir que las tres "
            "cohortes casi no concuerdan en la magnitud del efecto. El panel derecho "
            "descarta la explicación fácil: al quitar la cohorte de marcaje isobárico "
            "y dejar solo las dos de LFQ, la heterogeneidad apenas baja. No es la "
            "tecnología. El histograma naranja además muestra por qué ese cálculo es "
            "solo un diagnóstico: con dos cohortes I² se va a cero o a uno y casi "
            "nada queda en medio, así que su mediana es un número inestable."
        ),
        png=png,
        svg=svg,
        tsv=tsv,
    )


# --------------------------------------------------------------------- C11
def fig_c11_independencia(
    meta: dict, correlaciones: dict[str, pd.DataFrame], salida: Path
) -> Figura:
    """Las correlaciones entre muestras de depósitos distintos, con el umbral."""
    ind = meta["independencia_de_cohortes"]
    umbral = float(ind["umbral_de_sospecha"])

    filas = []
    series: dict[str, pd.Series] = {}
    for par, tabla in correlaciones.items():
        col = next((c for c in ("r", "correlacion", "pearson") if c in tabla.columns), None)
        if col is None:
            continue
        r = pd.to_numeric(tabla[col], errors="coerce").dropna()
        if r.empty:
            continue
        series[par] = r
        filas.append(
            {
                "par_de_depositos": par,
                "n_pares_de_muestras": int(len(r)),
                "r_maximo": float(r.max()),
                "r_mediano": float(r.median()),
                "umbral_de_sospecha": umbral,
                "n_por_encima_del_umbral": int((r >= umbral).sum()),
            }
        )
    df = pd.DataFrame(filas) if filas else pd.DataFrame(
        [
            {
                "par_de_depositos": "resumen de la auditoría",
                "n_pares_de_muestras": np.nan,
                "r_maximo": ind["r_maximo_entre_muestras"],
                "r_mediano": ind["r_mediano_entre_muestras"],
                "umbral_de_sospecha": umbral,
                "n_por_encima_del_umbral": 0,
            }
        ]
    )
    tsv = _tabla(df, salida, "c11_independencia_de_cohortes")

    fig, ax = plt.subplots(figsize=(7.4, 3.2))
    if series:
        bins = np.linspace(
            min(float(s.min()) for s in series.values()), 1.0, 40
        )
        for (par, r), color in zip(series.items(), (AZUL, NARANJA, VERDE, GRIS)):
            ax.hist(r, bins=bins, alpha=0.6, color=color,
                    label=f"{par} (n={len(r)})")
    _linea_vertical(ax, umbral, f"sospecha ≥ {_coma(umbral)}", color=ROJO, arriba=0.99)
    _linea_vertical(
        ax,
        float(df["r_maximo"].max()),
        f"máximo {_coma3(df['r_maximo'].max())}",
        color="#444444",
        estilo=":",
        arriba=0.80,
    )
    ax.set_xlabel("Correlación entre una muestra de un depósito y una del otro", fontsize=8.5)
    ax.set_ylabel("Número de pares de muestras", fontsize=8.5)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    _estilo(ax, formato_x=_coma)
    ax.set_title(
        "C11. Ninguna muestra se repite entre depósitos: son cohortes distintas",
        fontsize=10, loc="left",
    )
    png, svg = _guardar(fig, salida, "c11_independencia_de_cohortes")
    return Figura(
        nombre="c11_independencia_de_cohortes",
        titulo="Dos niveles para comprobar que las cohortes son independientes",
        pie=(
            "Comparar los checksums de PRIDE no basta: dos depósitos pueden ser "
            "ficheros distintos del mismo paciente. El segundo nivel correlaciona "
            "cada muestra de un depósito con cada muestra del otro sobre las "
            "proteínas comunes. Ningún par llega al umbral de sospecha, así que "
            "PXD019504 entra como tercera cohorte. PXD031705 sí compartía ficheros "
            "con PXD019504 y quedó fuera."
        ),
        png=png,
        svg=svg,
        tsv=tsv,
    )


# --------------------------------------------------------------------- índice R
NOMBRES_DEL_CHANGELOG = {
    "R1_concordancia_de_signo": "R1, concordancia de signo",
    "R2_detectabilidad": "R2, detectabilidad",
    "R4_composicion": "R4, sobrevive al ajuste",
    "R5_p_empirico": "R5, p del nulo",
    "concordancia_de_via": "Concordancia de vías",
    "polipo_magnitud_intermedia": "Pólipo a medio camino",
    "concordancia_tejido_plasma": "Concordancia tejido-plasma",
}


def fig_indice_r(changelog: pd.DataFrame, indice: pd.DataFrame, salida: Path) -> Figura:
    """Los cinco subíndices de v2 y, al lado, qué números cambiaron frente a v1.

    La primera versión de esta figura superponía v1 y v2 en el mismo eje. No
    servía: cuatro de los cinco subíndices casi no se mueven, los marcadores se
    tapaban y el cambio importante, que R3 pasara de no calculable a tener un
    valor, no se veía porque un valor ausente no se dibuja. Lo que de verdad
    cambia está en el changelog, así que el changelog es el panel.
    """
    df = indice.copy()
    orden = ["R1", "R2", "R3", "R4", "R5"]
    filas = []
    for pref in orden:
        fila = df.loc[df["codigo"].astype(str).str.startswith(pref)]
        if fila.empty:
            continue
        r = fila.iloc[0]
        filas.append(
            {
                "subindice": str(r["codigo"]).replace("_", " "),
                "valor": pd.to_numeric(r.get("valor"), errors="coerce"),
                "numerador": pd.to_numeric(r.get("numerador"), errors="coerce"),
                "denominador": pd.to_numeric(r.get("denominador"), errors="coerce"),
                "ci_low": pd.to_numeric(r.get("ci_low"), errors="coerce"),
                "ci_high": pd.to_numeric(r.get("ci_high"), errors="coerce"),
                "estado": str(r.get("estado", "")),
                "calculable_en_v1": pref != "R3",
            }
        )
    tabla = pd.DataFrame(filas)

    cl = changelog.copy()
    if {"numero", "v1", "v2"} <= set(cl.columns):
        cl["v1_num"] = pd.to_numeric(cl["v1"], errors="coerce")
        cl["v2_num"] = pd.to_numeric(cl["v2"], errors="coerce")
        cl = cl.loc[
            cl["v1_num"].between(0, 1) & cl["v2_num"].between(0, 1)
            & (cl["v1_num"] != cl["v2_num"])
        ]
        cl["etiqueta"] = cl["numero"].map(
            lambda v: NOMBRES_DEL_CHANGELOG.get(str(v), str(v).replace("_", " "))
        )
        cl = cl[["numero", "etiqueta", "v1_num", "v2_num", "correccion"]]
    else:
        cl = pd.DataFrame(columns=["numero", "etiqueta", "v1_num", "v2_num", "correccion"])

    tsv = _tabla(
        pd.concat(
            [
                tabla.assign(panel="subindices_de_v2"),
                cl.rename(columns={"v1_num": "v1", "v2_num": "v2"}).assign(
                    panel="numeros_que_cambian"
                ),
            ],
            ignore_index=True,
        ),
        salida,
        "indice_r_y_changelog",
    )

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(10.6, 3.5),
        gridspec_kw={"width_ratios": [1, 1.15], "wspace": 0.1},
    )

    ys = np.arange(len(tabla))[::-1]
    for y_, (_, r) in zip(ys, tabla.iterrows()):
        tope = r["estado"] == "calculado_como_tope"
        color = NARANJA if tope else AZUL
        _barras_de_error(ax1, y_, r["valor"], r["ci_low"], r["ci_high"], color=color)
        detalle = _coma_fino(r["valor"])
        if np.isfinite(r["numerador"]) and np.isfinite(r["denominador"]):
            detalle += f"  ({int(r['numerador'])}/{int(r['denominador'])})"
        elif np.isfinite(r["denominador"]):
            detalle += f"  (n={int(r['denominador'])})"
        if tope:
            detalle += "  tope"
        ax1.text(r["valor"], y_ + 0.22, detalle, fontsize=7.5, ha="center", color=color)
        if not r["calculable_en_v1"]:
            ax1.annotate(
                "no calculable en v1",
                xy=(r["valor"], y_),
                xytext=(14, -11),
                textcoords="offset points",
                fontsize=7,
                color=VERDE,
                style="italic",
            )
    ax1.set_yticks(ys)
    ax1.set_yticklabels(tabla["subindice"], fontsize=8.5)
    ax1.set_xlim(-0.05, 1.05)
    _holgura(ax1, len(tabla))
    ax1.set_xlabel("Valor en v2, con su intervalo", fontsize=8.5)
    _estilo(ax1, formato_x=_coma)
    ax1.set_title("Los cinco subíndices no se promedian", fontsize=9.5, loc="left")

    ys2 = np.arange(len(cl))[::-1]
    for y_, (_, r) in zip(ys2, cl.iterrows()):
        v1, v2 = float(r["v1_num"]), float(r["v2_num"])
        ax2.annotate(
            "",
            xy=(v2, y_),
            xytext=(v1, y_),
            arrowprops=dict(arrowstyle="-|>", color="#c9c9c9", lw=1.4,
                            shrinkA=4, shrinkB=4),
        )
        ax2.plot(v1, y_, "s", color=GRIS, markersize=5.5)
        ax2.plot(v2, y_, "o", color=ROJO if v2 < v1 else VERDE, markersize=6.5)
        ax2.annotate(
            f"{_coma_fino(v1)} → {_coma_fino(v2)}  ({r['correccion']})",
            xy=(min(v1, v2), y_),
            xytext=(0, 8),
            textcoords="offset points",
            fontsize=7,
            color="#444444",
        )
    ax2.plot([], [], "s", color=GRIS, label="v1")
    ax2.plot([], [], "o", color=ROJO, label="v2, más bajo")
    ax2.plot([], [], "o", color=VERDE, label="v2, más alto")
    ax2.set_yticks(ys2)
    ax2.set_yticklabels(cl["etiqueta"], fontsize=8)
    _filas_a_la_derecha(ax2)
    ax2.set_xlim(-0.05, 1.12)
    _holgura(ax2, len(cl))
    ax2.set_xlabel("Valor del número", fontsize=8.5)
    ax2.legend(fontsize=7, frameon=False, loc="upper left")
    _estilo(ax2, formato_x=_coma)
    ax2.set_title("Los números que cambian al corregir el método", fontsize=9.5, loc="left")

    fig.suptitle(
        "El índice R de v2 y el changelog frente a v1",
        fontsize=10, x=0.02, y=1.04, ha="left",
    )
    png, svg = _guardar(fig, salida, "indice_r_y_changelog")
    return Figura(
        nombre="indice_r_y_changelog",
        titulo="El índice R de v2 y qué números cambian frente a v1",
        pie=(
            "Izquierda: cada subíndice responde a una pregunta distinta, así que "
            "colapsarlos en un número daría un valor sin interpretación. R3 era no "
            "calculable en v1 porque hacían falta tres cohortes de tejido y solo "
            "había dos. R5 figura como tope y no como estimación: el nulo no produjo "
            "ninguna permutación que igualara lo observado, así que el valor es el "
            "mínimo alcanzable con ese número de permutaciones. Derecha: los números "
            "que se mueven al corregir el método, con el código de la corrección que "
            "los movió. El changelog completo incluye además dos números que no "
            "caben en este eje porque no son fracciones: los términos de vía pasan "
            "de 4 a 51, y el veredicto de H1 cambia de etiqueta."
        ),
        png=png,
        svg=svg,
        tsv=tsv,
    )


# --------------------------------------------------------------------- C7 y C9
def fig_c7_c9_composicion(certified: dict, salida: Path) -> Figura:
    """Qué sobrevive al ajuste por composición y qué le pasa a R1."""
    ajustes = certified["c7_composicion"]["ajustes"]
    bonito = {
        "ecm_leave_one_out": "Matriz extracelular\n(sin la proteína)",
        "ecm_sin_exclusion": "Matriz extracelular\n(circular, v1)",
        "musculo_liso": "Músculo liso",
        "inmune": "Inmune",
        "epitelial": "Epitelial",
        "conjunto_multiindice": "Las cuatro juntas",
    }
    filas = []
    for a in ajustes:
        f = a["fraccion"]
        filas.append(
            {
                "ajuste": bonito.get(a["indice"], a["indice"]),
                "clave": a["indice"],
                "sobreviven": f["numerador"],
                "certificadas": f["denominador"],
                "fraccion": f["valor"],
                "ci_low": f["ci_low"],
                "ci_high": f["ci_high"],
                "es_principal": bool(a.get("es_principal", False)),
            }
        )
    df = pd.DataFrame(filas)

    tras = certified["c15_r1_direccion"]["tras_el_ajuste_de_composicion"]
    base = certified["c15_r1_direccion"]["global"]
    filas_r1 = [
        {
            "calculo": "R1 sin ajustar",
            "valor": base["valor"],
            "ci_low": base["ci_low"],
            "ci_high": base["ci_high"],
            "numerador": base["numerador"],
            "denominador": base["denominador"],
        }
    ]
    if tras.get("estado") == "calculado":
        f = tras["fraccion"]
        filas_r1.append(
            {
                "calculo": "R1 con el efecto ajustado",
                "valor": f["valor"],
                "ci_low": f["ci_low"],
                "ci_high": f["ci_high"],
                "numerador": f["numerador"],
                "denominador": f["denominador"],
            }
        )
    df_r1 = pd.DataFrame(filas_r1)
    tsv = _tabla(
        pd.concat(
            [df.assign(panel="sobrevive_el_ajuste"), df_r1.assign(panel="r1_tras_el_ajuste")],
            ignore_index=True,
        ),
        salida,
        "c7_c9_composicion",
    )

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(10.4, 3.6),
        gridspec_kw={"width_ratios": [1.3, 1], "wspace": 0.1},
    )

    ys = np.arange(len(df))[::-1]
    for y_, (_, r) in zip(ys, df.iterrows()):
        color = ROJO if r["es_principal"] else GRIS
        _barras_de_error(ax1, y_, r["fraccion"], r["ci_low"], r["ci_high"], color=color)
        ax1.text(r["fraccion"], y_ + 0.26,
                 f"{int(r['sobreviven'])}/{int(r['certificadas'])}",
                 fontsize=7.5, ha="center", color=color)
    ax1.set_yticks(ys)
    ax1.set_yticklabels(df["ajuste"], fontsize=7.5)
    ax1.set_xlim(0.0, 1.02)
    _holgura(ax1, len(df))
    ax1.set_xlabel("Fracción de las certificadas que sigue significativa", fontsize=8.5)
    _estilo(ax1, formato_x=_coma)
    ax1.set_title("C7. En rojo, el método principal preregistrado", fontsize=9.5, loc="left")

    ys2 = np.arange(len(df_r1))[::-1]
    for y_, (_, r) in zip(ys2, df_r1.iterrows()):
        _barras_de_error(ax2, y_, r["valor"], r["ci_low"], r["ci_high"],
                         color=GRIS if y_ == 1 else AZUL)
        ax2.text(r["valor"], y_ + 0.2,
                 f"{_coma3(r['valor'])}  ({int(r['numerador'])}/{int(r['denominador'])})",
                 fontsize=8, ha="center", color=GRIS if y_ == 1 else AZUL)
    _linea_vertical(ax2, 0.70, "umbral de H1", color=ROJO, arriba=0.12)
    ax2.set_yticks(ys2)
    ax2.set_yticklabels(df_r1["calculo"], fontsize=8.5)
    _filas_a_la_derecha(ax2)
    ax2.set_xlim(0.55, 1.03)
    ax2.set_ylim(-0.7, len(df_r1) - 0.25)
    ax2.set_xlabel("Concordancia de dirección entre cohortes", fontsize=8.5)
    _estilo(ax2, formato_x=_coma)
    ax2.set_title("C9. Ajustar no baja la concordancia", fontsize=9.5, loc="left")

    fig.suptitle(
        "C7 y C9. Composición del tejido como covariable, sin circularidad",
        fontsize=10, x=0.02, y=1.04, ha="left",
    )
    png, svg = _guardar(fig, salida, "c7_c9_composicion")
    return Figura(
        nombre="c7_c9_composicion",
        titulo="Composición del tejido: qué sobrevive y qué le pasa a R1",
        pie=(
            "Izquierda: al ajustar por el índice estromal, recalculado sin la propia "
            "proteína, menos de la mitad de los hallazgos certificados sigue "
            "significativa. Eso es una cota del efecto de composición, no una pureza "
            "medida. Derecha: la concordancia de dirección entre cohortes no baja al "
            "usar el efecto ajustado, aunque el denominador se reduce. La cláusula de "
            "proteína de H1 exige bajar de 0,70 y no lo hace ni antes ni después."
        ),
        png=png,
        svg=svg,
        tsv=tsv,
    )
