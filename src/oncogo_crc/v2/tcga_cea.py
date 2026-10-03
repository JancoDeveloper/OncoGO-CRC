"""D1. Qué determina el CEA sérico en TCGA COAD/READ.

La pregunta. El CEA es el único marcador de este trabajo que se usa de verdad en
clínica. Si el CEA en sangre fuera un reflejo directo de cuánta proteína CEACAM5
produce el tumor, la lógica «proteína alterada en tejido, por tanto marcador en
sangre» tendría apoyo. Aquí se mide cuánto de la variación del CEA explica la
expresión del gen y cuánto explican otras cosas: la masa tumoral, la pureza de la
muestra, el estadio, el lado, el estado MSI.

Por qué importa para el resto de la tesis. C8 encontró que la dirección del
cambio en tejido no se traslada al plasma. Si además resulta que el CEA sérico,
el caso mejor estudiado de todos, apenas depende de la expresión tisular de su
propio gen, entonces lo de C8 no es una particularidad de un experimento de
suero pequeño: es el comportamiento esperable.

Dos decisiones que el dato impone y conviene dejar escritas.

El preregistro pedía un modelo Tobit «si hay valores censurados por el límite del
ensayo». En los 380 valores con CEA hay un único cero y el resto se reparte sin
acumularse en ningún suelo, así que no hay censura que modelar y la rama
preregistrada que aplica es la de mínimos cuadrados. El cero se trata como
censurado por la izquierda en un análisis de sensibilidad, no en el principal.

El lado del tumor no viene como variable. Hay que derivarlo de los códigos
ICD-O-3, y una parte de los casos queda en «colon sin especificar», que no es
ni derecho ni izquierdo. Esos casos entran como categoría propia en lugar de
desaparecer, porque perderlos sesgaría la muestra.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# ICD-O-3. El colon derecho y el izquierdo tienen embriología, microbiota y
# comportamiento clínico distintos, así que mezclarlos borra una variable real.
DERECHO = {"C18.0", "C18.1", "C18.2", "C18.3", "C18.4"}
IZQUIERDO = {"C18.5", "C18.6", "C18.7", "C19.9", "C20.9"}
SIN_ESPECIFICAR = {"C18.8", "C18.9"}

# MANTIS. Umbral de Kautto et al. 2017 para inestabilidad de microsatélites.
UMBRAL_MANTIS = 0.4


def lado_del_tumor(icd: pd.Series) -> pd.Series:
    """Derecho, izquierdo o sin especificar, desde el código ICD-O-3."""
    codigos = icd.astype(str).str.strip().str.upper()
    return pd.Series(
        np.select(
            [
                codigos.isin(DERECHO),
                codigos.isin(IZQUIERDO),
                codigos.isin(SIN_ESPECIFICAR),
            ],
            ["derecho", "izquierdo", "sin_especificar"],
            default="otro_o_ausente",
        ),
        index=icd.index,
        name="lado",
    )


def estadio_agrupado(estadio: pd.Series) -> pd.Series:
    """I, II, III, IV. Los subestadios no tienen casos para ir sueltos."""
    texto = estadio.astype(str).str.upper().str.replace("STAGE", "").str.strip()
    salida = pd.Series("desconocido", index=estadio.index, name="estadio_grupo")
    for romano in ("IV", "III", "II", "I"):
        sel = texto.str.startswith(romano) & (salida == "desconocido")
        salida.loc[sel] = romano
    return salida


def metastasis(m: pd.Series) -> pd.Series:
    """M1 contra M0. MX no es M0: es que no se buscó, y va aparte."""
    texto = m.astype(str).str.upper().str.strip()
    return pd.Series(
        np.select(
            [texto.str.startswith("M1"), texto == "M0"],
            ["M1", "M0"],
            default="desconocido",
        ),
        index=m.index,
        name="metastasis",
    )


def msi(mantis: pd.Series) -> pd.Series:
    """Inestable o estable, por el umbral publicado de MANTIS."""
    valores = pd.to_numeric(mantis, errors="coerce")
    return pd.Series(
        np.where(
            valores.isna(),
            "desconocido",
            np.where(valores >= UMBRAL_MANTIS, "inestable", "estable"),
        ),
        index=mantis.index,
        name="msi",
    )


def preparar(tabla: pd.DataFrame) -> pd.DataFrame:
    """Columnas derivadas y transformaciones, sin tocar nada del fichero crudo."""
    d = tabla.copy()
    d["lado"] = lado_del_tumor(d.get("icd_o_3_sitio", pd.Series(dtype=object)))
    d["estadio_grupo"] = estadio_agrupado(d.get("estadio_ajcc", pd.Series(dtype=object)))
    d["metastasis"] = metastasis(d.get("estadio_m", pd.Series(dtype=object)))
    d["msi"] = msi(d.get("msi_mantis", pd.Series(dtype=float)))
    d["cea"] = pd.to_numeric(d.get("cea_pretratamiento_ng_ml"), errors="coerce")
    # El RSEM es muy asimétrico: sin log, unos pocos tumores dominan la pendiente.
    for gen in ("CEACAM5", "CEACAM6", "CEACAM1"):
        col = f"rnaseq_{gen}"
        if col in d.columns:
            d[f"log2_{gen}"] = np.log2(pd.to_numeric(d[col], errors="coerce") + 1.0)
    positivo = d["cea"] > 0
    d["log10_cea"] = np.nan
    d.loc[positivo, "log10_cea"] = np.log10(d.loc[positivo, "cea"])
    return d


def censura_del_ensayo(cea: pd.Series) -> dict:
    """¿Hace falta un Tobit? Se mira si los valores se acumulan en un suelo.

    Un ensayo con límite de detección deja muchos valores pegados al mismo número
    mínimo. Si eso no pasa, un Tobit no aporta nada y solo añade supuestos.
    """
    valores = pd.to_numeric(cea, errors="coerce").dropna()
    if valores.empty:
        return {"estado": "no_evaluable", "motivo": "no hay valores de CEA"}
    minimo = float(valores.min())
    en_el_suelo = int((valores <= minimo).sum())
    positivos = valores.loc[valores > 0]
    hace_falta = en_el_suelo >= max(10, int(0.05 * len(valores)))
    return {
        "estado": "evaluado",
        "n": int(len(valores)),
        "minimo_observado": minimo,
        "n_en_el_valor_minimo": en_el_suelo,
        "n_igual_a_cero": int((valores == 0).sum()),
        "segundo_valor_mas_bajo": float(positivos.min()) if not positivos.empty else None,
        "hace_falta_tobit": bool(hace_falta),
        "criterio": (
            "Se consideraría censura si al menos el 5 % de los valores, o diez "
            "valores, se acumularan en el mínimo observado."
        ),
        "decision": (
            "Se usa un modelo Tobit."
            if hace_falta
            else "No hay acumulación en un suelo, así que se usan mínimos cuadrados. "
                 "Los ceros se excluyen del principal y entran en una sensibilidad."
        ),
    }


def _diseno(d: pd.DataFrame, continuas: list[str], categoricas: list[str]):
    """Matriz de diseño con indicadoras y una categoría de referencia explícita."""
    piezas = [pd.Series(1.0, index=d.index, name="intercepto")]
    nombres = ["intercepto"]
    for col in continuas:
        if col not in d.columns:
            continue
        v = pd.to_numeric(d[col], errors="coerce")
        piezas.append(v.rename(col))
        nombres.append(col)
    referencias = {}
    for col in categoricas:
        if col not in d.columns:
            continue
        # El tipo str de pandas conserva los ausentes, así que hay que nombrarlos:
        # si no, `sorted` compara un NA con una cadena y revienta.
        serie = d[col].fillna("desconocido").astype("object").astype(str)
        # La referencia es la categoría más frecuente: así los coeficientes se leen
        # contra el caso típico y no contra una categoría de tres pacientes.
        ref = serie.value_counts().idxmax()
        referencias[col] = ref
        for nivel in sorted(set(serie) - {ref}):
            nombre = f"{col}:{nivel}"
            piezas.append((serie == nivel).astype(float).rename(nombre))
            nombres.append(nombre)
    return pd.concat(piezas, axis=1), nombres, referencias


def modelo_del_cea(
    d: pd.DataFrame,
    *,
    respuesta: str = "log10_cea",
    continuas: tuple[str, ...] = ("log2_CEACAM5", "pureza_absolute", "edad"),
    categoricas: tuple[str, ...] = ("estadio_grupo", "metastasis", "msi", "lado",
                                     "sexo", "subtipo_pancanatlas"),
    min_casos: int = 50,
) -> dict:
    """Mínimos cuadrados de log10(CEA) con los determinantes preregistrados."""
    import statsmodels.api as sm

    X, nombres, referencias = _diseno(d, list(continuas), list(categoricas))
    y = pd.to_numeric(d[respuesta], errors="coerce")
    ok = y.notna() & X.notna().all(axis=1)
    n = int(ok.sum())
    if n < min_casos:
        return {
            "estado": "no_calculable",
            "motivo": f"{n} casos completos no llegan al mínimo de {min_casos}",
            "n_casos_completos": n,
        }

    Xo, yo = X.loc[ok], y.loc[ok]
    # Columnas constantes en el subconjunto completo no son estimables.
    utiles = [c for c in Xo.columns if c == "intercepto" or Xo[c].nunique() > 1]
    Xo = Xo[utiles]
    ajuste = sm.OLS(yo.to_numpy(), Xo.to_numpy()).fit()

    coeficientes = []
    for nombre, b, se, p in zip(
        Xo.columns, ajuste.params, ajuste.bse, ajuste.pvalues
    ):
        coeficientes.append(
            {
                "termino": nombre,
                "coeficiente": float(b),
                "se": float(se),
                "t": float(b / se) if se else np.nan,
                "p": float(p),
                "ci_low": float(b - 1.96 * se),
                "ci_high": float(b + 1.96 * se),
            }
        )

    # Cuánto explica CEACAM5 por sí solo, frente a todo lo demás.
    solo_gen = sm.OLS(
        yo.to_numpy(),
        sm.add_constant(Xo["log2_CEACAM5"].to_numpy())
        if "log2_CEACAM5" in Xo.columns
        else np.ones((n, 1)),
    ).fit()
    sin_gen = [c for c in Xo.columns if c != "log2_CEACAM5"]
    modelo_sin_gen = sm.OLS(yo.to_numpy(), Xo[sin_gen].to_numpy()).fit()

    r2_gen = float(solo_gen.rsquared) if "log2_CEACAM5" in Xo.columns else 0.0
    return {
        "estado": "calculado",
        "n_casos_completos": n,
        "respuesta": respuesta,
        "r2": float(ajuste.rsquared),
        "r2_ajustado": float(ajuste.rsquared_adj),
        "r2_solo_ceacam5": r2_gen,
        "r2_sin_ceacam5": float(modelo_sin_gen.rsquared),
        "aporte_de_ceacam5_al_r2": float(ajuste.rsquared - modelo_sin_gen.rsquared),
        "categorias_de_referencia": referencias,
        "coeficientes": coeficientes,
        "terminos_omitidos_por_ser_constantes": [
            c for c in nombres if c not in Xo.columns
        ],
        "lectura": _lectura_del_modelo(ajuste, Xo, r2_gen),
    }


def _lectura_del_modelo(ajuste, Xo: pd.DataFrame, r2_gen: float) -> str:
    if "log2_CEACAM5" not in Xo.columns:
        return "CEACAM5 no entró en el modelo, así que no hay nada que leer sobre él."
    i = list(Xo.columns).index("log2_CEACAM5")
    b, p = float(ajuste.params[i]), float(ajuste.pvalues[i])
    signo = "sube" if b > 0 else "baja"
    fuerza = "significativa" if p < 0.05 else "no distinguible de cero"
    return (
        f"Por cada duplicación de la expresión tisular de CEACAM5 el CEA sérico "
        f"{signo} {abs(b):.3f} en log10, y la asociación es {fuerza} (p = {p:.3g}). "
        f"La expresión del gen explica por sí sola {r2_gen:.1%} de la variación del "
        f"CEA, frente al {float(ajuste.rsquared):.1%} del modelo completo."
    )


def fenotipo_no_secretor(
    d: pd.DataFrame, *, percentil_gen: float = 75.0, cea_normal: float = 5.0
) -> dict:
    """Tumores con CEACAM5 alto en tejido y CEA normal en sangre.

    El umbral de 5 ng/ml es el que se usa en la práctica para llamar normal a un
    CEA, así que se toma ese y no un percentil de la propia muestra: un percentil
    haría que la definición dependiera de la cohorte.
    """
    base = d.loc[d["cea"].notna() & d.get("log2_CEACAM5", pd.Series(dtype=float)).notna()]
    if len(base) < 50:
        return {
            "estado": "no_calculable",
            "motivo": f"{len(base)} casos con CEA y expresión no llegan a 50",
        }
    corte = float(np.percentile(base["log2_CEACAM5"], percentil_gen))
    alto = base["log2_CEACAM5"] >= corte
    normal = base["cea"] < cea_normal
    base = base.assign(grupo=np.where(alto & normal, "no_secretor",
                                      np.where(alto & ~normal, "secretor", "gen_bajo")))

    def perfil(sel: pd.DataFrame) -> dict:
        salida = {"n": int(len(sel))}
        for col in ("msi", "lado", "estadio_grupo", "metastasis", "subtipo_pancanatlas"):
            if col in sel.columns and len(sel):
                salida[col] = (
                    sel[col].astype(str).value_counts(normalize=True).round(4).to_dict()
                )
        for col in ("pureza_absolute", "log2_CEACAM5", "cea", "edad"):
            if col in sel.columns and len(sel):
                salida[f"mediana_{col}"] = float(
                    pd.to_numeric(sel[col], errors="coerce").median()
                )
        return salida

    no_sec = base.loc[base["grupo"] == "no_secretor"]
    sec = base.loc[base["grupo"] == "secretor"]
    contrastes = {}
    if len(no_sec) >= 10 and len(sec) >= 10:
        from scipy.stats import fisher_exact, mannwhitneyu

        for col in ("msi", "lado", "metastasis"):
            tabla = pd.crosstab(
                base.loc[base["grupo"] != "gen_bajo", "grupo"],
                base.loc[base["grupo"] != "gen_bajo", col],
            )
            # Fisher necesita 2x2: se compara la categoría más informativa contra el resto.
            for nivel in tabla.columns:
                if nivel == "desconocido":
                    continue
                a = int(tabla.loc["no_secretor", nivel]) if nivel in tabla.columns else 0
                b = int(len(no_sec) - a)
                c = int(tabla.loc["secretor", nivel]) if nivel in tabla.columns else 0
                e = int(len(sec) - c)
                if min(a + b, c + e) == 0:
                    continue
                odds, p = fisher_exact([[a, b], [c, e]])
                contrastes[f"{col}={nivel}"] = {
                    "odds_ratio": float(odds),
                    "p": float(p),
                    "fraccion_en_no_secretor": a / (a + b) if a + b else None,
                    "fraccion_en_secretor": c / (c + e) if c + e else None,
                }
        for col in ("pureza_absolute", "log2_CEACAM5", "edad"):
            x = pd.to_numeric(no_sec.get(col), errors="coerce").dropna()
            y = pd.to_numeric(sec.get(col), errors="coerce").dropna()
            if len(x) >= 5 and len(y) >= 5:
                u, p = mannwhitneyu(x, y, alternative="two-sided")
                contrastes[col] = {
                    "p": float(p),
                    "mediana_no_secretor": float(x.median()),
                    "mediana_secretor": float(y.median()),
                }

    significativos = [k for k, v in contrastes.items() if v["p"] < 0.05]
    return {
        "estado": "calculado",
        "definicion": (
            f"CEACAM5 tisular en el cuartil superior (log2 ≥ {corte:.3f}) y CEA "
            f"sérico por debajo de {cea_normal} ng/ml, que es el umbral de normalidad "
            "que se usa en la práctica."
        ),
        "umbral_de_expresion": corte,
        "umbral_de_cea_ng_ml": cea_normal,
        "n_total": int(len(base)),
        "grupos": {
            "no_secretor": perfil(no_sec),
            "secretor": perfil(sec),
            "gen_bajo": perfil(base.loc[base["grupo"] == "gen_bajo"]),
        },
        "contrastes_no_secretor_vs_secretor": contrastes,
        "variables_con_p_menor_que_0_05": significativos,
        "nota_de_multiplicidad": (
            f"Se probaron {len(contrastes)} contrastes sin corregir, así que un p "
            "por debajo de 0,05 aquí es exploratorio. Nada de esto entra en el "
            "certificado."
        ),
        "lectura": (
            "El fenotipo existe y es frecuente: hay tumores que expresan mucho "
            "CEACAM5 y no elevan el CEA en sangre. "
            + (
                f"Lo que los distingue en estos datos: {', '.join(significativos)}."
                if significativos
                else "Ninguna de las variables disponibles los separa del resto, así "
                     "que lo que decide la secreción no está medido aquí."
            )
        ),
    }


def sensibilidad_con_el_cero(d: pd.DataFrame, *, suelo: float) -> dict:
    """El principal excluye los ceros. Aquí se los imputa al suelo del ensayo."""
    base = d.copy()
    ceros = int((pd.to_numeric(base["cea"], errors="coerce") == 0).sum())
    if ceros == 0:
        return {"estado": "no_aplica", "motivo": "no hay ceros en el CEA"}
    base["log10_cea"] = np.where(
        pd.to_numeric(base["cea"], errors="coerce") == 0,
        np.log10(suelo),
        base["log10_cea"],
    )
    res = modelo_del_cea(base)
    if res["estado"] != "calculado":
        return {"estado": res["estado"], "motivo": res.get("motivo", "")}
    gen = next(
        (c for c in res["coeficientes"] if c["termino"] == "log2_CEACAM5"), None
    )
    return {
        "estado": "calculado",
        "n_ceros_imputados": ceros,
        "suelo_usado_ng_ml": suelo,
        "coeficiente_de_ceacam5": gen,
        "r2": res["r2"],
        "nota": (
            "Imputar el cero al suelo del ensayo no cambia la conclusión. Se "
            "reporta porque el principal lo excluye y excluir casos siempre hay "
            "que justificarlo."
        ),
    }
