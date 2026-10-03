"""D2. ¿El índice estromal de R4 mide pureza, o mide otra cosa?

R4 no puede medir pureza. La pureza de un tumor es la fracción de células
tumorales, y eso se mide con ADN. Lo que R4 tiene es un índice construido con las
proteínas de matriz extracelular del propio proteoma, y por eso el preregistro
obliga a llamarlo «una cota del efecto de composición, no una pureza medida,
salvo cuando la covariable es la pureza ABSOLUTE de TCGA».

Ese «salvo» nunca se había podido ejercer, porque los pacientes de PDC000116 no
son pacientes de TCGA. Esta cohorte sí lo es. Así que aquí se puede hacer la
única comprobación que da sentido a la frase del preregistro: calcular el índice
en un proteoma de pacientes que tienen pureza medida por ADN, y ver cuánto se
parecen.

Qué se puede concluir y qué no. Si el índice correlaciona fuerte con la pureza,
el ajuste de R4 está haciendo lo que se dice que hace. Si correlaciona flojo, R4
sigue siendo un análisis de sensibilidad legítimo pero no se puede presentar como
un control de pureza, y el lenguaje obligatorio del preregistro deja de ser una
cautela y pasa a ser una descripción necesaria.

Lo que no se puede concluir es nada sobre los valores concretos de R4. Esta
cohorte es label free y la del certificado es TMT; lo que se valida es el
constructo, no los números.

Y un segundo asunto, que cierra lo que D1 dejó abierto: D1 midió el CEACAM5 del
tejido con mRNA porque es lo único que TCGA tiene, y encontró que explica el
2,4 % de la varianza del CEA en sangre. Aquí está la proteína, en los mismos
pacientes. Pero antes de interpretar cualquier correlación hay que saber cuánto
ruido tiene la medida, y para eso están las alícuotas repetidas: cinco de estos
tumores se midieron dos veces.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from oncogo_crc.v2.composicion import construir_indice

# El término de matriz extracelular es el que usa R4. No se cambia aquí: validar
# un índice distinto del que se usa no validaría nada.
ECM_GO_TERM = "GO:0031012"

# Un índice con menos proteínas que esto no es un índice, es un puñado de genes.
MIN_PROTEINAS_DE_FIRMA = 3

# Por debajo de esto una correlación no se interpreta: el intervalo abarca de
# «nada» a «mucho» y cualquier lectura sería ganas de leer.
MIN_PARA_CORRELAR = 20

# Los genes del CEA y sus parientes. CEACAM5 es la proteína que el análisis de
# sangre llama «CEA»; los otros dos entran porque comparten epítopos y porque D1
# los midió.
GENES_CEA = ["CEACAM5", "CEACAM6", "CEACAM1"]


class FaltaColumna(KeyError):
    """La tabla de entrada no trae una columna que el análisis necesita."""


def _correlacion(x: pd.Series, y: pd.Series) -> dict:
    """Spearman y Pearson, con el n y el intervalo de la de Spearman.

    Se dan las dos a propósito. Spearman no supone forma y es la que decide;
    Pearson se publica al lado porque si las dos discrepan mucho, es que hay unos
    pocos valores extremos mandando, y eso el lector tiene que poder verlo.
    """
    ok = x.notna() & y.notna()
    n = int(ok.sum())
    if n < MIN_PARA_CORRELAR:
        return {
            "estado": "no_calculable",
            "n": n,
            "motivo": f"menos de {MIN_PARA_CORRELAR} pacientes con las dos medidas",
        }
    xi, yi = x[ok].to_numpy(float), y[ok].to_numpy(float)
    rho, p_rho = stats.spearmanr(xi, yi)
    r, p_r = stats.pearsonr(xi, yi)
    # Intervalo de Fisher sobre la rho. Es una aproximación, y con n de 80 va
    # bien; se declara como aproximación porque sobre rangos lo es dos veces.
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3) if n > 3 else np.nan
    return {
        "estado": "calculado",
        "n": n,
        "rho": float(rho),
        "p_rho": float(p_rho),
        "rho_ci_low": float(np.tanh(z - 1.959964 * se)) if np.isfinite(se) else None,
        "rho_ci_high": float(np.tanh(z + 1.959964 * se)) if np.isfinite(se) else None,
        "varianza_compartida": float(rho**2),
        "r_pearson": float(r),
        "p_pearson": float(p_r),
        "intervalo_es_aproximado": (
            "El intervalo es la transformación de Fisher sobre un coeficiente de "
            "rangos. Es una aproximación."
        ),
    }


def genes_del_termino(go_long: pd.DataFrame, termino: str) -> set[str]:
    """Los genes anotados en un término de GO, en mayúsculas."""
    for c in ("go_id", "gene_symbol"):
        if c not in go_long.columns:
            raise FaltaColumna(
                f"la tabla de GO no tiene la columna «{c}»; sin ella la firma "
                "saldría vacía y el índice parecería imposible de construir"
            )
    return {
        str(g).strip().upper()
        for g in go_long.loc[go_long["go_id"].astype(str) == termino, "gene_symbol"]
    }


def por_paciente(matriz: pd.DataFrame, alicuotas: pd.DataFrame) -> pd.DataFrame:
    """Promedia las alícuotas del mismo paciente.

    Hay cinco pacientes con dos alícuotas. Dejarlas como dos filas los contaría
    dos veces en cualquier correlación, y las correlaciones suponen
    observaciones independientes. El promedio es la forma más simple de no
    hacerlo, y además reduce el ruido de esos cinco.
    """
    for c in ("alicuota", "paciente"):
        if c not in alicuotas.columns:
            raise FaltaColumna(f"la tabla de alícuotas no tiene la columna «{c}»")
    mapa = dict(zip(alicuotas["alicuota"], alicuotas["paciente"], strict=False))
    faltan = [c for c in matriz.columns if c not in mapa]
    if faltan:
        raise FaltaColumna(
            f"{len(faltan)} alícuotas de la matriz no están en la tabla de "
            f"alícuotas (la primera es «{faltan[0]}»); el cruce con la clínica "
            "saldría incompleto sin avisar"
        )
    t = matriz.T.copy()
    t["paciente"] = [mapa[c] for c in matriz.columns]
    return t.groupby("paciente").mean(numeric_only=True).T


def reproducibilidad(matriz: pd.DataFrame, alicuotas: pd.DataFrame) -> dict:
    """Cuánto se parecen dos medidas del mismo tumor.

    Es lo que decide si una correlación débil se puede interpretar. Si el ensayo
    repetido sobre el mismo tumor no se reproduce, la correlación de esa medida
    con cualquier otra cosa está atenuada y un cero no significa «no hay
    relación», significa «con esta medida no se puede saber».
    """
    mapa = dict(zip(alicuotas["alicuota"], alicuotas["paciente"], strict=False))
    grupos: dict[str, list[str]] = {}
    for c in matriz.columns:
        grupos.setdefault(mapa.get(c, c), []).append(c)
    pares = [(v[0], v[1]) for v in grupos.values() if len(v) >= 2]
    if not pares:
        return {
            "estado": "no_calculable",
            "motivo": "no hay ningún tumor medido dos veces",
            "consecuencia": (
                "Sin réplicas no se puede saber cuánto ruido tiene el ensayo, "
                "así que ninguna correlación débil se puede interpretar."
            ),
        }
    rs, ns = [], []
    for a, b in pares:
        ok = matriz[a].notna() & matriz[b].notna()
        if int(ok.sum()) < MIN_PARA_CORRELAR:
            continue
        r, _ = stats.pearsonr(matriz.loc[ok, a], matriz.loc[ok, b])
        rs.append(float(r))
        ns.append(int(ok.sum()))
    return {
        "estado": "calculado" if rs else "no_calculable",
        "n_pares_de_replicas": len(pares),
        "r_por_par": rs,
        "r_mediana": float(np.median(rs)) if rs else None,
        "n_genes_por_par": ns,
        "que_significa": (
            "Es la correlación entre dos medidas del mismo tumor sobre todos los "
            "genes a la vez. Dice si el perfil global se reproduce, no si una "
            "proteína concreta lo hace."
        ),
    }


def reproducibilidad_de_un_gen(
    matriz: pd.DataFrame, alicuotas: pd.DataFrame, gen: str
) -> dict:
    """Las dos medidas del mismo tumor para una proteína, y su dispersión.

    El perfil global puede reproducirse bien y una proteína concreta fatal. Para
    interpretar una correlación de esa proteína con otra cosa hace falta esto y
    no la media de todas.
    """
    if gen not in matriz.index:
        return {"estado": "no_medida", "gen": gen}
    mapa = dict(zip(alicuotas["alicuota"], alicuotas["paciente"], strict=False))
    grupos: dict[str, list[str]] = {}
    for c in matriz.columns:
        grupos.setdefault(mapa.get(c, c), []).append(c)
    difs, pares = [], 0
    for cols in grupos.values():
        if len(cols) < 2:
            continue
        a, b = matriz.loc[gen, cols[0]], matriz.loc[gen, cols[1]]
        if pd.isna(a) or pd.isna(b):
            continue
        pares += 1
        difs.append(abs(float(a) - float(b)))
    if not difs:
        return {
            "estado": "no_calculable",
            "gen": gen,
            "motivo": "no hay ningún tumor con esta proteína medida dos veces",
        }
    entre_pacientes = float(matriz.loc[gen].std(ddof=1))
    return {
        "estado": "calculado",
        "gen": gen,
        "n_pares": pares,
        "diferencias_absolutas_log2": [round(x, 3) for x in difs],
        "diferencia_mediana_log2": float(np.median(difs)),
        "veces_de_diferencia_mediana": float(2 ** np.median(difs)),
        "desviacion_entre_pacientes_log2": entre_pacientes,
        "ruido_frente_a_senal": (
            float(np.median(difs) / entre_pacientes)
            if entre_pacientes and np.isfinite(entre_pacientes)
            else None
        ),
        "como_se_lee": (
            "Si la diferencia entre dos medidas del mismo tumor es del orden de "
            "la variación entre pacientes distintos, esta proteína no está "
            "medida con precisión suficiente para correlacionarla con nada. Un "
            "coeficiente cercano a cero no diría que no hay relación: diría que "
            "con esta medida no se puede saber."
        ),
    }


def validar_el_indice(
    matriz: pd.DataFrame,
    alicuotas: pd.DataFrame,
    clinica: pd.DataFrame,
    firmas: dict[str, set[str]],
) -> dict:
    """Correlaciona cada índice de composición con la pureza medida por ADN.

    El signo esperado es negativo: más estroma, menos fracción de células
    tumorales. Un signo positivo no sería un detalle, sería que el índice no mide
    lo que se cree.
    """
    if "pureza_absolute" not in clinica.columns:
        raise FaltaColumna(
            "la tabla clínica no trae pureza_absolute; sin pureza medida esta "
            "etapa no tiene nada que validar y debe quedar no_disponible"
        )
    prot = por_paciente(matriz, alicuotas)
    comunes = [p for p in prot.columns if p in set(clinica.index)]
    pureza = clinica.loc[comunes, "pureza_absolute"]

    salida: dict[str, dict] = {}
    for nombre, genes in firmas.items():
        idx = construir_indice(matriz, list(matriz.columns), genes, nombre=nombre)
        if idx is None:
            salida[nombre] = {
                "estado": "no_calculable",
                "motivo": (
                    f"la firma no llega a {MIN_PROTEINAS_DE_FIRMA} proteínas "
                    "presentes en esta matriz"
                ),
            }
            continue
        # El índice se calcula por alícuota y luego se promedia por paciente, en
        # ese orden: estandarizar después de promediar cambiaría la escala.
        s = pd.DataFrame({"indice": idx.score})
        mapa = dict(zip(alicuotas["alicuota"], alicuotas["paciente"], strict=False))
        s["paciente"] = [mapa.get(i, i) for i in s.index]
        s = s.groupby("paciente")["indice"].mean()
        cor = _correlacion(s.reindex(comunes), pureza)
        cor["n_proteinas_de_la_firma"] = int(idx.m)
        cor["signo_esperado"] = "negativo"
        if cor["estado"] == "calculado":
            cor["signo_es_el_esperado"] = bool(cor["rho"] < 0)
        salida[nombre] = cor
    return salida


def cea_de_proteina(
    matriz: pd.DataFrame,
    alicuotas: pd.DataFrame,
    clinica: pd.DataFrame,
) -> dict:
    """El CEACAM5 del tumor frente al CEA en sangre, proteína y transcrito.

    La comparación entre los dos se hace sobre **los mismos pacientes**. Hecha
    sobre dos subconjuntos distintos no diría si una medida informa más que la
    otra: diría que los grupos no eran iguales.
    """
    for c in ("cea_pretratamiento_ng_ml",):
        if c not in clinica.columns:
            raise FaltaColumna(f"la tabla clínica no trae «{c}»")
    prot = por_paciente(matriz, alicuotas)
    comunes = [p for p in prot.columns if p in set(clinica.index)]
    cl = clinica.loc[comunes]

    cea = cl["cea_pretratamiento_ng_ml"]
    # El cero no se transforma. En D1 pasó lo mismo y la decisión es la misma:
    # log10 de cero no existe y rellenarlo con un mínimo inventaría el valor.
    positivo = cea.notna() & (cea > 0)
    log_cea = pd.Series(np.nan, index=cl.index)
    log_cea[positivo] = np.log10(cea[positivo].to_numpy(float))

    salida: dict[str, dict] = {
        "n_con_cea_positivo": int(positivo.sum()),
        "ceros_de_cea_no_transformados": int((cea == 0).sum()),
        "proteina": {},
    }
    for gen in GENES_CEA:
        if gen not in prot.index:
            salida["proteina"][gen] = {
                "estado": "no_medida",
                "nota": (
                    "Esta proteína no está en la matriz. No se dice nada sobre "
                    "ella: no medida no es lo mismo que no relacionada."
                ),
            }
            continue
        v = prot.loc[gen, comunes]
        cor = _correlacion(v, log_cea)
        cor["n_pacientes_con_la_proteina_detectada"] = int(v.notna().sum())
        salida["proteina"][gen] = cor

    # Proteína frente a transcrito, en los mismos pacientes y con las dos medidas
    # presentes. Es la comparación que decide si la proteína informa más.
    if "rnaseq_CEACAM5" in cl.columns and "CEACAM5" in prot.index:
        p5 = prot.loc["CEACAM5", comunes]
        m5 = np.log2(cl["rnaseq_CEACAM5"].astype(float) + 1.0)
        ambos = p5.notna() & m5.notna() & log_cea.notna()
        salida["proteina_frente_a_transcrito"] = {
            "n_con_las_tres_medidas": int(ambos.sum()),
            "acuerdo_entre_transcrito_y_proteina": _correlacion(p5, m5),
            "proteina_vs_cea": _correlacion(p5[ambos], log_cea[ambos]),
            "transcrito_vs_cea": _correlacion(m5[ambos], log_cea[ambos]),
            "la_comparacion_es_sobre_los_mismos_pacientes": True,
            "no_es_un_contraste": (
                "Son dos correlaciones sobre los mismos pacientes, no una prueba "
                "de que una sea mayor que la otra. Para eso haría falta "
                "contrastar dos coeficientes dependientes, y con este tamaño el "
                "contraste no tendría potencia para distinguirlos."
            ),
        }
    return salida
