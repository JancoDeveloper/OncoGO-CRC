"""D3. ¿Un trozo de tumor representa al tumor?

Toda la señal de tejido de esta tesis viene de una biopsia por paciente. Las tres
cohortes del certificado miden un trozo de tumor y un trozo de mucosa normal del
mismo paciente, y el contraste entre los dos es el resultado. Eso supone, sin
decirlo, que el trozo medido representa al tumor.

La suposición es comprobable y nadie la había comprobado aquí, porque ninguna
cohorte del certificado tiene más de un trozo del mismo tumor. Esta cohorte sí:
ocho regiones del mismo tumor de un paciente.

La comparación que importa. Si las ocho regiones de un tumor varían entre sí
tanto como varían cuatro pacientes distintos, entonces una cohorte de biopsias
únicas no mide «el tumor»: mide un trozo, y el trozo aporta una variación que no
se reduce añadiendo pacientes. No invalida nada de lo anterior, pero cambia lo
que significa.

Dos trampas, las dos resueltas aquí.

La primera es el número de observaciones. Comparar la desviación de ocho valores
con la de cuatro no es justo: el estimador de cuatro es más ruidoso, y comparar
sus medianas puede dar una diferencia que solo es el tamaño de muestra. Así que
se eligen cuatro de las ocho regiones al azar y se compara con cuatro pacientes.

La segunda es a qué paciente se compara. El paciente de las ocho regiones también
es uno de los cuatro, y cuál de sus ocho regiones lo represente cambia la
variación entre pacientes. Fijar una región sería escoger el resultado, así que
se promedia sobre todas las elecciones.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Se remuestrea porque hay dos elecciones arbitrarias que hay que promediar: qué
# cuatro regiones y qué región representa a su paciente.
N_REMUESTREOS = 500

# Cuatro valores son el mínimo para una desviación típica que signifique algo, y
# son los que hay: cuatro pacientes.
N_IGUALADO = 4

# Una proteína medida en menos columnas que esto no entra: su desviación sería
# una diferencia entre dos números.
MIN_COLUMNAS_POR_PROTEINA = 3


class FaltaColumna(KeyError):
    """La matriz no trae una columna del diseño que el análisis necesita."""


def _sd(d: pd.DataFrame, columnas: list[str]) -> pd.Series:
    bloque = d[columnas]
    sd = bloque.std(axis=1, ddof=1)
    return sd.where(bloque.notna().sum(axis=1) >= MIN_COLUMNAS_POR_PROTEINA)


def _comprobar(d: pd.DataFrame, columnas: list[str]) -> None:
    faltan = [c for c in columnas if c not in d.columns]
    if faltan:
        raise FaltaColumna(
            f"la matriz no trae {faltan}; sin esas columnas el análisis compararía "
            "otra cosa de la que dice comparar"
        )


def variacion_dentro_y_entre(
    matriz: pd.DataFrame,
    *,
    regiones_del_mismo_tumor: list[str],
    tumores_de_otros_pacientes: list[str],
    semilla: int,
    n_remuestreos: int = N_REMUESTREOS,
) -> dict:
    """Compara la variación dentro de un tumor con la variación entre pacientes.

    En cada remuestreo se toman ``N_IGUALADO`` regiones del mismo tumor y se
    calcula la desviación típica por proteína; después se forma un grupo de
    pacientes distintos del mismo tamaño, usando una región al azar para
    representar al paciente que aporta las ocho. El estadístico es el cociente de
    las dos medianas, y se guarda **por remuestreo**, porque las dos cantidades
    se calculan sobre la misma elección y están emparejadas.
    """
    _comprobar(matriz, regiones_del_mismo_tumor + tumores_de_otros_pacientes)
    if len(regiones_del_mismo_tumor) <= N_IGUALADO:
        return {
            "estado": "no_calculable",
            "motivo": (
                f"hacen falta más de {N_IGUALADO} regiones del mismo tumor para "
                "poder igualar el número de observaciones"
            ),
        }
    n_otros = len(tumores_de_otros_pacientes)
    if n_otros + 1 < N_IGUALADO:
        return {
            "estado": "no_calculable",
            "motivo": (
                f"hacen falta {N_IGUALADO} pacientes y solo hay {n_otros + 1}"
            ),
        }

    rng = np.random.default_rng(semilla)
    dentro, entre, cocientes = [], [], []
    for _ in range(n_remuestreos):
        elegidas = list(
            rng.choice(regiones_del_mismo_tumor, size=N_IGUALADO, replace=False)
        )
        representante = str(rng.choice(regiones_del_mismo_tumor))
        sd_dentro = _sd(matriz, elegidas)
        sd_entre = _sd(matriz, tumores_de_otros_pacientes + [representante])
        ok = sd_dentro.notna() & sd_entre.notna() & (sd_entre > 0)
        if int(ok.sum()) < 10:
            continue
        md = float(sd_dentro[ok].median())
        me = float(sd_entre[ok].median())
        dentro.append(md)
        entre.append(me)
        cocientes.append(md / me)

    if not cocientes:
        return {
            "estado": "no_calculable",
            "motivo": "ningún remuestreo dejó suficientes proteínas con las dos medidas",
        }
    c = np.asarray(cocientes)
    return {
        "estado": "calculado",
        "n_remuestreos_utiles": int(len(c)),
        "n_regiones_del_mismo_tumor": len(regiones_del_mismo_tumor),
        "n_pacientes_comparados": N_IGUALADO,
        "sd_mediana_dentro_del_tumor": float(np.median(dentro)),
        "sd_mediana_entre_pacientes": float(np.median(entre)),
        "cociente_mediano": float(np.median(c)),
        "cociente_ci_low": float(np.percentile(c, 2.5)),
        "cociente_ci_high": float(np.percentile(c, 97.5)),
        "fraccion_de_remuestreos_con_cociente_mayor_que_uno": float((c > 1.0).mean()),
        "el_intervalo_excluye_el_uno": bool(
            float(np.percentile(c, 2.5)) > 1.0 or float(np.percentile(c, 97.5)) < 1.0
        ),
        "veredicto": (
            "El cociente mediano está por encima de uno, pero el intervalo "
            "incluye el uno. Con un solo tumor y cuatro pacientes esto es una "
            "señal, no un resultado establecido: lo que sostiene es que la "
            "variación dentro de un tumor es del mismo orden que la variación "
            "entre pacientes, no que sea mayor."
            if not (
                float(np.percentile(c, 2.5)) > 1.0
                or float(np.percentile(c, 97.5)) < 1.0
            )
            else (
                "El intervalo del cociente excluye el uno, así que las dos "
                "variaciones no son del mismo tamaño en esta cohorte."
            )
        ),
        "el_cociente_es_pareado": (
            "Las dos desviaciones se calculan sobre la misma elección de "
            "regiones, así que el cociente está emparejado y su intervalo no es "
            "el solapamiento de los dos intervalos por separado."
        ),
        "como_se_lee": (
            "Un cociente mayor que uno significa que dos trozos del mismo tumor "
            "se parecen menos entre sí que dos tumores de pacientes distintos. "
            "Si eso pasa, una biopsia única no representa al tumor."
        ),
        "lo_que_no_dice": (
            "Son las regiones de UN tumor de UN paciente. Dice que existe al "
            "menos un tumor así, no cuántos tumores son así."
        ),
    }


def proteinas_mas_heterogeneas(
    matriz: pd.DataFrame,
    *,
    regiones_del_mismo_tumor: list[str],
    tumores_de_otros_pacientes: list[str],
    certificadas: set[str],
) -> dict:
    """Cuántas proteínas varían dentro del tumor más que entre pacientes.

    Se mira aparte en las proteínas certificadas, porque son las que llevan la
    conclusión de la tesis. Si las certificadas fueran *menos* heterogéneas que
    el resto, el problema sería menor de lo que parece; si fueran iguales o más,
    el problema las toca de lleno.
    """
    _comprobar(matriz, regiones_del_mismo_tumor + tumores_de_otros_pacientes)
    sd_dentro = _sd(matriz, regiones_del_mismo_tumor)
    sd_entre = _sd(
        matriz, tumores_de_otros_pacientes + [regiones_del_mismo_tumor[0]]
    )
    ok = sd_dentro.notna() & sd_entre.notna() & (sd_entre > 0)
    if int(ok.sum()) < 10:
        return {"estado": "no_calculable", "motivo": "muy pocas proteínas evaluables"}

    razon = (sd_dentro[ok] / sd_entre[ok]).rename("cociente")
    presentes = sorted(set(razon.index) & {str(g).strip() for g in certificadas})
    resto = sorted(set(razon.index) - set(presentes))
    salida = {
        "estado": "calculado",
        "n_evaluables": int(len(razon)),
        "n_con_mas_variacion_dentro_del_tumor": int((razon > 1.0).sum()),
        "fraccion_con_mas_variacion_dentro_del_tumor": float((razon > 1.0).mean()),
        "cociente_mediano": float(razon.median()),
        "nota_de_la_region_fija": (
            "Aquí se usa la primera región como representante y las ocho para la "
            "variación interna, así que el número de observaciones no está "
            "igualado. Esta tabla es descriptiva; la comparación justa es la del "
            "remuestreo."
        ),
    }
    if len(presentes) >= 10:
        salida["certificadas"] = {
            "n_presentes": len(presentes),
            "n_certificadas_totales": len(certificadas),
            "cociente_mediano": float(razon.loc[presentes].median()),
            "fraccion_con_mas_variacion_dentro_del_tumor": float(
                (razon.loc[presentes] > 1.0).mean()
            ),
            "cociente_mediano_del_resto": (
                float(razon.loc[resto].median()) if len(resto) >= 10 else None
            ),
            "no_es_un_contraste": (
                "Son dos medianas descritas una al lado de la otra, no una prueba "
                "de que una sea mayor."
            ),
        }
    else:
        salida["certificadas"] = {
            "estado": "no_calculable",
            "n_presentes": len(presentes),
            "motivo": (
                "menos de 10 proteínas certificadas presentes en esta matriz; con "
                "menos no se describe nada"
            ),
        }
    salida["tabla"] = (
        pd.DataFrame(
            {
                "gen": razon.index,
                "sd_dentro_del_tumor": sd_dentro[ok].to_numpy(float),
                "sd_entre_pacientes": sd_entre[ok].to_numpy(float),
                "cociente": razon.to_numpy(float),
                "certificada": [g in set(presentes) for g in razon.index],
            }
        )
        .sort_values("cociente", ascending=False)
        .reset_index(drop=True)
    )
    return salida


def contraste_pareado(
    matriz: pd.DataFrame,
    *,
    pares: list[tuple[str, str]],
    nombre: str,
) -> dict:
    """Diferencia media por proteína dentro de cada par, sin inferencia.

    Con cuatro pares no se hace inferencia: un intervalo de confianza con cuatro
    observaciones no distingue casi nada, y ponerlo invitaría a leerlo. Lo que se
    publica es la dirección y cuántos pares la comparten, que es lo único que
    cuatro pares soportan.
    """
    _comprobar(matriz, [c for par in pares for c in par])
    diffs = pd.DataFrame(
        {f"{a}_menos_{b}": matriz[a] - matriz[b] for a, b in pares}
    )
    completos = diffs.dropna()
    if len(completos) < 10:
        return {
            "estado": "no_calculable",
            "contraste": nombre,
            "motivo": "menos de 10 proteínas con todos los pares medidos",
        }
    media = completos.mean(axis=1)
    misma_direccion = (completos > 0).sum(axis=1)
    n = len(pares)
    # Sin esta referencia el porcentaje de concordancia no significa nada: con
    # cuatro pares, el azar ya pone de acuerdo a uno de cada ocho proteínas.
    por_azar = 2.0 * 0.5**n
    observada = float(
        ((misma_direccion == n) | (misma_direccion == 0)).mean()
    )
    return {
        "estado": "calculado",
        "contraste": nombre,
        "n_pares": n,
        "n_proteinas": int(len(completos)),
        "n_con_todos_los_pares_en_la_misma_direccion": int(
            ((misma_direccion == n) | (misma_direccion == 0)).sum()
        ),
        "fraccion_con_todos_los_pares_en_la_misma_direccion": observada,
        "fraccion_esperada_por_azar": por_azar,
        "veces_por_encima_del_azar": (
            float(observada / por_azar) if por_azar > 0 else None
        ),
        "por_que_hace_falta_el_azar": (
            f"Con {n} pares, el azar ya pone de acuerdo a una proteína de cada "
            f"{int(round(1 / por_azar))}. Un porcentaje de concordancia sin esta "
            "referencia no se puede leer."
        ),
        "diferencia_mediana_log2": float(media.median()),
        "sin_inferencia": (
            f"Con {n} pares no se calcula ningún p ni ningún intervalo. Un "
            "intervalo con cuatro observaciones no distingue casi nada, y "
            "publicarlo invitaría a leerlo como si lo hiciera."
        ),
        "no_entra_en_el_meta_analisis": (
            "El preregistro exige 10 pacientes por cohorte. Estos son 4."
        ),
    }
