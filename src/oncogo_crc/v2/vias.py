"""Cláusula de vías de H1 con n mínimo y prueba competitiva (C4, C13).

Qué falló en v1. La concordancia de vía se declaró 1,0 y la cláusula se dio por
satisfecha. Detrás había **cuatro** términos y un intervalo de confianza que iba
de 0,398 a 1,0. Con cuatro términos, un intervalo así no distingue una vía
perfectamente concordante de una moneda. Tres cosas se corrigen aquí.

1. **Mínimo de términos**. Por debajo de diez términos evaluables la cláusula
   queda ``no_evaluable``. No queda satisfecha ni rechazada: queda sin evaluar,
   y eso se propaga a la tabla de verdad de H1.
2. **Prueba competitiva**. v1 usaba sobrerrepresentación de Fisher sobre la lista
   de significativas y después tomaba la mediana del log2FC de los miembros como
   dirección. La sobrerrepresentación ignora la correlación entre proteínas de la
   misma vía y la mediana de los miembros no es un estadístico de conjunto. Aquí
   se usa una prueba competitiva al estilo de ``camera``, con el conjunto contra
   el resto y la varianza inflada por la correlación intergénica.
3. **Tres colecciones**. GO Biological Process, Reactome y MSigDB Hallmark. Una
   conclusión que solo aparece en una de las tres no es una conclusión robusta.

El nulo de C13 permuta qué proteína pertenece a qué término, que es la etapa que
el nulo de inversión de signos de R5 no cubría.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.v2.conjuntos import (
    ResultadoConjunto,
    camera,
    nulo_etiquetas_de_gen,
)
from oncogo_crc.v2.metricas import Fraccion


@dataclass(frozen=True)
class Coleccion:
    nombre: str
    terminos: dict[str, set[str]]
    estado: str
    motivo: str = ""
    fuente: str = ""

    @property
    def n_terminos(self) -> int:
        return len(self.terminos)


def leer_gmt(path: Path, *, nombre: str, fuente: str = "") -> Coleccion:
    """Lee un GMT. El identificador es el campo 1 si parece un identificador."""
    terminos: dict[str, set[str]] = {}
    for linea in path.read_text(encoding="utf-8", errors="replace").splitlines():
        campos = [c.strip() for c in linea.rstrip("\n").split("\t")]
        if len(campos) < 3:
            continue
        etiqueta, segundo = campos[0], campos[1]
        # Reactome pone el identificador R-HSA en el segundo campo; MSigDB pone
        # una URL. Se usa el identificador cuando existe y el nombre si no.
        clave = segundo if segundo.startswith("R-HSA-") else etiqueta
        genes = {g.upper() for g in campos[2:] if g}
        if genes:
            terminos.setdefault(clave, set()).update(genes)
    return Coleccion(
        nombre=nombre,
        terminos=terminos,
        estado="usada" if terminos else "no_disponible",
        motivo="" if terminos else f"{path} no produjo ningún término",
        fuente=fuente or str(path),
    )


def cargar_colecciones(
    *,
    go_long: pd.DataFrame | None,
    base_externa: Path,
) -> dict[str, Coleccion]:
    """GO BP del repositorio y las dos colecciones descargadas, si están."""
    salida: dict[str, Coleccion] = {}

    columna_aspecto = None
    if go_long is not None and not go_long.empty:
        for candidata in ("namespace", "aspect"):
            if candidata in go_long.columns:
                columna_aspecto = candidata
                break
    if columna_aspecto is not None:
        bp = go_long.loc[go_long[columna_aspecto].astype(str).str.upper().isin({"BP", "P"})]
        terminos: dict[str, set[str]] = {}
        for termino, grupo in bp.groupby(bp["go_id"].astype(str)):
            terminos[termino] = {str(g).upper() for g in grupo["gene_symbol"]}
        salida["go_bp"] = Coleccion(
            nombre="go_bp",
            terminos=terminos,
            estado="usada" if terminos else "no_disponible",
            motivo="" if terminos else "la anotación GO no trae aspecto BP",
            fuente="anotación GO del repositorio, la misma de v1",
        )
    else:
        salida["go_bp"] = Coleccion(
            nombre="go_bp", terminos={}, estado="no_disponible",
            motivo="no hay anotación GO con columna de aspecto (namespace)",
        )

    ficha_path = base_externa / "estado.json"
    ficha = json.loads(ficha_path.read_text(encoding="utf-8")) if ficha_path.exists() else {}
    for clave, fichero in (("reactome", "reactome_pathways.gmt"), ("hallmark", "msigdb_hallmark.gmt")):
        nombre = "msigdb_hallmark" if clave == "hallmark" else clave
        path = base_externa / fichero
        info = ficha.get(clave, {})
        if info.get("estado") == "no_disponible" or not path.exists():
            salida[nombre] = Coleccion(
                nombre=nombre, terminos={}, estado="no_disponible",
                motivo=info.get("motivo", f"no existe {path}"),
            )
            continue
        salida[nombre] = leer_gmt(
            path, nombre=nombre, fuente=f"{info.get('url', path)} sha256={info.get('sha256', '')[:12]}"
        )
    return salida


def _pertenencia(
    genes: np.ndarray, terminos: dict[str, set[str]], *, min_proteinas: int
) -> dict[str, np.ndarray]:
    """Máscara booleana por término, solo para los que alcanzan el mínimo."""
    arriba = np.array([str(g).upper() for g in genes])
    salida: dict[str, np.ndarray] = {}
    for termino, miembros in terminos.items():
        sel = np.isin(arriba, list(miembros))
        if int(sel.sum()) >= min_proteinas:
            salida[termino] = sel
    return salida


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p
    orden = np.argsort(p, kind="stable")
    ajust = np.empty(n)
    previo = 1.0
    for rango in range(n - 1, -1, -1):
        i = orden[rango]
        previo = min(previo, p[i] * n / (rango + 1))
        ajust[i] = previo
    return np.clip(ajust, 0.0, 1.0)


def prueba_por_cohorte(
    estadisticos: pd.Series,
    terminos: dict[str, set[str]],
    *,
    rho: float,
    df: float,
    min_proteinas: int,
    alpha: float,
) -> pd.DataFrame:
    """Prueba competitiva de todos los términos en una cohorte."""
    genes = estadisticos.index.to_numpy()
    pertenencia = _pertenencia(genes, terminos, min_proteinas=min_proteinas)
    resultados: list[ResultadoConjunto] = camera(
        estadisticos.to_numpy(dtype=float),
        pertenencia,
        rho=rho,
        df=df,
        min_proteinas=min_proteinas,
    )
    if not resultados:
        return pd.DataFrame(
            columns=[
                "termino", "n_en_conjunto", "n_fuera", "estadistico", "direccion",
                "p_value", "padj", "significativo", "vif", "rho", "metodo",
            ]
        )
    tabla = pd.DataFrame([r.__dict__ for r in resultados])
    tabla["padj"] = _bh(tabla["p_value"].to_numpy())
    tabla["significativo"] = tabla["padj"] < alpha
    return tabla.sort_values("p_value", kind="stable").reset_index(drop=True)


def concordancia_entre_cohortes(
    descubrimiento: pd.DataFrame,
    replica: pd.DataFrame,
    *,
    coleccion: str,
    n_min_terminos: int,
    umbral: float,
) -> dict:
    """Términos significativos en descubrimiento: ¿coincide el signo en la réplica?

    El término entra en el denominador solo si la réplica pudo evaluarlo, lo que
    exige que tuviera suficientes proteínas medidas allí. La selección del
    término no mira nunca la réplica.
    """
    if descubrimiento.empty or replica.empty:
        return {
            "coleccion": coleccion,
            "estado": "no_evaluable",
            "motivo": "una de las dos cohortes no produjo ningún término evaluable",
            "n_terminos_evaluados": 0,
        }
    sig = descubrimiento.loc[descubrimiento["significativo"], ["termino", "direccion", "estadistico"]]
    rep = replica.set_index("termino")
    filas = []
    for _, fila in sig.iterrows():
        termino = str(fila["termino"])
        if termino not in rep.index:
            continue
        d_rep = int(rep.loc[termino, "direccion"])
        d_des = int(fila["direccion"])
        if d_des == 0 or d_rep == 0:
            continue
        filas.append(
            {
                "termino": termino,
                "direccion_descubrimiento": d_des,
                "direccion_replica": d_rep,
                "estadistico_descubrimiento": float(fila["estadistico"]),
                "estadistico_replica": float(rep.loc[termino, "estadistico"]),
                "padj_replica": float(rep.loc[termino, "padj"]),
                "mismo_signo": d_des == d_rep,
            }
        )
    detalle = pd.DataFrame(filas)
    n = int(len(detalle))
    n_sig_descubrimiento = int(descubrimiento["significativo"].sum())
    base = {
        "coleccion": coleccion,
        "n_terminos_probados_descubrimiento": int(len(descubrimiento)),
        "n_terminos_significativos_descubrimiento": n_sig_descubrimiento,
        "n_terminos_evaluados": n,
        "n_min_terminos": int(n_min_terminos),
        "umbral": float(umbral),
    }
    if n < n_min_terminos:
        base |= {
            "estado": "no_evaluable",
            "motivo": (
                f"{n} términos evaluables no llegan al mínimo preregistrado de "
                f"{n_min_terminos}. Con tan pocos, el intervalo no distingue "
                "concordancia perfecta de azar."
            ),
            "concordancia": None,
            "detalle": detalle,
        }
        return base
    fraccion = Fraccion(
        nombre=f"concordancia_via_{coleccion}",
        numerador=int(detalle["mismo_signo"].sum()),
        denominador=n,
        unidad_denominador="terminos",
        nota="Términos significativos en descubrimiento y evaluables en la réplica.",
    )
    d = fraccion.to_dict()
    # La cláusula se sostiene solo si el intervalo entero supera el umbral: con
    # el punto estimado basta para engañarse, como pasó en v1.
    sostiene = d["ci_low"] is not None and d["ci_low"] > umbral
    margen = None if d["ci_low"] is None else float(d["ci_low"]) - float(umbral)
    base |= {
        "estado": "sostiene" if sostiene else "rechaza",
        "concordancia": d,
        "margen_del_ic_inferior_sobre_el_umbral": margen,
        "holgado": None if margen is None else bool(margen > 0.05),
        "criterio": (
            "El límite inferior del intervalo del 95 % tiene que superar el umbral. "
            "El punto estimado por sí solo no basta."
        ),
        "detalle": detalle,
    }
    return base


def tabla_verdad_h1(clausula_proteina: str, clausula_via: str) -> dict:
    """La tabla preregistrada, aplicada sin margen de interpretación."""
    if clausula_proteina == "no_evaluable":
        veredicto = "no_evaluable"
    elif clausula_proteina == "rechaza":
        veredicto = "rechaza_h1"
    elif clausula_via == "rechaza":
        veredicto = "rechaza_h1"
    elif clausula_via == "sostiene":
        veredicto = "sostiene_h1"
    else:
        veredicto = "parcial_solo_proteina"
    return {
        "clausula_proteina": clausula_proteina,
        "clausula_via": clausula_via,
        "veredicto": veredicto,
        "nota": (
            "H1 es la conjunción de las dos cláusulas. La de proteína basta para "
            "rechazar. Una cláusula de vía no evaluable no salva H1: deja el "
            "resultado parcial."
        ),
    }


def nulo_de_etiquetas(
    estadisticos: pd.Series,
    terminos: dict[str, set[str]],
    *,
    rho: float,
    df: float,
    min_proteinas: int,
    n_permutaciones: int,
    semilla: int,
) -> dict:
    pertenencia = _pertenencia(
        estadisticos.index.to_numpy(), terminos, min_proteinas=min_proteinas
    )
    p = nulo_etiquetas_de_gen(
        estadisticos.to_numpy(dtype=float),
        pertenencia,
        rho=rho,
        df=df,
        min_proteinas=min_proteinas,
        n_permutaciones=n_permutaciones,
        semilla=semilla,
    )
    if not p:
        return {"estado": "no_calculable", "n_terminos": 0}
    valores = np.array(list(p.values()))
    padj = _bh(valores)
    return {
        "estado": "calculado",
        "n_terminos": int(valores.size),
        "n_permutaciones": int(n_permutaciones),
        "semilla": int(semilla),
        "p_minimo_alcanzable": 1.0 / (n_permutaciones + 1),
        "n_con_p_en_el_tope": int((valores <= 1.0 / (n_permutaciones + 1)).sum()),
        "n_significativos_bh": int((padj < 0.05).sum()),
        "por_termino": {k: float(v) for k, v in p.items()},
        "nota": (
            "Permuta qué proteína pertenece a qué término, conservando el tamaño de "
            "cada conjunto. Es la etapa que el nulo de inversión de signos de R5 no "
            "cubría: R5 dice si el contraste por proteína es real, no si el "
            "enriquecimiento lo es."
        ),
    }
