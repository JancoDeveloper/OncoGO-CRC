"""Lee un report de DIA-NN sin modificarlo.

El análisis principal de v3.2 cuantifica PG.MaxLFQ por Protein.Group. Una
casilla entra si el precursor pasa Q.Value y el grupo pasa PG.Q.Value en esa
corrida, y si el Global.PG.Q.Value del grupo pasa el corte. Un grupo entra en
la matriz principal si, además, algún miembro lleva el sufijo _HUMAN.

No imputa. Una corrida sin precursor que pase el corte queda vacía.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

SHA256_REPORT = "A8A5AADD74391AFEF68CF4D6C3A7D3860023101E5C92590E2E87C152A39E59E8"
BYTES_REPORT = 213_313_241
_ESPECIE = re.compile(r"_([A-Z0-9]{2,6})$")

COLUMNAS = [
    "Run",
    "Protein.Group",
    "Protein.Ids",
    "Protein.Names",
    "Genes",
    "PG.MaxLFQ",
    "Genes.MaxLFQ",
    "Q.Value",
    "PG.Q.Value",
    "Global.PG.Q.Value",
    "GG.Q.Value",
    "Proteotypic",
    "Precursor.Quantity",
    "Stripped.Sequence",
]


def especies_de_nombres(texto: object) -> list[str]:
    vistos: list[str] = []
    for parte in str(texto).split(";"):
        token = parte.strip().split()[0] if parte.strip() else ""
        hallado = _ESPECIE.search(token)
        if hallado and hallado.group(1) not in vistos:
            vistos.append(hallado.group(1))
    return vistos


def leer_report(path: Path) -> pd.DataFrame:
    cabecera = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    uso = [columna for columna in COLUMNAS if columna in cabecera]
    faltan = {"Run", "Protein.Group", "Q.Value", "PG.Q.Value", "Global.PG.Q.Value", "PG.MaxLFQ"} - set(uso)
    if faltan:
        raise ValueError(f"DIA-NN sin columnas {sorted(faltan)}")
    tabla = pd.read_csv(path, sep="\t", usecols=uso, low_memory=False)
    for columna in ("Q.Value", "PG.Q.Value", "Global.PG.Q.Value", "GG.Q.Value", "PG.MaxLFQ", "Genes.MaxLFQ", "Proteotypic"):
        if columna in tabla.columns:
            tabla[columna] = pd.to_numeric(tabla[columna], errors="coerce")
    return tabla


def controles_del_report(tabla: pd.DataFrame) -> dict:
    """Cifras de la sección 3.1. Sirven para comprobar que el fichero es el mismo."""
    grupos = tabla.groupby("Protein.Group", sort=False)
    global_q = grupos["Global.PG.Q.Value"].min()
    pg_min = grupos["PG.Q.Value"].min()
    nombres = grupos["Protein.Names"].first() if "Protein.Names" in tabla.columns else pd.Series("", index=global_q.index)
    humanos = nombres.map(lambda texto: "HUMAN" in especies_de_nombres(texto))
    if "Genes" in tabla.columns:
        genes = grupos["Genes"].first().fillna("").astype(str).str.strip()
        tiene_gen = ~genes.isin(["", "nan"])
    else:
        tiene_gen = pd.Series(False, index=global_q.index)
    pasa_global = global_q <= 0.01
    pasa_pg = pg_min <= 0.01
    ambos = pasa_global & pasa_pg & humanos
    return {
        "filas": int(len(tabla)),
        "corridas": int(tabla["Run"].nunique()),
        "q_value_max": float(tabla["Q.Value"].max()),
        "pg_q_value_max": float(tabla["PG.Q.Value"].max()),
        "global_pg_q_value_max": float(tabla["Global.PG.Q.Value"].max()),
        "grupos": int(global_q.size),
        "grupos_global_pg_q_0_01": int(pasa_global.sum()),
        "grupos_sin_pg_q_0_01_en_ninguna_corrida": int((~pasa_pg).sum()),
        "humanos_pg_y_global": int(ambos.sum()),
        "humanos_pg_y_global_con_gen": int((ambos & tiene_gen).sum()),
    }


_ISOFORMA = re.compile(r"-\d+$")
_OX = re.compile(r"\bOX=(\d+)\b")
_GN = re.compile(r"\bGN=(\S+)")
HUMANO_OX = "9606"
QUERATINAS_EPITELIALES = frozenset({"KRT7", "KRT8", "KRT18", "KRT19", "KRT20"})
PLASMA_CLASICO = frozenset({"ALB", "APOA1", "C3", "TF", "SERPINA1", "A2M", "APOE", "FGA", "FGB", "FGG", "LRG1"})


def acceso_base(identificador: object) -> str:
    """Cont_P02768, P02768-2 y sp|Cont_P02768|ALBU_BOVIN dan P02768."""
    texto = str(identificador).strip()
    if "|" in texto:
        partes = texto.split("|")
        texto = partes[1] if len(partes) >= 2 else partes[0]
    if texto.startswith("Cont_"):
        texto = texto[len("Cont_"):]
    return _ISOFORMA.sub("", texto)


@dataclass(frozen=True)
class LibreriaContaminantes:
    accesos: frozenset[str]
    accesos_humanos: frozenset[str]
    genes_humanos: frozenset[str]
    n_entradas: int
    por_taxon: dict[str, int] = field(default_factory=dict)

    @property
    def accesos_no_humanos(self) -> frozenset[str]:
        return self.accesos - self.accesos_humanos


LIBRERIA_VACIA = LibreriaContaminantes(frozenset(), frozenset(), frozenset(), 0)


def leer_libreria(path: Path) -> LibreriaContaminantes:
    """Accesos de la librería, separados por especie. El GN solo se guarda para OX=9606."""
    accesos: set[str] = set()
    humanos: set[str] = set()
    genes: set[str] = set()
    por_taxon: dict[str, int] = {}
    n = 0
    for linea in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not linea.startswith(">"):
            continue
        n += 1
        acceso = acceso_base(linea[1:].split()[0])
        hallado = _OX.search(linea)
        taxon = hallado.group(1) if hallado else "sin_OX"
        por_taxon[taxon] = por_taxon.get(taxon, 0) + 1
        accesos.add(acceso)
        if taxon == HUMANO_OX:
            humanos.add(acceso)
            gen = _GN.search(linea)
            if gen:
                genes.add(gen.group(1).upper())
    return LibreriaContaminantes(frozenset(accesos), frozenset(humanos), frozenset(genes), n, por_taxon)


def trazabilidad(
    tabla: pd.DataFrame,
    *,
    libreria: LibreriaContaminantes | None = None,
    genes_piel: set[str] | None = None,
    paneles: dict[str, set[str]] | None = None,
) -> pd.DataFrame:
    libreria = libreria or LIBRERIA_VACIA
    genes_piel = {g.upper() for g in (genes_piel or set())}
    paneles = {nombre: {g.upper() for g in miembros} for nombre, miembros in (paneles or {}).items()}
    filas = []
    for grupo, bloque in tabla.groupby("Protein.Group", sort=False):
        nombres = str(bloque["Protein.Names"].iloc[0]) if "Protein.Names" in bloque.columns else ""
        genes = str(bloque["Genes"].iloc[0]) if "Genes" in bloque.columns else ""
        especies = especies_de_nombres(nombres)
        humano = "HUMAN" in especies
        simbolos = {parte.strip().upper() for parte in genes.split(";") if parte.strip() and parte.strip() != "nan"}
        identificadores: set[str] = set()
        if "Protein.Ids" in bloque.columns:
            for texto in bloque["Protein.Ids"].astype(str):
                identificadores.update(acceso_base(p) for p in texto.split(";") if p.strip() and p.strip() != "nan")
        if humano:
            if identificadores & libreria.accesos_humanos:
                motivo = "acceso_humano_en_libreria"
            elif not identificadores and simbolos & libreria.genes_humanos:
                motivo = "simbolo_humano_sin_acceso"
            else:
                motivo = ""
        else:
            motivo = "acceso_no_humano_en_libreria" if identificadores & libreria.accesos else ""
        comparte = humano and bool(identificadores & libreria.accesos_no_humanos)
        panel = ";".join(nombre for nombre, miembros in paneles.items() if simbolos & miembros)
        proteotipicos = float("nan")
        if "Proteotypic" in bloque.columns and bloque["Proteotypic"].notna().any():
            proteotipicos = float(bloque["Proteotypic"].mean())
        filas.append({
            "protein_group": grupo,
            "genes": genes,
            "especies": ";".join(especies),
            "humano": humano,
            "pg_q_min": float(bloque["PG.Q.Value"].min()),
            "global_pg_q": float(bloque["Global.PG.Q.Value"].min()),
            "n_precursores": int(len(bloque)),
            "fraccion_proteotipicos": proteotipicos,
            "contaminante": bool(motivo),
            "contaminante_motivo": motivo,
            "comparte_con_contaminante": comparte,
            "queratina_epitelial": bool(simbolos & QUERATINAS_EPITELIALES),
            "contaminante_piel": bool(simbolos & genes_piel),
            "panel": panel,
        })
    return pd.DataFrame(filas)


def matriz_filtrada(
    tabla: pd.DataFrame,
    *,
    q_value_max: float = 0.01,
    pg_q_value_max: float = 0.01,
    global_pg_q_value_max: float = 0.01,
    exigir_humano: bool = True,
    valor: str = "PG.MaxLFQ",
    indice: str = "Protein.Group",
    q_de_grupo: str = "PG.Q.Value",
) -> pd.DataFrame:
    """Matriz ancha. El índice es el grupo o el gen. No escribe nada."""
    if valor not in tabla.columns or indice not in tabla.columns:
        raise ValueError(f"faltan {valor} o {indice}")
    trabajo = tabla.dropna(subset=[indice, "Run"]).copy()
    global_q = trabajo.groupby(indice)["Global.PG.Q.Value"].transform("min")
    trabajo = trabajo.loc[global_q <= global_pg_q_value_max]
    if exigir_humano and "Protein.Names" in trabajo.columns:
        humanos = {
            grupo
            for grupo, nombres in trabajo.groupby(indice)["Protein.Names"].first().items()
            if "HUMAN" in especies_de_nombres(nombres)
        }
        trabajo = trabajo.loc[trabajo[indice].isin(humanos)]
    casilla = trabajo.loc[
        (trabajo["Q.Value"] <= q_value_max) & (trabajo[q_de_grupo] <= pg_q_value_max)
    ]
    casilla = casilla.dropna(subset=[valor])
    if casilla.empty:
        return pd.DataFrame(columns=["protein_id"])
    ancha = casilla.pivot_table(index=indice, columns="Run", values=valor, aggfunc="median")
    corridas = list(dict.fromkeys(trabajo["Run"].astype(str)))
    ancha = ancha.reindex(columns=corridas)
    return ancha.reset_index().rename(columns={indice: "protein_id"})


def postprocess_dia_nn_plasma_download(matrix_path: Path, metadata_path: Path) -> None:
    """Escribe la metadata de corridas. No toca el report ni la matriz ancha."""
    cabecera = pd.read_csv(matrix_path, sep="\t", nrows=0).columns
    if "Run" not in cabecera:
        return
    corridas = pd.read_csv(matrix_path, sep="\t", usecols=["Run"])["Run"].astype(str).unique().tolist()
    metadata = build_plasma_metadata_for_pxd055218_runs(corridas)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(metadata_path, sep="\t", index=False)


def build_plasma_metadata_for_pxd055218_runs(run_ids: list[str]) -> pd.DataFrame:
    """PXD055218 trae tres grupos en el nombre de la corrida.

    CP_* es cáncer colorrectal, CN_* es colonoscopia negativa y Normal_* es
    control sano. Los tres se conservan.
    """
    rows: list[dict[str, str]] = []
    for identificador in sorted(set(run_ids)):
        if not identificador or identificador == "Run":
            continue
        nombre = str(identificador).replace("\\", "/").split("/")[-1]
        if nombre.startswith("CP_"):
            condicion = "crc"
        elif nombre.startswith("CN_"):
            condicion = "colonoscopy_negative"
        elif nombre.startswith("Normal_"):
            condicion = "healthy"
        else:
            continue
        rows.append({"sample_id": identificador, "condition": condicion, "batch": "PXD055218"})
    return pd.DataFrame(rows)
