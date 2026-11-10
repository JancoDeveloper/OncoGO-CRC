"""Triangulación exploratoria. Rellena columnas con ficheros que ya están. No hace pruebas nuevas."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.dia_nn_plasma import acceso_base
from oncogo_crc.v3.salida import corrida_limpia, volcar_json, volcar_tsv

NO_MEDIDA = "no_medida"


def tabla_vacia(genes: list[str], columnas: list[str]) -> pd.DataFrame:
    tabla = pd.DataFrame({"gen": genes})
    for columna in columnas:
        tabla[columna] = NO_MEDIDA
    return tabla


def _poner(tabla: pd.DataFrame, gen: str, columna: str, valor) -> None:
    if valor is None or (isinstance(valor, float) and pd.isna(valor)):
        return
    tabla.loc[tabla["gen"].eq(gen), columna] = str(valor)


def _clase_por_gen(ruta: Path, columna_valor: str) -> dict[str, object]:
    if not ruta.is_file():
        return {}
    datos = pd.read_csv(ruta, sep="\t", dtype=str)
    if "gen" not in datos.columns and "protein_id" in datos.columns:
        datos = datos.rename(columns={"protein_id": "gen"})
    if "gen" not in datos.columns or columna_valor not in datos.columns:
        return {}
    datos["gen"] = datos["gen"].str.upper()
    return dict(zip(datos["gen"], datos[columna_valor], strict=False))


def _de_esta_corrida(salida: Path, nombre: str, anterior: Path) -> Path:
    local = salida / nombre
    if local.is_file() or corrida_limpia():
        return local
    return anterior


def ejecutar(settings: Settings, salida: Path, sello: dict, prereg) -> dict:
    genes = [str(g) for g in prereg.get("enmienda_v3_4.triangulacion.genes")]
    columnas = [str(c) for c in prereg.get("enmienda_v3_4.triangulacion.columnas")]
    tabla = tabla_vacia(genes, columnas)
    buscados = set(genes)

    tanaka = _de_esta_corrida(salida, "k1_primario_frente_a_mucosa.tsv", settings.resolve("results/v3_3/k1_primario_frente_a_mucosa.tsv"))
    uni = settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv")
    if tanaka.is_file() and uni.is_file():
        mapa = pd.read_csv(uni, sep="\t", usecols=["gene_symbol", "uniprot_accession"], dtype=str).dropna()
        mapa["gene_symbol"] = mapa["gene_symbol"].str.upper()
        acceso_a_gen = {}
        for gen, acceso in zip(mapa["gene_symbol"], mapa["uniprot_accession"], strict=True):
            acceso_a_gen.setdefault(acceso_base(acceso), set()).add(gen)
        unico = {a: next(iter(g)) for a, g in acceso_a_gen.items() if len(g) == 1}
        bloque = pd.read_csv(tanaka, sep="\t")
        bloque["gen"] = bloque["protein_id"].map(lambda pid: unico.get(acceso_base(str(pid).split(";")[0])))
        bloque = bloque.dropna(subset=["gen"]).drop_duplicates("gen")
        for _, fila in bloque.loc[bloque["gen"].isin(buscados)].iterrows():
            _poner(tabla, fila["gen"], "tanaka_log2fc", fila.get("log2fc"))
            _poner(tabla, fila["gen"], "tanaka_padj", fila.get("padj"))
            _poner(tabla, fila["gen"], "tanaka_n_obs", fila.get("n_obs"))

    cptac = settings.resolve("results/certified/cptac_unshared_treat.tsv")
    if cptac.is_file():
        bloque = pd.read_csv(cptac, sep="\t")
        bloque["gen"] = bloque["protein_id"].astype(str).str.upper()
        for _, fila in bloque.loc[bloque["gen"].isin(buscados)].iterrows():
            _poner(tabla, fila["gen"], "cptac_coad_log2fc", fila.get("log2fc"))
            _poner(tabla, fila["gen"], "cptac_coad_treat", fila.get("significant"))

    for columna, ruta, campo in (
        ("clase_h5", salida / "h5_especificidad_organo.tsv", "clase"),
        ("hpa_logfc_sanos", _de_esta_corrida(salida, "h10_especificidad.tsv", settings.resolve("results/v3_4/h10_especificidad.tsv")), "logfc_sanos"),
        ("hpa_padj_sanos", _de_esta_corrida(salida, "h10_especificidad.tsv", settings.resolve("results/v3_4/h10_especificidad.tsv")), "padj_sanos"),
        ("clase_h10", _de_esta_corrida(salida, "h10_especificidad.tsv", settings.resolve("results/v3_4/h10_especificidad.tsv")), "clase"),
        ("pxd055218_k4_crc_sanos", _de_esta_corrida(salida, "k4_genes_maxlfq_crc_frente_a_sanos.tsv", settings.resolve("results/v3_3/k4_genes_maxlfq_crc_frente_a_sanos.tsv")), "clase"),
        ("pxd055218_k4_crc_cn", _de_esta_corrida(salida, "k4_genes_maxlfq_crc_frente_a_colonoscopia_negativa.tsv", settings.resolve("results/v3_3/k4_genes_maxlfq_crc_frente_a_colonoscopia_negativa.tsv")), "clase"),
    ):
        if ruta == salida / "h5_especificidad_organo.tsv" and not ruta.is_file() and not corrida_limpia():
            ruta = settings.resolve("results/v3_4/h5_especificidad_organo.tsv")
        for gen, valor in _clase_por_gen(ruta, campo).items():
            if gen in buscados:
                _poner(tabla, gen, columna, valor)

    wannberg = settings.resolve("external/wannberg_2026/wannberg_2026_tabla_s5.tsv")
    if wannberg.is_file() and "wannberg" in tabla.columns:
        bloque = pd.read_csv(wannberg, sep="\t", dtype=str)
        bloque["proteina"] = bloque["proteina"].str.upper()
        for _, fila in bloque.loc[bloque["proteina"].isin(buscados)].drop_duplicates("proteina").iterrows():
            _poner(tabla, fila["proteina"], "wannberg", fila.get("odds_ratio"))

    h8 = salida / "h8_carga_tumoral.tsv"
    if h8.is_file():
        bloque = pd.read_csv(h8, sep="\t", dtype=str)
        presentes = set(bloque["gen"].str.upper())
        for gen in genes:
            _poner(tabla, gen, "pxd064655_presente", "si" if gen in presentes else "no")
        for _, fila in bloque.loc[bloque["gen"].str.upper().isin(buscados)].iterrows():
            _poner(tabla, str(fila["gen"]).upper(), "pxd064655_iv_frente_i_ii", fila.get("clase"))

    papier = settings.resolve("external/papier_2024/41467_2024_48017_MOESM4_ESM.xlsx")
    if papier.is_file() and "papier_hr" in tabla.columns:
        bloque = pd.read_excel(papier, sheet_name="Supp Data 3", header=3)
        bloque = bloque.loc[bloque["Cancer site"].eq("Colorectal Cancer")].copy()
        bloque["gen"] = bloque["Protein"].astype(str).str.upper()
        for _, fila in bloque.loc[bloque["gen"].isin(buscados)].drop_duplicates("gen").iterrows():
            _poner(tabla, fila["gen"], "papier_hr", fila.get("Hazard Ratio"))

    for columna in columnas:
        tabla[columna] = tabla[columna].fillna(NO_MEDIDA).replace("", NO_MEDIDA)
    vacias = int(tabla[columnas].isna().sum().sum()) + int((tabla[columnas] == "").sum().sum())
    volcar_tsv(tabla, salida / "triangulacion_v3_5.tsv", sello)
    payload = {
        "tipo": "exploratoria_post_hoc",
        "sin_pruebas_nuevas": True,
        "n_filas": int(len(tabla)),
        "n_columnas": int(len(columnas)),
        "celdas_vacias": vacias,
        "wannberg": "external/wannberg_2026/wannberg_2026_tabla_s5.tsv" if wannberg.is_file() else "no_disponible",
        "nota": "No es una validación. Cada celda sale de un fichero ya calculado o dice no_medida.",
    }
    volcar_json(payload, salida / "triangulacion_v3_5.json", sello)
    return payload
