"""H9. Cociente de Wald con un cis-pQTL líder. No decide la tesis."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from oncogo_crc.settings import Settings
from oncogo_crc.v3.salida import volcar_json, volcar_tsv

PROTEINAS = ("SAA1", "LBP", "LRG1", "CEACAM5")


def wald(beta_x: float, se_x: float, beta_y: float, se_y: float) -> dict:
    """Cociente beta_y/beta_x y su error por el método delta. Un solo instrumento."""
    if not all(np.isfinite(v) for v in (beta_x, se_x, beta_y, se_y)) or beta_x == 0:
        return {"estado": "no_evaluable"}
    ratio = beta_y / beta_x
    varianza = se_y**2 / beta_x**2 + beta_y**2 * se_x**2 / beta_x**4
    if varianza <= 0:
        return {"estado": "no_evaluable"}
    se = float(np.sqrt(varianza))
    z = ratio / se
    return {
        "estado": "exploratorio",
        "beta": float(ratio),
        "se": se,
        "z": float(z),
        "p": float(2 * norm.sf(abs(z))),
    }


def _md5(path: Path) -> str:
    resumen = hashlib.md5()
    with path.open("rb") as entrada:
        while True:
            bloque = entrada.read(1 << 20)
            if not bloque:
                break
            resumen.update(bloque)
    return resumen.hexdigest()


def _lideres(path: Path) -> pd.DataFrame:
    tabla = pd.read_excel(path, sheet_name="ST9", header=4)
    cis = tabla.loc[tabla["cis/trans"].astype(str).str.lower().eq("cis")].copy()
    cis["proteina"] = cis["Assay Target"].astype(str).str.upper()
    cis = cis.loc[cis["proteina"].isin(PROTEINAS)]
    cis["logp"] = pd.to_numeric(cis["log10(p) (discovery)"], errors="coerce")
    cis = cis.sort_values("logp", ascending=False, kind="stable").drop_duplicates("proteina")
    return cis


def _gwas(path: Path, rsids: set[str]) -> dict[str, dict]:
    hallados: dict[str, dict] = {}
    with path.open(encoding="utf-8", errors="replace") as entrada:
        columnas = entrada.readline().rstrip("\n").split("\t")
        indice = {nombre: i for i, nombre in enumerate(columnas)}
        for linea in entrada:
            rs = linea.split("\t", 1)[0]
            if rs not in rsids:
                continue
            partes = linea.rstrip("\n").split("\t")
            hallados[rs] = {nombre: partes[i] for nombre, i in indice.items()}
            if len(hallados) == len(rsids):
                break
    return hallados


def _alinear(efecto: str, otro: str, a1: str, a0: str, beta: float, se: float) -> tuple[float, float] | None:
    efecto, otro, a1, a0 = (x.upper() for x in (efecto, otro, a1, a0))
    if efecto == a1 and otro == a0:
        return beta, se
    if efecto == a0 and otro == a1:
        return -beta, se
    return None


def ejecutar(settings: Settings, salida: Path, sello: dict, prereg) -> dict:
    bloque = prereg.get("enmienda_v3_4.h9_gwas_verificado") or {}
    gwas = settings.resolve(str(bloque.get("archivo", "")))
    sun = settings.resolve("external/sun_2023/41586_2023_6592_MOESM3_ESM.xlsx")
    esperado = str(bloque.get("md5", "")).lower()
    if not bloque or not gwas.is_file() or not sun.is_file() or not esperado:
        payload = {"estado": "requiere_acceso", "motivo": "la enmienda no registra un GWAS verificado"}
        volcar_json(payload, salida / "h9_wald.json", sello)
        return payload
    observado = _md5(gwas)
    if observado != esperado:
        payload = {"estado": "no_evaluable", "motivo": f"el md5 del GWAS es {observado} y el registrado es {esperado}"}
        volcar_json(payload, salida / "h9_wald.json", sello)
        return payload
    lideres = _lideres(sun)
    if lideres.empty:
        payload = {"estado": "no_evaluable", "motivo": "ST9 no trae cis-pQTL de descubrimiento para estas proteínas"}
        volcar_json(payload, salida / "h9_wald.json", sello)
        return payload
    rsids = set(lideres["rsID"].astype(str))
    resultado_gwas = _gwas(gwas, rsids)
    filas = []
    for proteina in PROTEINAS:
        fila = lideres.loc[lideres["proteina"].eq(proteina)]
        if fila.empty:
            filas.append({"proteina": proteina, "estado": "no_medido", "motivo": "no está en el descubrimiento de ST9"})
            continue
        lider = fila.iloc[0]
        variante = str(lider["Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)"])
        partes = variante.split(":")
        if len(partes) < 4:
            filas.append({"proteina": proteina, "estado": "no_evaluable", "motivo": f"identificador no alineable: {variante}"})
            continue
        a0, a1 = partes[2], partes[3]
        rs = str(lider["rsID"])
        crudo = resultado_gwas.get(rs)
        if crudo is None:
            filas.append({"proteina": proteina, "estado": "no_evaluable", "motivo": f"{rs} no está en el GWAS", "rsid": rs})
            continue
        alineado = _alinear(
            crudo["effect_allele"], crudo["other_allele"], a1, a0,
            float(crudo["beta"]), float(crudo["standard_error"]),
        )
        if alineado is None:
            filas.append({"proteina": proteina, "estado": "no_evaluable", "motivo": "los alelos del GWAS no coinciden con A0/A1", "rsid": rs})
            continue
        beta_y, se_y = alineado
        estimado = wald(
            float(lider["BETA (discovery, wrt. A1)"]), float(lider["SE (discovery)"]), beta_y, se_y,
        )
        filas.append({
            "proteina": proteina,
            "rsid": rs,
            "variante": variante,
            "beta_exposicion": float(lider["BETA (discovery, wrt. A1)"]),
            "se_exposicion": float(lider["SE (discovery)"]),
            "beta_desenlace_alineado": beta_y,
            "se_desenlace": se_y,
            **estimado,
        })
    tabla = pd.DataFrame(filas)
    volcar_tsv(tabla, salida / "h9_wald.tsv", sello)
    medidos = [fila for fila in filas if fila.get("estado") == "exploratorio"]
    payload = {
        "estado": "exploratorio" if medidos else "no_evaluable",
        "motivo": "Wald con un instrumento por proteína. No decide la tesis. IVW, Egger y la mediana piden al menos tres.",
        "gwas": "GCST90255676",
        "md5": observado,
        "n_exploratorios": len(medidos),
        "proteinas": filas,
    }
    volcar_json(payload, salida / "h9_wald.json", sello)
    return payload
