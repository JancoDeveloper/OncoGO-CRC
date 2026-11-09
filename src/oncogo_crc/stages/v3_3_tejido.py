"""T1 y T2. Filtro MaxQuant y mínimo de 8 pares. No reetiqueta muestras."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.v3.genes_tanaka import genes_de, gen_unico, mapa_genes
from oncogo_crc.utils.sello import exigir_arbol_limpio, hacer_sello
from oncogo_crc.v2.metricas import Fraccion
from oncogo_crc.v3.concordancia import alinear_por_clave, concordancia_de_lin, concordancia_de_signo, spearman
from oncogo_crc.v3.etiquetas import centrar_por_mediana, contraste_no_pareado, contraste_pareado, filas_desde_columnas, matriz_log2_lfq
from oncogo_crc.v3.k2k3 import pendiente_deming
from oncogo_crc.v3.prereg import PreregistroV3
from oncogo_crc.v3.k2k3 import escribir_k2, escribir_k3
from oncogo_crc.v3.salida import carpeta, exigir_vigente, manifiesto_congelado, volcar_json, volcar_tsv

HIGADO = ("ALB", "APOH", "CPS1", "ARG1", "CYP2E1", "FGB", "HP")
MUCOSA = ("CA1", "CA2", "AQP8", "SLC26A3", "CLCA4", "GUCA2A", "GUCA2B")
SEXO = ("RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D")


def ids_descartados(protein_ids: pd.Series, protein_groups: pd.DataFrame | None) -> set[str]:
    """CON__, REV__ y las tres marcas de MaxQuant. No se elige proteína a proteína."""
    texto = protein_ids.astype(str)
    malos = set(texto[texto.str.startswith("CON__") | texto.str.startswith("REV__")])
    if protein_groups is None or protein_groups.empty:
        return malos
    grupos = protein_groups.copy()
    grupos["pid"] = grupos["Protein IDs"].astype(str).str.split(";").str[0]
    for columna in ("Reverse", "Potential contaminant", "Only identified by site"):
        if columna not in grupos.columns:
            continue
        marca = grupos[columna].astype(str).str.strip().eq("+")
        malos.update(grupos.loc[marca, "pid"].astype(str))
    return malos


def _mapa_gen(settings: Settings) -> dict[str, frozenset[str]]:
    return mapa_genes(
        settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt"),
        settings.resolve("metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv"),
    )


def _filas_gen(log2mat: pd.DataFrame, mapa: dict[str, frozenset[str]], genes: tuple[str, ...]) -> list[str]:
    buscados = {g.upper() for g in genes}
    return [pid for pid in log2mat.index.astype(str) if genes_de(pid, mapa) & buscados]


def _por_gen(tabla: pd.DataFrame, mapa: dict[str, frozenset[str]]) -> pd.DataFrame:
    """Una fila por gen. Si hay varias, se queda la de más pares observados."""
    replica = tabla.copy()
    replica["gen"] = [gen_unico(pid, mapa) for pid in replica["protein_id"].astype(str)]
    replica = replica.dropna(subset=["gen"])
    columna_n = next((c for c in ("n_obs", "n_pares_observados", "n_pares") if c in replica.columns), None)
    if columna_n is not None:
        replica = replica.sort_values(columna_n, ascending=False, kind="stable")
    return replica.drop_duplicates("gen", keep="first")


def _mediana_diferencia(log2mat: pd.DataFrame, filas: list[str], col_a: str, col_b: str) -> float:
    if not filas or col_a not in log2mat.columns or col_b not in log2mat.columns:
        return float("nan")
    dif = log2mat.loc[filas, col_a].to_numpy(dtype=float) - log2mat.loc[filas, col_b].to_numpy(dtype=float)
    dif = dif[np.isfinite(dif)]
    if dif.size < 2:
        return float("nan")
    return float(np.median(dif))


def _qc(log2mat: pd.DataFrame, meta: pd.DataFrame, uni_a_gen: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    higado = _filas_gen(log2mat, uni_a_gen, HIGADO)
    mucosa = _filas_gen(log2mat, uni_a_gen, MUCOSA)
    sexo = _filas_gen(log2mat, uni_a_gen, SEXO)
    filas_h, filas_m, filas_s = [], [], []
    fallos: dict[str, int] = {}
    for participante, bloque in meta.groupby("participant_id"):
        cols = {fila.tissue: fila.sample_id for fila in bloque.itertuples()}
        p, m, n = cols.get("primary"), cols.get("liver_metastasis"), cols.get("mucosa")
        d_mp = _mediana_diferencia(log2mat, higado, m, p) if m and p else float("nan")
        d_mn = _mediana_diferencia(log2mat, higado, m, n) if m and n else float("nan")
        falla_h = bool(np.isfinite(d_mp) and np.isfinite(d_mn) and (d_mp <= 0 or d_mn <= 0))
        if not higado:
            falla_h = False
            estado_h = "no_detectado"
        elif not np.isfinite(d_mp):
            estado_h = "no_evaluable"
            falla_h = False
        else:
            estado_h = "falla" if falla_h else "pasa"
        d_np = _mediana_diferencia(log2mat, mucosa, n, p) if n and p else float("nan")
        if not mucosa:
            falla_mu, estado_mu = False, "no_detectado"
        elif not np.isfinite(d_np):
            falla_mu, estado_mu = False, "no_evaluable"
        else:
            falla_mu = d_np <= 0
            estado_mu = "falla" if falla_mu else "pasa"
        presentes = []
        for col in (p, m, n):
            if col is None or not sexo:
                presentes.append(float("nan"))
            else:
                presentes.append(float(np.isfinite(log2mat.loc[sexo, col].to_numpy(dtype=float)).sum()))
        if not sexo:
            falla_s, estado_s = False, "no_detectado"
        elif not any(np.isfinite(presentes)):
            falla_s, estado_s = False, "no_evaluable"
        else:
            finitos = [v for v in presentes if np.isfinite(v)]
            falla_s = max(finitos) >= 2 and min(finitos) == 0
            estado_s = "falla" if falla_s else "pasa"
        fallos[str(participante)] = int(falla_h) + int(falla_mu) + int(falla_s)
        filas_h.append({"participant_id": participante, "mediana_M_menos_P": d_mp, "mediana_M_menos_N": d_mn, "n_genes": len(higado), "estado": estado_h})
        filas_m.append({"participant_id": participante, "mediana_N_menos_P": d_np, "n_genes": len(mucosa), "estado": estado_mu})
        filas_s.append({"participant_id": participante, "genes_Y_en_P": presentes[0], "genes_Y_en_M": presentes[1], "genes_Y_en_N": presentes[2], "n_genes_en_matriz": len(sexo), "estado": estado_s})
    excluir = [pid for pid, n in fallos.items() if n >= 2]
    return pd.DataFrame(filas_h), pd.DataFrame(filas_m), pd.DataFrame(filas_s), excluir


def _resumen(certificadas, tabla, mapa, n_cert, nombre):
    replica = _por_gen(tabla, mapa)
    cruzado = alinear_por_clave(
        certificadas,
        replica,
        clave_descubrimiento=certificadas["protein_id"].astype(str).str.upper(),
        clave_replica=replica["gen"],
    )
    fraccion = concordancia_de_signo(
        cruzado["log2fc_cptac"].to_numpy(),
        cruzado["log2fc_replica"].to_numpy(),
        nombre=f"R1_{nombre}",
        universo=n_cert,
        universo_nombre="proteinas_certificadas_cptac",
    )
    bajo, alto = fraccion.intervalo()
    return {
        "numerador": fraccion.numerador,
        "denominador": fraccion.denominador,
        "estimacion": fraccion.valor,
        "ci_low": bajo,
        "ci_high": alto,
        "spearman": spearman(cruzado["log2fc_cptac"].to_numpy(), cruzado["log2fc_replica"].to_numpy()),
        "lin": concordancia_de_lin(cruzado["log2fc_cptac"].to_numpy(), cruzado["log2fc_replica"].to_numpy()),
    }


def run(settings: Settings) -> None:
    prereg = PreregistroV3(settings.root)
    salida = carpeta(settings)
    exigir_vigente(salida / "marca.txt")
    salida.mkdir(parents=True, exist_ok=True)
    manifest = manifiesto_congelado(settings)
    sello = hacer_sello(settings.root, etapa="v3_3_tejido", manifest=manifest)
    exigir_arbol_limpio(sello)
    min_pares = int(prereg.get("enmienda_v3_2.tejido.min_pares_observados"))
    alpha = float(prereg.get("inferencia.alpha"))
    lfc = float(prereg.get("inferencia.treat_lfc"))
    ds = settings.datasets["datasets"]["pride_crc_meta"]
    matriz_path = settings.resolve(ds["matrix_path"])
    if not matriz_path.is_file():
        volcar_json({"estado": "no_disponible", "motivo": str(matriz_path)}, salida / "k1_etiquetas.json", sello)
        return
    cruda = pd.read_csv(matriz_path, sep="\t", low_memory=False)
    grupos_path = settings.resolve("data/raw/pride/Cohort1_proteinGroups.txt")
    flags = None
    if grupos_path.is_file():
        flags = pd.read_csv(
            grupos_path,
            sep="\t",
            usecols=lambda c: c in {"Protein IDs", "Reverse", "Potential contaminant", "Only identified by site"},
            low_memory=False,
        )
    malos = ids_descartados(cruda["protein_id"], flags)
    cruda = cruda.loc[~cruda["protein_id"].astype(str).isin(malos)].copy()
    meta = filas_desde_columnas([c for c in cruda.columns if c != "protein_id"])
    log2mat = matriz_log2_lfq(cruda)
    principal, n_pares = contraste_pareado(log2mat, meta, grupo_a="primary", grupo_b="mucosa", nombre="primario_frente_a_mucosa", lfc=lfc, alpha=alpha, min_pares=min_pares)
    metastasis, n_met = contraste_pareado(log2mat, meta, grupo_a="liver_metastasis", grupo_b="mucosa", nombre="metastasis_hepatica_frente_a_mucosa", lfc=lfc, alpha=alpha, min_pares=min_pares)
    frente_primario, n_mp = contraste_pareado(log2mat, meta, grupo_a="liver_metastasis", grupo_b="primary", nombre="metastasis_hepatica_frente_a_primario", lfc=lfc, alpha=alpha, min_pares=min_pares)
    no_pareado = contraste_no_pareado(log2mat, meta, grupo_a="primary", grupo_b="mucosa", nombre="primario_frente_a_mucosa_no_pareado", lfc=lfc, alpha=alpha, min_n=min_pares)
    for nombre, tabla in (
        ("primario_frente_a_mucosa", principal),
        ("metastasis_hepatica_frente_a_mucosa", metastasis),
        ("metastasis_hepatica_frente_a_primario", frente_primario),
        ("primario_frente_a_mucosa_no_pareado", no_pareado),
    ):
        volcar_tsv(tabla, salida / f"k1_{nombre}.tsv", sello)
    cptac = pd.read_csv(settings.resolve("results/certified/cptac_unshared_treat.tsv"), sep="\t")
    certificadas = cptac.loc[cptac["significant"].astype(bool)].copy()
    n_cert = int(certificadas["protein_id"].nunique())
    mapa = _mapa_gen(settings)
    resumenes = {
        nombre: _resumen(certificadas, tabla, mapa, n_cert, nombre)
        for nombre, tabla in (
            ("primario_frente_a_mucosa", principal),
            ("metastasis_hepatica_frente_a_mucosa", metastasis),
            ("primario_frente_a_mucosa_no_pareado", no_pareado),
        )
    }
    qc_h, qc_m, qc_s, excluir = _qc(log2mat, meta, mapa)
    volcar_tsv(qc_h, salida / "qc_etiquetas_higado.tsv", sello)
    volcar_tsv(qc_m, salida / "qc_etiquetas_mucosa.tsv", sello)
    volcar_tsv(qc_s, salida / "qc_sexo_por_participante.tsv", sello)
    if excluir:
        meta_s = meta.loc[~meta["participant_id"].astype(str).isin(excluir)]
        sens, n_sens = contraste_pareado(log2mat, meta_s, grupo_a="primary", grupo_b="mucosa", nombre="sensibilidad_etiquetas", lfc=lfc, alpha=alpha, min_pares=min_pares)
        volcar_tsv(sens, salida / "k1_sensibilidad_etiquetas.tsv", sello)
    else:
        n_sens = n_pares
    r1 = resumenes["primario_frente_a_mucosa"]
    centrada, desplazamiento, regla_ref = centrar_por_mediana(log2mat)
    centrado, _ = contraste_pareado(
        centrada, meta, grupo_a="primary", grupo_b="mucosa",
        nombre="primario_frente_a_mucosa_centrado", lfc=lfc, alpha=alpha, min_pares=min_pares,
    )
    volcar_tsv(centrado, salida / "k1_primario_frente_a_mucosa_centrado.tsv", sello)
    r1_centrado = _resumen(certificadas, centrado, mapa, n_cert, "primario_frente_a_mucosa_centrado")
    presencia = log2mat.notna().mean(axis=1)
    n_referencia = int((presencia >= (1.0 if regla_ref == "100 %" else 0.9)).sum())

    def _lectura(tabla, resumen):
        fc = pd.to_numeric(tabla["log2fc"], errors="coerce")
        fc = fc[np.isfinite(fc)]
        replica = _por_gen(tabla, mapa)
        cruzado = alinear_por_clave(
            certificadas, replica,
            clave_descubrimiento=certificadas["protein_id"].astype(str).str.upper(),
            clave_replica=replica["gen"],
        )
        ajuste = pendiente_deming(cruzado["log2fc_cptac"].to_numpy(), cruzado["log2fc_replica"].to_numpy())
        return {
            "r1_numerador": resumen["numerador"],
            "r1_denominador": resumen["denominador"],
            "r1_estimacion": resumen["estimacion"],
            "r1_ci_low": resumen["ci_low"],
            "r1_ci_high": resumen["ci_high"],
            "lin": resumen["lin"],
            "deming_pendiente": ajuste["pendiente"],
            "deming_intercepto": ajuste["intercepto"],
            "mediana_log2fc": float(fc.median()) if len(fc) else None,
            "porcentaje_cambios_positivos": float((fc > 0).mean() * 100) if len(fc) else None,
        }

    umbral = float(prereg.get("hipotesis.h1.umbral_concordancia_proteina"))
    estimacion_centrada = r1_centrado.get("estimacion")
    fragil = (
        r1.get("estimacion") is not None
        and estimacion_centrada is not None
        and (float(r1["estimacion"]) >= umbral) != (float(estimacion_centrada) >= umbral)
    )
    volcar_json({
        "regla_proteinas_de_referencia": regla_ref,
        "n_proteinas_de_referencia": n_referencia,
        "desplazamiento_por_muestra": {str(k): float(v) for k, v in desplazamiento.items()},
        "antes": _lectura(principal, r1),
        "despues": _lectura(centrado, r1_centrado),
        "h1_fragil_a_normalizacion": bool(fragil),
        "umbral_h1": umbral,
    }, salida / "k1_sensibilidad_normalizacion.json", sello)
    if r1["ci_high"] is None:
        clausula = "no_evaluable"
    elif float(r1["ci_high"]) < umbral:
        clausula = "sostiene"
    else:
        clausula = "rechaza"
    payload = {
        "estado": "calculado",
        "min_pares_observados": min_pares,
        "n_filas_excluidas_contaminante_o_decoy": int(len(malos)),
        "n_pares_primario_mucosa": n_pares,
        "n_pares_metastasis_mucosa": n_met,
        "n_pares_metastasis_primario": n_mp,
        "n_pares_sensibilidad_etiquetas": n_sens,
        "participantes_excluidos_solo_en_sensibilidad": excluir,
        "n_certificadas": n_cert,
        "contrastes": resumenes,
        "h1_clausula_proteina": clausula,
        "h1_veredicto": "rechaza_h1" if clausula == "rechaza" else ("parcial_solo_proteina" if clausula == "sostiene" else "no_evaluable"),
        "h1_fragil_a_normalizacion": bool(fragil),
        "cobertura": Fraccion("R1_cobertura", int(r1["denominador"] or 0), n_cert, "proteinas_unicas").to_dict(),
    }
    volcar_json(payload, salida / "k1_etiquetas.json", sello)
    anterior = settings.resolve("results/v3/k1_etiquetas.json")
    valor_v3 = "results/v3/k1_etiquetas.json"
    if anterior.is_file():
        import json
        viejo = json.loads(anterior.read_text(encoding="utf-8"))
        bloque = viejo.get("contrastes", {}).get("primario_frente_a_mucosa", {})
        valor_v3 = f"{bloque.get('numerador')}/{bloque.get('denominador')}"
    volcar_tsv(pd.DataFrame([{
        "metrica": "R1_concordancia_de_signo_primario_frente_a_mucosa",
        "valor_v3": valor_v3,
        "valor_v3_2": f"{r1['numerador']}/{r1['denominador']}",
        "numerador": r1["numerador"],
        "denominador": r1["denominador"],
        "motivo": "Se quitaron CON__, REV__ y las marcas de MaxQuant. El mínimo de pares vuelve a 8, el del cierre, por la enmienda y no por el recuento.",
        "correccion": "T1",
    }]), salida / "changelog_numeros.tsv", sello)
    escribir_k2(settings, salida, sello)
    replica = _por_gen(principal, mapa)
    cruzado = alinear_por_clave(
        certificadas,
        replica,
        clave_descubrimiento=certificadas["protein_id"].astype(str).str.upper(),
        clave_replica=replica["gen"],
    )
    escribir_k3(cruzado, salida, sello)
    print(f"[v3.2] tejido R1 {r1['numerador']}/{r1['denominador']} con {len(malos)} filas fuera")
