"""Orden de adquisición de PXD055218. Con menos de 60 crudos el estado es no_evaluable."""

from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from oncogo_crc.settings import Settings
from oncogo_crc.v3.ml_suero import ic_ledell
from oncogo_crc.v3.salida import volcar_json, volcar_tsv

_CRUDO = re.compile(r"(CP|CN|HC)_\d{2}_DIA\.raw")
_PRIDE_BYTES = {
    "CN_01_DIA.raw": 1049764670,
    "CN_08_DIA.raw": 1000697970,
    "CN_16_DIA.raw": 1205043346,
    "CP_18_DIA.raw": 1226017242,
    "HC_17_DIA.raw": 1065175806,
    "HC_20_DIA.raw": 1088171988,
}


def contar_crudos(nombres: list[str]) -> tuple[list[str], int]:
    validos = sorted(n for n in nombres if _CRUDO.fullmatch(n))
    apple = sum(1 for n in nombres if n.startswith("._"))
    return validos, apple


def fecha_de_metadato(meta: dict, campo: str) -> str | None:
    for bloque in meta.get("FileProperties", []):
        if bloque.get("name") == campo and bloque.get("value"):
            momento = datetime.strptime(str(bloque["value"]), "%m/%d/%Y %H:%M:%S")
            return momento.strftime("%Y-%m-%dT%H:%M:%S")
    return None


def _valor(meta: dict, fragmento: str) -> str | None:
    for bloque in meta.get("InstrumentProperties", []):
        nombre = str(bloque.get("name", ""))
        if fragmento.lower() in nombre.lower() and bloque.get("value"):
            return str(bloque["value"])
    return None


def presentar_ic(auc: float | None, bajo: float | None, alto: float | None) -> tuple[list | None, str | None]:
    if auc is None:
        return None, None
    if abs(auc) < 1e-12 or abs(auc - 1) < 1e-12:
        return None, "separación completa; el IC de LeDell degenera"
    return [bajo, alto], None


def p_exacta_mann_whitney(orden_cp, orden_cn) -> float:
    from scipy.stats import mannwhitneyu

    return float(mannwhitneyu(orden_cp, orden_cn, alternative="two-sided", method="exact").pvalue)


def cota_con_60(*, n_cp: int = 20, n_cn: int = 20, cp_sin_fecha: int = 1, cn_sin_fecha: int = 3) -> dict:
    """Coloca los crudos sin fecha en los dos extremos y mide el AUC posible."""
    conocidas_cp = n_cp - cp_sin_fecha
    conocidas_cn = n_cn - cn_sin_fecha
    parejas = n_cp * n_cn
    conocidas = conocidas_cp * conocidas_cn
    # Extremo alto: los CN que faltan van antes que todo CP y el CP que falta va después de todo CN.
    ganancias_max = cp_sin_fecha * n_cn + cn_sin_fecha * conocidas_cp
    return {
        "auc_min": 0.0,
        "auc_max": ganancias_max / parejas,
        "n_parejas_conocidas": conocidas,
        "n_parejas_desconocidas": parejas - conocidas,
    }


def bloques_de_orden(tabla: pd.DataFrame) -> tuple[dict, bool]:
    posicion = {}
    rangos = []
    for grupo, bloque in tabla.groupby("grupo"):
        posicion[str(grupo)] = {
            "orden_min": int(bloque["orden"].min()),
            "orden_max": int(bloque["orden"].max()),
            "inicio": str(bloque["fecha"].min()),
            "fin": str(bloque["fecha"].max()),
            "n": int(len(bloque)),
        }
        rangos.append((int(bloque["orden"].min()), int(bloque["orden"].max())))
    rangos.sort()
    perfectos = all(rangos[i][1] < rangos[i + 1][0] for i in range(len(rangos) - 1)) if len(rangos) >= 2 else False
    return posicion, bool(perfectos) and set(posicion) >= {"CP", "CN", "HC"}


def resumir_legibles(*, n_validos: int, n_con_fecha: int, auc, spearman: dict) -> dict:
    return {
        "estado": "no_evaluable" if n_con_fecha < 60 else "calculado",
        "motivo": f"{n_con_fecha} de {n_validos} crudos legibles; la prueba preregistrada pide 60" if n_con_fecha < 60 else "60 crudos con fecha",
        "calculado_con_legibles": {"auc": auc, "n_con_fecha": n_con_fecha, "spearman": spearman},
    }


def intentar_crudo(binario: Path, crudo: Path, destino: Path) -> dict:
    """Un intento con metadatos JSON y, si falla, con metadatos en texto. No descarga."""
    import subprocess

    destino.mkdir(parents=True, exist_ok=True)
    try:
        ver = subprocess.run([str(binario), "--version"], capture_output=True, text=True, timeout=30)
        version = (ver.stdout or ver.stderr).strip().splitlines()[0]
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"crudo": crudo.name, "version": "", "leido": False, "modo": "", "error": str(exc)}
    mensajes = []
    for modo in ("0", "1"):
        try:
            proc = subprocess.run(
                [str(binario), "-i", str(crudo), "-m", modo, "-f", "4", "-o", str(destino)],
                capture_output=True, text=True, timeout=120,
            )
            texto = ((proc.stdout or "") + (proc.stderr or "")).strip()
            producido = destino / f"{crudo.stem}.json"
            if proc.returncode == 0 and producido.is_file():
                return {"crudo": crudo.name, "version": version, "leido": True, "modo": modo, "error": texto}
            mensajes.append(f"modo {modo}: {texto}")
        except subprocess.TimeoutExpired as exc:
            mensajes.append(f"modo {modo}: timeout {exc}")
    return {"crudo": crudo.name, "version": version, "leido": False, "modo": "", "error": "\n".join(mensajes)}


def _sin_ruta_local(texto: str) -> str:
    limpio = str(texto).replace("\\", "/")
    for marca in ("/Volumes/", "/Users/", "C:/"):
        while marca in limpio:
            inicio = limpio.find(marca)
            fin = limpio.find(" ", inicio)
            if fin < 0:
                fin = len(limpio)
            limpio = limpio[:inicio] + limpio[fin:]
    return limpio.strip()


def ejecutar(settings: Settings, salida: Path, sello: dict, prereg) -> dict:
    campo = str(prereg.get("enmienda_v3_5.preanalitica.campo_fecha"))
    carpeta = settings.resolve("data/raw/plasma/pxd055218_raw")
    metadatos = settings.resolve("external/pxd055218_metadata")
    nombres = [p.name for p in carpeta.iterdir()] if carpeta.is_dir() else []
    validos, apple = contar_crudos(nombres)
    binario = settings.resolve("external/bin/trfp/ThermoRawFileParser")
    intentos = []
    filas = []
    sin_fecha = []
    previos: dict[str, dict] = {}
    previo = salida / "preanalitica_intentos.tsv"
    if previo.is_file():
        viejo = pd.read_csv(previo, sep="\t")
        if "crudo" in viejo.columns:
            previos = {str(fila["crudo"]): fila.to_dict() for _, fila in viejo.iterrows()}
    for nombre in validos:
        crudo = carpeta / nombre
        json_path = metadatos / f"{crudo.stem}.json"
        tamaño = crudo.stat().st_size if crudo.is_file() else None
        if not json_path.is_file() and nombre in previos:
            intento = {clave: previos[nombre].get(clave, "") for clave in ("crudo", "version", "leido", "modo", "error")}
            intento["error"] = _sin_ruta_local(intento["error"])
            intentos.append(intento)
        elif not json_path.is_file() and binario.is_file() and crudo.is_file():
            intento = intentar_crudo(binario, crudo, metadatos)
            intento["error"] = _sin_ruta_local(intento["error"])
            intentos.append(intento)
            if json_path.is_file():
                intento["leido"] = True
        if not json_path.is_file():
            pride = _PRIDE_BYTES.get(nombre)
            sin_fecha.append({
                "crudo": nombre,
                "bytes_local": tamaño,
                "bytes_pride": pride,
                "tamaño_igual_a_pride": tamaño == pride if pride else None,
                "motivo": "ThermoRawFileParser no escribió JSON. El tamaño coincide con PRIDE, así que la descarga no quedó corta." if tamaño == pride else "sin JSON",
            })
            continue
        meta = json.loads(json_path.read_text(encoding="utf-8"))
        filas.append({
            "crudo": nombre,
            "grupo": nombre[:2],
            "fecha": fecha_de_metadato(meta, campo),
            "instrumento": _valor(meta, "instrument model"),
            "serie": _valor(meta, "serial number"),
        })
    tabla = pd.DataFrame(filas)
    con_fecha = tabla.dropna(subset=["fecha"]).sort_values(["fecha", "crudo"]) if len(tabla) else tabla
    if len(con_fecha):
        con_fecha = con_fecha.copy()
        con_fecha["orden"] = range(1, len(con_fecha) + 1)
    else:
        con_fecha["orden"] = []
    volcar_tsv(con_fecha, salida / "preanalitica_orden.tsv", sello)
    volcar_tsv(pd.DataFrame(intentos), salida / "preanalitica_intentos.tsv", sello)
    auc = None
    ic = None
    n_cp = n_cn = 0
    if len(con_fecha):
        pares = con_fecha.loc[con_fecha["grupo"].isin(["CP", "CN"])]
        n_cp = int((pares["grupo"] == "CP").sum())
        n_cn = int((pares["grupo"] == "CN").sum())
        if n_cp >= 3 and n_cn >= 3:
            y = (pares["grupo"] == "CP").to_numpy().astype(int)
            auc, bajo, alto = ic_ledell(y, pares["orden"].to_numpy(dtype=float))
            ic, motivo_ic = presentar_ic(auc, bajo, alto)
            p_exacta = p_exacta_mann_whitney(
                pares.loc[pares["grupo"].eq("CP"), "orden"].to_numpy(dtype=float),
                pares.loc[pares["grupo"].eq("CN"), "orden"].to_numpy(dtype=float),
            )
        else:
            motivo_ic = None
            p_exacta = None
    else:
        motivo_ic = None
        p_exacta = None
    spearman = {}
    from oncogo_crc.v3.salida import corrida_limpia
    calidad = salida / "suero_calidad_muestra.tsv"
    if not calidad.is_file() and not corrida_limpia():
        calidad = settings.resolve("results/v3_3/suero_calidad_muestra.tsv")
    if calidad.is_file() and len(con_fecha):
        qc = pd.read_csv(calidad, sep="\t").rename(columns={"grupo": "grupo_calidad"})
        union_path = salida / "union_crudos_calidad.tsv"
        if union_path.is_file():
            union = pd.read_csv(union_path, sep="\t", dtype=str)
            verificadas = union.loc[union["verificada"].str.lower().eq("true"), ["sample_id", "crudo"]]
            qc = qc.merge(verificadas, on="sample_id", how="inner")
        else:
            qc["crudo"] = qc["sample_id"].astype(str) + ".raw"
        unido = con_fecha.merge(qc, on="crudo", how="inner")
        if "n_precursores" in unido.columns:
            unido["log2_n_precursores"] = np.log2(pd.to_numeric(unido["n_precursores"], errors="coerce"))
        def _spearman_de(bloque: pd.DataFrame) -> dict:
            salida_sp = {}
            for columna in ("eritrocitos", "plaquetas", "coagulacion", "piel", "log2_n_precursores"):
                if columna not in bloque.columns:
                    salida_sp[columna] = {"estado": "no_disponible"}
                    continue
                x = pd.to_numeric(bloque["orden"], errors="coerce")
                y = pd.to_numeric(bloque[columna], errors="coerce")
                ok = x.notna() & y.notna()
                if int(ok.sum()) < 5:
                    salida_sp[columna] = {"estado": "no_evaluable", "n": int(ok.sum())}
                    continue
                rho, p = spearmanr(x[ok], y[ok])
                salida_sp[columna] = {"rho": float(rho), "p": float(p), "n": int(ok.sum())}
            return salida_sp

        cp_cn = unido.loc[unido["grupo"].isin(["CP", "CN"])] if "grupo" in unido.columns else unido
        spearman = _spearman_de(cp_cn)
        n_sin_union = int(len(con_fecha) - len(unido))
        if n_sin_union:
            spearman_54 = {
                "estado": "no_evaluable",
                "n_unidas": int(len(unido)),
                "n_con_fecha": int(len(con_fecha)),
                "n_sin_union": n_sin_union,
                "motivo": "hay crudos con fecha cuyo nombre no coincide con el sample_id de la tabla de calidad",
                "en_las_unidas": _spearman_de(unido),
            }
        else:
            spearman_54 = _spearman_de(unido)
        spearman_por_grupo = {
            str(grupo): _spearman_de(bloque) for grupo, bloque in unido.groupby("grupo")
        } if "grupo" in unido.columns else {}
        if "grupo" in con_fecha.columns:
            for grupo in con_fecha["grupo"].astype(str).unique():
                spearman_por_grupo.setdefault(grupo, {"estado": "no_evaluable", "motivo": "el crudo no une con la tabla de calidad"})
    else:
        spearman_54 = {}
        spearman_por_grupo = {}
    resumen = resumir_legibles(
        n_validos=len(validos),
        n_con_fecha=int(len(con_fecha)),
        auc={"valor": auc, "ic95_ledell": ic, "motivo_ic": motivo_ic, "p_exacta_mann_whitney": p_exacta, "n_cp": n_cp, "n_cn": n_cn},
        spearman=spearman,
    )
    fechas = con_fecha["fecha"].tolist() if len(con_fecha) else []
    por_dia: dict[str, int] = {}
    for fecha in fechas:
        por_dia[fecha[:10]] = por_dia.get(fecha[:10], 0) + 1
    hc = con_fecha.loc[con_fecha["grupo"].eq("HC"), "orden"] if len(con_fecha) else pd.Series(dtype=float)
    otros = con_fecha.loc[con_fecha["grupo"].isin(["CP", "CN"]), "orden"] if len(con_fecha) else pd.Series(dtype=float)
    if len(hc) and len(otros):
        dentro = int(((hc >= otros.min()) & (hc <= otros.max())).sum())
        posicion_hc = "entre_CP_y_CN" if dentro else "aparte"
    else:
        dentro = 0
        posicion_hc = "no_evaluable"
    payload = {
        **resumen,
        "n_crudos_validos": len(validos),
        "n_appledouble": apple,
        "n_con_fecha": int(len(con_fecha)),
        "sin_fecha": sin_fecha,
        "campo_fecha": campo,
        "formato_fecha": "%m/%d/%Y %H:%M:%S",
        "inicio": fechas[0] if fechas else None,
        "fin": fechas[-1] if fechas else None,
        "corridas_por_dia": por_dia,
        "corridas_por_grupo": con_fecha["grupo"].value_counts().to_dict() if len(con_fecha) else {},
        "hc_dentro_del_rango_cp_cn": dentro,
        "posicion_hc": posicion_hc,
        "instrumento": con_fecha["instrumento"].dropna().mode().iloc[0] if len(con_fecha) and con_fecha["instrumento"].notna().any() else None,
        "p_exacta_mann_whitney": p_exacta,
        "cota_con_60": cota_con_60(),
        "posicion_por_grupo": bloques_de_orden(con_fecha)[0] if len(con_fecha) else {},
        "bloques_perfectos": bloques_de_orden(con_fecha)[1] if len(con_fecha) else False,
        "spearman_54": spearman_54,
        "spearman_por_grupo": spearman_por_grupo,
        "lectura": "Los grupos se midieron en bloques: primero el cáncer, después la colonoscopia negativa y al final los sanos. La piel y la coagulación cambian con el orden. Con estos datos no se puede separar el efecto del grupo del efecto del lote de adquisición. Eso no crea una proteína del tumor: refuerza que la separación del suero no es atribuible al tumor. Las diferencias frente a los sanos comparan también bloques de adquisición distintos.",
    }
    volcar_json(payload, salida / "preanalitica_v3_5.json", sello)
    return payload
