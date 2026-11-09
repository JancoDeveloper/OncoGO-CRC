"""Cifras que cambian entre results/v3_2 y results/v3_3. Lee los JSON; no calcula de nuevo."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def _json(path: Path) -> dict:
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _fila(metrica: str, viejo, nuevo, causa: str, fichero: str, nota: str = "") -> dict | None:
    if viejo == nuevo:
        return None
    return {
        "metrica": metrica,
        "valor_v3_2": viejo if viejo is not None else "no_calculado",
        "valor_v3_3": nuevo if nuevo is not None else "no_calculado",
        "causa": causa,
        "fichero": fichero,
        "nota": nota,
    }


def _clase(resumen: dict, version: str, contraste: str, campo: str):
    return (
        resumen.get("clases", {})
        .get(version, {})
        .get(contraste, {})
        .get(campo)
    )


def _contraste(k1: dict, nombre: str) -> str | None:
    bloque = k1.get("contrastes", {}).get(nombre)
    if not bloque:
        return None
    return f"{bloque.get('numerador')}/{bloque.get('denominador')}"


def filas(viejo: Path, nuevo: Path) -> pd.DataFrame:
    """Una fila por cifra distinta. La causa es la corrección 2.1 a 2.7 que toca ese fichero."""
    suero_v = _json(viejo / "suero_resumen.json")
    suero_n = _json(nuevo / "suero_resumen.json")
    k1_v = _json(viejo / "k1_etiquetas.json")
    k1_n = _json(nuevo / "k1_etiquetas.json")
    norm = _json(nuevo / "k1_sensibilidad_normalizacion.json")
    ml_v = _json(viejo / "ml_suero_permutacion.json")
    ml_n = _json(nuevo / "ml_suero_permutacion.json")
    h5_v = _json(viejo / "h5_especificidad_organo.json")
    h5_n = _json(nuevo / "h5_especificidad_organo.json")
    h10_v = _json(viejo / "h10_especificidad.json")
    h10_n = _json(nuevo / "h10_especificidad.json")
    h7_n = _json(nuevo / "h7_visibilidad_en_sangre.json")
    salida: list[dict] = []

    def anotar(metrica, viejo_valor, nuevo_valor, causa, fichero, nota=""):
        fila = _fila(metrica, viejo_valor, nuevo_valor, causa, fichero, nota)
        if fila is not None:
            salida.append(fila)

    anotar("n_principal", suero_v.get("n_principal"), suero_n.get("n_principal"), "2.1", "suero_resumen.json")
    anotar("frase_suero", suero_v.get("frase"), suero_n.get("frase"), "2.1", "suero_resumen.json")
    for contraste in (
        "crc_frente_a_colonoscopia_negativa",
        "crc_frente_a_sanos",
        "colonoscopia_negativa_frente_a_sanos",
    ):
        for campo in ("n", "diferente", "equivalente", "indeterminada", "conflicto"):
            anotar(
                f"principal.{contraste}.{campo}",
                _clase(suero_v, "principal", contraste, campo),
                _clase(suero_n, "principal", contraste, campo),
                "2.1",
                "suero_resumen.json",
            )
    for nombre in (
        "primario_frente_a_mucosa",
        "metastasis_hepatica_frente_a_mucosa",
        "primario_frente_a_mucosa_no_pareado",
    ):
        anotar(nombre, _contraste(k1_v, nombre), _contraste(k1_n, nombre), "2.3", "k1_etiquetas.json")
    anotar(
        "h1_veredicto",
        k1_v.get("h1_veredicto"),
        k1_n.get("h1_veredicto"),
        "2.3",
        "k1_etiquetas.json",
    )
    anotar(
        "h1_fragil_a_normalizacion",
        k1_v.get("h1_fragil_a_normalizacion"),
        k1_n.get("h1_fragil_a_normalizacion"),
        "2.4",
        "k1_etiquetas.json",
    )
    antes = norm.get("antes") or {}
    despues = norm.get("despues") or {}
    anotar(
        "mediana_log2fc_pareado",
        None,
        antes.get("mediana_log2fc"),
        "2.4",
        "k1_sensibilidad_normalizacion.json",
        "v3.2 no escribió esta sensibilidad. El valor nuevo es el de antes de centrar.",
    )
    anotar(
        "mediana_log2fc_pareado_centrado",
        antes.get("mediana_log2fc"),
        despues.get("mediana_log2fc"),
        "2.4",
        "k1_sensibilidad_normalizacion.json",
    )
    anotar(
        "deming_intercepto",
        antes.get("deming_intercepto"),
        despues.get("deming_intercepto"),
        "2.4",
        "k1_sensibilidad_normalizacion.json",
        "Antes y después del centrado por mediana, en la misma corrida v3.3.",
    )
    pares_ml = (
        ("filtrada_humana", "M1_proteinas_filtrada_humana", "auc_media", "auc_oof_media"),
        ("sin_contaminantes_ni_paneles", "M2_proteinas_principal", "auc_media", "auc_oof_media"),
        ("solo_coagulacion", "M4_solo_coagulacion", "auc_media", "auc_oof_media"),
        ("filtrada_humana", "M1_proteinas_filtrada_humana", "p", "p_permutacion"),
        ("sin_contaminantes_ni_paneles", "M2_proteinas_principal", "p", "p_permutacion"),
        ("solo_coagulacion", "M4_solo_coagulacion", "p", "p_permutacion"),
    )
    versiones_v = ml_v.get("versiones") or {}
    versiones_n = ml_n.get("versiones") or {}
    for clave_v, clave_n, campo_v, campo_n in pares_ml:
        anotar(
            f"ml.{clave_n}.{campo_n}",
            (versiones_v.get(clave_v) or {}).get(campo_v),
            (versiones_n.get(clave_n) or {}).get(campo_n),
            "2.2",
            "ml_suero_permutacion.json",
            "v3.2 informa el AUC medio de pliegues con muestras duplicadas. v3.3 informa el AUC fuera de pliegue con 40 filas.",
        )
    anotar(
        "titular_ml",
        (ml_v.get("titular") or {}).get("clave"),
        (ml_n.get("titular") or {}).get("clave"),
        "2.2",
        "ml_suero_permutacion.json",
    )
    anotar("h5.estado", h5_v.get("estado"), h5_n.get("estado"), "2.5", "h5_especificidad_organo.json")
    anotar(
        "h5.fraccion_pan_carcinoma",
        h5_v.get("fraccion_pan_carcinoma"),
        h5_n.get("fraccion_pan_carcinoma"),
        "2.5",
        "h5_especificidad_organo.json",
    )
    anotar(
        "h10.estado",
        h10_v.get("estado"),
        h10_n.get("estado"),
        "2.6",
        "h10_especificidad.json",
    )
    anotar(
        "h10.n_especificas_crc",
        h10_v.get("n_especificas", h10_v.get("n_especificas_crc")),
        h10_n.get("n_especificas_crc"),
        "2.6",
        "h10_especificidad.json",
        "v3.2 guardó el recuento en n_especificas. v3.3 lo guarda en n_especificas_crc.",
    )
    hip_v = _hipotesis(viejo / "hipotesis_v3_2.tsv")
    anotar(
        "h7.estado",
        hip_v.get("H7"),
        h7_n.get("estado"),
        "2.7",
        "h7_visibilidad_en_sangre.json",
        "En v3.2 el estado estaba escrito en la tabla de hipótesis, sin modelo.",
    )
    for control, pasa_v, pasa_n in _controles(viejo / "controles.tsv", nuevo / "controles.tsv"):
        anotar(f"control.{control}.pasa", pasa_v, pasa_n, "2.3", "controles.tsv")
    return pd.DataFrame(salida)


def _hipotesis(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    tabla = pd.read_csv(path, sep="\t")
    return dict(zip(tabla["hipotesis"].astype(str), tabla["estado"].astype(str), strict=False))


def _controles(viejo: Path, nuevo: Path) -> list[tuple[str, str, str]]:
    if not viejo.is_file() or not nuevo.is_file():
        return []
    a = pd.read_csv(viejo, sep="\t")
    b = pd.read_csv(nuevo, sep="\t")
    pares = []
    for control in ("CA1", "AQP8", "SLC26A3", "GUCA2A"):
        fila_a = a.loc[a["control"].astype(str).eq(control)]
        fila_b = b.loc[b["control"].astype(str).eq(control)]
        if fila_a.empty or fila_b.empty:
            continue
        pares.append((control, str(fila_a["pasa"].iloc[0]), str(fila_b["pasa"].iloc[0])))
    return pares
