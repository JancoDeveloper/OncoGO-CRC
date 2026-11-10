#!/usr/bin/env python3
"""Lee solo results/v3_6/ y escribe las macros de la tesis y el JSON de Word."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from oncogo_crc.formato import formatear

REPO = Path(__file__).resolve().parents[2]
V36 = REPO / "results" / "v3_6"
TEX = REPO / "tesis-ucaldas" / "latex" / "cifras.tex"
WORD = REPO / "tesis-ucaldas" / "word" / "cifras.json"
MD = REPO / "tesis-ucaldas" / "paquete-redaccion-v3_6" / "cifras-v3_6.md"


def _json(nombre: str) -> dict:
    return json.loads((V36 / nombre).read_text(encoding="utf-8"))


def _tex(texto: str) -> str:
    return (
        texto.replace("\\", "\\textbackslash{}")
        .replace("&", "\\&")
        .replace("%", "\\%")
        .replace("#", "\\#")
        .replace("_", "\\_")
    )


def main() -> None:
    embudo = _json("suero_embudo.json")
    controles = embudo["controles"]
    celda = embudo["celda_preregistrada"]
    ml = _json("ml_suero_permutacion.json")["versiones"]
    multi = _json("multiverso_resumen.json")
    tejido = _json("k1_etiquetas.json")["contrastes"]
    h5 = _json("h5_especificidad_organo.json")
    h6 = _json("h6_origen_celular.json")
    h7 = _json("h7_visibilidad_en_sangre.json")
    h8 = _json("h8_carga_tumoral.json")
    h9 = _json("h9_wald.json")
    h10 = _json("h10_especificidad.json")
    pre = _json("preanalitica.json")
    calidad = _json("k7_control_tecnico.json")["confusion_preanalitica"]
    certificado = _json("certificado.json")
    piel = {fila["indicador"]: fila for fila in calidad}
    auc = pre["calculado_con_legibles"]["auc"]
    spearman = pre["calculado_con_legibles"]["spearman"]["piel"]
    primario = tejido["primario_frente_a_mucosa"]
    metastasis = tejido["metastasis_hepatica_frente_a_mucosa"]
    genes = [
        fila["proteina"]
        for fila in h9.get("proteinas", [])
        if fila.get("estado") == "exploratorio"
    ]
    especificas = []
    tabla_h10 = pd.read_csv(V36 / "h10_especificidad.tsv", sep="\t")
    if "clase" in tabla_h10.columns:
        especificas = tabla_h10.loc[tabla_h10["clase"].eq("especifica_crc"), "gen"].astype(str).tolist()
    entradas = [
        ("SueroN", celda["n"], "entero", "results/v3_6/suero_embudo.json#celda_preregistrada.n"),
        ("SueroDif", celda["diferente"], "entero", "results/v3_6/suero_embudo.json#celda_preregistrada.diferente"),
        ("SueroEq", celda["equivalente"], "entero", "results/v3_6/suero_embudo.json#celda_preregistrada.equivalente"),
        ("SueroInd", celda["indeterminada"], "entero", "results/v3_6/suero_embudo.json#celda_preregistrada.indeterminada"),
        ("EmbudoGrupos", controles["grupos"], "entero", "results/v3_6/suero_embudo.json#controles.grupos"),
        ("EmbudoQ", controles["grupos_global_pg_q_0_01"], "entero", "results/v3_6/suero_embudo.json#controles.grupos_global_pg_q_0_01"),
        ("EmbudoHumanos", controles["humanos_pg_y_global"], "entero", "results/v3_6/suero_embudo.json#controles.humanos_pg_y_global"),
        ("MlTodas", ml["M2_proteinas_principal"]["auc_oof_media"], "auc", "results/v3_6/ml_suero_permutacion.json#M2_proteinas_principal.auc_oof_media"),
        ("MlCalidad", ml["M3_solo_qc"]["auc_oof_media"], "auc", "results/v3_6/ml_suero_permutacion.json#M3_solo_qc.auc_oof_media"),
        ("MlSinCalidad", ml["M5_principal_residualizada_qc"]["auc_oof_media"], "auc", "results/v3_6/ml_suero_permutacion.json#M5_principal_residualizada_qc.auc_oof_media"),
        ("MlSinCalidadP", ml["M5_principal_residualizada_qc"]["p_permutacion"], "p", "results/v3_6/ml_suero_permutacion.json#M5_principal_residualizada_qc.p_permutacion"),
        ("MlHumana", ml["M1_proteinas_filtrada_humana"]["auc_oof_media"], "auc", "results/v3_6/ml_suero_permutacion.json#M1_proteinas_filtrada_humana.auc_oof_media"),
        ("MlCoagulacion", ml["M4_solo_coagulacion"]["auc_oof_media"], "auc", "results/v3_6/ml_suero_permutacion.json#M4_solo_coagulacion.auc_oof_media"),
        ("MlPca", ml["M6_principal_pca5"]["auc_oof_media"], "auc", "results/v3_6/ml_suero_permutacion.json#M6_principal_pca5.auc_oof_media"),
        ("MlContaminantes", float(pd.read_csv(V36 / "multiverso_ml.tsv", sep="\t").loc[lambda t: t["contaminantes"].eq("sin_quitar"), "auc_oof_media"].median()), "auc", "results/v3_6/multiverso_ml.tsv#sin_quitar.auc_oof_media"),
        ("TejidoNum", primario["numerador"], "entero", "results/v3_6/k1_etiquetas.json#primario_frente_a_mucosa.numerador"),
        ("TejidoDen", primario["denominador"], "entero", "results/v3_6/k1_etiquetas.json#primario_frente_a_mucosa.denominador"),
        ("TejidoMetaNum", metastasis["numerador"], "entero", "results/v3_6/k1_etiquetas.json#metastasis_hepatica_frente_a_mucosa.numerador"),
        ("TejidoMetaDen", metastasis["denominador"], "entero", "results/v3_6/k1_etiquetas.json#metastasis_hepatica_frente_a_mucosa.denominador"),
        ("HcincoEstado", h5["estado"], "estado", "results/v3_6/h5_especificidad_organo.json#estado"),
        ("HcincoUniverso", h5["universo_principal_n"], "entero", "results/v3_6/h5_especificidad_organo.json#universo_principal_n"),
        ("HcincoDif", round(1 - float(h5["maxima_concordancia_otro_organo"]), 3), "auc", "results/v3_6/h5_especificidad_organo.json#maxima_concordancia_otro_organo"),
        ("HseisEstado", h6["estado"], "estado", "results/v3_6/h6_origen_celular.json#estado"),
        ("HsieteEstado", h7["estado"], "estado", "results/v3_6/h7_visibilidad_en_sangre.json#estado"),
        ("HsieteAlvez", h7["por_recurso"]["alvez_2023"]["auc"], "auc", "results/v3_6/h7_visibilidad_en_sangre.json#por_recurso.alvez_2023.auc"),
        ("HsieteHpa", h7["por_recurso"]["hpa_disease_blood_atlas"]["auc"], "auc", "results/v3_6/h7_visibilidad_en_sangre.json#por_recurso.hpa_disease_blood_atlas.auc"),
        ("HochoEstado", h8["estado"], "estado", "results/v3_6/h8_carga_tumoral.json#estado"),
        ("HochoUniverso", h8["n_universo"], "entero", "results/v3_6/h8_carga_tumoral.json#n_universo"),
        ("HochoTumorales", h8["n_tumorales_epiteliales"], "entero", "results/v3_6/h8_carga_tumoral.json#n_tumorales_epiteliales"),
        ("HochoPlaquetasP", h8["qc"]["plaquetas"]["p_kruskal"], "p", "results/v3_6/h8_carga_tumoral.json#qc.plaquetas.p_kruskal"),
        ("HnueveEstado", h9["estado"], "estado", "results/v3_6/h9_wald.json#estado"),
        ("HdiezEstado", h10["estado"], "estado", "results/v3_6/h10_especificidad.json#estado"),
        ("HdiezN", h10["n_especificas_crc"], "entero", "results/v3_6/h10_especificidad.json#n_especificas_crc"),
        ("PreEstado", pre["estado"], "estado", "results/v3_6/preanalitica.json#estado"),
        ("PreAuc", auc["valor"], "auc", "results/v3_6/preanalitica.json#calculado_con_legibles.auc.valor"),
        ("PreP", auc["p_exacta_mann_whitney"], "p", "results/v3_6/preanalitica.json#calculado_con_legibles.auc.p_exacta_mann_whitney"),
        ("PreCotaMax", pre["cota_con_60"]["auc_max"], "auc", "results/v3_6/preanalitica.json#cota_con_60.auc_max"),
        ("PreRhoPiel", spearman["rho"], "rho", "results/v3_6/preanalitica.json#spearman.piel.rho"),
        ("PreFecha", pre["n_con_fecha"], "entero", "results/v3_6/preanalitica.json#n_con_fecha"),
        ("PreCrudos", pre["n_crudos_validos"], "entero", "results/v3_6/preanalitica.json#n_crudos_validos"),
        ("CalidadCoagP", piel["coagulacion"]["p_kruskal"], "p", "results/v3_6/k7_control_tecnico.json#coagulacion.p_kruskal"),
        ("CalidadPielP", piel["piel"]["p_kruskal"], "p", "results/v3_6/k7_control_tecnico.json#piel.p_kruskal"),
        ("CertificadoEstado", certificado["estado"], "estado", "results/v3_6/certificado.json#estado"),
        ("MultiSuero", multi["n_celdas_suero"], "entero", "results/v3_6/multiverso_resumen.json#n_celdas_suero"),
        ("MultiMl", multi["n_celdas_ml"], "entero", "results/v3_6/multiverso_resumen.json#n_celdas_ml"),
        ("MultiDif", multi["n_con_diferente"], "entero", "results/v3_6/multiverso_resumen.json#n_con_diferente"),
    ]
    lineas = ["% Generado por scripts/v3/exportar_cifras_tesis.py. No editar a mano.", "% fuente: results/v3_6/"]
    word = {}
    for nombre, valor, tipo, fuente in entradas:
        if valor is None:
            raise SystemExit(f"falta {nombre}")
        texto = formatear(valor, tipo)
        lineas.append(f"% fuente: {fuente}")
        lineas.append(f"\\newcommand{{\\{nombre}}}{{{_tex(texto)}}}")
        word[nombre] = {"valor": texto, "fuente": fuente, "crudo": valor}
    lineas.append("% fuente: results/v3_6/h10_especificidad.tsv")
    lineas.append("\\newcommand{\\HdiezGenes}{" + _tex(", ".join(especificas)) + "}")
    word["HdiezGenes"] = {"valor": ", ".join(especificas), "fuente": "results/v3_6/h10_especificidad.tsv", "crudo": especificas}
    lineas.append("% fuente: results/v3_6/h9_wald.json")
    lineas.append("\\newcommand{\\HnueveGenes}{" + _tex(", ".join(genes)) + "}")
    word["HnueveGenes"] = {"valor": ", ".join(genes), "fuente": "results/v3_6/h9_wald.json", "crudo": genes}
    lineas.append("\\newcommand{\\Advertencia}{Ningún resultado de esta tesis sirve para decidir sobre un paciente.}")
    word["Advertencia"] = {"valor": "Ningún resultado de esta tesis sirve para decidir sobre un paciente.", "fuente": "results/v3_6/certificado.json", "crudo": None}
    TEX.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    WORD.write_text(json.dumps(word, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    MD.parent.mkdir(parents=True, exist_ok=True)
    filas_md = ["| Macro | Valor | Fichero |", "| --- | --- | --- |"]
    for nombre, bloque in word.items():
        filas_md.append(f"| \\{nombre} | {bloque['valor']} | {bloque['fuente']} |")
    MD.write_text("# Cifras v3.6\n\nSalen de `results/v3_6/` por `scripts/v3/exportar_cifras_tesis.py`.\n\n" + "\n".join(filas_md) + "\n", encoding="utf-8")
    print(f"{len(word)} macros")


if __name__ == "__main__":
    main()
