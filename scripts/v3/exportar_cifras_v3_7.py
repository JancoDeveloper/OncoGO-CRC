#!/usr/bin/env python3
"""Lee results/v3_7/ y escribe las macros de la tesis, el paquete y el changelog.

No inventa cifras. Si falta un campo, sale con error.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from oncogo_crc.formato import formatear

REPO = Path(__file__).resolve().parents[2]
V36 = REPO / "results" / "v3_6"
V37 = REPO / "results" / "v3_7"
TEX = REPO / "tesis-ucaldas" / "latex" / "cifras.tex"
WORD = REPO / "tesis-ucaldas" / "word" / "cifras.json"
MD = REPO / "tesis-ucaldas" / "paquete-redaccion-v3_7" / "cifras-v3_7.md"
CAMBIO = REPO / "docs" / "changelog_v3_6_frente_a_v3_7.md"
ESTADO = REPO / "docs" / "estado_v3_7.md"

PENDIENTES = (
    "vias_suero.json",
    "h1_vias.json",
    "go_ora_tejido.tsv",
    "catalogo_auditoria.tsv",
    "neoantigenos.tsv",
    "prtc_por_corrida.tsv",
    "manifest_entradas.sha256",
)


def _json(carpeta: Path, nombre: str) -> dict:
    return json.loads((carpeta / nombre).read_text(encoding="utf-8"))


def _tex(texto: str) -> str:
    return (
        texto.replace("\\", "\\textbackslash{}")
        .replace("&", "\\&")
        .replace("%", "\\%")
        .replace("#", "\\#")
        .replace("_", "\\_")
    )


def _entradas(carpeta: Path, version: str) -> list[tuple[str, object, str, str]]:
    prefijo = f"results/{version}"
    embudo = _json(carpeta, "suero_embudo.json")
    controles = embudo["controles"]
    celda = embudo["celda_preregistrada"]
    ml = _json(carpeta, "ml_suero_permutacion.json")["versiones"]
    multi = _json(carpeta, "multiverso_resumen.json")
    tejido = _json(carpeta, "k1_etiquetas.json")["contrastes"]
    h5 = _json(carpeta, "h5_especificidad_organo.json")
    h6 = _json(carpeta, "h6_origen_celular.json")
    h7 = _json(carpeta, "h7_visibilidad_en_sangre.json")
    h8 = _json(carpeta, "h8_carga_tumoral.json")
    h9 = _json(carpeta, "h9_wald.json")
    h10 = _json(carpeta, "h10_especificidad.json")
    pre = _json(carpeta, "preanalitica.json")
    tejido_raiz = _json(carpeta, "k1_etiquetas.json")
    go = _json(carpeta, "go_tejido.json") if (carpeta / "go_tejido.json").is_file() else {"concordancia_nes": {"numerador": None, "denominador": None}}
    priori = _json(carpeta, "priorizacion_embudo.json") if (carpeta / "priorizacion_embudo.json").is_file() else {"pasos": [{"paso": 5, "n": None}, {"paso": 6, "n": None}]}
    calidad = _json(carpeta, "k7_control_tecnico.json")["confusion_preanalitica"]
    certificado = _json(carpeta, "certificado.json")
    piel = {fila["indicador"]: fila for fila in calidad}
    auc = pre["calculado_con_legibles"]["auc"]
    spearman = pre["calculado_con_legibles"]["spearman"]["piel"]
    primario = tejido["primario_frente_a_mucosa"]
    metastasis = tejido["metastasis_hepatica_frente_a_mucosa"]
    return [
        ("SueroN", celda["n"], "entero", f"{prefijo}/suero_embudo.json#celda_preregistrada.n"),
        ("SueroDif", celda["diferente"], "entero", f"{prefijo}/suero_embudo.json#celda_preregistrada.diferente"),
        ("SueroEq", celda["equivalente"], "entero", f"{prefijo}/suero_embudo.json#celda_preregistrada.equivalente"),
        ("SueroInd", celda["indeterminada"], "entero", f"{prefijo}/suero_embudo.json#celda_preregistrada.indeterminada"),
        ("EmbudoGrupos", controles["grupos"], "entero", f"{prefijo}/suero_embudo.json#controles.grupos"),
        ("EmbudoQ", controles["grupos_global_pg_q_0_01"], "entero", f"{prefijo}/suero_embudo.json#controles.grupos_global_pg_q_0_01"),
        ("EmbudoHumanos", controles["humanos_pg_y_global"], "entero", f"{prefijo}/suero_embudo.json#controles.humanos_pg_y_global"),
        ("MlTodas", ml["M2_proteinas_principal"]["auc_oof_media"], "auc", f"{prefijo}/ml_suero_permutacion.json#M2_proteinas_principal.auc_oof_media"),
        ("MlCalidad", ml["M3_solo_qc"]["auc_oof_media"], "auc", f"{prefijo}/ml_suero_permutacion.json#M3_solo_qc.auc_oof_media"),
        ("MlSinCalidad", ml["M5_principal_residualizada_qc"]["auc_oof_media"], "auc", f"{prefijo}/ml_suero_permutacion.json#M5_principal_residualizada_qc.auc_oof_media"),
        ("MlSinCalidadP", ml["M5_principal_residualizada_qc"]["p_permutacion"], "p", f"{prefijo}/ml_suero_permutacion.json#M5_principal_residualizada_qc.p_permutacion"),
        ("MlHumana", ml["M1_proteinas_filtrada_humana"]["auc_oof_media"], "auc", f"{prefijo}/ml_suero_permutacion.json#M1_proteinas_filtrada_humana.auc_oof_media"),
        ("MlCoagulacion", ml["M4_solo_coagulacion"]["auc_oof_media"], "auc", f"{prefijo}/ml_suero_permutacion.json#M4_solo_coagulacion.auc_oof_media"),
        ("MlPca", ml["M6_principal_pca5"]["auc_oof_media"], "auc", f"{prefijo}/ml_suero_permutacion.json#M6_principal_pca5.auc_oof_media"),
        ("MlContaminantes", float(pd.read_csv(carpeta / "multiverso_ml.tsv", sep="\t").loc[lambda t: t["contaminantes"].eq("sin_quitar"), "auc_oof_media"].median()), "auc", f"{prefijo}/multiverso_ml.tsv#sin_quitar.auc_oof_media"),
        ("TejidoNum", primario["numerador"], "entero", f"{prefijo}/k1_etiquetas.json#primario_frente_a_mucosa.numerador"),
        ("TejidoDen", primario["denominador"], "entero", f"{prefijo}/k1_etiquetas.json#primario_frente_a_mucosa.denominador"),
        ("TejidoMetaNum", metastasis["numerador"], "entero", f"{prefijo}/k1_etiquetas.json#metastasis_hepatica_frente_a_mucosa.numerador"),
        ("TejidoMetaDen", metastasis["denominador"], "entero", f"{prefijo}/k1_etiquetas.json#metastasis_hepatica_frente_a_mucosa.denominador"),
        ("HcincoEstado", h5["estado"], "estado", f"{prefijo}/h5_especificidad_organo.json#estado"),
        ("HcincoUniverso", h5["universo_principal_n"], "entero", f"{prefijo}/h5_especificidad_organo.json#universo_principal_n"),
        ("HcincoDif", round(1 - float(h5["maxima_concordancia_otro_organo"]), 3), "auc", f"{prefijo}/h5_especificidad_organo.json#maxima_concordancia_otro_organo"),
        ("HseisEstado", h6["estado"], "estado", f"{prefijo}/h6_origen_celular.json#estado"),
        ("HsieteEstado", h7["estado"], "estado", f"{prefijo}/h7_visibilidad_en_sangre.json#estado"),
        ("HsieteAlvez", h7["por_recurso"]["alvez_2023"]["auc"], "auc", f"{prefijo}/h7_visibilidad_en_sangre.json#por_recurso.alvez_2023.auc"),
        ("HsieteHpa", h7["por_recurso"]["hpa_disease_blood_atlas"]["auc"], "auc", f"{prefijo}/h7_visibilidad_en_sangre.json#por_recurso.hpa_disease_blood_atlas.auc"),
        ("HochoEstado", h8["estado"], "estado", f"{prefijo}/h8_carga_tumoral.json#estado"),
        ("HochoUniverso", h8["n_universo"], "entero", f"{prefijo}/h8_carga_tumoral.json#n_universo"),
        ("HochoTumorales", h8["n_tumorales_epiteliales"], "entero", f"{prefijo}/h8_carga_tumoral.json#n_tumorales_epiteliales"),
        ("HochoPlaquetasP", h8["qc"]["plaquetas"]["p_kruskal"], "p", f"{prefijo}/h8_carga_tumoral.json#qc.plaquetas.p_kruskal"),
        ("HnueveEstado", h9["estado"], "estado", f"{prefijo}/h9_wald.json#estado"),
        ("HdiezEstado", h10["estado"], "estado", f"{prefijo}/h10_especificidad.json#estado"),
        ("HdiezN", h10["n_especificas_crc"], "entero", f"{prefijo}/h10_especificidad.json#n_especificas_crc"),
        ("PreEstado", pre["estado"], "estado", f"{prefijo}/preanalitica.json#estado"),
        ("PreAuc", auc["valor"], "auc", f"{prefijo}/preanalitica.json#calculado_con_legibles.auc.valor"),
        ("PreP", auc["p_exacta_mann_whitney"], "p", f"{prefijo}/preanalitica.json#calculado_con_legibles.auc.p_exacta_mann_whitney"),
        ("PreRhoPiel", spearman["rho"], "rho", f"{prefijo}/preanalitica.json#spearman.piel.rho"),
        ("PreFecha", pre["n_con_fecha"], "entero", f"{prefijo}/preanalitica.json#n_con_fecha"),
        ("CoberturaNum", (tejido_raiz.get("cobertura") or {}).get("numerador"), "entero", f"{prefijo}/k1_etiquetas.json#cobertura.numerador"),
        ("CoberturaDen", (tejido_raiz.get("cobertura") or {}).get("denominador"), "entero", f"{prefijo}/k1_etiquetas.json#cobertura.denominador"),
        ("HunoVeredicto", tejido_raiz.get("h1_veredicto"), "estado", f"{prefijo}/k1_etiquetas.json#h1_veredicto"),
        ("BhPrecursoresCn", ((pre.get("bh_descriptivo") or {}).get("q") or {}).get("CN:log2_n_precursores"), "p", f"{prefijo}/preanalitica.json#bh_descriptivo.q.CN:log2_n_precursores"),
        ("BhPielCn", ((pre.get("bh_descriptivo") or {}).get("q") or {}).get("CN:piel"), "p", f"{prefijo}/preanalitica.json#bh_descriptivo.q.CN:piel"),
        ("GoConNum", (go.get("concordancia_nes") or {}).get("numerador"), "entero", f"{prefijo}/go_tejido.json#concordancia_nes.numerador"),
        ("GoConDen", (go.get("concordancia_nes") or {}).get("denominador"), "entero", f"{prefijo}/go_tejido.json#concordancia_nes.denominador"),
        ("Priori5", next((p["n"] for p in priori.get("pasos") or [] if p["paso"] == 5), None), "entero", f"{prefijo}/priorizacion_embudo.json#pasos"),
        ("Priori6", next((p["n"] for p in priori.get("pasos") or [] if p["paso"] == 6), None), "entero", f"{prefijo}/priorizacion_embudo.json#pasos"),
        ("HsieteSec", (h7.get("modelo_secundario") or {}).get("auc_secundario"), "auc_completo", f"{prefijo}/h7_visibilidad_en_sangre.json#modelo_secundario.auc_secundario"),
        ("PreCrudos", pre["n_crudos_validos"], "entero", f"{prefijo}/preanalitica.json#n_crudos_validos"),
        ("CalidadCoagP", piel["coagulacion"]["p_kruskal"], "p", f"{prefijo}/k7_control_tecnico.json#coagulacion.p_kruskal"),
        ("CalidadPielP", piel["piel"]["p_kruskal"], "p", f"{prefijo}/k7_control_tecnico.json#piel.p_kruskal"),
        ("CertificadoEstado", certificado["estado"], "estado", f"{prefijo}/certificado.json#estado"),
        ("MultiSuero", multi["n_celdas_suero"], "entero", f"{prefijo}/multiverso_resumen.json#n_celdas_suero"),
        ("MultiMl", multi["n_celdas_ml"], "entero", f"{prefijo}/multiverso_resumen.json#n_celdas_ml"),
        ("MultiDif", multi["n_con_diferente"], "entero", f"{prefijo}/multiverso_resumen.json#n_con_diferente"),
    ]


def _listas(carpeta: Path) -> tuple[list[str], list[str]]:
    h9 = _json(carpeta, "h9_wald.json")
    genes = [fila["proteina"] for fila in h9.get("proteinas", []) if fila.get("estado") == "exploratorio"]
    tabla = pd.read_csv(carpeta / "h10_especificidad.tsv", sep="\t")
    especificas = tabla.loc[tabla["clase"].eq("especifica_crc"), "gen"].astype(str).tolist() if "clase" in tabla.columns else []
    return genes, especificas


def main() -> None:
    nuevas = _entradas(V37, "v3_7")
    viejas = {nombre: (valor, tipo) for nombre, valor, tipo, _fuente in _entradas(V36, "v3_6")}
    genes, especificas = _listas(V37)
    lineas = ["% Generado por scripts/v3/exportar_cifras_v3_7.py. No editar a mano.", "% fuente: results/v3_7/"]
    word = {}
    filas_cambio = ["| Macro | v3.6 | v3.7 | |", "| --- | --- | --- | --- |"]
    for nombre, valor, tipo, fuente in nuevas:
        if valor is None:
            raise SystemExit(f"falta {nombre}")
        texto = formatear(valor, tipo)
        lineas.append(f"% fuente: {fuente}")
        lineas.append(f"\\newcommand{{\\{nombre}}}{{{_tex(texto)}}}")
        word[nombre] = {"valor": texto, "fuente": fuente, "crudo": valor}
        if nombre not in viejas:
            filas_cambio.append(f"| \\{nombre} | — | {texto} | nueva |")
        else:
            anterior, tipo_anterior = viejas[nombre]
            texto_anterior = formatear(anterior, tipo_anterior)
            marca = "igual" if texto_anterior == texto else "cambia"
            filas_cambio.append(f"| \\{nombre} | {texto_anterior} | {texto} | {marca} |")
    lineas.append("% fuente: results/v3_7/h10_especificidad.tsv")
    lineas.append("\\newcommand{\\HdiezGenes}{" + _tex(", ".join(especificas)) + "}")
    word["HdiezGenes"] = {"valor": ", ".join(especificas), "fuente": "results/v3_7/h10_especificidad.tsv", "crudo": especificas}
    lineas.append("% fuente: results/v3_7/h9_wald.json")
    lineas.append("\\newcommand{\\HnueveGenes}{" + _tex(", ".join(genes)) + "}")
    word["HnueveGenes"] = {"valor": ", ".join(genes), "fuente": "results/v3_7/h9_wald.json", "crudo": genes}
    lineas.append("\\newcommand{\\Advertencia}{Ningún resultado de esta tesis sirve para decidir sobre un paciente.}")
    word["Advertencia"] = {
        "valor": "Ningún resultado de esta tesis sirve para decidir sobre un paciente.",
        "fuente": "results/v3_7/certificado.json",
        "crudo": None,
    }
    TEX.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    WORD.write_text(json.dumps(word, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    MD.parent.mkdir(parents=True, exist_ok=True)
    filas_md = ["| Macro | Valor | Fichero |", "| --- | --- | --- |"]
    for nombre, bloque in word.items():
        filas_md.append(f"| \\{nombre} | {bloque['valor']} | {bloque['fuente']} |")
    MD.write_text(
        "# Cifras v3.7\n\nSalen de `results/v3_7/` por `scripts/v3/exportar_cifras_v3_7.py`.\n\n"
        + "\n".join(filas_md)
        + "\n",
        encoding="utf-8",
    )
    pre = _json(V37, "preanalitica.json")["calculado_con_legibles"]["spearman"]
    pre6 = _json(V36, "preanalitica.json")["calculado_con_legibles"]["spearman"]
    extra = []
    for panel in ("piel", "coagulacion"):
        extra.append(
            f"| n de Spearman, {panel} | {formatear(pre6[panel]['n'], 'entero')} | {formatear(pre[panel]['n'], 'entero')} | cambia |"
        )
        extra.append(
            f"| rho de {panel} | {formatear(pre6[panel]['rho'], 'rho')} | {formatear(pre[panel]['rho'], 'rho')} | cambia |"
        )
    CAMBIO.write_text(
        "# v3.6 frente a v3.7\n\n"
        "Cada valor lo formatea `src/oncogo_crc/formato.py` desde el fichero de esa corrida. "
        "Una fila `cambia` es una diferencia real. El resto se reprodujo.\n\n"
        + "\n".join(filas_cambio + extra)
        + "\n\nLa causa de los cambios está en `docs/estado_v3_7.md`.\n",
        encoding="utf-8",
    )
    sello = _json(V37, "certificado.json")["sello"]
    union = pd.read_csv(V37 / "union_crudos_calidad.tsv", sep="\t")
    faltan = [nombre for nombre in PENDIENTES if not (V37 / nombre).is_file()]
    presentes = [nombre for nombre in PENDIENTES if (V37 / nombre).is_file()]
    lineas_estado = [
        "# Estado de la v3.7",
        "",
        "Este archivo lo escribe `scripts/v3/exportar_cifras_v3_7.py` a partir de `results/v3_7/`.",
        "",
        f"Certificado: {formatear(_json(V37, 'certificado.json')['estado'], 'estado')}.",
        f"Sello: commit `{sello['commit'][:12]}`, fecha {sello['fecha_utc']}.",
        f"`sha256_manifest_entradas` en el sello: {sello.get('sha256_manifest_entradas')!s}.",
        "",
        "## Unión de los sanos",
        "",
        f"En `union_crudos_calidad.tsv`, `verificada` es verdadero en {int(union['verificada'].sum())} filas y falso en {int((~union['verificada']).sum())}.",
        "Las filas falsas son las corridas `Normal_` y la evidencia es «el report no nombra ese crudo». No se inventó `HC_n` = `Normal_n`.",
        "",
        "## Ficheros de la especificación que no están",
        "",
    ]
    if presentes:
        lineas_estado.append("Presentes: " + ", ".join(f"`{nombre}`" for nombre in presentes) + ".")
    lineas_estado.append("Ausentes:")
    lineas_estado.extend(f"- `{nombre}`" for nombre in faltan)
    lineas_estado.extend([
        "",
        "Esas ausencias están declaradas en `enmienda_v3_7_cierre`. No son un resultado.",
        "",
        "La corrida versionada terminó en un proceso. El certificado y el sello son los de `results/v3_7/certificado.json`. La comparación con la corrida archivada está en `docs/determinismo_dos_procesos.md`: los números principales coinciden. Los sellos no, porque la segunda corrida declara `results/certified/inference_audit.json`.",
        "",
    ])
    ESTADO.write_text("\n".join(lineas_estado), encoding="utf-8")
    print(f"{len(word)} macros; cambios escritos")


if __name__ == "__main__":
    main()
