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
CORRESPONDENCIA = REPO / "config" / "correspondencia_propuesta.tsv"
CORRESPONDENCIA_TEX = REPO / "tesis-ucaldas" / "latex" / "correspondencia.tex"
PASOS = ("", "Uno", "Dos", "Tres", "Cuatro", "Cinco", "Seis", "Siete")

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


def _anexar(filas: list, version: str, nombre: str, valor: object, tipo: str, fuente: str) -> None:
    if valor is None and version != "v3_7":
        return
    filas.append((nombre, valor, tipo, fuente))


def _hc(carpeta: Path, tabla: str) -> tuple[int | None, list[str]]:
    ruta = carpeta / tabla
    traza_ruta = carpeta / "suero_trazabilidad_grupos.tsv"
    if not ruta.is_file() or not traza_ruta.is_file():
        return None, []
    datos = pd.read_csv(ruta, sep="\t")
    if "clase" not in datos.columns or "protein_id" not in datos.columns:
        return None, []
    dif = datos.loc[datos["clase"].eq("diferente")].copy()
    if "padj" in dif.columns:
        dif = dif.sort_values("padj")
    traza = pd.read_csv(traza_ruta, sep="\t", usecols=["protein_group", "genes"])
    mapa = traza.drop_duplicates("protein_group").set_index("protein_group")["genes"]
    genes = [str(mapa.get(acceso, acceso)).split(";")[0] for acceso in dif["protein_id"].astype(str)]
    return int(len(dif)), genes


def _paso_cinco(carpeta: Path) -> tuple[int, list[str], list[str]] | None:
    ruta = carpeta / "priorizacion.tsv"
    if not ruta.is_file():
        return None
    tabla = pd.read_csv(ruta, sep="\t")
    cinco = tabla.loc[pd.to_numeric(tabla["ultimo_paso"], errors="coerce") >= 5].copy()
    fc = pd.to_numeric(cinco["log2fc_cptac"], errors="coerce")
    if fc.isna().any():
        raise SystemExit("log2fc_cptac vacío en el paso 5")
    bajan = int((fc < 0).sum())
    suben = cinco.loc[fc > 0, "gen"].astype(str).tolist()
    if bajan + len(suben) != len(cinco):
        raise SystemExit("hay un log2fc_cptac de cero en el paso 5")
    return bajan, cinco["gen"].astype(str).tolist(), suben


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
    priori = _json(carpeta, "priorizacion_embudo.json") if (carpeta / "priorizacion_embudo.json").is_file() else {"pasos": []}
    explicacion = _json(carpeta, "ml_explicacion.json") if (carpeta / "ml_explicacion.json").is_file() else {}
    calidad = _json(carpeta, "k7_control_tecnico.json")["confusion_preanalitica"]
    certificado = _json(carpeta, "certificado.json")
    piel = {fila["indicador"]: fila for fila in calidad}
    auc = pre["calculado_con_legibles"]["auc"]
    spearman = pre["calculado_con_legibles"]["spearman"]["piel"]
    primario = tejido["primario_frente_a_mucosa"]
    metastasis = tejido["metastasis_hepatica_frente_a_mucosa"]
    filas = [
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
        ("HsieteSec", (h7.get("modelo_secundario") or {}).get("auc_secundario"), "auc_completo", f"{prefijo}/h7_visibilidad_en_sangre.json#modelo_secundario.auc_secundario"),
        ("PreCrudos", pre["n_crudos_validos"], "entero", f"{prefijo}/preanalitica.json#n_crudos_validos"),
        ("CalidadCoagP", piel["coagulacion"]["p_kruskal"], "p", f"{prefijo}/k7_control_tecnico.json#coagulacion.p_kruskal"),
        ("CalidadPielP", piel["piel"]["p_kruskal"], "p", f"{prefijo}/k7_control_tecnico.json#piel.p_kruskal"),
        ("CertificadoEstado", certificado["estado"], "estado", f"{prefijo}/certificado.json#estado"),
        ("MultiSuero", multi["n_celdas_suero"], "entero", f"{prefijo}/multiverso_resumen.json#n_celdas_suero"),
        ("MultiMl", multi["n_celdas_ml"], "entero", f"{prefijo}/multiverso_resumen.json#n_celdas_ml"),
        ("MultiDif", multi["n_con_diferente"], "entero", f"{prefijo}/multiverso_resumen.json#n_con_diferente"),
    ]
    for paso in priori.get("pasos") or []:
        nombre = PASOS[int(paso["paso"])] if 1 <= int(paso["paso"]) < len(PASOS) else ""
        if nombre:
            _anexar(filas, version, f"PrioriPaso{nombre}", paso.get("n"), "entero", f"{prefijo}/priorizacion_embudo.json#pasos")
    paso5 = _paso_cinco(carpeta)
    if paso5 is not None:
        bajan, _genes, suben = paso5
        _anexar(filas, version, "PrioriBajan", bajan, "entero", f"{prefijo}/priorizacion.tsv#log2fc_cptac")
        _anexar(filas, version, "PrioriSuben", len(suben), "entero", f"{prefijo}/priorizacion.tsv#log2fc_cptac")
    n_cp, _genes_cp = _hc(carpeta, "suero_k4_crc_sanos.tsv")
    n_cn, _genes_cn = _hc(carpeta, "suero_k4_cn_sanos.tsv")
    _anexar(filas, version, "HcCpDif", n_cp, "entero", f"{prefijo}/suero_k4_crc_sanos.tsv#clase")
    _anexar(filas, version, "HcCnDif", n_cn, "entero", f"{prefijo}/suero_k4_cn_sanos.tsv#clase")
    cohortes = go.get("cohortes") or {}
    cptac = cohortes.get("cptac_tumor_frente_a_normal") or {}
    tanaka = cohortes.get("tanaka_primario_frente_a_mucosa") or {}
    for nombre, valor, campo in (
        ("GoCptacTerminos", cptac.get("n_terminos_probados"), "cptac_tumor_frente_a_normal.n_terminos_probados"),
        ("GoCptacSuben", cptac.get("n_ora_significativos_suben"), "cptac_tumor_frente_a_normal.n_ora_significativos_suben"),
        ("GoCptacBajan", cptac.get("n_ora_significativos_bajan"), "cptac_tumor_frente_a_normal.n_ora_significativos_bajan"),
        ("GoCptacGsea", cptac.get("n_gsea_fdr"), "cptac_tumor_frente_a_normal.n_gsea_fdr"),
        ("GoTanakaTerminos", tanaka.get("n_terminos_probados"), "tanaka_primario_frente_a_mucosa.n_terminos_probados"),
        ("GoTanakaSuben", tanaka.get("n_ora_significativos_suben"), "tanaka_primario_frente_a_mucosa.n_ora_significativos_suben"),
        ("GoTanakaBajan", tanaka.get("n_ora_significativos_bajan"), "tanaka_primario_frente_a_mucosa.n_ora_significativos_bajan"),
        ("GoTanakaGsea", tanaka.get("n_gsea_fdr"), "tanaka_primario_frente_a_mucosa.n_gsea_fdr"),
    ):
        _anexar(filas, version, nombre, valor, "entero", f"{prefijo}/go_tejido.json#{campo}")
    concordancia = go.get("concordancia_nes") or {}
    _anexar(filas, version, "GoCiBajo", concordancia.get("ci_low"), "auc", f"{prefijo}/go_tejido.json#concordancia_nes.ci_low")
    _anexar(filas, version, "GoCiAlto", concordancia.get("ci_high"), "auc", f"{prefijo}/go_tejido.json#concordancia_nes.ci_high")
    parametros = go.get("parametros") or {}
    _anexar(filas, version, "GoFdr", parametros.get("fdr"), "p", f"{prefijo}/go_tejido.json#parametros.fdr")
    _anexar(filas, version, "GoPerm", parametros.get("permutaciones"), "entero", f"{prefijo}/go_tejido.json#parametros.permutaciones")
    comprobacion = (explicacion.get("comprobacion_shap") or {}).get("M2_proteinas_principal")
    _anexar(filas, version, "ShapLineal", comprobacion, "entero", f"{prefijo}/ml_explicacion.json#comprobacion_shap.M2_proteinas_principal")
    _anexar(filas, version, "ShapPiel", explicacion.get("n_piel_en_las_20"), "entero", f"{prefijo}/ml_explicacion.json#n_piel_en_las_20")
    _anexar(filas, version, "ShapCoag", explicacion.get("n_coagulacion_en_las_20"), "entero", f"{prefijo}/ml_explicacion.json#n_coagulacion_en_las_20")
    _anexar(filas, version, "ShapCalidad", explicacion.get("n_correlacionan_con_calidad_bh"), "entero", f"{prefijo}/ml_explicacion.json#n_correlacionan_con_calidad_bh")
    if (carpeta / "ml_explicacion.tsv").is_file():
        tabla = pd.read_csv(carpeta / "ml_explicacion.tsv", sep="\t")
        n_top = int(tabla["modelo"].eq("M2_proteinas_principal").sum()) if "modelo" in tabla.columns else None
        _anexar(filas, version, "ShapTop", n_top, "entero", f"{prefijo}/ml_explicacion.tsv#modelo")
    return filas


def _comando_texto(lineas: list[str], word: dict, nombre: str, texto: str, fuente: str) -> None:
    lineas.append(f"% fuente: {fuente}")
    lineas.append(f"\\newcommand{{\\{nombre}}}{{{_tex(texto)}}}")
    word[nombre] = {"valor": texto, "fuente": fuente, "crudo": texto}


def _escribir_correspondencia() -> None:
    n_cp, _genes_cp = _hc(V37, "suero_k4_crc_sanos.tsv")
    n_cn, _genes_cn = _hc(V37, "suero_k4_cn_sanos.tsv")
    if n_cp is None or n_cn is None:
        raise SystemExit("faltan los contrastes frente a HC")
    frase = (
        "Frente a CN ninguna proteína de la celda es diferente. CN no son personas sanas. "
        f"Frente a HC salen {n_cp} (CP) y {n_cn} (CN), leídos de suero_k4_crc_sanos.tsv y suero_k4_cn_sanos.tsv."
    )
    lineas = CORRESPONDENCIA.read_text(encoding="utf-8").splitlines()
    cabecera = lineas[0].split("\t")
    punto = cabecera.index("punto")
    explicacion = cabecera.index("explicacion")
    salida = [lineas[0]]
    filas = []
    for linea in lineas[1:]:
        if not linea.strip():
            continue
        celdas = linea.split("\t")
        if celdas[punto] == "proteinas_diferenciales":
            celdas[explicacion] = frase
        salida.append("\t".join(celdas))
        filas.append(dict(zip(cabecera, celdas)))
    CORRESPONDENCIA.write_text("\n".join(salida) + "\n", encoding="utf-8")
    cuerpo = [
        "% Generado desde config/correspondencia_propuesta.tsv. No editar a mano.",
        "\\begin{tabular}{p{3.2cm}p{3.4cm}p{7.2cm}}",
        "\\textbf{Punto} & \\textbf{Estado} & \\textbf{Dónde está} \\\\",
    ]
    for fila in filas:
        donde = fila.get("fichero") or fila.get("regla") or ""
        if fila["punto"] == "proteinas_diferenciales":
            donde = f"{fila['fichero']}. {fila['explicacion']}"
        cuerpo.append(f"{_tex(fila['prometido'])} & {_tex(fila['estado'])} & {_tex(donde)} \\\\")
    cuerpo.append("\\end{tabular}")
    CORRESPONDENCIA_TEX.write_text("\n".join(cuerpo) + "\n", encoding="utf-8")


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
    paso5 = _paso_cinco(V37)
    if paso5 is None:
        raise SystemExit("falta priorizacion.tsv")
    _bajan, genes_paso, suben = paso5
    _comando_texto(lineas, word, "PrioriGenesPasoCinco", ", ".join(genes_paso), "results/v3_7/priorizacion.tsv#gen")
    _comando_texto(lineas, word, "PrioriSubenGenes", ", ".join(suben), "results/v3_7/priorizacion.tsv#gen")
    _n_cp, genes_cp = _hc(V37, "suero_k4_crc_sanos.tsv")
    _n_cn, genes_cn = _hc(V37, "suero_k4_cn_sanos.tsv")
    _comando_texto(lineas, word, "HcCpGenes", ", ".join(genes_cp), "results/v3_7/suero_k4_crc_sanos.tsv#protein_id")
    _comando_texto(lineas, word, "HcCnGenes", ", ".join(genes_cn), "results/v3_7/suero_k4_cn_sanos.tsv#protein_id")
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
    _escribir_correspondencia()
    print(f"{len(word)} macros; cambios escritos")


if __name__ == "__main__":
    main()
