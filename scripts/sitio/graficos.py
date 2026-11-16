"""Gráficos de results/v3_7. Reusa el lector del sitio v3.7 y cambia el color a la tinta."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import re
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SITIO = REPO / "sitio"


def _exportador():
    spec = importlib.util.spec_from_file_location(
        "exportar_sitio_v3_7", REPO / "scripts" / "v3" / "exportar_sitio_v3_7.py"
    )
    modulo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(modulo)
    return modulo


def _etiquetas() -> dict[str, str]:
    import yaml
    crudo = yaml.safe_load((REPO / "sitio_fuente" / "etiquetas.yaml").read_text(encoding="utf-8")) or {}
    salida = {}
    for clave, valor in crudo.items():
        salida[str(clave)] = str(valor)
        salida[str(clave).casefold()] = str(valor)
    return salida


def _tinta(texto: str) -> str:
    texto = texto.replace("#0b4f55", "currentColor").replace("#6b2d12", "currentColor")
    texto = texto.replace("595118ae09ea6006c744c86b28cb0a10d42e0442", "v3.7-final")
    texto = re.sub(r">origen [^<]+<", ">Origen en el repositorio<", texto)
    etiquetas = _etiquetas()

    def etiquetar(coincidencia):
        interior = coincidencia.group(1)

        def palabra(pieza):
            token = pieza.group(0)
            if token in etiquetas and (
                "_" in token or token in {"true", "false", "True", "False"}
            ):
                return etiquetas[token]
            return token

        return ">" + re.sub(r"\b[A-Za-z][A-Za-z0-9_]*\b", palabra, interior) + "<"

    texto = re.sub(r">([^<]*)<", etiquetar, texto)

    def envolver(coincidencia):
        plano = re.sub(r"<[^>]+>", " ", coincidencia.group(0))
        if re.search(r"results/|config/|\.tsv\b|\.json\b|\.yaml\b", plano):
            return (
                '<details class="especialistas"><summary>Para especialistas</summary>'
                + coincidencia.group(0) + "</details>"
            )
        return coincidencia.group(0)

    texto = re.sub(r"<p\b[^>]*>.*?</p>", envolver, texto, flags=re.S)
    if "<desc>" in texto and "como-leer" not in texto:
        desc = re.search(r"<desc>(.*?)</desc>", texto)
        if desc:
            bloque = (
                '<details class="como-leer"><summary>Cómo leer este gráfico</summary>'
                f"<p>{desc.group(1)}</p></details>"
            )
            texto = texto.replace("</svg>", "</svg>" + bloque, 1)
    return texto


def _cerrar_sello(fid: str) -> None:
    ruta = SITIO / "figuras" / f"{fid}.svg"
    texto = _tinta(ruta.read_text(encoding="utf-8"))
    ruta.write_text(texto, encoding="utf-8")
    sello_ruta = SITIO / "figuras" / f"{fid}.sello.json"
    sello = json.loads(sello_ruta.read_text(encoding="utf-8"))
    sello["svg_sha256"] = hashlib.sha256(texto.encode()).hexdigest()
    sello["generador"] = "scripts/sitio/graficos.py"
    sello_ruta.write_text(json.dumps(sello, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _volcan(mod, fid: str, titulo: str, tabla: pd.DataFrame, fuente: str) -> str:
    frame = tabla.copy()
    frame["log2fc"] = pd.to_numeric(frame["log2fc"], errors="coerce")
    frame["padj"] = pd.to_numeric(frame["padj"], errors="coerce")
    frame["ic95_bajo"] = pd.to_numeric(frame["ic95_bajo"], errors="coerce")
    frame["ic95_alto"] = pd.to_numeric(frame["ic95_alto"], errors="coerce")
    frame = frame.dropna(subset=["log2fc", "padj", "ic95_bajo", "ic95_alto"])
    positivo = frame.loc[frame["padj"] > 0, "padj"]
    piso = float(positivo.min()) if len(positivo) else 1e-300
    frame["menos_log10"] = [-math.log10(max(float(p), piso)) for p in frame["padj"]]
    xs = frame["ic95_bajo"].tolist() + frame["ic95_alto"].tolist()
    sx, *_ = mod.escala(xs, 50, 600)
    sy, *_ = mod.escala(frame["menos_log10"].tolist() + [0], 300, 40)
    rayas = {"diferente": "0", "equivalente": "4 2", "indeterminada": "1 2"}
    cuerpo, filas = [], []
    for _, fila in frame.iterrows():
        x, y = float(fila["log2fc"]), float(fila["menos_log10"])
        raya = rayas.get(str(fila["clase"]), "2 2")
        proteina = str(fila["protein_id"])
        cuerpo.append(
            f'<line x1="{sx(float(fila["ic95_bajo"])):.1f}" y1="{sy(y):.1f}" '
            f'x2="{sx(float(fila["ic95_alto"])):.1f}" y2="{sy(y):.1f}" stroke="currentColor" stroke-width="1"/>'
            f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3" fill="none" stroke="currentColor" '
            f'stroke-dasharray="{raya}" data-value="{format(x, ".16g")}" '
            f'data-fuente="{fuente}#protein_id={mod.e(proteina)}.log2fc" tabindex="0" data-tip="{mod.e(proteina)}"/>'
        )
        filas.append({
            "protein_id": proteina,
            "log2fc": x,
            "ic95_bajo": float(fila["ic95_bajo"]),
            "ic95_alto": float(fila["ic95_alto"]),
            "padj": float(fila["padj"]),
            "clase": fila["clase"],
        })
    return mod.figura(
        fid, titulo,
        "Cada punto es una proteína. La raya horizontal es el intervalo de 95 %. El trazo distingue la clase.",
        "".join(cuerpo), filas, [fuente],
        "¿Dónde cae cada proteína?",
    )


def _go_puntos(mod, gsea: pd.DataFrame, go_json) -> str:
    fdr = float(go_json["parametros"]["fdr"])
    cptac = gsea[gsea["cohorte"] == "cptac_tumor_frente_a_normal"]
    bajos = cptac[pd.to_numeric(cptac["FDR q-val"], errors="coerce") <= fdr]
    puntos = []
    for _, fila in bajos.iterrows():
        nes = float(fila["NES"])
        if nes == 0:
            continue
        puntos.append((str(fila["Term"]), nes, float(fila["FDR q-val"])))
    if not puntos:
        return mod.figura("r04", "Puntos de GO", "No hubo términos dentro del FDR.", "", [], ["gsea_tejido.tsv"], "¿Qué términos quedan?")
    sx, *_ = mod.escala([p[1] for p in puntos] + [0], 160, 600)
    cuerpo, filas = [], []
    alto = 24 + len(puntos) * 18
    for i, (termino, nes, fdr_valor) in enumerate(puntos):
        y = 20 + i * 18
        cuerpo.append(
            f'<circle cx="{sx(nes):.1f}" cy="{y}" r="4" fill="none" stroke="currentColor" '
            f'data-value="{format(nes, ".16g")}" '
            f'data-fuente="results/v3_7/gsea_tejido.tsv#Term={mod.e(termino)}|cohorte=cptac_tumor_frente_a_normal.NES" '
            f'tabindex="0" data-tip="{mod.e(termino)}"/>'
        )
        filas.append({"termino": termino, "nes": nes, "fdr": fdr_valor})
    return mod.figura(
        "r04", "Puntos de GO en CPTAC",
        "Cada punto es un término de CPTAC con FDR dentro del umbral del fichero y NES distinto de cero.",
        "".join(cuerpo), filas, ["gsea_tejido.tsv", "go_tejido.json"],
        "¿Hacia dónde apunta cada término?",
        viewbox=f"0 0 640 {max(360, alto)}",
    )


def _shap_calor(mod, tabla: pd.DataFrame) -> str:
    columnas = ["rho_eritrocitos", "rho_plaquetas", "rho_coagulacion", "rho_piel", "rho_n_precursores", "rho_orden"]
    bloque = tabla[tabla["modelo"] == "M2_proteinas_principal"].copy()
    cuerpo, filas, numeros = [], [], []
    for _, fila in bloque.iterrows():
        for col in columnas:
            valor = pd.to_numeric(pd.Series([fila[col]]), errors="coerce").iloc[0]
            if pd.notna(valor):
                numeros.append(abs(float(valor)))
    limite = max(numeros) if numeros else 1
    for i, (_, fila) in enumerate(bloque.iterrows()):
        gen = str(fila["gen"]) if pd.notna(fila["gen"]) else str(fila["variable"])
        cuerpo.append(f'<text x="4" y="{28 + i * 16}" fill="currentColor" font-size="10">{mod.e(gen)}</text>')
        for j, col in enumerate(columnas):
            valor = pd.to_numeric(pd.Series([fila[col]]), errors="coerce").iloc[0]
            if pd.isna(valor):
                continue
            tono = int(255 * (1 - min(1, abs(float(valor)) / limite)))
            raya = "0" if float(valor) >= 0 else "3 2"
            fuente = f"results/v3_7/ml_explicacion.tsv#modelo=M2_proteinas_principal|variable={fila['variable']}.{col}"
            cuerpo.append(
                f'<rect x="{120 + j * 70}" y="{16 + i * 16}" width="64" height="14" '
                f'fill="rgb({tono},{tono},{tono})" stroke="#4e5d6c" stroke-width="1" stroke-dasharray="{raya}" '
                f'data-value="{format(float(valor), ".16g")}" data-fuente="{mod.e(fuente)}" tabindex="0" '
                f'data-tip="{mod.e(gen)} {mod.e(col)}"/>'
            )
            filas.append({"gen": gen, "columna": col, "rho": float(valor)})
    alto = 40 + max(1, len(bloque)) * 16
    return mod.figura(
        "r05", "SHAP de M2 frente a la calidad",
        "Cada celda es la correlación de esa proteína con una variable de calidad. El trazo discontinuo es una correlación negativa.",
        "".join(cuerpo), filas, ["ml_explicacion.tsv"],
        "¿La explicación se confunde con la calidad?",
        viewbox=f"0 0 640 {max(360, alto)}",
    )


def _bosque_paso5(mod, pri: pd.DataFrame, suero: pd.DataFrame, mapa: pd.DataFrame) -> str:
    genes = pri.loc[pri["ultimo_paso"] == 5, "gen"].astype(str)
    ids = set(suero["protein_id"].astype(str))
    cruce = mapa[mapa["gene_symbol"].isin(genes) & mapa["uniprot_accession"].isin(ids)]
    conteo = cruce.groupby("gene_symbol")["uniprot_accession"].nunique()
    unicos = conteo[conteo == 1].index
    cruce = cruce[cruce["gene_symbol"].isin(unicos)].drop_duplicates("gene_symbol")
    cruce = cruce[["gene_symbol", "uniprot_accession"]]
    tabla = suero.merge(cruce, left_on="protein_id", right_on="uniprot_accession", how="inner")
    tabla["log2fc"] = pd.to_numeric(tabla["log2fc"], errors="coerce")
    tabla["ic95_bajo"] = pd.to_numeric(tabla["ic95_bajo"], errors="coerce")
    tabla["ic95_alto"] = pd.to_numeric(tabla["ic95_alto"], errors="coerce")
    tabla = tabla.dropna(subset=["log2fc", "ic95_bajo", "ic95_alto"]).sort_values("log2fc")
    xs = tabla["ic95_bajo"].tolist() + tabla["ic95_alto"].tolist() + [0]
    sx, *_ = mod.escala(xs, 120, 600)
    cuerpo, filas = [], []
    for i, (_, fila) in enumerate(tabla.iterrows()):
        y = 24 + i * 22
        gen = str(fila["gene_symbol"])
        cuerpo.append(
            f'<line x1="{sx(float(fila["ic95_bajo"])):.1f}" y1="{y}" x2="{sx(float(fila["ic95_alto"])):.1f}" y2="{y}" stroke="currentColor"/>'
            f'<circle cx="{sx(float(fila["log2fc"])):.1f}" cy="{y}" r="4" fill="none" stroke="currentColor" '
            f'data-value="{format(float(fila["log2fc"]), ".16g")}" '
            f'data-fuente="results/v3_7/suero_k4_crc_cn.tsv#protein_id={mod.e(str(fila["protein_id"]))}.log2fc" '
            f'tabindex="0" data-tip="{mod.e(gen)}"/>'
            f'<text x="8" y="{y + 4}" fill="currentColor" font-size="12">{mod.e(gen)}</text>'
        )
        filas.append({
            "gen": gen,
            "protein_id": fila["protein_id"],
            "log2fc": float(fila["log2fc"]),
            "ic95_bajo": float(fila["ic95_bajo"]),
            "ic95_alto": float(fila["ic95_alto"]),
            "clase": fila["clase"],
        })
    alto = 40 + len(tabla) * 22
    return mod.figura(
        "r06", "Las 20 del paso 5 en CP frente a CN",
        "Cada punto es el log2FC del suero. La raya es el intervalo de 95 %. El gen se une al acceso solo cuando el mapa revisado da uno.",
        "".join(cuerpo), filas,
        ["suero_k4_crc_cn.tsv", "priorizacion.tsv", "metadata/protein_id_maps/gene_uniprot_human_reviewed_bulk.tsv"],
        "¿Qué efecto tienen en el contraste que decide?",
        viewbox=f"0 0 640 {max(360, alto)}",
    )


def _tablero(mod, k1, h5, h6, h7, h8, h9, h10, pre) -> str:
    umbral = mod.caminar(pre, "hipotesis.h4.umbral")
    umbral_visible = str(umbral).replace(".", ",")
    filas_src = [
        ("H1", "No se sostiene", k1["h1_veredicto"], "results/v3_7/k1_etiquetas.json#h1_veredicto"),
        ("H2", "No evaluada", mod.caminar(pre, "hipotesis.h2.estado_heredado"), "config/preregistro_v3.yaml#hipotesis.h2.estado_heredado"),
        ("H3", "Adoptada, enmienda de cierre v3.7", "adoptada", "config/preregistro_v3.yaml#enmienda_v3_7_cierre.fecha"),
        ("H4", f"Umbral {umbral_visible} · sin veredicto", umbral, "config/preregistro_v3.yaml#hipotesis.h4.umbral"),
        ("H5", "Se sostiene", h5["estado"], "results/v3_7/h5_especificidad_organo.json#estado"),
        ("H6", "No se sostiene", h6["estado"], "results/v3_7/h6_origen_celular.json#estado"),
        ("H7", "Refutada", h7["estado"], "results/v3_7/h7_visibilidad_en_sangre.json#estado"),
        ("H8", "No se pudo evaluar", h8["estado"], "results/v3_7/h8_carga_tumoral.json#estado"),
        ("H9", "Exploratoria", h9["estado"], "results/v3_7/h9_wald.json#estado"),
        ("H10", "Descriptiva", h10["estado"], "results/v3_7/h10_especificidad.json#estado"),
    ]
    cuerpo, filas = [], []
    for i, (nombre, visible, bruto, fuente) in enumerate(filas_src):
        col, ren = i % 5, i // 5
        valor = format(bruto, ".16g") if isinstance(bruto, float) else str(bruto)
        cuerpo.append(
            f'<rect x="{20 + col * 122}" y="{40 + ren * 120}" width="112" height="90" fill="none" stroke="currentColor" '
            f'stroke-width="2" data-value="{mod.e(valor)}" data-fuente="{mod.e(fuente)}" tabindex="0" '
            f'data-tip="{mod.e(nombre)} {mod.e(visible)}"/>'
            f'<text x="{28 + col * 122}" y="{72 + ren * 120}" fill="currentColor">{mod.e(nombre)}</text>'
            f'<text x="{28 + col * 122}" y="{96 + ren * 120}" fill="currentColor" font-size="11">{mod.e(visible[:22])}</text>'
        )
        filas.append({"hipotesis": nombre, "valor": visible, "fuente": fuente})
    return mod.figura(
        "n11", "Tablero de H1 a H10",
        "Cada celda muestra el estado en lenguaje llano. H3 es la priorización sin pesos, adoptada en la enmienda de cierre v3.7. H4 muestra el umbral, no un veredicto.",
        "".join(cuerpo), filas,
        ["k1_etiquetas.json", "h5_especificidad_organo.json", "h6_origen_celular.json",
         "h7_visibilidad_en_sangre.json", "h8_carga_tumoral.json", "h9_wald.json", "h10_especificidad.json"],
        "¿Cuál es el estado de cada hipótesis?",
    )


def generar() -> dict[str, str]:
    mod = _exportador()
    v37 = REPO / "results" / "v3_7"
    cargar = mod.cargar
    tsv = mod.tabla_tsv
    emb = cargar("suero_embudo.json")
    k1 = cargar("k1_etiquetas.json")
    go_json = cargar("go_tejido.json")
    ml = cargar("ml_suero_permutacion.json")
    k5 = cargar("k5_auc.json")
    pre = yaml_preregistro()
    pri = cargar("priorizacion_embudo.json")
    cert = cargar("certificado.json")
    crc_cn = tsv("suero_k4_crc_cn.tsv")
    crc_hc = tsv("suero_k4_crc_sanos.tsv")
    cn_hc = tsv("suero_k4_cn_sanos.tsv")
    corr = pd.read_csv(REPO / "config" / "correspondencia_propuesta.tsv", sep="\t")
    estado = {fila.punto: fila.estado for fila in corr.itertuples(index=False)}
    piezas = {
        "f01": mod.mapa(),
        "f02": mod.embudo(emb),
        "f03": mod.celda_svg(emb, crc_hc, cn_hc),
        "f04": mod.multiverso(cargar("multiverso_resumen.json"), tsv("multiverso_suero.tsv")),
        "f05": mod.r1(k1, tsv("h7_tabla_proteinas.tsv"), tsv("triangulacion.tsv")),
        "f06": mod.controles(tsv("controles.tsv")),
        "f07": mod.go(go_json, tsv("gsea_tejido.tsv")),
        "f08": mod.bosque(ml, k5),
        "f09": mod.shap(tsv("ml_explicacion.tsv")),
        "n10": mod.muestras(tsv("preanalitica_orden.tsv"), tsv("suero_calidad_muestra.tsv")),
        "n11": _tablero(
            mod, k1, cargar("h5_especificidad_organo.json"), cargar("h6_origen_celular.json"),
            cargar("h7_visibilidad_en_sangre.json"), cargar("h8_carga_tumoral.json"),
            cargar("h9_wald.json"), cargar("h10_especificidad.json"), pre,
        ),
        "n12": mod.embudo_priori(pri),
        "n13": mod.calor(tsv("priorizacion.tsv"), tsv("triangulacion.tsv")),
        "n14": mod.certificado_svg(cert),
        "p01": mod.mini("p01", "celda", "Proteínas de la celda", emb["celda_preregistrada"]["n"], "entero", "results/v3_7/suero_embudo.json#celda_preregistrada.n"),
        "p02": mod.mini("p02", "diferente", "Diferentes frente a CN", emb["celda_preregistrada"]["diferente"], "entero", "results/v3_7/suero_embudo.json#celda_preregistrada.diferente"),
        "p03": mod.mini("p03", "cobertura", "Cobertura de R1", k1["cobertura"]["numerador"], "entero", "results/v3_7/k1_etiquetas.json#cobertura.numerador"),
        "p04": mod.mini("p04", "paso5", "Paso 5", pri["pasos"][4]["n"], "entero", "results/v3_7/priorizacion_embudo.json#pasos.4.n"),
        "p05": mod.mini("p05", "nes", "Concordancia de NES", go_json["concordancia_nes"]["numerador"], "entero", "results/v3_7/go_tejido.json#concordancia_nes.numerador"),
        "p06": mod.mini("p06", "redes", "Redes", estado["redes"], "texto", "config/correspondencia_propuesta.tsv#punto=redes.estado"),
        "p07": mod.mini("p07", "m3", "AUC de M3", ml["versiones"]["M3_solo_qc"]["auc_oof_media"], "auc", "results/v3_7/ml_suero_permutacion.json#versiones.M3_solo_qc.auc_oof_media"),
        "p08": mod.mini("p08", "paso6", "Paso 6", pri["pasos"][5]["n"], "entero", "results/v3_7/priorizacion_embudo.json#pasos.5.n"),
        "p09": mod.mini("p09", "cohortes", "Cohortes nuevas", estado["ampliar_cohortes"], "texto", "config/correspondencia_propuesta.tsv#punto=ampliar_cohortes.estado"),
        "p10": mod.mini("p10", "certificado", "Certificado", cert["estado"], "estado", "results/v3_7/certificado.json#estado"),
        "r01": _volcan(mod, "r01", "Volcán CP frente a CN", crc_cn, "results/v3_7/suero_k4_crc_cn.tsv"),
        "r02": _volcan(mod, "r02", "Volcán CP frente a HC", crc_hc, "results/v3_7/suero_k4_crc_sanos.tsv"),
        "r03": _volcan(mod, "r03", "Volcán CN frente a HC", cn_hc, "results/v3_7/suero_k4_cn_sanos.tsv"),
        "r04": _go_puntos(mod, tsv("gsea_tejido.tsv"), go_json),
        "r05": _shap_calor(mod, tsv("ml_explicacion.tsv")),
        "r06": _bosque_paso5(
            mod, tsv("priorizacion.tsv"), crc_cn,
            pd.read_csv(REPO / "metadata" / "protein_id_maps" / "gene_uniprot_human_reviewed_bulk.tsv", sep="\t"),
        ),
    }
    for fid in piezas:
        _cerrar_sello(fid)
        piezas[fid] = _tinta(piezas[fid])
    return piezas


def yaml_preregistro():
    import yaml
    return yaml.safe_load((REPO / "config" / "preregistro_v3.yaml").read_text(encoding="utf-8"))
