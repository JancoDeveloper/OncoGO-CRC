"""Figuras del sitio, leídas de la corrida sellada.

No modifica results/. Escribe SVG en web/figuras/. Las etiquetas usan coma
decimal y tres decimales cuando redondean un AUC; el valor original queda
en el elemento desc.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
V37 = REPO / "results" / "v3_7"
FUENTES = REPO / "docs" / "fuentes_de_datos.md"
META_CPTAC = REPO / "metadata" / "sample_metadata" / "cptac_metadata.tsv"

FUENTE = "Segoe UI, Roboto, Helvetica, Arial, sans-serif"
TINTA = "#1a1a1a"
MUTED = "#3d3a37"
PAPEL = "#f7f6f3"
LINEA = "#d4d0cb"
AZUL = "#0072B2"
CIELO = "#56B4E9"
VERDE = "#009E73"
NARANJA = "#E69F00"
VERMILLON = "#D55E00"
MORADO = "#CC79A7"

NOMBRES_MODELO = {
    "M1_proteinas_filtrada_humana": "M1  Proteínas humanas filtradas",
    "M2_proteinas_principal": "M2  Proteínas de la matriz principal",
    "M3_solo_qc": "M3  Solo calidad de la muestra",
    "M4_solo_coagulacion": "M4  Solo coagulación",
    "M5_principal_residualizada_qc": "M5  Proteínas menos la calidad",
    "M6_principal_pca5": "M6  Proteínas en 5 componentes",
}
COLOR_MODELO = {
    "M1_proteinas_filtrada_humana": AZUL,
    "M2_proteinas_principal": CIELO,
    "M3_solo_qc": VERMILLON,
    "M4_solo_coagulacion": NARANJA,
    "M5_principal_residualizada_qc": VERDE,
    "M6_principal_pca5": MORADO,
}


def esc(texto: object) -> str:
    return (
        str(texto)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def coma_tres(texto: str) -> str:
    """Tres decimales con coma, a partir del texto de la tabla."""
    cuantizado = Decimal(texto).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)
    return format(cuantizado, "f").replace(".", ",")


def coma_literal(texto: str) -> str:
    return texto.replace(".", ",")


def lienzo(ancho: int, alto: int, titulo: str, desc: str, cuerpo: str, pie: str) -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {ancho} {alto}" '
        'role="img" aria-labelledby="titulo desc">\n'
        f"  <title id=\"titulo\">{esc(titulo)}</title>\n"
        f"  <desc id=\"desc\">{esc(desc)}</desc>\n"
        f'  <rect width="{ancho}" height="{alto}" fill="{PAPEL}"/>\n'
        f'  <text x="32" y="42" font-family="{FUENTE}" font-size="26" '
        f'font-weight="700" fill="{TINTA}">{esc(titulo)}</text>\n'
        f"{cuerpo}\n"
        f'  <text x="32" y="{alto - 24}" font-family="{FUENTE}" font-size="15" '
        f'fill="{MUTED}">{esc(pie)}</text>\n'
        "</svg>\n"
    )


def leer_json(nombre: str) -> dict:
    return json.loads((V37 / nombre).read_text(encoding="utf-8"))


def leer_modelos() -> list[dict[str, str]]:
    with (V37 / "ml_modelos.tsv").open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def contar_cptac() -> dict[str, int]:
    with META_CPTAC.open(encoding="utf-8", newline="") as fh:
        filas = list(csv.DictReader(fh, delimiter="\t"))
    condiciones = Counter(fila["condition"] for fila in filas)
    participantes = {fila["participant_id"] for fila in filas if fila["participant_id"]}
    return {
        "filas": len(filas),
        "tumor": condiciones.get("tumor", 0),
        "normal": condiciones.get("normal", 0),
        "unknown": condiciones.get("unknown", 0),
        "otras": sum(n for clave, n in condiciones.items() if clave not in {"tumor", "normal", "unknown"}),
        "participantes": len(participantes),
        "sin_participante": sum(1 for fila in filas if not fila["participant_id"]),
    }


def matriz_cptac_documentada() -> tuple[int, int, int]:
    texto = FUENTES.read_text(encoding="utf-8")
    hallado = re.search(r"(\d+) tumor, (\d+) normal, (\d+) sin tipo", texto)
    if hallado is None:
        raise SystemExit("docs/fuentes_de_datos.md no trae el conteo de la matriz CPTAC.")
    return (int(hallado.group(1)), int(hallado.group(2)), int(hallado.group(3)))


def figura_grupos_suero(pre: dict) -> str:
    grupos = pre["corridas_por_grupo"]
    orden = [("CP", "Cáncer", AZUL), ("CN", "Colonoscopia negativa", VERMILLON), ("HC", "Personas sanas", VERDE)]
    cuerpo = []
    y = 100
    for clave, etiqueta, color in orden:
        n = int(grupos[clave])
        ancho = n * 18
        cuerpo.append(
            f'  <rect x="300" y="{y}" width="{ancho}" height="48" fill="{color}"/>'
        )
        cuerpo.append(
            f'  <text x="280" y="{y + 32}" text-anchor="end" font-family="{FUENTE}" '
            f'font-size="20" fill="{TINTA}">{esc(etiqueta)}</text>'
        )
        cuerpo.append(
            f'  <text x="{320 + ancho}" y="{y + 32}" font-family="{FUENTE}" '
            f'font-size="22" font-weight="700" fill="{TINTA}">{n}</text>'
        )
        y += 78
    desc = (
        "PXD055218. "
        + ", ".join(f"{clave}={grupos[clave]}" for clave, _, _ in orden)
        + ". Fuente: results/v3_7/preanalitica.json, corridas_por_grupo."
    )
    return lienzo(
        820,
        400,
        "Suero: tres grupos de 20",
        desc,
        "\n".join(cuerpo),
        "Cada barra es un grupo. Fuente: results/v3_7/preanalitica.json.",
    )


def figura_estadios(h8: dict) -> str:
    conteo = h8["n_por_grupo"]
    orden = [("I_II", "Estadios I y II", AZUL), ("III", "Estadio III", NARANJA), ("IV", "Estadio IV", VERMILLON)]
    cuerpo = []
    y = 100
    for clave, etiqueta, color in orden:
        n = int(conteo[clave])
        cuerpo.append(f'  <rect x="280" y="{y}" width="{n * 8}" height="48" fill="{color}"/>')
        cuerpo.append(
            f'  <text x="260" y="{y + 32}" text-anchor="end" font-family="{FUENTE}" '
            f'font-size="20" fill="{TINTA}">{esc(etiqueta)}</text>'
        )
        cuerpo.append(
            f'  <text x="{296 + n * 8}" y="{y + 32}" font-family="{FUENTE}" '
            f'font-size="22" font-weight="700" fill="{TINTA}">{n}</text>'
        )
        y += 78
    desc = "PXD064655 n_por_grupo " + json.dumps(conteo, ensure_ascii=False)
    return lienzo(
        860,
        400,
        "Suero por estadio: 100 muestras",
        desc,
        "\n".join(cuerpo),
        "Fuente: results/v3_7/h8_carga_tumoral.json, campo n_por_grupo.",
    )


def figura_concordancia(k1: dict) -> str:
    pares = [
        ("primario_frente_a_mucosa", "Tumor primario y mucosa"),
        ("metastasis_hepatica_frente_a_mucosa", "Metástasis hepática y mucosa"),
    ]
    cuerpo = [
        f'  <text x="32" y="84" font-family="{FUENTE}" font-size="18" fill="{MUTED}">'
        "Proteínas con el mismo sentido del cambio, entre las comparables.</text>"
    ]
    y = 120
    for clave, etiqueta in pares:
        bloque = k1["contrastes"][clave]
        num = int(bloque["numerador"])
        den = int(bloque["denominador"])
        ancho = int(round(520 * num / den))
        cuerpo.append(f'  <text x="32" y="{y}" font-family="{FUENTE}" font-size="20" fill="{TINTA}">{esc(etiqueta)}</text>')
        cuerpo.append(f'  <rect x="32" y="{y + 14}" width="520" height="36" fill="#fff" stroke="{LINEA}"/>')
        cuerpo.append(f'  <rect x="32" y="{y + 14}" width="{ancho}" height="36" fill="{AZUL}"/>')
        cuerpo.append(
            f'  <text x="568" y="{y + 40}" font-family="{FUENTE}" font-size="22" '
            f'font-weight="700" fill="{TINTA}">{num} de {den}</text>'
        )
        y += 110
    cuerpo.append(
        f'  <text x="32" y="{y + 10}" font-family="{FUENTE}" font-size="18" fill="{MUTED}">'
        f'Cobertura del ancla: {k1["cobertura"]["numerador"]} de '
        f'{k1["cobertura"]["denominador"]} proteínas con cambio en CPTAC.</text>'
    )
    desc = (
        "primario "
        f'{k1["contrastes"]["primario_frente_a_mucosa"]["numerador"]}/'
        f'{k1["contrastes"]["primario_frente_a_mucosa"]["denominador"]}; '
        "metastasis "
        f'{k1["contrastes"]["metastasis_hepatica_frente_a_mucosa"]["numerador"]}/'
        f'{k1["contrastes"]["metastasis_hepatica_frente_a_mucosa"]["denominador"]}; '
        f'cobertura {k1["cobertura"]["numerador"]}/{k1["cobertura"]["denominador"]}'
    )
    return lienzo(
        860,
        460,
        "El cambio en el tejido coincide",
        desc,
        "\n".join(cuerpo),
        "Fuente: results/v3_7/k1_etiquetas.json.",
    )


def figura_clases(embudo: dict) -> str:
    celda = embudo["celda_preregistrada"]
    partes = [
        ("diferente", "Diferentes", VERMILLON),
        ("equivalente", "Equivalentes", VERDE),
        ("indeterminada", "Sin clasificar", NARANJA),
    ]
    total = int(celda["n"])
    cuerpo = [
        f'  <text x="32" y="86" font-family="{FUENTE}" font-size="18" fill="{MUTED}">'
        f"Contraste principal del suero: cáncer frente a colonoscopia negativa. n = {total}.</text>"
    ]
    x = 32
    y = 140
    for clave, etiqueta, color in partes:
        n = int(celda[clave])
        if n == 0:
            continue
        ancho = int(round(860 * n / total))
        cuerpo.append(f'  <rect x="{x}" y="{y}" width="{ancho}" height="64" fill="{color}"/>')
        x += ancho
    y_texto = 250
    for clave, etiqueta, color in partes:
        n = int(celda[clave])
        cuerpo.append(f'  <rect x="32" y="{y_texto}" width="22" height="22" fill="{color}"/>')
        cuerpo.append(
            f'  <text x="64" y="{y_texto + 18}" font-family="{FUENTE}" font-size="22" fill="{TINTA}">'
            f"{esc(etiqueta)}: {n}</text>"
        )
        y_texto += 40
    desc = "celda_preregistrada " + json.dumps(celda, ensure_ascii=False)
    return lienzo(
        960,
        430,
        "Suero: ninguna proteína diferente",
        desc,
        "\n".join(cuerpo),
        "Margen de 1,5 veces. Fuente: results/v3_7/suero_embudo.json.",
    )


def figura_auc(modelos: list[dict[str, str]]) -> str:
    cuerpo = [
        f'  <text x="32" y="80" font-family="{FUENTE}" font-size="18" fill="{MUTED}">'
        "AUC interno, fuera del pliegue. La línea marca el azar (0,50). El bigote es el intervalo de LeDell.</text>"
    ]
    plot_x = 430
    plot_w = 420
    y = 110

    def px(valor: float) -> int:
        return plot_x + int(round(plot_w * valor))

    cuerpo.append(f'  <line x1="{px(0.5)}" y1="100" x2="{px(0.5)}" y2="520" stroke="{TINTA}" stroke-dasharray="4 4"/>')
    cuerpo.append(
        f'  <text x="{px(0.5)}" y="98" text-anchor="middle" font-family="{FUENTE}" '
        f'font-size="14" fill="{MUTED}">azar</text>'
    )
    exactos = []
    for fila in modelos:
        nombre = NOMBRES_MODELO[fila["modelo"]]
        color = COLOR_MODELO[fila["modelo"]]
        auc = float(fila["auc_oof_media"])
        bajo = float(fila["ic95_ledell_bajo"])
        alto = float(fila["ic95_ledell_alto"])
        etiqueta = coma_tres(fila["auc_oof_media"])
        exactos.append(f'{fila["modelo"]}={fila["auc_oof_media"]};p={fila["p_permutacion"]}')
        cuerpo.append(
            f'  <text x="24" y="{y + 28}" font-family="{FUENTE}" font-size="18" fill="{TINTA}">{esc(nombre)}</text>'
        )
        cuerpo.append(f'  <rect x="{plot_x}" y="{y + 8}" width="{max(px(auc) - plot_x, 1)}" height="22" fill="{color}"/>')
        cuerpo.append(f'  <line x1="{px(bajo)}" y1="{y + 19}" x2="{px(alto)}" y2="{y + 19}" stroke="#ffffff" stroke-width="6"/>')
        cuerpo.append(f'  <line x1="{px(bajo)}" y1="{y + 19}" x2="{px(alto)}" y2="{y + 19}" stroke="{TINTA}" stroke-width="2"/>')
        cuerpo.append(f'  <line x1="{px(bajo)}" y1="{y + 11}" x2="{px(bajo)}" y2="{y + 27}" stroke="{TINTA}" stroke-width="2"/>')
        cuerpo.append(f'  <line x1="{px(alto)}" y1="{y + 11}" x2="{px(alto)}" y2="{y + 27}" stroke="{TINTA}" stroke-width="2"/>')
        cuerpo.append(
            f'  <text x="{plot_x + plot_w + 16}" y="{y + 28}" font-family="{FUENTE}" '
            f'font-size="20" font-weight="700" fill="{TINTA}">{etiqueta}</text>'
        )
        y += 64
    desc = "auc_oof_media exacto. " + " | ".join(exactos)
    return lienzo(
        980,
        580,
        "Capacidad de separar los grupos",
        desc,
        "\n".join(cuerpo),
        "Tres decimales del AUC. Fuente: results/v3_7/ml_modelos.tsv.",
    )


def figura_orden(pre: dict) -> str:
    grupos = pre["posicion_por_grupo"]
    colores = {"CP": AZUL, "CN": VERMILLON, "HC": VERDE}
    etiquetas = {"CP": "Cáncer", "CN": "Colonoscopia negativa", "HC": "Personas sanas"}
    n = int(pre["n_con_fecha"])
    cuerpo = [
        f'  <text x="32" y="84" font-family="{FUENTE}" font-size="18" fill="{MUTED}">'
        f"{n} mediciones con fecha, en el orden en que salieron del instrumento.</text>"
    ]
    ancho = 14
    x0 = 32
    for i in range(1, n + 1):
        color = "#ffffff"
        for clave, info in grupos.items():
            if int(info["orden_min"]) <= i <= int(info["orden_max"]):
                color = colores[clave]
                break
        cuerpo.append(f'  <rect x="{x0 + (i - 1) * ancho}" y="120" width="{ancho - 2}" height="48" fill="{color}"/>')
    leyenda_y = 210
    for clave in ("CP", "CN", "HC"):
        info = grupos[clave]
        cuerpo.append(f'  <rect x="32" y="{leyenda_y}" width="22" height="22" fill="{colores[clave]}"/>')
        cuerpo.append(
            f'  <text x="64" y="{leyenda_y + 18}" font-family="{FUENTE}" font-size="20" fill="{TINTA}">'
            f'{esc(etiquetas[clave])}: posiciones {info["orden_min"]} a {info["orden_max"]} '
            f'({info["n"]})</text>'
        )
        leyenda_y += 36
    desc = "bloques_perfectos=" + str(pre["bloques_perfectos"]).lower() + " " + json.dumps(grupos, ensure_ascii=False)
    return lienzo(
        920,
        400,
        "El orden de medición coincide con el grupo",
        desc,
        "\n".join(cuerpo),
        "Fuente: results/v3_7/preanalitica.json, posicion_por_grupo.",
    )


def figura_multiverso(resumen: dict) -> str:
    total = int(resumen["n_celdas_suero"])
    con = int(resumen["n_con_diferente"])
    sin = total - con
    ancho_sin = int(round(800 * sin / total))
    ancho_con = 800 - ancho_sin
    cuerpo = [
        f'  <text x="32" y="88" font-family="{FUENTE}" font-size="18" fill="{MUTED}">'
        "Combinaciones de decisiones de análisis en el suero.</text>",
        f'  <rect x="32" y="130" width="{ancho_sin}" height="64" fill="{VERDE}"/>',
        f'  <rect x="{32 + ancho_sin}" y="130" width="{ancho_con}" height="64" fill="{VERMILLON}"/>',
        f'  <text x="32" y="240" font-family="{FUENTE}" font-size="22" fill="{TINTA}">'
        f"{sin} de {total}: ninguna proteína diferente</text>",
        f'  <text x="32" y="278" font-family="{FUENTE}" font-size="22" fill="{TINTA}">'
        f"{con} de {total}: alguna proteína diferente</text>",
        f'  <text x="32" y="316" font-family="{FUENTE}" font-size="18" fill="{MUTED}">'
        "Esas "
        f"{con} dejan los contaminantes sin quitar "
        f'({", ".join(resumen["contaminantes_de_las_celdas_con_diferente"])}).</text>',
    ]
    desc = (
        f"n_celdas_suero={total}; n_con_diferente={con}; sin={sin}; "
        f"n_celdas_ml={resumen['n_celdas_ml']}"
    )
    return lienzo(
        920,
        400,
        "Casi ninguna decisión cambia el suero",
        desc,
        "\n".join(cuerpo),
        "Fuente: results/v3_7/multiverso_resumen.json. El verde es total menos las 13.",
    )


def figura_cptac(meta: dict[str, int], matriz: tuple[int, int, int]) -> str:
    tumor, normal, sin_tipo = matriz
    cuerpo = [
        f'  <text x="32" y="90" font-family="{FUENTE}" font-size="20" font-weight="700" fill="{TINTA}">'
        "Matriz que analiza la tesis</text>",
        f'  <text x="32" y="126" font-family="{FUENTE}" font-size="22" fill="{TINTA}">'
        f"{tumor} tumor · {normal} normal · {sin_tipo} sin tipo</text>",
        f'  <text x="32" y="158" font-family="{FUENTE}" font-size="16" fill="{MUTED}">'
        "docs/fuentes_de_datos.md, la misma cifra que la tabla 7 de la tesis.</text>",
        f'  <text x="32" y="214" font-family="{FUENTE}" font-size="20" font-weight="700" fill="{TINTA}">'
        "Metadato guardado en el repositorio</text>",
        f'  <text x="32" y="250" font-family="{FUENTE}" font-size="22" fill="{TINTA}">'
        f'{meta["tumor"]} tumor · {meta["normal"]} normal · {meta["unknown"]} sin tipo</text>',
        f'  <text x="32" y="286" font-family="{FUENTE}" font-size="22" fill="{TINTA}">'
        f'{meta["filas"]} filas · {meta["participantes"]} participantes con identificador'
        f' · {meta["sin_participante"]} sin él</text>',
        f'  <text x="32" y="318" font-family="{FUENTE}" font-size="16" fill="{MUTED}">'
        "metadata/sample_metadata/cptac_metadata.tsv. No es el tamaño de la matriz.</text>",
        f'  <text x="32" y="368" font-family="{FUENTE}" font-size="18" fill="{TINTA}">'
        "La tesis dice que el número de participantes no quedó confirmado con el metadato.</text>",
    ]
    if meta["otras"]:
        cuerpo.append(
            f'  <text x="32" y="396" font-family="{FUENTE}" font-size="16" fill="{MUTED}">'
            f'Otras condiciones en el metadato: {meta["otras"]}.</text>'
        )
    desc = (
        f"matriz {tumor}/{normal}/{sin_tipo}; metadato filas={meta['filas']} "
        f"tumor={meta['tumor']} normal={meta['normal']} unknown={meta['unknown']} "
        f"participantes={meta['participantes']}"
    )
    return lienzo(
        920,
        460,
        "CPTAC: dos conteos, no los sume",
        desc,
        "\n".join(cuerpo),
        "La matriz es la de la tesis. El metadato es un conteo de este repositorio.",
    )


def construir() -> dict[str, str]:
    pre = leer_json("preanalitica.json")
    h8 = leer_json("h8_carga_tumoral.json")
    k1 = leer_json("k1_etiquetas.json")
    embudo = leer_json("suero_embudo.json")
    modelos = leer_modelos()
    resumen = leer_json("multiverso_resumen.json")
    return {
        "grupos_suero.svg": figura_grupos_suero(pre),
        "estadios_pxd064655.svg": figura_estadios(h8),
        "concordancia_tejido.svg": figura_concordancia(k1),
        "clases_suero.svg": figura_clases(embudo),
        "auc_modelos.svg": figura_auc(modelos),
        "orden_adquisicion.svg": figura_orden(pre),
        "multiverso_suero.svg": figura_multiverso(resumen),
        "cptac_dos_conteos.svg": figura_cptac(contar_cptac(), matriz_cptac_documentada()),
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Escribe las figuras SVG del sitio.")
    parser.add_argument("--salida", type=Path, default=REPO / "web" / "figuras")
    args = parser.parse_args(argv)
    args.salida.mkdir(parents=True, exist_ok=True)
    for nombre, contenido in construir().items():
        (args.salida / nombre).write_text(contenido, encoding="utf-8")


if __name__ == "__main__":
    main()
