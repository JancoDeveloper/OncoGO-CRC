"""Tabla y matriz GO de las proteínas que llegan al suero en la corrida v3.7.

No modifica results/. Lee la lista de priorizacion_embudo.json, los cambios de
priorizacion.tsv, triangulacion.tsv y los contrastes K4, el SHAP de
ml_explicacion.tsv y las tres variables GO de h7_tabla_proteinas.tsv.

La matriz de procesos biológicos no es un enriquecimiento de esas proteínas.
La corrida no lo calculó. La celda dice si el gen está en
metadata/gene_ontology/go_gene_long.tsv. Un término entra si lo comparten al
menos tres genes, o si lo comparten dos y, además, es significativo en el ORA
de CPTAC en la misma dirección que el log2FC de esos genes. Si dos términos
tienen el mismo conjunto de genes, se queda el de menor q en ese ORA.
"""

from __future__ import annotations

import argparse
import csv
import json
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
V37 = REPO / "results" / "v3_7"
GO_LARGO = REPO / "metadata" / "gene_ontology" / "go_gene_long.tsv"
MAPA = REPO / "metadata" / "protein_id_maps" / "gene_uniprot_human_reviewed_bulk.tsv"

FUENTE = "Segoe UI, Roboto, Helvetica, Arial, sans-serif"
TINTA = "#1a1a1a"
MUTED = "#3d3a37"
PAPEL = "#f7f6f3"
LINEA = "#d4d0cb"
AZUL = "#084c7a"
VERDE = "#009E73"
NARANJA = "#E69F00"
VERMILLON = "#D55E00"
MORADO = "#CC79A7"
CIELO = "#56B4E9"

# Nombres leídos el 20 de noviembre de 2026 en go.obo, data-version releases/2026-07-26.
# No son un resultado de la corrida: son el nombre oficial del identificador.
NOMBRES = {
    "GO:0005576": "Región extracelular",
    "GO:0005886": "Membrana plasmática",
    "GO:0070062": "Exosoma extracelular",
    "GO:0006869": "Transporte de lípidos",
    "GO:0033344": "Salida de colesterol",
    "GO:0033700": "Salida de fosfolípidos",
    "GO:0034375": "Remodelado de HDL",
    "GO:0042157": "Metabolismo de lipoproteínas",
    "GO:0042632": "Homeostasis del colesterol",
    "GO:0043691": "Transporte reverso de colesterol",
    "GO:0030301": "Transporte de colesterol",
    "GO:0006879": "Homeostasis del hierro",
    "GO:0006954": "Respuesta inflamatoria",
    "GO:0006953": "Respuesta de fase aguda",
    "GO:0030593": "Quimiotaxis de neutrófilos",
    "GO:0042742": "Defensa contra bacterias",
    "GO:0043542": "Migración endotelial",
    "GO:0061844": "Péptido antimicrobiano",
    "GO:0030199": "Fibrillas de colágeno",
    "GO:0050766": "Fagocitosis",
    "GO:0050821": "Estabilización de proteínas",
    "GO:0051496": "Fibras de estrés",
}
INGLES = {
    "GO:0005576": "extracellular region",
    "GO:0005886": "plasma membrane",
    "GO:0070062": "extracellular exosome",
    "GO:0006869": "lipid transport",
    "GO:0033344": "cholesterol efflux",
    "GO:0033700": "phospholipid efflux",
    "GO:0034375": "high-density lipoprotein particle remodeling",
    "GO:0042157": "lipoprotein metabolic process",
    "GO:0042632": "cholesterol homeostasis",
    "GO:0043691": "reverse cholesterol transport",
    "GO:0030301": "cholesterol transport",
    "GO:0006879": "intracellular iron ion homeostasis",
    "GO:0006954": "inflammatory response",
    "GO:0006953": "acute-phase response",
    "GO:0030593": "neutrophil chemotaxis",
    "GO:0042742": "defense response to bacterium",
    "GO:0043542": "endothelial cell migration",
    "GO:0061844": "antimicrobial humoral immune response mediated by antimicrobial peptide",
    "GO:0030199": "collagen fibril organization",
    "GO:0050766": "positive regulation of phagocytosis",
    "GO:0050821": "protein stabilization",
    "GO:0051496": "positive regulation of stress fiber assembly",
}
GRUPO = {
    "GO:0005576": "Lugar en la célula",
    "GO:0005886": "Lugar en la célula",
    "GO:0070062": "Lugar en la célula",
    "GO:0006869": "Lípidos",
    "GO:0033344": "Lípidos",
    "GO:0033700": "Lípidos",
    "GO:0034375": "Lípidos",
    "GO:0042157": "Lípidos",
    "GO:0042632": "Lípidos",
    "GO:0043691": "Lípidos",
    "GO:0030301": "Lípidos",
    "GO:0006879": "Hierro",
    "GO:0006954": "Inflamación",
    "GO:0006953": "Inflamación",
    "GO:0030593": "Inflamación",
    "GO:0042742": "Inflamación",
    "GO:0043542": "Inflamación",
    "GO:0061844": "Inflamación",
    "GO:0030199": "Matriz",
    "GO:0050766": "Otras anotaciones",
    "GO:0050821": "Otras anotaciones",
    "GO:0051496": "Otras anotaciones",
}
ORDEN_GRUPO = (
    "Lugar en la célula",
    "Lípidos",
    "Hierro",
    "Inflamación",
    "Matriz",
    "Otras anotaciones",
)
COMPONENTES = (
    ("GO:0005576", "go_extracelular"),
    ("GO:0005886", "go_membrana"),
    ("GO:0070062", "go_exosoma"),
)
# La tesis las nombra como cambios grandes del suero que no cumplen el criterio.
# No son las veinte de la integración tejido–suero.
OTRAS = ("SAA1", "CRP", "LBP", "LRG1")

SUPER = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def esc(texto: object) -> str:
    return (
        str(texto)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _decimal(texto: str) -> Decimal:
    return Decimal(str(texto).strip())


def fmt_fc(texto: str) -> str:
    cuantizado = _decimal(texto).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP)
    return format(cuantizado, "f").replace(".", ",")


def fmt_q(texto: str) -> str:
    numero = _decimal(texto)
    if numero == 0:
        return "0"
    if abs(numero) >= Decimal("0.01"):
        return format(numero.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP), "f").replace(".", ",")
    exponente = numero.adjusted()
    mantisa = numero.scaleb(-exponente).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP)
    if mantisa >= Decimal("10"):
        mantisa = Decimal("1.0")
        exponente += 1
    return f"{format(mantisa, 'f').replace('.', ',')} × 10{str(exponente).translate(SUPER)}"


def fmt_shap(texto: str) -> str:
    cuantizado = _decimal(texto).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)
    return format(cuantizado, "f").replace(".", ",")


def fmt_veces(texto: str) -> str:
    veces = (Decimal(2) ** _decimal(texto)).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP)
    return format(veces, "f").replace(".", ",")


def direccion(texto: str) -> str:
    numero = _decimal(texto)
    if numero > 0:
        return "Sube"
    if numero < 0:
        return "Baja"
    return "Sin cambio"


def _leer_tsv(ruta: Path) -> list[dict[str, str]]:
    with ruta.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def genes_lista() -> list[str]:
    embudo = json.loads((V37 / "priorizacion_embudo.json").read_text(encoding="utf-8"))
    genes = list(embudo["hipotesis_para_cohortes_futuras"])
    if embudo.get("nadie_pasa_el_paso_6") is not True or len(genes) != 20:
        raise SystemExit("priorizacion_embudo.json no trae las 20 proteínas que se quedan en el paso 5.")
    return genes


def accesos(genes: list[str]) -> dict[str, str]:
    hallados: dict[str, set[str]] = {gen: set() for gen in genes}
    with MAPA.open(encoding="utf-8", newline="") as fh:
        for fila in csv.DictReader(fh, delimiter="\t"):
            gen = fila["gene_symbol"].upper()
            if gen in hallados and fila["reviewed_status"] == "reviewed":
                hallados[gen].add(fila["uniprot_base"])
    salida = {}
    for gen, bases in hallados.items():
        if len(bases) != 1:
            raise SystemExit(f"{gen} no tiene un único acceso revisado: {bases}")
        salida[gen] = next(iter(bases))
    return salida


def _k4(nombre: str) -> dict[str, dict[str, str]]:
    filas = _leer_tsv(V37 / nombre)
    return {fila["protein_id"].upper(): fila for fila in filas}


def _shap_m2() -> dict[str, str]:
    mejor: dict[str, Decimal] = {}
    for fila in _leer_tsv(V37 / "ml_explicacion.tsv"):
        if fila["modelo"] != "M2_proteinas_principal":
            continue
        gen = fila["gen"].upper()
        valor = _decimal(fila["media_abs_shap"])
        if gen not in mejor or valor > mejor[gen]:
            mejor[gen] = valor
    return {gen: format(valor, "f") for gen, valor in mejor.items()}


def tabla(genes: list[str] | None = None) -> list[dict[str, str]]:
    """Una fila por gen, con el texto de la tabla y el valor crudo para ordenar."""
    genes = genes_lista() if genes is None else genes
    prioridad = {fila["gen"]: fila for fila in _leer_tsv(V37 / "priorizacion.tsv")}
    triangulo = {fila["gen"]: fila for fila in _leer_tsv(V37 / "triangulacion.tsv")}
    cn = _k4("k4_genes_maxlfq_crc_frente_a_colonoscopia_negativa.tsv")
    sanos = _k4("k4_genes_maxlfq_crc_frente_a_sanos.tsv")
    shap = _shap_m2()
    otro = _k4("k4_genes_maxlfq_colonoscopia_negativa_frente_a_sanos.tsv")
    uni = accesos(genes)
    filas = []
    for gen in genes:
        base = prioridad[gen]
        tri = triangulo[gen]
        suero = cn[gen]
        sano = sanos[gen]
        frente = otro[gen]
        if tri["pxd055218_k4_crc_cn"] != suero["clase"] or tri["pxd055218_k4_crc_sanos"] != sano["clase"]:
            raise SystemExit(f"{gen}: la clase de triangulacion.tsv no coincide con K4.")
        if _decimal(tri["cptac_coad_log2fc"]) != _decimal(base["log2fc_cptac"]):
            raise SystemExit(f"{gen}: el log2FC de CPTAC no coincide entre archivos.")
        filas.append(
            {
                "gen": gen,
                "uniprot": uni[gen],
                "log2fc_cptac": base["log2fc_cptac"],
                "q_cptac": base["padj_cptac"],
                "direccion": direccion(base["log2fc_cptac"]),
                "log2fc_tanaka": tri["tanaka_log2fc"],
                "q_tanaka": tri["tanaka_padj"],
                "log2fc_cn": suero["log2fc"],
                "q_cn": suero["padj"],
                "clase_cn": suero["clase"],
                "log2fc_sanos": sano["log2fc"],
                "q_sanos": sano["padj"],
                "clase_sanos": sano["clase"],
                "log2fc_cn_sanos": frente["log2fc"],
                "clase_cn_sanos": frente["clase"],
                "fc_cn_sanos": fmt_fc(frente["log2fc"]),
                "clase_h5": base["clase_h5"],
                "shap": shap.get(gen, ""),
                "fc_cptac": fmt_fc(base["log2fc_cptac"]),
                "q_cptac_txt": fmt_q(base["padj_cptac"]),
                "fc_tanaka": fmt_fc(tri["tanaka_log2fc"]),
                "q_tanaka_txt": fmt_q(tri["tanaka_padj"]),
                "fc_cn": fmt_fc(suero["log2fc"]),
                "q_cn_txt": fmt_q(suero["padj"]),
                "fc_sanos": fmt_fc(sano["log2fc"]),
                "veces_cn": fmt_veces(suero["log2fc"]),
                "shap_txt": fmt_shap(shap[gen]) if gen in shap else "—",
            }
        )
    return filas


def _anotacion_bp(genes: set[str]) -> dict[str, set[str]]:
    terminos: dict[str, set[str]] = {}
    with GO_LARGO.open(encoding="utf-8", newline="") as fh:
        for fila in csv.DictReader(fh, delimiter="\t"):
            if fila["namespace"] != "BP":
                continue
            gen = fila["gene_symbol"].upper()
            if gen not in genes:
                continue
            terminos.setdefault(fila["go_id"], set()).add(gen)
    return terminos


def _ora_cptac() -> dict[tuple[str, str], Decimal]:
    """q del ORA de CPTAC solo donde el archivo marca el término como significativo."""
    salida: dict[tuple[str, str], Decimal] = {}
    for fila in _leer_tsv(V37 / "go_ora_tejido.tsv"):
        if fila["cohorte"] != "cptac_tumor_frente_a_normal" or fila["significant"] != "True":
            continue
        clave = (fila["go_id"], fila["direccion"])
        q = _decimal(fila["padj"])
        if clave not in salida or q < salida[clave]:
            salida[clave] = q
    return salida


def terminos_proceso(genes: list[str], log2fc: dict[str, str]) -> list[dict[str, object]]:
    anotados = _anotacion_bp(set(genes))
    ora = _ora_cptac()
    elegibles: list[tuple[str, frozenset[str]]] = []
    for go_id, miembros in anotados.items():
        if len(miembros) >= 3:
            elegibles.append((go_id, frozenset(miembros)))
            continue
        if len(miembros) < 2:
            continue
        for direccion_ora in ("suben", "bajan"):
            if (go_id, direccion_ora) not in ora:
                continue
            signos = {_decimal(log2fc[gen]) > 0 for gen in miembros}
            espera_subir = direccion_ora == "suben"
            if signos == {espera_subir}:
                elegibles.append((go_id, frozenset(miembros)))
                break
    por_conjunto: dict[frozenset[str], list[str]] = {}
    for go_id, miembros in elegibles:
        por_conjunto.setdefault(miembros, []).append(go_id)

    def clave(go_id: str, miembros: frozenset[str]) -> tuple:
        mejores = []
        for direccion_ora, sube in (("suben", True), ("bajan", False)):
            if (go_id, direccion_ora) not in ora:
                continue
            if {_decimal(log2fc[gen]) > 0 for gen in miembros} == {sube}:
                mejores.append(ora[(go_id, direccion_ora)])
        q = min(mejores) if mejores else Decimal("1")
        return (0 if mejores else 1, q, go_id)

    elegidos = []
    for miembros, ids in por_conjunto.items():
        elegido = min(ids, key=lambda go_id: clave(go_id, miembros))
        elegidos.append((elegido, miembros, sorted(go for go in ids if go != elegido)))
    columnas = []
    for go_id, _miembros, tambien in elegidos:
        if go_id not in NOMBRES or go_id not in GRUPO:
            raise SystemExit(f"Falta el nombre público de {go_id}.")
        columnas.append(
            {
                "go_id": go_id,
                "nombre": NOMBRES[go_id],
                "ingles": INGLES[go_id],
                "grupo": GRUPO[go_id],
                "genes": sorted(anotados[go_id]),
                "n": len(anotados[go_id]),
                "fuente": "metadata/gene_ontology/go_gene_long.tsv",
                "tambien": tambien,
            }
        )
    columnas.sort(key=lambda col: (ORDEN_GRUPO.index(str(col["grupo"])), -int(col["n"]), str(col["go_id"])))
    return columnas


def componentes_h7(genes: list[str]) -> list[dict[str, object]]:
    filas = {fila["gen"]: fila for fila in _leer_tsv(V37 / "h7_tabla_proteinas.tsv")}
    columnas = []
    for go_id, campo in COMPONENTES:
        presentes = []
        for gen in genes:
            valor = filas[gen][campo]
            if valor not in {"0.0", "1.0"}:
                raise SystemExit(f"{gen} {campo} no es 0.0 ni 1.0: {valor}")
            if valor == "1.0":
                presentes.append(gen)
        columnas.append(
            {
                "go_id": go_id,
                "nombre": NOMBRES[go_id],
                "ingles": INGLES[go_id],
                "grupo": "Lugar en la célula",
                "genes": presentes,
                "n": len(presentes),
                "fuente": "results/v3_7/h7_tabla_proteinas.tsv",
                "tambien": [],
            }
        )
    return columnas


def senaladas() -> list[dict[str, str]]:
    """SAA1, CRP, LBP y LRG1: el suero las deja indeterminadas frente a colonoscopia negativa."""
    cn = _k4("k4_genes_maxlfq_crc_frente_a_colonoscopia_negativa.tsv")
    sanos = _k4("k4_genes_maxlfq_crc_frente_a_sanos.tsv")
    uni = accesos(list(OTRAS))
    filas = []
    for gen in OTRAS:
        suero = cn[gen]
        sano = sanos[gen]
        filas.append(
            {
                "gen": gen,
                "uniprot": uni[gen],
                "log2fc_cn": suero["log2fc"],
                "q_cn": suero["padj"],
                "clase_cn": suero["clase"],
                "fc_cn": fmt_fc(suero["log2fc"]),
                "q_cn_txt": fmt_q(suero["padj"]),
                "veces_cn": fmt_veces(suero["log2fc"]),
                "clase_sanos": sano["clase"],
                "fc_sanos": fmt_fc(sano["log2fc"]),
                "q_sanos_txt": fmt_q(sano["padj"]),
            }
        )
    return filas


def matriz(genes: list[str] | None = None) -> tuple[list[str], list[dict[str, object]]]:
    genes = genes_lista() if genes is None else genes
    log2fc = {fila["gen"]: fila["log2fc_cptac"] for fila in _leer_tsv(V37 / "priorizacion.tsv")}
    columnas = componentes_h7(genes) + terminos_proceso(genes, log2fc)
    return genes, columnas


def figura(genes: list[str], columnas: list[dict[str, object]]) -> str:
    ancho_celda = 36
    alto_celda = 24
    izquierda = 118
    techo = 248
    ancho = max(izquierda + ancho_celda * len(columnas) + 36, 980)
    alto = techo + alto_celda * len(genes) + 96
    partes = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {ancho} {alto}" role="img" aria-labelledby="titulo desc">',
        '<title id="titulo">Procesos que comparten las veinte proteínas y procesos que las distinguen</title>',
    ]
    desc = []
    for col in columnas:
        extras = ""
        if col["tambien"]:
            extras = " tambien=" + ",".join(col["tambien"])
        desc.append(f'{col["go_id"]} {col["ingles"]}: {",".join(col["genes"])}{extras} ({col["fuente"]})')
    partes.append(f'<desc id="desc">{esc(" | ".join(desc))}</desc>')
    partes.append(f'<rect width="{ancho}" height="{alto}" fill="{PAPEL}"/>')
    partes.append(
        f'<text x="16" y="34" font-family="{FUENTE}" font-size="20" font-weight="700" fill="{TINTA}">'
        "Qué comparten y qué las distingue</text>"
    )
    partes.append(
        f'<text x="16" y="58" font-family="{FUENTE}" font-size="14" fill="{MUTED}">'
        "Celda llena: el gen está anotado a ese término. No es un enriquecimiento de esta lista.</text>"
    )
    # Barras de grupo.
    indice = 0
    while indice < len(columnas):
        grupo = columnas[indice]["grupo"]
        fin = indice
        while fin < len(columnas) and columnas[fin]["grupo"] == grupo:
            fin += 1
        x = izquierda + indice * ancho_celda
        w = (fin - indice) * ancho_celda - 4
        color = {
            "Lugar en la célula": AZUL,
            "Lípidos": NARANJA,
            "Hierro": VERDE,
            "Inflamación": VERMILLON,
            "Matriz": MORADO,
            "Otras anotaciones": CIELO,
        }[str(grupo)]
        partes.append(f'<rect x="{x}" y="74" width="{w}" height="8" fill="{color}"/>')
        indice = fin
    for i, col in enumerate(columnas):
        x = izquierda + i * ancho_celda + 18
        partes.append(
            f'<text x="{x}" y="230" font-family="{FUENTE}" font-size="11" fill="{TINTA}" '
            f'transform="rotate(-65 {x} 230)">{esc(col["nombre"])}</text>'
        )
    presentes = {str(col["go_id"]): set(col["genes"]) for col in columnas}
    for r, gen in enumerate(genes):
        y = techo + r * alto_celda
        partes.append(
            f'<text x="{izquierda - 8}" y="{y + 16}" text-anchor="end" font-family="{FUENTE}" '
            f'font-size="13" fill="{TINTA}">{esc(gen)}</text>'
        )
        for c, col in enumerate(columnas):
            x = izquierda + c * ancho_celda
            si = gen in presentes[str(col["go_id"])]
            relleno = AZUL if si else "#ffffff"
            partes.append(
                f'<rect x="{x}" y="{y}" width="30" height="20" fill="{relleno}" stroke="{LINEA}">'
                f'<title>{esc(gen)} {col["go_id"]} {"sí" if si else "no"}</title></rect>'
            )
    y_leyenda = techo + len(genes) * alto_celda + 22
    partes.append(f'<rect x="16" y="{y_leyenda}" width="18" height="14" fill="{AZUL}" stroke="{LINEA}"/>')
    partes.append(
        f'<text x="40" y="{y_leyenda + 12}" font-family="{FUENTE}" font-size="13" fill="{TINTA}">Anotada</text>'
    )
    partes.append(f'<rect x="120" y="{y_leyenda}" width="18" height="14" fill="#ffffff" stroke="{LINEA}"/>')
    partes.append(
        f'<text x="144" y="{y_leyenda + 12}" font-family="{FUENTE}" font-size="13" fill="{TINTA}">No anotada</text>'
    )
    leyenda = (
        ("Lugar en la célula", AZUL),
        ("Lípidos", NARANJA),
        ("Hierro", VERDE),
        ("Inflamación", VERMILLON),
        ("Matriz", MORADO),
        ("Otras anotaciones", CIELO),
    )
    for i, (nombre, color) in enumerate(leyenda):
        x_ley = 250 + (i % 3) * 230
        y_item = y_leyenda + (i // 3) * 22
        partes.append(f'<rect x="{x_ley}" y="{y_item}" width="14" height="14" fill="{color}"/>')
        partes.append(
            f'<text x="{x_ley + 18}" y="{y_item + 12}" font-family="{FUENTE}" font-size="12" fill="{TINTA}">'
            f"{esc(nombre)}</text>"
        )
    partes.append(
        f'<text x="16" y="{alto - 18}" font-family="{FUENTE}" font-size="13" fill="{MUTED}">'
        "Lugar en la célula: h7_tabla_proteinas.tsv. Procesos: go_gene_long.tsv. Nombres: Gene Ontology 2026-07-26.</text>"
    )
    partes.append("</svg>\n")
    return "\n".join(partes) + "\n"


def construir() -> str:
    genes, columnas = matriz()
    return figura(genes, columnas)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Escribe la matriz GO de las proteínas candidatas.")
    parser.add_argument("--salida", type=Path, default=REPO / "web" / "figuras" / "go_candidatas.svg")
    args = parser.parse_args(argv)
    args.salida.parent.mkdir(parents=True, exist_ok=True)
    args.salida.write_text(construir(), encoding="utf-8")


if __name__ == "__main__":
    main()
