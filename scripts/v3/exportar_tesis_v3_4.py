"""Sección v3.4 de la tesis. Las cifras salen de results/v3_4/."""

from __future__ import annotations

import json
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[2]
BASE = RAIZ / "results" / "v3_4"
CAPITULO = RAIZ / "tesis-ucaldas" / "latex" / "chapters" / "09-v34.tex"


def _json(nombre: str) -> dict:
    return json.loads((BASE / nombre).read_text(encoding="utf-8"))


def _tex(texto: object) -> str:
    return (
        str(texto)
        .replace("\\", "\\textbackslash{}")
        .replace("&", "\\&")
        .replace("%", "\\%")
        .replace("_", "\\_")
        .replace("#", "\\#")
    )


def _oracion(prefijo: str, estado: object, motivo: object) -> str:
    texto = f"{prefijo} {_tex(estado)}."
    if motivo:
        resto = str(motivo).strip()
        if not resto.endswith("."):
            resto += "."
        texto += " " + _tex(resto)
    return texto


def _coma(valor: object, digitos: int = 3) -> str:
    if isinstance(valor, bool) or not isinstance(valor, (int, float)):
        return _tex(valor)
    return f"{float(valor):.{digitos}f}".replace(".", ",")


def main() -> None:
    certificado = _json("certificado_v3_4.json")
    h5 = _json("h5_especificidad_organo.json")
    h6 = _json("h6_origen_celular.json")
    h7 = _json("h7_visibilidad_en_sangre.json")
    h8 = _json("h8_carga_tumoral.json")
    h9 = _json("h9_wald.json")
    h10 = _json("h10_especificidad.json")
    multi = _json("multiverso_v3_4.json")
    pre = _json("preanalitica_v3_4.json")
    sec = _json("secuencia_v3_4.json")
    controles = _json("controles_v3_4.json")
    sello = certificado.get("sello") or {}
    filas = []
    for casilla in certificado["casillas"]:
        filas.append(
            f"{_tex(casilla['id'])} & {_tex(casilla['estado'])} & {_tex(casilla['motivo'])} \\\\"
        )
    recursos = []
    for nombre, bloque in (h7.get("por_recurso") or {}).items():
        if bloque.get("estado") == "calculado":
            recursos.append(
                f"{_tex(nombre)}: AUC {_coma(bloque.get('auc'))} "
                f"(intervalo {_coma(bloque.get('ic95_bajo'))} a {_coma(bloque.get('ic95_alto'))}, "
                f"{bloque.get('n_pos')} y {bloque.get('n_neg')})"
            )
        else:
            recursos.append(
                f"{_tex(nombre)}: {_tex(bloque.get('estado'))}, "
                f"{bloque.get('n_pos')} positivos y {bloque.get('n_neg')} negativos"
            )
    controles_txt = ", ".join(
        f"{_tex(fila['gen'])} {_tex(fila['estado'])}" for fila in controles.get("controles_fisicos") or []
    )
    secundario = h7.get("modelo_secundario") or {}
    lineas = [
        "% Generado por scripts/v3/exportar_tesis_v3_4.py. No editar a mano.",
        "\\section{Corrida v3.4}",
        "Esta sección sale de \\texttt{results/v3\\_4/}. No reutiliza las cifras de v3.3.",
        "",
        f"El certificado queda en {_tex(certificado['estado'])}. "
        f"La corrida firmada es el commit {_tex(sello.get('commit'))}, {_tex(sello.get('fecha_utc'))}.",
        "",
        "\\begin{tabular}{lll}",
        "\\toprule",
        "Casilla & Estado & Motivo \\\\",
        "\\midrule",
        *filas,
        "\\bottomrule",
        "\\end{tabular}",
        "",
        f"H5 queda en {_tex(h5.get('estado'))}. Clasifica {h5.get('n_clasificadas')} proteínas. "
        f"El universo principal tiene {h5.get('universo_principal_n')} genes. "
        f"La fracción pan-carcinoma es {_coma(h5.get('fraccion_pan_carcinoma'))}.",
        "",
        f"Los controles físicos del tejido dan veredicto {_tex(controles.get('veredicto_v3_4'))}. "
        f"{controles_txt}.",
        "",
        f"H6 queda en {_tex(h6.get('estado'))}. {_tex(h6.get('motivo'))}.",
        "",
        f"H7 queda en {_tex(h7.get('estado'))}. El modelo principal está completo: "
        f"{_tex(h7.get('modelo_completo'))}. "
        + " ".join(recursos)
        + ".",
        "",
        f"El modelo secundario queda en {_tex(secundario.get('estado'))}. "
        f"AUC principal {_coma(secundario.get('auc_principal'))}, "
        f"AUC secundaria {_coma(secundario.get('auc_secundario'))}, "
        f"incremento {_coma(secundario.get('incremento'), 4)}. "
        f"{_tex(secundario.get('nota'))} "
        f"Hay secuencia para {sec.get('n_con_secuencia')} genes y no la hay para {sec.get('n_sin_secuencia')}.",
        "",
        _oracion("H8 queda en", h8.get("estado"), h8.get("motivo")),
        "",
        _oracion("H9 queda en", h9.get("estado"), h9.get("motivo")),
        "",
        f"H10 queda en {_tex(h10.get('estado'))}, con {h10.get('n_especificas_crc')} proteínas específicas.",
        "",
        f"El multiverso tiene {multi.get('n_celdas_suero')} celdas de suero, "
        f"{multi.get('n_celdas_ml')} de aprendizaje y {multi.get('n_r1')} normalizaciones.",
        "",
        _oracion(
            f"La preanalítica queda en",
            f"{pre.get('estado')}, con {pre.get('n_crudos')} crudos contados",
            pre.get("motivo"),
        ),
        "",
        _tex(certificado.get("advertencia")),
        "",
    ]
    CAPITULO.write_text("\n".join(lineas), encoding="utf-8")
    print(f"[v3.4] capítulo escrito en {CAPITULO}")


if __name__ == "__main__":
    main()
