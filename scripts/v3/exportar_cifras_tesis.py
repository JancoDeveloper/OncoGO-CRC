"""Macros de la tesis y el resumen para el director. Las cifras salen de results/v3_3/."""

from __future__ import annotations

import json
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[2]
BASE = RAIZ / "results" / "v3_3"


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


def main() -> None:
    certificado = _json("certificado_v3_3.json")
    suero = _json("suero_resumen.json")
    principal = suero["clases"]["principal"]["crc_frente_a_colonoscopia_negativa"]
    k1 = _json("k1_etiquetas.json")
    r1 = k1["contrastes"]["primario_frente_a_mucosa"]
    macros = {
        "FraseSuero": certificado["frase_suero"],
        "TitularMl": certificado["titular_ml"]["texto"],
        "ClaveTitular": certificado["titular_ml"]["clave"],
        "SueroN": principal["n"],
        "SueroDiferente": principal["diferente"],
        "SueroEquivalente": principal["equivalente"],
        "SueroIndeterminada": principal["indeterminada"],
        "Huno": certificado["h1"],
        "Runo": f"{r1['numerador']}/{r1['denominador']}",
        "HunoFragil": certificado["h1_fragil_a_normalizacion"],
        "Hcinco": certificado["h5"]["estado"],
        "Hsiete": certificado["h7"]["estado"],
        "Hdiez": certificado["h10"]["estado"],
        "HdiezN": certificado["h10"]["n_especificas_crc"],
        "Advertencia": certificado["advertencia"],
    }
    lineas = ["% Generado por scripts/v3/exportar_cifras_tesis.py. No editar a mano."]
    for nombre, valor in macros.items():
        lineas.append(f"\\newcommand{{\\{nombre}}}{{{_tex(valor)}}}")
    destino = RAIZ / "tesis-ucaldas" / "latex" / "cifras.tex"
    destino.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    (RAIZ / "tesis-ucaldas" / "word" / "cifras.json").write_text(
        json.dumps(macros, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    capitulo = RAIZ / "tesis-ucaldas" / "latex" / "chapters" / "09-v33.tex"
    capitulo.write_text(
        "\n".join([
            "% Generado por scripts/v3/exportar_cifras_tesis.py.",
            "\\section{Corrida v3.3}",
            "Esta sección usa las macros de \\texttt{cifras.tex}. Esas macros salen de \\texttt{results/v3\\_3/}.",
            "",
            "\\FraseSuero\\ El contraste principal tiene \\SueroN\\ proteínas: \\SueroDiferente\\ diferentes, \\SueroEquivalente\\ equivalentes y \\SueroIndeterminada\\ indeterminadas.",
            "",
            "\\TitularMl",
            "",
            "H1 queda en \\Huno. La concordancia de signo del primario frente a la mucosa es \\Runo. Frágil a la normalización: \\HunoFragil.",
            "",
            "H5 queda en \\Hcinco. H7 queda en \\Hsiete. H10 queda en \\Hdiez, con \\HdiezN\\ proteínas específicas.",
            "",
            "CA1 pasa. AQP8, SLC26A3 y GUCA2A no tienen un log2FC calculable, así que las etiquetas de tumor y mucosa no quedan confirmadas. ESM-2 y TMbed no entraron en H7: el fichero los lista como variables faltantes.",
            "",
            "\\Advertencia",
            "",
        ]),
        encoding="utf-8",
    )
    resumen = RAIZ / "docs" / "resumen_director_v3_3.md"
    resumen.write_text(
        "\n".join([
            "# Resumen para el director, corrida v3.3",
            "",
            "Las cifras salen de `results/v3_3/certificado_v3_3.json` y de `suero_resumen.json`.",
            "",
            f"- {certificado['frase_suero']} La matriz principal tiene {principal['n']} proteínas: {principal['diferente']} diferentes, {principal['equivalente']} equivalentes y {principal['indeterminada']} indeterminadas.",
            f"- {certificado['titular_ml']['texto']} Clave: `{certificado['titular_ml']['clave']}`.",
            f"- H1: `{certificado['h1']}`. Concordancia {r1['numerador']}/{r1['denominador']}. Frágil a la normalización: {certificado['h1_fragil_a_normalizacion']}.",
            f"- H5: `{certificado['h5']['estado']}`. H7: `{certificado['h7']['estado']}`. H10: `{certificado['h10']['estado']}`, {certificado['h10']['n_especificas_crc']} específicas.",
            "- CA1 pasa. AQP8, SLC26A3 y GUCA2A fallan. Las etiquetas de tumor y mucosa no quedan confirmadas.",
            "- ESM-2, TMbed y Papier 2024 están en `variables_faltantes` de H7.",
            f"- {certificado['advertencia']}",
            "",
            "No está publicado. El checklist de la sección 16 no está en verde.",
            "",
        ]),
        encoding="utf-8",
    )
    print("[v3.3] cifras de la tesis escritas")
    escribir_v35(lineas, macros)


def _coma(valor: object) -> str:
    if isinstance(valor, float):
        return format(valor, ".12g").replace(".", ",")
    return str(valor).replace(".", ",")


def escribir_v35(lineas: list[str], macros: dict) -> None:
    """Añade las macros de v3.5. No pisa las de v3.3."""
    base = RAIZ / "results" / "v3_5"
    v33 = RAIZ / "results" / "v3_3"

    def leer(carpeta: Path, nombre: str) -> dict:
        return json.loads((carpeta / nombre).read_text(encoding="utf-8"))

    cert = leer(base, "certificado_v3_5.json")
    h5 = leer(base, "h5_especificidad_organo.json")
    h6 = leer(base, "h6_origen_celular.json")
    h7 = leer(base, "h7_visibilidad_en_sangre.json")
    h8 = leer(base, "h8_carga_tumoral.json")
    h9 = leer(base, "h9_wald.json")
    h10 = leer(base, "h10_especificidad.json")
    pre = leer(base, "preanalitica_v3_5.json")
    ml = leer(base, "ml_celda_preregistrada.json")
    multi = leer(base, "multiverso_v3_4.json")
    resumen = leer(base, "multiverso_resumen_v3_5.json")
    celda = multi["celda_preregistrada"]
    modelos = (v33 / "ml_suero_modelos.tsv").read_text(encoding="utf-8").splitlines()
    cab = modelos[0].split("\t")
    por_modelo = {fila.split("\t")[0]: dict(zip(cab, fila.split("\t"), strict=False)) for fila in modelos[1:]}
    nuevas = {
        "VcincoEstado": cert["estado"],
        "VcincoSueroN": celda["n"],
        "VcincoSueroDif": celda["diferente"],
        "VcincoSueroEq": celda["equivalente"],
        "VcincoSueroInd": celda["indeterminada"],
        "VcincoMdos": _coma(ml["auc"]),
        "VcincoMtres": _coma(por_modelo["M3_solo_qc"]["auc_oof_media"]),
        "VcincoMcinco": _coma(por_modelo["M5_principal_residualizada_qc"]["auc_oof_media"]),
        "VcincoMcincoP": _coma(por_modelo["M5_principal_residualizada_qc"]["p_permutacion"]),
        "VcincoHcinco": h5["estado"],
        "VcincoHseis": h6["estado"],
        "VcincoHsiete": h7["estado"],
        "VcincoHocho": h8["estado"],
        "VcincoHnueve": h9["estado"],
        "VcincoHdiez": h10["estado"],
        "VcincoHdiezN": h10["n_especificas_crc"],
        "VcincoPreCrudos": pre["n_crudos_validos"],
        "VcincoPreFecha": pre["n_con_fecha"],
        "VcincoPreEstado": pre["estado"],
        "VcincoTejido": macros["Runo"],
        "VcincoAdvertencia": cert["advertencia"],
        "VcincoAucContaminantes": _coma(resumen["auc_ml_sin_quitar_mediana"]),
    }
    for nombre, valor in nuevas.items():
        lineas.append(f"\\newcommand{{\\{nombre}}}{{{_tex(valor)}}}")
    destino = RAIZ / "tesis-ucaldas" / "latex" / "cifras.tex"
    destino.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    macros.update({k: str(v) for k, v in nuevas.items()})
    (RAIZ / "tesis-ucaldas" / "word" / "cifras.json").write_text(
        json.dumps(macros, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("[v3.5] cifras de la tesis escritas")


if __name__ == "__main__":
    main()
