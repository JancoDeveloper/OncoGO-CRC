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


if __name__ == "__main__":
    main()
