"""La tesis v3.6 solo usa macros de results/v3_6/ y no arrastra cifras viejas."""

from __future__ import annotations

import re
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
LATEX = RAIZ / "tesis-ucaldas" / "latex"
PROHIBIDAS = (
    "205 proteínas",
    "105 equivalentes",
    "0,981",
    "el orden de corrida no está",
    "no cambia el titular",
    "89 de 90",
)


def _tex() -> list[Path]:
    return [p for p in LATEX.rglob("*.tex") if p.name != "91-anexos.tex"]


def test_cada_macro_definida_se_usa_y_cada_uso_esta_definido():
    cifras = (LATEX / "cifras.tex").read_text(encoding="utf-8")
    definidas = set(re.findall(r"\\newcommand\{\\(\w+)\}", cifras))
    usadas = set()
    for path in LATEX.rglob("*.tex"):
        if path.name == "cifras.tex":
            continue
        usadas.update(re.findall(r"\\([A-Za-z]+)", path.read_text(encoding="utf-8")))
    assert definidas
    assert definidas <= usadas
    assert not (usadas & {nombre for nombre in usadas if nombre.startswith("Vcinco")})


def test_las_fuentes_de_las_macros_salen_de_v3_6():
    cifras = (LATEX / "cifras.tex").read_text(encoding="utf-8")
    fuentes = re.findall(r"% fuente: (\S+)", cifras)
    assert fuentes
    assert all(fuente.startswith("results/v3_6/") for fuente in fuentes)


def test_no_hay_flotantes_sin_redondear_ni_frases_viejas():
    for path in _tex():
        texto = path.read_text(encoding="utf-8")
        # 0,00099 es el formato de p con dos cifras significativas. El fallo que se vigila es el flotante sin redondear.
        assert not re.search(r"\d,\d{6,}", texto), path.name
        for frase in PROHIBIDAS:
            assert frase not in texto, f"{frase} en {path.name}"
