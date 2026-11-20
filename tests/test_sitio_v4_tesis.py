"""La web de v4 habla de la tesis y no de la pregunta anterior."""

from __future__ import annotations

from pathlib import Path

WEB = Path(__file__).resolve().parents[1] / "web"
PAGINAS = {
    "index.html",
    "flujo.html",
    "metodos.html",
    "resultados.html",
    "robustez.html",
    "conclusiones.html",
    "datos.html",
    "descargas.html",
    "glosario.html",
    "citar.html",
}
PROHIBIDO = ("hipótesis", "hipotesis", "historia de la pregunta", "Aprende")


def test_estan_las_paginas_de_la_tesis():
    assert {path.name for path in WEB.glob("*.html")} == PAGINAS


def test_no_hay_lenguaje_de_la_pregunta_anterior():
    for path in WEB.glob("*.html"):
        texto = path.read_text(encoding="utf-8")
        for palabra in PROHIBIDO:
            assert palabra not in texto, f"{path.name} contiene {palabra}"


def test_no_hay_fichas_de_proteinas_ni_api_de_firebase():
    assert not (WEB / "proteinas").exists()
    texto = "\n".join(path.read_text(encoding="utf-8") for path in WEB.rglob("*") if path.is_file())
    assert "apiKey" not in texto
    assert "firebase" not in texto.lower()
