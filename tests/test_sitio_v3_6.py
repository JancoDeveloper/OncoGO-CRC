"""El sitio v3.6 sale de results/v3_6/ y no arrastra cifras de corridas anteriores."""

from __future__ import annotations

import re
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
SITIO = RAIZ / "sitio" / "archivo" / "v3_6"
VIGENTE = RAIZ / "sitio"
AVISO = "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
PAGINAS = [
    "index.html", "como-se-hizo.html", "resultados.html", "aprende.html",
    "datos.html", "glosario.html", "preguntas.html", "reproducir.html",
    "historial.html", "404.html",
]
PROHIBIDAS = ("205 proteínas", "105 equivalentes", "0,981", "89 de 90", "no cambia el titular")


def test_paginas_aviso_y_figuras():
    for nombre in PAGINAS:
        texto = (SITIO / nombre).read_text(encoding="utf-8")
        assert texto.count(AVISO) >= 2, nombre
        assert "/Volumes/" not in texto and "/Users/" not in texto and "C:\\" not in texto
        for frase in PROHIBIDAS:
            assert frase not in texto
    resultados = (SITIO / "resultados.html").read_text(encoding="utf-8")
    for i in range(1, 21):
        assert f'id="g{i}"' in resultados
    assert len(re.findall(r"<dt[ >]", (SITIO / "glosario.html").read_text(encoding="utf-8"))) >= 40


def test_el_estado_v3_7_coincide_con_sus_ficheros():
    import json

    pagina = (VIGENTE / "index.html").read_text(encoding="utf-8")
    pares = re.findall(r'<data value="([^"]*)" data-fuente="(results/v3_7/[^"]+)"', pagina)
    assert pares
    for valor, fuente in pares:
        ruta, _, campo = fuente.partition("#")
        datos = json.loads((RAIZ / ruta).read_text(encoding="utf-8"))
        actual = datos
        for parte in campo.split("."):
            actual = actual[int(parte)] if isinstance(actual, list) else actual[parte]
        esperado = str(actual).lower() if isinstance(actual, bool) else str(actual)
        assert esperado == valor, fuente


def test_cada_cifra_tiene_fuente_de_v3_6():
    for path in SITIO.glob("*.html"):
        for valor, fuente in re.findall(r'<data value="([^"]*)" data-fuente="([^"]+)"', path.read_text(encoding="utf-8")):
            assert fuente.startswith("results/v3_6/") or fuente.startswith("results/v3_7/"), fuente
            assert valor
