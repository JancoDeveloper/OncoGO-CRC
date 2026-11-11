"""La página v3.5 muestra el valor del fichero y no las frases ya invalidadas."""

from __future__ import annotations

import json
import re
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
SITIO = RAIZ / "sitio" / "archivo" / "v3_6"
AVISO = "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
PAGINAS = [
    "index.html",
    "como-se-hizo.html",
    "resultados.html",
    "aprende.html",
    "datos.html",
    "glosario.html",
    "preguntas.html",
    "reproducir.html",
    "historial.html",
    "404.html",
]
PROHIBIDAS = (
    "directLFQ no está instalado",
    "120 crudos",
    "H5 refuta",
    "las etiquetas no quedan confirmadas",
    "89/90",
    "205 proteínas",
)


def test_cada_cifra_sale_de_su_fichero():
    for nombre in PAGINAS:
        html = (SITIO / nombre).read_text(encoding="utf-8")
        assert html.count(AVISO) >= 2
        for valor, fuente in re.findall(r'<data value="([^"]*)" data-fuente="([^"]+)"', html):
            ruta, _, campo = fuente.partition("#")
            assert (RAIZ / ruta).is_file(), fuente
            if not campo:
                continue
            datos = json.loads((RAIZ / ruta).read_text(encoding="utf-8"))
            for parte in campo.split("."):
                datos = datos[parte]
            assert str(datos) == valor, fuente


def test_no_hay_frases_invalidadas_fuera_del_historial():
    for nombre in PAGINAS:
        if nombre == "historial.html":
            continue
        html = (SITIO / nombre).read_text(encoding="utf-8")
        for frase in PROHIBIDAS:
            assert frase not in html


def test_glosario_tiene_al_menos_treinta_terminos():
    html = (SITIO / "glosario.html").read_text(encoding="utf-8")
    assert len(re.findall(r"<dt ", html)) >= 30
