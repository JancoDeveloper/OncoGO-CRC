"""La página v3.5 muestra el valor del fichero y no las frases ya invalidadas."""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

_ruta = Path(__file__).resolve().parents[1] / "scripts" / "v3" / "exportar_sitio_v3_5.py"
_spec = importlib.util.spec_from_file_location("exportar_sitio_v3_5", _ruta)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
AVISO, leer, texto = _mod.AVISO, _mod.leer, _mod.texto

RAIZ = Path(__file__).resolve().parents[1]
SITIO = RAIZ / "sitio"
PAGINAS = [
    "index.html",
    "metodos.html",
    "resultados.html",
    "datos.html",
    "glosario.html",
    "preguntas.html",
    "reproducir.html",
    "historial.html",
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
        for fuente, valor in re.findall(r'data-fuente="([^"]+)">([^<]*)</data>', html):
            ruta, campo = fuente.split("#", 1)
            assert texto(leer(ruta, campo)) == valor, fuente


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
