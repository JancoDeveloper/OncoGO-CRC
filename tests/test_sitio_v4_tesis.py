"""La web habla de la tesis, en lenguaje llano, y no de la pregunta anterior."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

WEB = Path(__file__).resolve().parents[1] / "web"
REPO = WEB.parent
PAGINAS = {
    "index.html",
    "problema.html",
    "flujo.html",
    "datos.html",
    "resultados.html",
    "robustez.html",
    "reproducibilidad.html",
    "conclusiones.html",
    "autor.html",
    "referencias.html",
    "glosario.html",
    "aprender.html",
    "aprender-marco.html",
    "aprender-metodos.html",
    "creditos-imagenes.html",
}
# «Aprender» es la sección nueva de fondo. Sigue prohibido el lenguaje de la pregunta anterior.
PROHIBIDO = ("hipótesis", "hipotesis", "historia de la pregunta")
ENLACE = re.compile(r'(?:href|src)="([^"]+)"')
TEXTO = {".html", ".css", ".svg", ".md", ".js", ".txt", ".json"}


def _archivos_texto():
    for path in WEB.rglob("*"):
        if path.is_file() and path.suffix.lower() in TEXTO:
            yield path


def test_estan_las_paginas_de_la_tesis():
    assert {path.name for path in WEB.glob("*.html")} == PAGINAS


def test_no_hay_lenguaje_de_la_pregunta_anterior():
    for path in _archivos_texto():
        texto = path.read_text(encoding="utf-8")
        for palabra in PROHIBIDO:
            assert palabra not in texto, f"{path.relative_to(WEB)} contiene {palabra}"


def test_no_hay_fichas_de_proteinas_ni_api_de_firebase():
    assert not (WEB / "proteinas").exists()
    texto = "\n".join(path.read_text(encoding="utf-8") for path in _archivos_texto())
    assert "apiKey" not in texto
    assert "firebase" not in texto.lower()


def test_los_enlaces_locales_existen():
    paginas = {path.resolve(): path.read_text(encoding="utf-8") for path in WEB.rglob("*.html")}
    for path, texto in paginas.items():
        for url in ENLACE.findall(texto):
            if url.startswith(("http://", "https://", "mailto:")):
                continue
            if url.startswith("#"):
                destino, ancla = path, url[1:]
            else:
                base, _, ancla = url.partition("#")
                destino = (path.parent / base).resolve()
                assert destino.is_file(), f"{path.name} apunta a {url}"
            if not ancla:
                continue
            cuerpo = paginas.get(destino)
            if cuerpo is None:
                cuerpo = destino.read_text(encoding="utf-8")
            assert f'id="{ancla}"' in cuerpo, f"{path.name} apunta a #{ancla} y no existe"


def test_bonilla_solo_aparece_en_la_ficha():
    texto = "\n".join(path.read_text(encoding="utf-8") for path in _archivos_texto())
    assert "en castellano" not in texto.lower()
    assert "Carlos Eduardo" not in texto
    assert "coautor" not in texto.lower()
    assert "revisión colombiana" not in texto.lower()
    assert texto.count("Bonilla") == 1
    assert "y Bonilla, C. E. (2026)" in (WEB / "referencias.html").read_text(encoding="utf-8")


def test_las_cifras_visibles_salen_de_v3_7():
    resultados = (WEB / "resultados.html").read_text(encoding="utf-8")
    robustez = (WEB / "robustez.html").read_text(encoding="utf-8")
    k1 = json.loads((REPO / "results/v3_7/k1_etiquetas.json").read_text(encoding="utf-8"))
    embudo = json.loads((REPO / "results/v3_7/suero_embudo.json").read_text(encoding="utf-8"))
    resumen = json.loads((REPO / "results/v3_7/multiverso_resumen.json").read_text(encoding="utf-8"))
    primario = k1["contrastes"]["primario_frente_a_mucosa"]
    metastasis = k1["contrastes"]["metastasis_hepatica_frente_a_mucosa"]
    celda = embudo["celda_preregistrada"]
    assert f'{primario["numerador"]} de {primario["denominador"]}' in resultados
    assert f'{metastasis["numerador"]} de {metastasis["denominador"]}' in resultados
    assert str(celda["n"]) in resultados
    assert str(celda["equivalente"]) in resultados
    assert str(celda["indeterminada"]) in resultados
    sin = int(resumen["n_celdas_suero"]) - int(resumen["n_con_diferente"])
    assert f'{sin} de {resumen["n_celdas_suero"]}' in robustez


def test_las_figuras_se_regeneran_igual(tmp_path):
    script = REPO / "scripts/web/figuras_sitio.py"
    subprocess.run([sys.executable, str(script), "--salida", str(tmp_path)], check=True)
    generadas = sorted(path.name for path in tmp_path.glob("*.svg"))
    publicadas = sorted(path.name for path in (WEB / "figuras").glob("*.svg"))
    assert generadas == publicadas
    for nombre in generadas:
        assert (tmp_path / nombre).read_text(encoding="utf-8") == (WEB / "figuras" / nombre).read_text(encoding="utf-8")
