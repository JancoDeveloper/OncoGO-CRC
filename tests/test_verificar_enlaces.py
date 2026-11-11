"""El verificador no pide la página del editor. No sale a la red."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TITULO_RAI = (
    "Data-Independent Acquisition Analysis to Identify Serum Biomarkers and "
    "Neoantigens for Differentiating Colorectal Cancer from Inflammatory Bowel Diseases"
)


def modulo():
    spec = importlib.util.spec_from_file_location(
        "verificar_enlaces", REPO / "scripts" / "v3" / "verificar_enlaces.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_un_doi_solo_pide_el_handle_y_crossref():
    mod = modulo()
    pedidas = []

    def obtener(url: str):
        pedidas.append(url)
        if "/api/handles/" in url:
            cuerpo = json.dumps({
                "responseCode": 1,
                "values": [{"type": "URL", "data": {"value": "https://pubs.acs.org/doi/10.1021/acs.jproteome.5c00285"}}],
            })
            return 200, url, cuerpo
        cuerpo = json.dumps({"message": {"title": [TITULO_RAI]}})
        return 200, url, cuerpo

    fila = {
        "id": "doi_rai",
        "url": "https://doi.org/10.1021/acs.jproteome.5c00285",
        "metodo": "doi",
        "titulo_esperado": TITULO_RAI,
        "doi_esperado": "",
    }
    resultado = mod.comprobar_fila(fila, obtener)
    assert resultado["estado_http"] == "200"
    assert pedidas == [
        "https://doi.org/api/handles/10.1021/acs.jproteome.5c00285",
        "https://api.crossref.org/works/10.1021/acs.jproteome.5c00285",
    ]
    assert not any("pubs.acs.org" in url or "linkinghub.elsevier.com" in url or "europepmc.org" in url for url in pedidas)


def test_europepmc_solo_pide_la_api_de_ebi():
    mod = modulo()
    pedidas = []

    def obtener(url: str):
        pedidas.append(url)
        cuerpo = json.dumps({
            "hitCount": 1,
            "resultList": {"result": [{"title": TITULO_RAI + ".", "doi": "10.1021/acs.jproteome.5c00285"}]},
        })
        return 200, url, cuerpo

    fila = {
        "id": "rai_europepmc",
        "url": "https://europepmc.org/article/MED/40652393",
        "metodo": "europepmc",
        "titulo_esperado": TITULO_RAI,
        "doi_esperado": "10.1021/acs.jproteome.5c00285",
    }
    resultado = mod.comprobar_fila(fila, obtener)
    assert resultado["estado_http"] == "200"
    assert pedidas == [
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:40652393%20AND%20SRC:MED&format=json"
    ]


def test_doi_datacite_solo_pide_la_api_de_datacite():
    mod = modulo()
    pedidas = []

    def obtener(url: str):
        pedidas.append(url)
        cuerpo = json.dumps({
            "data": {"attributes": {"titles": [{"title": "A Unified Approach to Interpreting Model Predictions"}]}},
        })
        return 200, url, cuerpo

    fila = {
        "id": "doi_lundberg2017",
        "url": "https://doi.org/10.48550/arXiv.1705.07874",
        "metodo": "doi_datacite",
        "titulo_esperado": "A Unified Approach to Interpreting Model Predictions",
        "doi_esperado": "",
    }
    resultado = mod.comprobar_fila(fila, obtener)
    assert resultado["estado_http"] == "200"
    assert pedidas == ["https://api.datacite.org/dois/10.48550/arXiv.1705.07874"]
    assert not any("arxiv.org/pdf" in url or "europepmc.org" in url for url in pedidas)


def test_el_user_agent_lleva_el_correo(monkeypatch):
    mod = modulo()
    visto = {}

    class Respuesta:
        def __init__(self, url):
            self._url = url

        def getcode(self):
            return 200

        def geturl(self):
            return self._url

        def read(self, _n):
            return b"{}"

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

    def falso(peticion, timeout=45):
        visto["agente"] = peticion.get_header("User-agent")
        return Respuesta(peticion.full_url)

    monkeypatch.setattr(mod.urllib.request, "urlopen", falso)
    mod.obtener_red("https://api.crossref.org/works/10.1021/acs.jproteome.5c00285")
    assert visto["agente"] == "OncoGO-CRC/3.7 (mailto:jancod@jancod.co)"
