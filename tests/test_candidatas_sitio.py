"""La página de candidatas y sus figuras salen de results/v3_7, sin inventar un enriquecimiento."""

from __future__ import annotations

import csv
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
WEB = REPO / "web"
V37 = REPO / "results" / "v3_7"
GITHUB = re.compile(
    r"https://github\.com/JancoDeveloper/OncoGO-CRC(?:/(blob|tree)/([A-Za-z0-9._-]+)(?:/([^\"'\s#]+))?)?"
)


def _modulo():
    spec = importlib.util.spec_from_file_location(
        "candidatas_sitio", REPO / "scripts" / "web" / "candidatas_sitio.py"
    )
    modulo = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(modulo)
    return modulo


def test_la_figura_se_regenera_igual(tmp_path):
    destino = tmp_path / "go_candidatas.svg"
    subprocess.run(
        [sys.executable, str(REPO / "scripts/web/candidatas_sitio.py"), "--salida", str(destino)],
        check=True,
    )
    assert destino.read_text(encoding="utf-8") == (WEB / "figuras" / "go_candidatas.svg").read_text(encoding="utf-8")


def test_la_tabla_usa_los_archivos_de_la_corrida():
    modulo = _modulo()
    pagina = (WEB / "candidatas.html").read_text(encoding="utf-8")
    embudo = __import__("json").loads((V37 / "priorizacion_embudo.json").read_text(encoding="utf-8"))
    filas = modulo.tabla()
    assert [fila["gen"] for fila in filas] == embudo["hipotesis_para_cohortes_futuras"]
    assert len(filas) == 20
    assert embudo["nadie_pasa_el_paso_6"] is True
    for fila in filas:
        assert f'>{fila["gen"]}</a>' in pagina
        assert f'uniprotkb/{fila["uniprot"]}' in pagina
        assert fila["fc_cptac"] in pagina
        assert fila["q_cptac_txt"] in pagina
        assert fila["clase_cn"] in pagina
    assert pagina.count(">HP</a>") == 1
    assert "HS</a>" not in pagina


def test_las_cuatro_del_suero_no_estan_en_las_veinte():
    modulo = _modulo()
    pagina = (WEB / "candidatas.html").read_text(encoding="utf-8")
    veinte = {fila["gen"] for fila in modulo.tabla()}
    for fila in modulo.senaladas():
        assert fila["gen"] not in veinte
        assert fila["clase_cn"] == "indeterminada"
        assert fila["veces_cn"] in pagina
        assert fila["fc_cn"] in pagina
    clases = {fila["gen"]: fila["clase_sanos"] for fila in modulo.senaladas()}
    assert clases == {"SAA1": "diferente", "CRP": "indeterminada", "LBP": "diferente", "LRG1": "diferente"}


def test_la_matriz_no_inventa_el_enriquecimiento_de_la_lista():
    modulo = _modulo()
    genes, columnas = modulo.matriz()
    assert len(columnas) == 12
    por_id = {col["go_id"]: col for col in columnas}
    assert set(por_id["GO:0005576"]["genes"]) == set(genes)
    assert "LUM" not in por_id["GO:0005886"]["genes"]
    assert set(por_id["GO:0042632"]["genes"]) == {"APOA1", "APOA2", "APOA4"}
    assert "GO:0006869" in por_id["GO:0042632"]["tambien"]
    svg = (WEB / "figuras" / "go_candidatas.svg").read_text(encoding="utf-8")
    assert "LUM GO:0005576 sí" in svg
    assert "LUM GO:0005886 no" in svg
    assert "enriquecimiento de esta lista" in svg
    lideres = set()
    with (V37 / "gsea_tejido.tsv").open(encoding="utf-8", newline="") as fh:
        for fila in csv.DictReader(fh, delimiter="\t"):
            if fila["cohorte"] != "cptac_tumor_frente_a_normal":
                continue
            if float(fila["FDR q-val"]) > 0.25:
                continue
            lideres.update(gen for gen in fila["Lead_genes"].split(";") if gen)
    assert set(genes).isdisjoint(lideres)


def test_los_enlaces_de_github_apuntan_a_v4_0_0_o_a_main():
    texto = "\n".join(path.read_text(encoding="utf-8") for path in WEB.glob("*.html"))
    assert "v3.7-final" not in texto
    assert "tree/v3.7" not in texto
    hallados = GITHUB.findall(texto)
    assert hallados, "el sitio no enlaza el repositorio"
    for tipo, ref, ruta in hallados:
        if not tipo:
            continue
        assert ref in {"v4.0.0", "main"}, ref
        if not ruta:
            continue
        subprocess.run(["git", "cat-file", "-e", f"v4.0.0:{ruta}"], cwd=REPO, check=True)


def test_la_pagina_de_lineas_abiertas_no_las_presenta_como_resultado():
    pagina = (WEB / "hacia-donde.html").read_text(encoding="utf-8")
    assert "Esta tesis no estudió neoantígenos ni péptidos presentados por HLA" in pagina
    for ancla in ("burgess2013", "chong2022", "xiang2026", "seo2026", "tesis2026", "collins2024"):
        assert f'referencias.html#{ancla}' in pagina
        assert f'id="{ancla}"' in (WEB / "referencias.html").read_text(encoding="utf-8")
    assert "h9_wald.json" in pagina
    assert "no decide la tesis" in pagina
