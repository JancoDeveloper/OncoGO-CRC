"""La comparación de citas no necesita red.

Cubre los tres errores que obligaron a la versión 3 y el autor colectivo de
Zhang 2014, que Crossref registra como the NCI CPTAC.
"""

from oncogo_crc.v3.citas import (
    autores_coinciden,
    cita_local_markdown,
    comparar,
    titulos_coinciden,
)


def test_proteograph_no_es_el_titulo_de_alvez():
    assert not titulos_coinciden(
        "Next generation pan-cancer blood proteome profiling using proteograph",
        "Next generation pan-cancer blood proteome profiling using proximity extension assay",
    )
    assert titulos_coinciden(
        "Next generation pan-cancer blood proteome profiling using proximity extension assay",
        "Next generation pan-cancer blood proteome profiling using proximity extension assay",
    )


def test_el_doi_de_cascianelli_no_sostiene_el_titulo_de_praus():
    estado = comparar(
        titulo_local="Panomics reveals patient individuality as the major driver of colorectal cancer progression",
        autor_local="Herzog",
        anio_local="2023",
        titulo_remoto="Multi-label transcriptional classification of colorectal cancer reflects tumor cell population heterogeneity",
        autor_remoto="Cascianelli",
        anio_remoto="2023",
        resuelve=True,
    )
    assert "titulo_no_coincide" in estado
    assert "autor_no_coincide" in estado


def test_wannberg_2026_no_es_el_anio_de_crossref():
    estado = comparar(
        titulo_local="Plasma protein profiling predicts cancer in patients with non-specific symptoms",
        autor_local="Wannberg",
        anio_local="2026",
        titulo_remoto="Plasma protein profiling predicts cancer in patients with non-specific symptoms",
        autor_remoto="Wannberg",
        anio_remoto="2025",
        resuelve=True,
    )
    assert estado == "anio_no_coincide"


def test_autor_colectivo_de_zhang_2014():
    assert autores_coinciden(
        "the NCI CPTAC and Zhang, Bing and Wang, Jing",
        "the NCI CPTAC",
    )
    assert not autores_coinciden("Zhang, Bing and Wang, Jing", "the NCI CPTAC")


def test_un_doi_que_no_resuelve_falla():
    assert comparar(
        titulo_local="",
        autor_local="",
        anio_local="",
        titulo_remoto="",
        autor_remoto="",
        anio_remoto="",
        resuelve=False,
    ) == "doi_no_resuelve"


def test_markdown_toma_la_cita_pegada_al_doi_y_no_el_encabezado():
    texto = (
        "### Wannberg et al. 2024 — un encabezado que no es la cita.\n"
        "- **Cita**: Wannberg F, Qvick A, et al. Plasma protein profiling "
        "predicts cancer in patients with non-specific symptoms. *Nature "
        "Communications* 2025;17:151. doi:10.1038/s41467-025-67688-3\n"
    )
    local = cita_local_markdown(texto, "10.1038/s41467-025-67688-3")
    assert autores_coinciden(local["autor"], "Wannberg")
    assert local["anio"] == "2025"
    assert "plasma protein profiling" in local["titulo"].lower()
    assert "encabezado" not in local["titulo"].lower()
