"""Tests de C4: n mínimo de términos, estado no_evaluable y tabla de verdad de H1."""

from __future__ import annotations

import numpy as np
import pandas as pd

from oncogo_crc.v2 import vias


def _tabla(terminos, direcciones, significativos, padj=None):
    n = len(terminos)
    return pd.DataFrame(
        {
            "termino": list(terminos),
            "direccion": list(direcciones),
            "estadistico": [float(d) * 3.0 for d in direcciones],
            "significativo": list(significativos),
            "padj": list(padj) if padj is not None else [0.001] * n,
        }
    )


def test_pocos_terminos_dejan_la_clausula_no_evaluable():
    """El fallo de v1: cuatro términos y la cláusula se dio por satisfecha."""
    terminos = [f"T{i}" for i in range(4)]
    des = _tabla(terminos, [1, 1, -1, 1], [True] * 4)
    rep = _tabla(terminos, [1, 1, -1, 1], [True] * 4)
    res = vias.concordancia_entre_cohortes(
        des, rep, coleccion="prueba", n_min_terminos=10, umbral=0.80
    )
    assert res["estado"] == "no_evaluable"
    assert res["concordancia"] is None
    assert res["n_terminos_evaluados"] == 4
    assert "mínimo preregistrado de 10" in res["motivo"]


def test_con_terminos_suficientes_y_concordancia_perfecta_sostiene():
    terminos = [f"T{i}" for i in range(20)]
    signos = [1] * 20
    des = _tabla(terminos, signos, [True] * 20)
    rep = _tabla(terminos, signos, [True] * 20)
    res = vias.concordancia_entre_cohortes(
        des, rep, coleccion="prueba", n_min_terminos=10, umbral=0.80
    )
    assert res["estado"] == "sostiene"
    assert res["concordancia"]["valor"] == 1.0
    assert res["margen_del_ic_inferior_sobre_el_umbral"] > 0


def test_el_punto_estimado_no_basta_hace_falta_el_intervalo():
    """20 de 22 da 0,909 de punto pero el IC baja de 0,80: no sostiene."""
    terminos = [f"T{i}" for i in range(22)]
    signos_des = [1] * 22
    signos_rep = [1] * 20 + [-1, -1]
    des = _tabla(terminos, signos_des, [True] * 22)
    rep = _tabla(terminos, signos_rep, [True] * 22)
    res = vias.concordancia_entre_cohortes(
        des, rep, coleccion="prueba", n_min_terminos=10, umbral=0.80
    )
    assert res["concordancia"]["valor"] > 0.90
    assert res["concordancia"]["ci_low"] < 0.80
    assert res["estado"] == "rechaza"


def test_la_seleccion_del_termino_no_mira_la_replica():
    """Solo los significativos en descubrimiento entran en el denominador."""
    terminos = [f"T{i}" for i in range(20)]
    des = _tabla(terminos, [1] * 20, [True] * 12 + [False] * 8)
    rep = _tabla(terminos, [1] * 12 + [-1] * 8, [True] * 20)
    res = vias.concordancia_entre_cohortes(
        des, rep, coleccion="prueba", n_min_terminos=10, umbral=0.80
    )
    # Los 8 discordantes no eran significativos en descubrimiento, así que no entran.
    assert res["n_terminos_evaluados"] == 12
    assert res["concordancia"]["valor"] == 1.0


def test_un_termino_ausente_en_la_replica_no_cuenta_como_discordante():
    terminos = [f"T{i}" for i in range(20)]
    des = _tabla(terminos, [1] * 20, [True] * 20)
    rep = _tabla(terminos[:14], [1] * 14, [True] * 14)
    res = vias.concordancia_entre_cohortes(
        des, rep, coleccion="prueba", n_min_terminos=10, umbral=0.80
    )
    assert res["n_terminos_evaluados"] == 14
    assert res["concordancia"]["denominador"] == 14


def test_la_tabla_de_verdad_de_h1_no_deja_resquicios():
    t = vias.tabla_verdad_h1
    assert t("rechaza", "sostiene")["veredicto"] == "rechaza_h1"
    assert t("rechaza", "no_evaluable")["veredicto"] == "rechaza_h1"
    assert t("rechaza", "rechaza")["veredicto"] == "rechaza_h1"
    assert t("sostiene", "rechaza")["veredicto"] == "rechaza_h1"
    assert t("sostiene", "sostiene")["veredicto"] == "sostiene_h1"
    # Una cláusula de vía sin evaluar no salva H1: la deja parcial.
    assert t("sostiene", "no_evaluable")["veredicto"] == "parcial_solo_proteina"
    assert t("no_evaluable", "sostiene")["veredicto"] == "no_evaluable"
    assert t("no_evaluable", "no_evaluable")["veredicto"] == "no_evaluable"


def test_leer_gmt_usa_el_identificador_de_reactome_y_el_nombre_de_hallmark(tmp_path):
    gmt = tmp_path / "SINTETICO_prueba.gmt"
    gmt.write_text(
        "Vía uno\tR-HSA-1234\tGENA\tgenb\n"
        "HALLMARK_ALGO\thttps://ejemplo\tGENC\tGEND\n"
        "linea mala\tsolo dos campos\n",
        encoding="utf-8",
    )
    col = vias.leer_gmt(gmt, nombre="prueba")
    assert col.estado == "usada"
    assert set(col.terminos) == {"R-HSA-1234", "HALLMARK_ALGO"}
    # Los símbolos se normalizan a mayúsculas.
    assert col.terminos["R-HSA-1234"] == {"GENA", "GENB"}


def test_un_gmt_vacio_queda_no_disponible(tmp_path):
    gmt = tmp_path / "SINTETICO_vacio.gmt"
    gmt.write_text("solo\tdos\n", encoding="utf-8")
    col = vias.leer_gmt(gmt, nombre="vacia")
    assert col.estado == "no_disponible"
    assert col.n_terminos == 0


def test_la_prueba_competitiva_encuentra_el_conjunto_desplazado():
    rng = np.random.default_rng(20260914)
    genes = np.array([f"G{i}" for i in range(400)])
    stats = rng.normal(size=400)
    stats[:30] += 3.0
    serie = pd.Series(stats, index=genes)
    terminos = {
        "desplazado": {f"G{i}" for i in range(30)},
        "neutro": {f"G{i}" for i in range(200, 240)},
        "pequeno": {"G0", "G1"},
    }
    tabla = vias.prueba_por_cohorte(
        serie, terminos, rho=0.01, df=90.0, min_proteinas=10, alpha=0.05
    )
    # El término de 2 proteínas no llega al mínimo y no aparece.
    assert set(tabla["termino"]) == {"desplazado", "neutro"}
    fila = tabla.set_index("termino").loc["desplazado"]
    assert bool(fila["significativo"]) is True
    assert int(fila["direccion"]) == 1
    assert bool(tabla.set_index("termino").loc["neutro", "significativo"]) is False


def test_bh_es_monotono_y_no_pasa_de_uno():
    p = np.array([0.001, 0.01, 0.2, 0.5, 0.9])
    q = vias._bh(p)
    assert np.all(np.diff(q) >= -1e-12)
    assert q.max() <= 1.0
    assert q[0] >= p[0]
