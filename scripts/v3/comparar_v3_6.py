#!/usr/bin/env python3
"""Compara las cifras de control de v3.5 con la corrida limpia. Corre fuera de la etapa."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
V36 = REPO / "results" / "v3_6"
V35 = REPO / "results" / "v3_5"
V33 = REPO / "results" / "v3_3"
SALIDA = V36 / "changelog_numeros.tsv"


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _fila(cifra: str, v35, v36, causa: str, seccion: str) -> dict:
    igual = v35 == v36
    if isinstance(v35, float) and isinstance(v36, float):
        igual = abs(v35 - v36) <= 1e-12
    return {
        "cifra": cifra,
        "valor_v3_5": v35,
        "valor_v3_6": v36,
        "relacion": "igual" if igual else "cambia",
        "causa": "" if igual else causa,
        "seccion": seccion,
    }


def main() -> None:
    embudo = _json(V36 / "suero_embudo.json")
    controles = embudo["controles"]
    celda = embudo["celda_preregistrada"]
    viejo_embudo = _json(V33 / "suero_controles.json")["controles"]
    ml = _json(V36 / "ml_suero_permutacion.json")["versiones"]
    ml_viejo = _json(V33 / "ml_suero_permutacion.json")["versiones"]
    h5 = _json(V36 / "h5_especificidad_organo.json")
    h5v = _json(V35 / "h5_especificidad_organo.json")
    h6 = _json(V36 / "h6_origen_celular.json")["fisher"]
    h6v = _json(V35 / "h6_origen_celular.json")["fisher"]
    h7 = _json(V36 / "h7_visibilidad_en_sangre.json")
    h7v = _json(V35 / "h7_visibilidad_en_sangre.json")
    h8 = _json(V36 / "h8_carga_tumoral.json")
    h8v = _json(V35 / "h8_carga_tumoral.json")
    pre = _json(V36 / "preanalitica.json")
    prev = _json(V35 / "preanalitica_v3_5.json")
    multi = _json(V36 / "multiverso_resumen.json")
    multiv = _json(V35 / "multiverso_resumen_v3_5.json")
    diferencia_h5 = round(1 - float(h5["maxima_concordancia_otro_organo"]), 3)
    diferencia_h5_v = round(1 - float(h5v["maxima_concordancia_otro_organo"]), 3)
    filas = [
        _fila("grupos_report", viejo_embudo["grupos"], controles["grupos"], "", "1.2"),
        _fila("grupos_q_global", viejo_embudo["grupos_global_pg_q_0_01"], controles["grupos_global_pg_q_0_01"], "", "1.2"),
        _fila("celda_n", 284, celda["n"], "", "1.2"),
        _fila("celda_diferente", 0, celda["diferente"], "", "1.2"),
        _fila("celda_equivalente", 166, celda["equivalente"], "", "1.2"),
        _fila("celda_indeterminada", 118, celda["indeterminada"], "", "1.2"),
        _fila("ml_m2", ml_viejo["M2_proteinas_principal"]["auc_oof_media"], ml["M2_proteinas_principal"]["auc_oof_media"], "", "1.2"),
        _fila("ml_m3", ml_viejo["M3_solo_qc"]["auc_oof_media"], ml["M3_solo_qc"]["auc_oof_media"], "", "1.2"),
        _fila("ml_m4", ml_viejo["M4_solo_coagulacion"]["auc_oof_media"], ml["M4_solo_coagulacion"]["auc_oof_media"], "", "1.2"),
        _fila("ml_m5", ml_viejo["M5_principal_residualizada_qc"]["auc_oof_media"], ml["M5_principal_residualizada_qc"]["auc_oof_media"], "", "1.2"),
        _fila("ml_m5_p", ml_viejo["M5_principal_residualizada_qc"]["p_permutacion"], ml["M5_principal_residualizada_qc"]["p_permutacion"], "", "1.2"),
        _fila("multiverso_suero", multiv.get("n_celdas"), multi["n_celdas_suero"], "", "1.2"),
        _fila("multiverso_ml", int(len(pd.read_csv(V35 / "multiverso_ml_v3_4.tsv", sep="\t"))), multi["n_celdas_ml"], "", "1.2"),
        _fila("multiverso_con_diferente", multiv.get("n_celdas_con_proteinas_diferentes"), multi["n_con_diferente"], "", "1.2"),
        _fila("h5_estado", h5v["estado"], h5["estado"], "", "1.2"),
        _fila("h5_universo", h5v["universo_principal_n"], h5["universo_principal_n"], "", "1.2"),
        _fila("h5_diferencia_con_el_otro_maximo", diferencia_h5_v, diferencia_h5, "", "1.2"),
        _fila("h6_estado", "no_sostiene", _json(V36 / "h6_origen_celular.json")["estado"], "", "1.2"),
        _fila("h6_tabla", json.dumps(h6v["tabla_mismo_signo_luego_distinto"]), json.dumps(h6["tabla_mismo_signo_luego_distinto"]), "", "1.2"),
        _fila(
            "h7_auc_secundario",
            h7v["modelo_secundario"]["auc_secundario"],
            h7["modelo_secundario"]["auc_secundario"],
            "La corrida limpia vuelve a calcular el modelo secundario, que añade ESM-2. El AUC principal de Álvez coincide. El veredicto sigue en refuta.",
            "4.1",
        ),
        _fila("h7_alvez", h7v["por_recurso"]["alvez_2023"]["auc"], h7["por_recurso"]["alvez_2023"]["auc"], "", "1.2"),
        _fila("h8_universo", h8v["n_universo"], h8["n_universo"], "", "1.2"),
        _fila("h8_estado", h8v["estado"], h8["estado"], "", "1.2"),
        _fila("pre_auc", prev["calculado_con_legibles"]["auc"]["valor"], pre["calculado_con_legibles"]["auc"]["valor"], "", "1.2"),
        _fila("pre_n_con_fecha", prev["n_con_fecha"], pre["n_con_fecha"], "", "1.2"),
        _fila("pre_bloques", None, pre["bloques_perfectos"], "Descriptivo nuevo de la enmienda v3.6.", "2.3"),
        _fila("pre_ic", prev["calculado_con_legibles"]["auc"]["ic95_ledell"], pre["calculado_con_legibles"]["auc"]["ic95_ledell"], "Con AUC 0 el IC de LeDell queda nulo. Es regla de presentación.", "2.3"),
    ]
    pd.DataFrame(filas).to_csv(SALIDA, sep="\t", index=False, na_rep="")
    cambios = [fila["cifra"] for fila in filas if fila["relacion"] == "cambia"]
    print("cambian:", ", ".join(cambios) if cambios else "ninguna")


if __name__ == "__main__":
    main()
