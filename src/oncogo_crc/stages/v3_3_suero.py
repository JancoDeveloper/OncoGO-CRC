"""S1. Matriz filtrada de PXD055218. No toca el report ni results/v3."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path


from oncogo_crc.settings import Settings
from oncogo_crc.v3.rutas_publicas import ruta_publica
from oncogo_crc.utils.dia_nn_plasma import (
    BYTES_REPORT,
    PLASMA_CLASICO,
    SHA256_REPORT,
    LIBRERIA_VACIA,
    controles_del_report,
    leer_libreria,
    leer_report,
    matriz_filtrada,
    trazabilidad,
)
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.utils.sello import escribir_sello, exigir_arbol_limpio, hacer_sello
from oncogo_crc.v3.prereg import PreregistroV3
from oncogo_crc.v3.salida import carpeta, exigir_vigente, manifiesto_congelado
from oncogo_crc.v3.suero_v32 import analizar_suero

ESPERADO = {
    "filas": 283_488,
    "corridas": 60,
    "grupos_global_pg_q_0_01": 329,
    "grupos_sin_pg_q_0_01_en_ninguna_corrida": 432,
    "humanos_pg_y_global": 312,
    "humanos_pg_y_global_con_gen": 307,
}


def _salida(path: Path) -> None:
    exigir_vigente(path)


def _sha256(path: Path) -> str:
    resumen = hashlib.sha256()
    with path.open("rb") as entrada:
        while True:
            bloque = entrada.read(1 << 20)
            if not bloque:
                break
            resumen.update(bloque)
    return resumen.hexdigest().upper()


def run(settings: Settings) -> None:
    prereg = PreregistroV3(settings.root)
    salida = carpeta(settings)
    _salida(salida / "marca.txt")
    salida.mkdir(parents=True, exist_ok=True)
    manifest = manifiesto_congelado(settings)
    sello = hacer_sello(settings.root, etapa="v3_3_suero", manifest=manifest)
    exigir_arbol_limpio(sello)
    report = settings.resolve(str(prereg.get("enmienda_v3_2.report_pxd055218.destino")))
    if not report.is_file():
        raise FileNotFoundError(report)
    if report.stat().st_size != BYTES_REPORT or _sha256(report) != SHA256_REPORT:
        raise RuntimeError(
            f"el report no coincide con la enmienda: {report.stat().st_size} bytes, sha {_sha256(report)}"
        )
    try:
        os.chmod(report, 0o444)
    except OSError:
        pass
    tabla = leer_report(report)
    cifras = controles_del_report(tabla)
    cifras["sha256"] = SHA256_REPORT
    cifras["bytes"] = BYTES_REPORT
    cifras["solo_lectura"] = oct(report.stat().st_mode & 0o777)
    distinto = {clave: cifras.get(clave) for clave in ESPERADO if cifras.get(clave) != ESPERADO[clave]}
    fasta = settings.resolve(str(prereg.get("enmienda_v3_2.contaminantes.fichero")))
    libreria = leer_libreria(fasta) if fasta.is_file() else LIBRERIA_VACIA
    piel = set(prereg.get("enmienda_v3_2.contaminante_piel.genes"))
    paneles = {
        "eritrocitos": set(prereg.get("enmienda_v3_2.paneles_geyer_2019.eritrocitos")),
        "plaquetas": set(prereg.get("enmienda_v3_2.paneles_geyer_2019.plaquetas")),
        "coagulacion": set(prereg.get("enmienda_v3_2.paneles_geyer_2019.coagulacion")),
    }
    traza = trazabilidad(tabla, libreria=libreria, genes_piel=piel, paneles=paneles)
    humanos_marcados = traza.loc[traza["contaminante"].astype(bool) & traza["humano"].astype(bool)]
    genes_marcados = set()
    for texto in humanos_marcados["genes"].astype(str):
        genes_marcados.update(parte.strip().upper() for parte in texto.split(";") if parte.strip())
    cruce = sorted(genes_marcados & PLASMA_CLASICO)
    if cruce:
        raise RuntimeError(f"proteína de plasma clásica marcada como contaminante: {cruce}")
    ancha = matriz_filtrada(tabla)
    destino_matriz = settings.resolve(str(prereg.get("enmienda_v3_2.report_pxd055218.matriz_filtrada")))
    destino_matriz.parent.mkdir(parents=True, exist_ok=True)
    ancha.to_csv(destino_matriz, sep="\t", index=False)
    _salida(salida / "suero_trazabilidad_grupos.tsv")
    write_tsv(traza, salida / "suero_trazabilidad_grupos.tsv")
    payload = {
        "sello": sello,
        "controles": cifras,
        "esperado": ESPERADO,
        "coincide_con_la_revision": distinto == {},
        "diferencias": distinto,
        "n_grupos_en_matriz_principal": int(ancha["protein_id"].nunique()) if "protein_id" in ancha.columns else 0,
        "n_contaminantes_en_la_traza": int(traza["contaminante"].sum()) if len(traza) else 0,
        "n_piel_en_la_traza": int(traza["contaminante_piel"].sum()) if len(traza) else 0,
        "libreria_n_entradas": libreria.n_entradas,
        "libreria_por_taxon": libreria.por_taxon,
        "n_contaminantes_humanos": int((traza["contaminante"] & traza["humano"]).sum()) if len(traza) else 0,
        "n_contaminantes_no_humanos": int((traza["contaminante"] & ~traza["humano"]).sum()) if len(traza) else 0,
        "n_grupos_mixtos": int(traza["comparte_con_contaminante"].sum()) if len(traza) else 0,
        "n_queratinas_epiteliales": int(traza["queratina_epitelial"].sum()) if len(traza) else 0,
        "grupos_humanos_contaminantes": [
            {"protein_group": str(fila.protein_group), "genes": str(fila.genes), "motivo": str(fila.contaminante_motivo)}
            for fila in humanos_marcados.itertuples(index=False)
        ],
        "contaminantes_no_marcan_plasma_clasico": "pasa",
        "matriz": ruta_publica(destino_matriz, settings)["ruta"],
    }
    destino = salida / "suero_controles.json"
    _salida(destino)
    destino.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    escribir_sello(salida / "suero_controles.json.sello.json", sello)
    escribir_sello(salida / "suero_trazabilidad_grupos.tsv.sello.json", sello)
    escribir_sello(salida / "sello_corrida.json", sello)
    if distinto or cifras["q_value_max"] > 0.01:
        raise RuntimeError(f"los controles no coinciden con la revisión: {distinto or cifras}")
    analizar_suero(tabla, traza, ancha, prereg, salida, sello, settings)
    print(f"[v3.2] matriz principal {payload['n_grupos_en_matriz_principal']} grupos -> {destino_matriz}")
