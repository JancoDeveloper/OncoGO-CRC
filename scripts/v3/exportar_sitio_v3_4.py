"""Escribe sitio/datos/v3_4.json y la página desde results/v3_4.

No toca el sitio si el certificado no dice publicable. Las cifras salen de
los JSON de esa corrida.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from oncogo_crc.stages.v3_4_report import puede_publicar

RAIZ = Path(__file__).resolve().parents[2]
SITIO = RAIZ / "sitio"
RESULTADOS = Path("/Volumes/JANCODATA/OncoGO-CRC/results/v3_4")


def _html(texto: object) -> str:
    return str(texto).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _leer(base: Path, nombre: str) -> dict | None:
    ruta = base / nombre
    if not ruta.is_file():
        return None
    return json.loads(ruta.read_text(encoding="utf-8"))


def datos_publicos(base: Path) -> dict:
    certificado = _leer(base, "certificado_v3_4.json")
    if certificado is None:
        raise FileNotFoundError(base / "certificado_v3_4.json")
    if not puede_publicar(certificado):
        raise RuntimeError(f"el certificado está en {certificado.get('estado')}; no se escribe el sitio")
    piezas = {
        "h5": _leer(base, "h5_especificidad_organo.json"),
        "h6": _leer(base, "h6_origen_celular.json"),
        "h7": _leer(base, "h7_visibilidad_en_sangre.json"),
        "h8": _leer(base, "h8_carga_tumoral.json"),
        "h9": _leer(base, "h9_wald.json"),
        "h10": _leer(base, "h10_especificidad.json"),
        "multiverso": _leer(base, "multiverso_v3_4.json"),
        "preanalitica": _leer(base, "preanalitica_v3_4.json"),
        "secuencia": _leer(base, "secuencia_v3_4.json"),
    }
    sello = certificado.get("sello") or {}
    return {
        "certificado": {clave: certificado[clave] for clave in certificado if clave != "sello"},
        "sello": {"commit": sello.get("commit"), "fecha_utc": sello.get("fecha_utc"), "etapa": sello.get("etapa")},
        "piezas": {
            nombre: {"estado": pieza.get("estado"), "motivo": pieza.get("motivo")}
            for nombre, pieza in piezas.items()
            if isinstance(pieza, dict)
        },
    }


def pagina(datos: dict) -> str:
    certificado = datos["certificado"]
    filas = "".join(
        "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            _html(c["id"]), _html(c["estado"]), _html(c["motivo"]), _html(c["fuente"])
        )
        for c in certificado["casillas"]
    )
    piezas = "".join(
        f"<li>{_html(nombre)}: {_html(pieza.get('estado'))}. {_html(pieza.get('motivo') or '')}</li>"
        for nombre, pieza in datos["piezas"].items()
    )
    sello = datos["sello"]
    return f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>OncoGO-CRC · corrida v3.4</title>
<style>
  :root {{ color-scheme: light; --ink:#142028; --paper:#f4f7f6; --card:#fff; --line:#c5d0cc; }}
  body {{ margin:0; font-family:"Avenir Next","Segoe UI",sans-serif; background:var(--paper); color:var(--ink); line-height:1.5; }}
  main {{ max-width:48rem; margin:0 auto; padding:1.5rem; }}
  table {{ border-collapse:collapse; width:100%; font-size:0.9rem; }}
  td, th {{ border-bottom:1px solid var(--line); text-align:left; padding:0.3rem; vertical-align:top; }}
  a {{ color:#083f3b; }}
</style>
</head>
<body>
<main>
<h1>OncoGO-CRC, corrida v3.4</h1>
<p>{_html(certificado["advertencia"])}</p>
<p>Certificado: {_html(certificado["estado"])}. Commit {_html(sello.get("commit"))}, {_html(sello.get("fecha_utc"))}.</p>
<h2>Checklist</h2>
<table>
<thead><tr><th>Casilla</th><th>Estado</th><th>Motivo</th><th>Fuente</th></tr></thead>
<tbody>{filas}</tbody>
</table>
<h2>Estados que escribió la corrida</h2>
<ul>{piezas}</ul>
<p><a href="v3_3.html">Corrida v3.3, página anterior.</a> <a href="datos/v3_4.json">JSON de esta página.</a></p>
<p>Los números los escribe <code>scripts/v3/exportar_sitio_v3_4.py</code> desde <code>results/v3_4/</code>. Esta página no decide sobre un paciente.</p>
</main>
</body>
</html>
"""


def escribir(base: Path, sitio: Path = SITIO) -> None:
    datos = datos_publicos(base)
    sitio.mkdir(parents=True, exist_ok=True)
    indice = sitio / "index.html"
    historico = sitio / "v3_3.html"
    if indice.is_file() and not historico.is_file() and "corrida v3.3" in indice.read_text(encoding="utf-8"):
        historico.write_text(indice.read_text(encoding="utf-8"), encoding="utf-8")
    (sitio / "datos").mkdir(parents=True, exist_ok=True)
    (sitio / "datos" / "v3_4.json").write_text(json.dumps(datos, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    indice.write_text(pagina(datos), encoding="utf-8")


def main() -> None:
    base = Path(sys.argv[1]) if len(sys.argv) > 1 else RESULTADOS
    try:
        escribir(base)
    except (RuntimeError, FileNotFoundError) as error:
        print(error)
        sys.exit(1)
    print(f"[v3.4] sitio escrito desde {base}")


if __name__ == "__main__":
    main()
