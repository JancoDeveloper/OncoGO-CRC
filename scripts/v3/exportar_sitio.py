"""Escribe sitio/datos/v3_3.json y las páginas. Las cifras salen de results/v3_3/."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parents[2]
BASE = RAIZ / "results" / "v3_3"
SITIO = RAIZ / "sitio"


def _json(nombre: str) -> dict:
    return json.loads((BASE / nombre).read_text(encoding="utf-8"))


def _html(texto: object) -> str:
    return (
        str(texto)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def main() -> None:
    certificado = _json("certificado_v3_3.json")
    controles = pd.read_csv(BASE / "controles.tsv", sep="\t")
    mucosa = controles.loc[controles["control"].isin(["CA1", "AQP8", "SLC26A3", "GUCA2A"])]
    changelog = pd.read_csv(BASE / "changelog_numeros.tsv", sep="\t")
    multiverso = _json("multiverso.json")
    datos = {
        "certificado": {clave: certificado[clave] for clave in certificado if clave != "sello"},
        "controles_mucosa": mucosa[["control", "observado", "pasa"]].to_dict(orient="records"),
        "multiverso": {clave: multiverso[clave] for clave in ("conclusion_preregistrada", "fraccion_de_celdas_con_esa_conclusion", "n_celdas", "nota")},
        "changelog": changelog.fillna("").to_dict(orient="records"),
    }
    (SITIO / "datos").mkdir(parents=True, exist_ok=True)
    (SITIO / "datos" / "v3_3.json").write_text(json.dumps(datos, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (SITIO / "index.html").write_text(_pagina(datos), encoding="utf-8")
    (SITIO / "simulador.html").write_text(_simulador(), encoding="utf-8")
    print("[v3.3] sitio escrito desde results/v3_3")


def _pagina(datos: dict) -> str:
    certificado = datos["certificado"]
    filas_control = "".join(
        f"<li>{_html(fila['control'])}: {_html(fila['pasa'])}, observado {_html(fila['observado'])}</li>"
        for fila in datos["controles_mucosa"]
    )
    filas_cambio = "".join(
        "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            _html(fila["metrica"]), _html(fila["valor_v3_2"]), _html(fila["valor_v3_3"]), _html(fila["causa"])
        )
        for fila in datos["changelog"]
    )
    figuras = [
        ("F1", "Resumen de la corrida"),
        ("F2", "Grupos que quedan en el contraste principal"),
        ("F3", "Indicadores de calidad por grupo"),
        ("F4", "log2FC del contraste principal"),
        ("F5", "AUC de v3.2 y de v3.3"),
        ("F6", "AUC de ruido con y sin filas duplicadas"),
        ("F7", "Controles del tejido"),
        ("F8", "log2FC antes y después del centrado"),
        ("F9", "De las proteínas certificadas a las específicas"),
        ("F10", "Clases de H5"),
        ("F11", "Rejilla del multiverso de K4"),
        ("F12", "Curvas de H7"),
        ("F13", "Proteínas de partida en el atlas"),
        ("F14", "Fase de Pepe de esta corrida"),
    ]
    bloques = "\n".join(
        f'<figure><img src="figuras/{nombre}.png" alt="{_html(alt)}"><figcaption>{_html(alt)}. Fuente: results/v3_3.</figcaption></figure>'
        for nombre, alt in figuras
    )
    faltan = ", ".join(certificado["h7"]["variables_faltantes"])
    return f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>OncoGO-CRC · corrida v3.3</title>
<style>
  :root {{ color-scheme: light; --ink:#142028; --paper:#f4f7f6; --card:#fff; --line:#c5d0cc; }}
  body {{ margin:0; font-family:"Avenir Next","Segoe UI",sans-serif; background:var(--paper); color:var(--ink); line-height:1.5; }}
  main {{ max-width:46rem; margin:0 auto; padding:1.5rem; }}
  figure {{ margin:1.5rem 0; background:var(--card); border:1px solid var(--line); padding:0.75rem; }}
  img {{ max-width:100%; height:auto; }}
  table {{ border-collapse:collapse; width:100%; font-size:0.9rem; }}
  td, th {{ border-bottom:1px solid var(--line); text-align:left; padding:0.3rem; vertical-align:top; }}
  a {{ color:#083f3b; }}
</style>
</head>
<body>
<main>
<p><a href="#contenido">Saltar al contenido</a></p>
<h1 id="contenido">OncoGO-CRC, corrida v3.3</h1>
<p>{_html(certificado["advertencia"])}</p>
<p>{_html(certificado["frase_suero"])}</p>
<p>{_html(certificado["titular_ml"]["texto"])}</p>
<p>H1: {_html(certificado["h1"])}. Frágil a la normalización: {_html(certificado["h1_fragil_a_normalizacion"])}.</p>
<p>H5: {_html(certificado["h5"]["estado"])}. H7: {_html(certificado["h7"]["estado"])}. H10: {_html(certificado["h10"]["estado"])}, {_html(certificado["h10"]["n_especificas_crc"])} proteínas específicas, fase {_html(certificado["h10"]["fase_pepe"])}.</p>
<h2>Controles del tejido</h2>
<ul>{filas_control}</ul>
<p>Con AQP8, SLC26A3 y GUCA2A en falla, las etiquetas de tumor y mucosa no quedan confirmadas.</p>
<h2>Lo que esta corrida no cierra</h2>
<ul>
<li>Variables que H7 lista como faltantes: {_html(faltan)}.</li>
<li>El multiverso calculado tiene { _html(datos["multiverso"]["n_celdas"]) } celdas. {_html(datos["multiverso"]["nota"])}</li>
</ul>
<h2>Figuras</h2>
{bloques}
<h2>Qué cambió frente a v3.2</h2>
<table>
<thead><tr><th>Cifra</th><th>v3.2</th><th>v3.3</th><th>Causa</th></tr></thead>
<tbody>{filas_cambio}</tbody>
</table>
<p><a href="lectura-31-oct.html">Lectura del 31 de octubre, histórica.</a> <a href="simulador.html">Simulador del intervalo de una AUC.</a></p>
<p>Los números de esta página los escribe <code>scripts/v3/exportar_sitio.py</code> desde <code>results/v3_3/</code>.</p>
</main>
</body>
</html>
"""


def _simulador() -> str:
    return """<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Simulador del intervalo de una AUC</title>
<style>
  body { font-family:"Avenir Next","Segoe UI",sans-serif; max-width:40rem; margin:1.5rem auto; line-height:1.5; color:#142028; }
  label { display:block; margin-top:0.8rem; }
  output { display:block; margin-top:1rem; font-size:1.1rem; }
</style>
</head>
<body>
<p><a href="index.html">Volver a la corrida v3.3</a></p>
<h1>Anchura del intervalo de una AUC</h1>
<p>Esta página no usa un resultado del estudio. Aplica la fórmula de Hanley y McNeil (1982) a los dos números que escriba el lector. Riley y colaboradores (2020) advierten que el tamaño de un modelo de predicción no se decide solo con este intervalo.</p>
<form>
<label>AUC esperada <input id="auc" type="number" min="0.5" max="0.99" step="0.01" value="0.70"></label>
<label>Muestras en el grupo 1 <input id="n1" type="number" min="2" step="1" value="20"></label>
<label>Muestras en el grupo 0 <input id="n0" type="number" min="2" step="1" value="20"></label>
</form>
<output id="salida"></output>
<noscript><p>Sin JavaScript, la fórmula queda escrita: el error estándar es la raíz de [AUC(1−AUC) + (n1−1)(Q1−AUC²) + (n0−1)(Q2−AUC²)] / (n1 n0), con Q1 = AUC/(2−AUC) y Q2 = 2 AUC²/(1+AUC). La anchura del intervalo del 95 % es 2 × 1,96 × ese error.</p></noscript>
<script>
function calcular() {
  const auc = Number(document.getElementById("auc").value);
  const n1 = Number(document.getElementById("n1").value);
  const n0 = Number(document.getElementById("n0").value);
  const q1 = auc / (2 - auc);
  const q2 = 2 * auc * auc / (1 + auc);
  const varianza = auc * (1 - auc) + (n1 - 1) * (q1 - auc * auc) + (n0 - 1) * (q2 - auc * auc);
  const se = Math.sqrt(varianza / (n1 * n0));
  const anchura = 2 * 1.96 * se;
  document.getElementById("salida").textContent = "Anchura aproximada del intervalo del 95 %: " + anchura.toFixed(3) + ". No es un resultado de OncoGO-CRC.";
}
document.querySelectorAll("input").forEach((campo) => campo.addEventListener("input", calcular));
calcular();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
