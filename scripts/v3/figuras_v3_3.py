"""Figuras F1 a F14. Cada número sale de results/v3_3/. No recalcula un contraste."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve

from oncogo_crc.utils.sello import escribir_sello, exigir_arbol_limpio, hacer_sello
from oncogo_crc.v3.salida import manifiesto_congelado
from oncogo_crc.settings import load_settings

RAIZ = Path(__file__).resolve().parents[2]
VERMELLON = "#D55E00"
AZUL = "#0072B2"
VERDE = "#009E73"
CIELO = "#56B4E9"
NARANJA = "#E69F00"
MORADO = "#CC79A7"
TINTA = "#000000"
CLASE = {
    "diferente": VERMELLON,
    "equivalente": VERDE,
    "indeterminada": CIELO,
    "conflicto": NARANJA,
}
GRUPO = {"crc": VERMELLON, "colonoscopy_negative": AZUL, "healthy": VERDE}


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _guardar(fig, nombre: str, titulo: str, sello: dict, destinos: list[Path]) -> None:
    fig.suptitle(titulo, fontsize=11)
    for destino in destinos:
        destino.mkdir(parents=True, exist_ok=True)
        fig.savefig(destino / f"{nombre}.png", dpi=300, bbox_inches="tight")
        fig.savefig(destino / f"{nombre}.svg", bbox_inches="tight")
        escribir_sello(destino / f"{nombre}.sello.json", {**sello, "figura": nombre, "titulo": titulo})
    plt.close(fig)


def _f1(base: Path, sello: dict, destinos: list[Path]) -> None:
    certificado = _json(base / "certificado_v3_3.json")
    controles = _tsv(base / "controles.tsv")
    mucosa = controles.loc[controles["analisis"].eq("tumor_frente_a_mucosa") & controles["control"].isin(["CA1", "AQP8", "SLC26A3", "GUCA2A"])]
    texto_controles = ", ".join(f"{fila.control} {fila.pasa}" for fila in mucosa.itertuples())
    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.axis("off")
    bloques = [
        "Pregunta\nCáncer frente a colonoscopia negativa",
        f"Suero\n{certificado['frase_suero']}",
        f"Modelo\n{certificado['titular_ml']['texto']}",
        f"Tejido\nH1 {certificado['h1']}. Controles: {texto_controles}.",
        f"Sangre desde el tejido\nH7 {certificado['h7']['estado']}. H10 {certificado['h10']['estado']}, {certificado['h10']['n_especificas_crc']} específicas.",
    ]
    for i, texto in enumerate(bloques):
        ax.text(0.02, 0.92 - i * 0.18, texto, va="top", fontsize=8, transform=ax.transAxes)
    _guardar(fig, "F1", certificado["titular_ml"]["texto"], sello, destinos)


def _f2(base: Path, sello: dict, destinos: list[Path]) -> None:
    controles = _json(base / "suero_controles.json")["controles"]
    resumen = _json(base / "suero_resumen.json")
    etiquetas = ["Grupos del report", "Global.PG.Q ≤ 0,01", "Humanos con PG.Q y Global.PG.Q", "Matriz del contraste principal"]
    valores = [controles["grupos"], controles["grupos_global_pg_q_0_01"], controles["humanos_pg_y_global"], resumen["n_principal"]]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.barh(etiquetas[::-1], valores[::-1], color=[CIELO, AZUL, NARANJA, VERMELLON])
    for y, valor in enumerate(valores[::-1]):
        ax.text(valor, y, f" {valor}", va="center", fontsize=8)
    ax.set_xlabel("Grupos de proteínas")
    _guardar(fig, "F2", f"De {controles['grupos']} grupos del report, el contraste principal usa {resumen['n_principal']}.", sello, destinos)


def _f3(base: Path, sello: dict, destinos: list[Path]) -> None:
    calidad = _tsv(base / "suero_calidad_muestra.tsv")
    pruebas = _tsv(base / "k7_pruebas.tsv").set_index("indicador")
    columnas = ["eritrocitos", "plaquetas", "coagulacion", "piel"]
    fig, ejes = plt.subplots(1, len(columnas), figsize=(10, 3.6), sharey=False)
    rng = np.random.default_rng(20260914)
    for ax, columna in zip(ejes, columnas, strict=True):
        for i, grupo in enumerate(("crc", "colonoscopy_negative", "healthy")):
            valores = calidad.loc[calidad["grupo"].eq(grupo), columna].to_numpy(dtype=float)
            ax.scatter(np.full(valores.size, i) + rng.uniform(-0.12, 0.12, valores.size), valores, s=12, color=GRUPO[grupo], label=grupo)
        p = float(pruebas.loc[columna, "p_kruskal"])
        ax.set_title(f"{columna}\np = {p:.3g}", fontsize=8)
        ax.set_xticks(range(3), ["CP", "CN", "Sanos"], fontsize=7)
    _guardar(fig, "F3", "Los indicadores de calidad no son iguales entre cáncer, colonoscopia negativa y personas sanas.", sello, destinos)


def _f4(base: Path, sello: dict, destinos: list[Path]) -> None:
    tabla = _tsv(base / "k4_crc_frente_a_colonoscopia_negativa.tsv")
    tabla = tabla.sort_values("log2fc", kind="stable")
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.axvspan(-0.585, 0.585, color="#F0E442", alpha=0.35)
    for clase, bloque in tabla.groupby("clase", sort=False):
        ax.scatter(bloque["log2fc"], np.arange(len(tabla))[tabla["clase"].eq(clase)], s=8, color=CLASE.get(clase, TINTA), label=f"{clase} ({len(bloque)})")
    ax.axvline(0, color=TINTA, lw=0.6)
    ax.set_xlabel("log2FC, cáncer frente a colonoscopia negativa")
    ax.set_yticks([])
    ax.legend(fontsize=7, frameon=False)
    conteo = tabla["clase"].value_counts()
    _guardar(
        fig,
        "F4",
        f"Ninguna proteína es diferente: {int(conteo.get('equivalente', 0))} equivalentes y {int(conteo.get('indeterminada', 0))} indeterminadas.",
        sello,
        destinos,
    )


def _f5(base: Path, sello: dict, destinos: list[Path]) -> None:
    viejo = _json(RAIZ / "results/v3_2/ml_suero_permutacion.json")["versiones"]
    nuevo = _tsv(base / "ml_suero_modelos.tsv").set_index("modelo")
    pares = [
        ("filtrada_humana", "M1_proteinas_filtrada_humana"),
        ("sin_contaminantes_ni_paneles", "M2_proteinas_principal"),
        ("solo_coagulacion", "M4_solo_coagulacion"),
    ]
    fig, ax = plt.subplots(figsize=(8, 3.8))
    for i, (clave_v, clave_n) in enumerate(pares):
        ax.plot([viejo[clave_v]["auc_media"], nuevo.loc[clave_n, "auc_oof_media"]], [i, i], color=TINTA, lw=1)
        ax.scatter([viejo[clave_v]["auc_media"]], [i], color=NARANJA, label="v3.2, pliegues con filas duplicadas" if i == 0 else None)
        bajo = float(nuevo.loc[clave_n, "ic95_ledell_bajo"])
        alto = float(nuevo.loc[clave_n, "ic95_ledell_alto"])
        ax.plot([bajo, alto], [i, i], color=AZUL, lw=4, alpha=0.35, solid_capstyle="butt")
        ax.scatter([nuevo.loc[clave_n, "auc_oof_media"]], [i], color=AZUL, label="v3.3, fuera de pliegue" if i == 0 else None)
    media_nulo = float(nuevo.loc["M2_proteinas_principal", "media_del_nulo"])
    ax.axvline(media_nulo, color=VERDE, ls="--", label=f"media del nulo de M2, {media_nulo:.3f}")
    ax.set_yticks(range(len(pares)), ["Matriz humana", "Principal", "Solo coagulación"])
    ax.set_xlabel("AUC")
    ax.set_xlim(0.45, 1.02)
    ax.legend(fontsize=7, frameon=False)
    _guardar(fig, "F5", "Al quitar las filas duplicadas, el AUC del contraste principal baja y el de la matriz humana también.", sello, destinos)


def _f6(base: Path, sello: dict, destinos: list[Path]) -> None:
    fuga = _json(base / "ml_demostracion_fuga.json")
    fig, ax = plt.subplots(figsize=(6, 3.2))
    nombres = ["Sin duplicar", "Filas duplicadas"]
    valores = [fuga["auc_sin_duplicar"], fuga["auc_con_filas_duplicadas"]]
    ax.scatter(valores, [0, 1], s=40, color=[VERDE, VERMELLON])
    ax.set_yticks([0, 1], nombres)
    ax.set_xlabel("AUC media de cinco pliegues, sobre ruido normal")
    ax.set_xlim(0.4, 1.0)
    _guardar(fig, "F6", f"En ruido, duplicar las filas sube el AUC de {fuga['auc_sin_duplicar']} a {fuga['auc_con_filas_duplicadas']}.", sello, destinos)


def _f7(base: Path, sello: dict, destinos: list[Path]) -> None:
    controles = _tsv(base / "controles.tsv")
    mucosa = controles.loc[controles["control"].isin(["CA1", "AQP8", "SLC26A3", "GUCA2A"])].copy()
    fig, ax = plt.subplots(figsize=(7, 3.6))
    for i, fila in enumerate(mucosa.itertuples()):
        valor = pd.to_numeric(pd.Series([fila.observado]), errors="coerce").iloc[0]
        if np.isfinite(valor):
            ax.scatter([valor], [i], color=VERDE, s=40)
            ax.text(valor, i, f"  {fila.control} log2FC {valor:.2f}", va="center", fontsize=8)
        else:
            ax.text(0.15, i, f"{fila.control}: sin log2FC. Estado: {fila.pasa}", va="center", fontsize=8, color=VERMELLON)
    ax.axvline(0, color=TINTA, lw=0.6)
    ax.set_xlim(-3.2, 2.4)
    ax.set_yticks([])
    ax.set_xlabel("log2FC, tumor frente a mucosa")
    _guardar(fig, "F7", "CA1 baja en el tumor. AQP8, SLC26A3 y GUCA2A no tienen un log2FC calculable.", sello, destinos)


def _f8(base: Path, sello: dict, destinos: list[Path]) -> None:
    sensibilidad = _json(base / "k1_sensibilidad_normalizacion.json")
    antes = _tsv(base / "k1_primario_frente_a_mucosa.tsv")["log2fc"].dropna()
    despues = _tsv(base / "k1_primario_frente_a_mucosa_centrado.tsv")["log2fc"].dropna()
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    ax.hist(antes, bins=40, density=True, histtype="step", color=NARANJA, label="Sin centrar")
    ax.hist(despues, bins=40, density=True, histtype="step", color=AZUL, label="Centrado por mediana")
    ax.axvline(sensibilidad["antes"]["mediana_log2fc"], color=NARANJA, ls="--", lw=0.8)
    ax.axvline(sensibilidad["despues"]["mediana_log2fc"], color=AZUL, ls="--", lw=0.8)
    ax.set_xlabel("log2FC pareado, primario frente a mucosa")
    ax.legend(fontsize=7, frameon=False)
    fragil = sensibilidad["h1_fragil_a_normalizacion"]
    _guardar(
        fig,
        "F8",
        f"La mediana del log2FC pasa de {sensibilidad['antes']['mediana_log2fc']:.3f} a {sensibilidad['despues']['mediana_log2fc']:.3f}. H1 frágil: {fragil}.",
        sello,
        destinos,
    )


def _f9(base: Path, sello: dict, destinos: list[Path]) -> None:
    k1 = _json(base / "k1_etiquetas.json")
    resumen = _json(base / "suero_resumen.json")
    h10 = _json(base / "h10_especificidad.json")
    cobertura = k1["cobertura"]
    etiquetas = [
        "Certificadas en CPTAC",
        "Con réplica de signo en Tanaka",
        "Matriz principal de suero",
        "Específicas de CRC en el atlas",
    ]
    valores = [k1["n_certificadas"], cobertura["numerador"], resumen["n_principal"], h10["n_especificas_crc"]]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.barh(etiquetas[::-1], valores[::-1], color=[AZUL, CIELO, NARANJA, VERMELLON])
    for y, valor in enumerate(valores[::-1]):
        ax.text(valor, y, f" {valor}", va="center", fontsize=8)
    ax.set_xlabel("Proteínas")
    _guardar(fig, "F9", f"De {k1['n_certificadas']} proteínas certificadas, {h10['n_especificas_crc']} salen específicas de cáncer colorrectal en el atlas.", sello, destinos)


def _f10(base: Path, sello: dict, destinos: list[Path]) -> None:
    h5 = _json(base / "h5_especificidad_organo.json")
    conteo = h5["conteo"]
    fig, ax = plt.subplots(figsize=(6, 3.4))
    izquierda = 0
    for clave, color in (("pan_carcinoma", VERMELLON), ("especifica_colon", AZUL), ("intermedia", CIELO)):
        ax.barh([0], [conteo[clave]], left=izquierda, color=color, label=f"{clave} ({conteo[clave]})")
        izquierda += conteo[clave]
    ax.axvline(0.20 * h5["n_clasificadas"], color=TINTA, ls="--", label="0,20 de las clasificadas")
    ax.set_yticks([])
    ax.legend(fontsize=7, frameon=False)
    _guardar(fig, "F10", f"La fracción pan-carcinoma es {h5['fraccion_pan_carcinoma']:.3f}. H5 {h5['estado']}.", sello, destinos)


def _f11(base: Path, sello: dict, destinos: list[Path]) -> None:
    tabla = _tsv(base / "multiverso.tsv").reset_index(drop=True)
    meta = _json(base / "multiverso.json")
    tabla = tabla.sort_values(["diferente", "equivalente"], kind="stable")
    fig, ax = plt.subplots(figsize=(8, 3.8))
    ax.plot(range(len(tabla)), tabla["diferente"], color=AZUL)
    marcada = tabla["es_la_celda_preregistrada"].astype(str).str.lower().isin(["true", "1"])
    if marcada.any():
        posicion = np.flatnonzero(marcada.to_numpy())
        ax.scatter(posicion, tabla["diferente"].to_numpy()[posicion], color=VERMELLON, zorder=3, label="Celda preregistrada")
    cuantificacion = ", ".join(sorted(tabla["cuantificacion"].astype(str).unique()))
    ax.set_xlabel("Celdas de la rejilla, ordenadas")
    ax.set_ylabel("Proteínas diferentes")
    ax.legend(fontsize=7, frameon=False)
    _guardar(
        fig,
        "F11",
        f"En la rejilla {cuantificacion}, {meta['fraccion_de_celdas_con_esa_conclusion']:.3f} de {meta['n_celdas']} celdas dan {meta['conclusion_preregistrada']}.",
        sello,
        destinos,
    )


def _f12(base: Path, sello: dict, destinos: list[Path]) -> None:
    predicciones = _tsv(base / "h7_predicciones.tsv")
    h7 = _json(base / "h7_visibilidad_en_sangre.json")
    fig, ax = plt.subplots(figsize=(6.2, 5))
    ax.plot([0, 1], [0, 1], color="#999999", lw=0.7)
    colores = {"alvez_2023": VERMELLON, "hpa_disease_blood_atlas": AZUL}
    for recurso, bloque in predicciones.groupby("recurso"):
        fpr, tpr, _ = roc_curve(bloque["y"], bloque["p"])
        informe = h7["por_recurso"][recurso]
        ax.plot(fpr, tpr, color=colores.get(recurso, TINTA), label=f"{recurso}: AUC {informe['auc']:.3f}")
    ax.set_xlabel("Falsos positivos")
    ax.set_ylabel("Verdaderos positivos")
    ax.legend(fontsize=7, frameon=False)
    incluye = all(info["ic95_bajo"] <= 0.5 <= info["ic95_alto"] for info in h7["por_recurso"].values())
    frase = "el intervalo incluye 0,50 en los dos recursos" if incluye else "el intervalo no incluye 0,50 en todos los recursos"
    _guardar(fig, "F12", f"H7 {h7['estado']}: {frase}.", sello, destinos)


def _f13(base: Path, sello: dict, destinos: list[Path]) -> None:
    tabla = _tsv(base / "h10_proteinas_de_partida.tsv")
    columnas = ["logfc_sanos", "logfc_misma_clase", "logfc_todas_las_demas"]
    matriz = tabla[columnas].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    imagen = ax.imshow(matriz, cmap="coolwarm", aspect="auto")
    ax.set_xticks(range(3), ["Sanos", "Misma clase", "Otras enfermedades"], fontsize=8)
    ax.set_yticks(range(len(tabla)), [f"{gen} ({clase})" for gen, clase in zip(tabla["gen"], tabla["clase"], strict=True)], fontsize=8)
    fig.colorbar(imagen, ax=ax, label="logFC")
    medidas = tabla.loc[tabla["clase"].ne("no_medida"), "gen"].astype(str).tolist()
    clases = tabla.loc[tabla["clase"].ne("no_medida"), "clase"].astype(str).tolist()
    detalle = ", ".join(f"{gen} {clase}" for gen, clase in zip(medidas, clases, strict=True)) or "ninguna medida"
    _guardar(fig, "F13", f"Proteínas de partida medidas en el atlas: {detalle}.", sello, destinos)


def _f14(base: Path, sello: dict, destinos: list[Path]) -> None:
    certificado = _json(base / "certificado_v3_3.json")
    ml = _json(base / "ml_suero_permutacion.json")
    h7 = _json(base / "h7_visibilidad_en_sangre.json")
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.axis("off")
    lineas = [
        "Fase 1 de Pepe. Hallazgos de esta corrida, no una prueba.",
        f"Suero: {certificado['frase_suero']}",
        f"Modelo, fase {ml['fase_pepe']}: {certificado['titular_ml']['clave']}",
        f"H7, fase {h7['fase_pepe']}: {certificado['h7']['estado']}",
        f"H10, fase {certificado['h10']['fase_pepe']}: {certificado['h10']['estado']}, {certificado['h10']['n_especificas_crc']} específicas",
        certificado["advertencia"],
    ]
    ax.text(0.02, 0.95, "\n\n".join(lineas), va="top", fontsize=9, transform=ax.transAxes)
    _guardar(fig, "F14", "Los hallazgos de v3.3 están en la fase 1 de Pepe.", sello, destinos)


def main() -> None:
    ajustes = load_settings()
    sello = hacer_sello(ajustes.root, etapa="v3_3_figuras", manifest=manifiesto_congelado(ajustes))
    exigir_arbol_limpio(sello)
    base = RAIZ / "results" / "v3_3"
    destinos = [
        base / "report" / "figuras",
        RAIZ / "tesis-ucaldas" / "latex" / "figures",
        RAIZ / "sitio" / "figuras",
    ]
    for funcion in (_f1, _f2, _f3, _f4, _f5, _f6, _f7, _f8, _f9, _f10, _f11, _f12, _f13, _f14):
        funcion(base, sello, destinos)
    volumen = ajustes.resolve("results/v3_3/report/figuras")
    if volumen.resolve() != (base / "report" / "figuras").resolve():
        volumen.mkdir(parents=True, exist_ok=True)
        for fichero in (base / "report" / "figuras").iterdir():
            if fichero.name.startswith("._"):
                continue
            shutil.copy2(fichero, volumen / fichero.name)
    print("[v3.3] 14 figuras")


if __name__ == "__main__":
    main()
