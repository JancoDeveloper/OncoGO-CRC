"""Resumen, abstract, 18 diapositivas y el banco de preguntas. Las cifras salen de results/v3_3/."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parents[2]
BASE = RAIZ / "results" / "v3_3"


def _json(nombre: str) -> dict:
    return json.loads((BASE / nombre).read_text(encoding="utf-8"))


def _coma(valor: float, digitos: int = 3) -> str:
    return f"{valor:.{digitos}f}".replace(".", ",")


def _palabras(texto: str) -> int:
    limpio = []
    for linea in texto.splitlines():
        if linea.startswith("\\") or linea.startswith("%") or linea.startswith("{"):
            continue
        limpio.append(linea.replace("\\,", " "))
    return len(" ".join(limpio).split())


def _datos() -> dict:
    certificado = _json("certificado_v3_3.json")
    suero = _json("suero_resumen.json")
    k1 = _json("k1_etiquetas.json")
    ml = _json("ml_suero_permutacion.json")["versiones"]
    fuga = _json("ml_demostracion_fuga.json")
    principal = suero["clases"]["principal"]["crc_frente_a_colonoscopia_negativa"]
    sanos = suero["clases"]["principal"]["crc_frente_a_sanos"]
    r1 = k1["contrastes"]["primario_frente_a_mucosa"]
    met = k1["contrastes"]["metastasis_hepatica_frente_a_mucosa"]
    k4 = pd.read_csv(BASE / "k4_crc_frente_a_sanos.tsv", sep="\t")
    traza = pd.read_csv(BASE / "suero_trazabilidad_grupos.tsv", sep="\t")
    genes = dict(zip(traza["protein_group"].astype(str), traza["genes"].astype(str), strict=False))
    distintas = k4.loc[k4["clase"].eq("diferente"), ["protein_id", "log2fc", "padj", "n_obs"]]
    nombres = [f"{genes.get(fila.protein_id, fila.protein_id)} ({fila.protein_id})" for fila in distintas.itertuples()]
    return {
        "frase": certificado["frase_suero"],
        "titular": certificado["titular_ml"]["texto"],
        "clave": certificado["titular_ml"]["clave"],
        "n": principal["n"],
        "diferente": principal["diferente"],
        "equivalente": principal["equivalente"],
        "indeterminada": principal["indeterminada"],
        "sanos_n": sanos["diferente"],
        "sanos_nombres": ", ".join(nombres),
        "r1": f"{r1['numerador']}/{r1['denominador']}",
        "ci": f"{_coma(r1['ci_low'])} a {_coma(r1['ci_high'], 0) if r1['ci_high'] == 1 else _coma(r1['ci_high'])}",
        "lin": _coma(r1["lin"]["ccc"]),
        "met": f"{met['numerador']}/{met['denominador']}",
        "lin_met": _coma(met["lin"]["ccc"]),
        "excluidas": k1["n_filas_excluidas_contaminante_o_decoy"],
        "certificadas": k1["n_certificadas"],
        "grupos": suero["n_grupos_report"],
        "humanos": suero["n_matriz_humana_fdr"],
        "h5": certificado["h5"]["estado"],
        "h7": certificado["h7"]["estado"],
        "h10": certificado["h10"]["estado"],
        "especificas": certificado["h10"]["n_especificas_crc"],
        "fragil": certificado["h1_fragil_a_normalizacion"],
        "m2": _coma(ml["M2_proteinas_principal"]["auc_oof_media"]),
        "p2": _coma(ml["M2_proteinas_principal"]["p_permutacion"], 6),
        "m5": _coma(ml["M5_principal_residualizada_qc"]["auc_oof_media"]),
        "p5": _coma(ml["M5_principal_residualizada_qc"]["p_permutacion"], 3),
        "nulo": _coma(ml["M2_proteinas_principal"]["media_del_nulo"]),
        "fuga_sin": _coma(fuga["auc_sin_duplicar"]),
        "fuga_con": _coma(fuga["auc_con_filas_duplicadas"]),
        "pan": _coma(certificado["h5"]["fraccion_pan_carcinoma"]),
        "auc_alvez": _coma(certificado["h7"]["por_recurso"]["alvez_2023"]["auc"]),
        "auc_hpa": _coma(certificado["h7"]["por_recurso"]["hpa_disease_blood_atlas"]["auc"]),
    }


def _resumen(d: dict) -> str:
    return f"""\\chapter*{{Resumen}}
\\addcontentsline{{toc}}{{chapter}}{{Resumen}}
% fuente: results/v3_3/certificado_v3_3.json
% fuente: results/v3_3/suero_resumen.json
% fuente: results/v3_3/k1_etiquetas.json
% fuente: results/v3_3/k4_crc_frente_a_sanos.tsv
El cambio proteico del tumor colorrectal ya está descrito en el tejido. La corrida v3.3 pregunta si ese cambio llega a la sangre y si el suero separa el cáncer de la colonoscopia negativa. No se fusionaron cohortes. En el contraste no se imputaron ausencias. El ancla son las {d["certificadas"]} proteínas certificadas de CPTAC. La réplica es PXD034575, con 16 pacientes. P es el tumor primario, M la metástasis hepática y N la mucosa.

Se excluyen {d["excluidas"]} filas de contaminante o señuelo y se piden al menos ocho pares. El primario frente a la mucosa coincide en signo en {d["r1"]}, con intervalo de {d["ci"]} y concordancia de Lin de {d["lin"]}. La metástasis frente a la mucosa coincide en {d["met"]}, con Lin de {d["lin_met"]}. CA1 baja en el tumor. AQP8, SLC26A3 y GUCA2A no tienen log2FC, así que las etiquetas no quedan confirmadas.

El suero es PXD055218, con 20 personas por grupo. El informe trae {d["grupos"]} grupos. Con la tasa de descubrimientos falsos del 1\\,\\% y la especie humana quedan {d["humanos"]}. El contraste principal tiene {d["n"]} proteínas: {d["diferente"]} diferentes, {d["equivalente"]} equivalentes y {d["indeterminada"]} indeterminadas frente a la colonoscopia negativa. Frente a personas sanas se separan {d["sanos_nombres"]}. {d["titular"]} H5 queda en {d["h5"]} y H7 en {d["h7"]}. H10 es {d["h10"]} y cuenta {d["especificas"]} proteínas específicas. El trabajo no propone un examen de sangre.

\\noindent\\textbf{{Palabras clave:}} cáncer colorrectal; proteómica; suero; especificidad; calidad de muestra.
\\clearpage
"""


def _abstract(d: dict) -> str:
    return f"""\\begin{{otherlanguage}}{{english}}
\\chapter*{{Abstract}}
\\addcontentsline{{toc}}{{chapter}}{{Abstract}}
% fuente: results/v3_3/certificado_v3_3.json
% fuente: results/v3_3/suero_resumen.json
% fuente: results/v3_3/k1_etiquetas.json
The protein change of the colorectal tumor is already described in tissue. Run v3.3 asks whether that change reaches the blood and whether serum separates cancer from a negative colonoscopy. Cohorts were not merged. Missing values were not imputed in the contrast. The anchor is the {d["certificadas"]} certified CPTAC proteins. The replicate is PXD034575, 16 patients. P is the primary tumor, M the liver metastasis and N the mucosa.

The run drops {d["excluidas"]} contaminant or decoy rows and requires at least eight pairs. Primary versus mucosa agrees in sign for {d["r1"]}, with an interval from {d["ci"]} and Lin concordance {d["lin"]}. Metastasis versus mucosa agrees for {d["met"]}, with Lin concordance {d["lin_met"]}. CA1 falls in the tumor. AQP8, SLC26A3 and GUCA2A have no log2 fold change, so the tissue labels are not confirmed.

Serum is PXD055218, 20 people per group. The report has {d["grupos"]} groups. A 1\\% protein false-discovery rate and the human species filter leave {d["humanos"]}. The principal contrast has {d["n"]} proteins: {d["diferente"]} different, {d["equivalente"]} equivalent and {d["indeterminada"]} indeterminate against the negative colonoscopy. Against healthy people the separated proteins are {d["sanos_nombres"]}. The separation seen by the model is explained by sample quality. H5 is {d["h5"]} and H7 is {d["h7"]}. H10 is {d["h10"]} and counts {d["especificas"]} specific proteins. The work does not propose a blood test.

\\noindent\\textbf{{Keywords:}} colorectal cancer; proteomics; serum; specificity; sample quality.
\\end{{otherlanguage}}
\\clearpage
"""


def _diapositivas(d: dict) -> str:
    marcos = [
        ("Pregunta", "¿El suero separa el cáncer colorrectal de la colonoscopia negativa, más allá de la calidad de la muestra?"),
        ("Diseño", "No se fusionan cohortes. No se imputan ausencias en el contraste. P es tumor primario, M metástasis hepática y N mucosa."),
        ("Ancla y réplica", f"CPTAC aporta {d['certificadas']} proteínas certificadas. La réplica es PXD034575, 16 pacientes. Se excluyen {d['excluidas']} filas."),
        ("Controles del tejido", "CA1 baja en el tumor y pasa. AQP8, SLC26A3 y GUCA2A no tienen log2FC. Las etiquetas no quedan confirmadas."),
        ("Concordancia", f"Primario frente a mucosa: {d['r1']}. Intervalo {d['ci']}. Lin {d['lin']}. H1: rechaza\\_h1. Frágil a la normalización: {d['fragil']}."),
        ("Metástasis", f"Metástasis hepática frente a mucosa: {d['met']}. Lin {d['lin_met']}."),
        ("Filtro del suero", f"PXD055218, 20 personas por grupo. {d['grupos']} grupos en el informe. {d['humanos']} humanos tras el FDR. Matriz principal: {d['n']}."),
        ("Contraste principal", f"{d['frase']} Son {d['diferente']} diferentes, {d['equivalente']} equivalentes y {d['indeterminada']} indeterminadas."),
        ("Frente a personas sanas", f"Se separan {d['sanos_n']} proteínas: {d['sanos_nombres']}. Esa no es la pregunta clínica."),
        ("Modelo", f"{d['titular']} Clave: {d['clave']}."),
        ("Las dos AUC", f"Con las proteínas principales, AUC fuera de pliegue {d['m2']} y p {d['p2']}. Al quitar la calidad, AUC {d['m5']} y p {d['p5']}. Media del nulo de M2: {d['nulo']}."),
        ("La fuga", f"En ruido, el AUC es {d['fuga_sin']} sin duplicar filas y {d['fuga_con']} al duplicarlas. Los AUC de v3.2 quedan retirados."),
        ("H5", f"Estado: {d['h5']}. Fracción pan-carcinoma: {d['pan']}. UCEC sigue fuera: no comparte identificadores."),
        ("H7", f"Estado: {d['h7']}. AUC {d['auc_alvez']} en Álvez 2023 y {d['auc_hpa']} en el atlas del HPA. El intervalo incluye 0,50. No es una prueba."),
        ("H10", f"Estado: {d['h10']}. Proteínas específicas de cáncer colorrectal: {d['especificas']}. Fase 1 de Pepe."),
        ("Lo que no cierra", "AQP8 y GUCA2A no están en MaxQuant. SLC26A3 no tiene log2FC. El multiverso calculado es la rejilla K4. ESM-2 y TMbed no se ejecutaron."),
        ("Advertencia", "Ningún resultado de esta tesis sirve para decidir sobre un paciente."),
        ("Cierre", "La corrida v3.3 está en results/v3\\_3. La publicación espera a que el checklist quede en verde."),
    ]
    cuerpo = []
    for titulo, texto in marcos:
        cuerpo.append("\\begin{frame}{" + titulo + "}\n" + texto + "\n\\end{frame}\n")
    return (
        "\\documentclass{beamer}\n\\usepackage[spanish]{babel}\n\\usepackage[utf8]{inputenc}\n"
        "\\title{OncoGO-CRC, corrida v3.3}\n\\author{Juan Diego Polanco Franco}\n\\date{9 de noviembre de 2026}\n"
        "\\begin{document}\n"
        + "\n".join(cuerpo)
        + "\\end{document}\n"
    )


def _preguntas(d: dict) -> str:
    pares = [
        ("¿Cuál es la pregunta de la corrida v3.3?", "Si el suero separa el cáncer de la colonoscopia negativa más allá de la calidad de muestra. Fuente: certificado, campo titular_ml."),
        ("¿Cuántas proteínas tiene el contraste principal?", f"{d['n']}. Fuente: suero_resumen.json, n_principal."),
        ("¿Cuántas son diferentes, equivalentes e indeterminadas?", f"{d['diferente']}, {d['equivalente']} y {d['indeterminada']}. {d['frase']}"),
        ("¿Qué dice el titular del modelo?", d["titular"]),
        ("¿Por qué el titular no es que el suero separa el cáncer?", f"M2 tiene AUC {d['m2']} y p {d['p2']}, pero M5, ya sin la calidad, tiene AUC {d['m5']} y p {d['p5']}. La regla pide las dos p menores que 0,05 y un AUC de M5 de al menos 0,70."),
        ("¿Qué pasó con los AUC 0,996, 0,981 y 0,819 de v3.2?", "Quedan retirados. Esas corridas duplicaban cada muestra. La demostración con ruido pasa de " + d["fuga_sin"] + " a " + d["fuga_con"] + "."),
        ("¿La media del nulo cae donde debía?", f"La de M2 es {d['nulo']}. La regla pedía entre 0,45 y 0,55."),
        ("¿Cuántas filas entran al modelo?", "40: 20 de cáncer y 20 de colonoscopia negativa. Fuente: ml_suero_permutacion.json, n_por_clase."),
        ("¿H1 se sostiene?", f"No. El estado es rechaza_h1. La concordancia es {d['r1']}, y la regla solo sostiene H1 si el límite alto del intervalo queda bajo 0,70."),
        ("¿Esa concordancia es frágil al centrado?", f"No. h1_fragil_a_normalizacion es {d['fragil']}."),
        ("¿Las etiquetas P y N quedaron confirmadas?", "No. CA1 pasa. AQP8, SLC26A3 y GUCA2A fallan. Fuente: controles.tsv."),
        ("¿Por qué fallan AQP8 y GUCA2A?", "No están en la salida de MaxQuant. No se inventa una fila."),
        ("¿Por qué falla SLC26A3?", "Tiene una fila y el log2FC no es finito. El control lo marca no evaluable y la etapa lo cuenta como falla."),
        ("¿Qué es la columna P?", "El tumor primario de PXD034575. No es un pólipo."),
        ("¿Cuántos pares tiene el tejido?", "16 de primario frente a mucosa y 16 de metástasis frente a mucosa. El mínimo sigue en 8."),
        ("¿Cuántas filas se excluyen del tejido?", f"{d['excluidas']}, por CON__, REV__ y las marcas de MaxQuant."),
        ("¿Cuál es la concordancia de la metástasis?", f"{d['met']}, con Lin {d['lin_met']}."),
        ("¿Cuántos grupos trae el informe de suero?", f"{d['grupos']}. Tras el FDR de proteína y la especie humana quedan {d['humanos']}."),
        ("¿Por qué la matriz principal no es 312?", f"312 es la matriz humana con FDR. El contraste principal quita contaminantes de piel y paneles y se queda en {d['n']}."),
        ("¿Qué proteínas separan el cáncer de las personas sanas?", d["sanos_nombres"] + ". Fuente: k4_crc_frente_a_sanos.tsv, clase diferente, cruzada con la traza."),
        ("¿Esa separación frente a sanos es la pregunta clínica?", "No. La colonoscopia negativa reúne otras enfermedades intestinales. Comparar solo con sanos infla la señal inflamatoria."),
        ("¿SAA1 o LBP son un examen?", "No. Se separan de las personas sanas y el contraste principal, frente a colonoscopia negativa, no tiene proteínas diferentes."),
        ("¿H5 qué dice?", f"Estado {d['h5']}. La fracción pan-carcinoma es {d['pan']}."),
        ("¿Por qué UCEC no entra en H5?", "Sus columnas normales se llaman NX y no comparten identificador con el tumor. El mínimo de ocho pares no se bajó."),
        ("¿H7 qué dice?", f"Estado {d['h7']}. AUC {d['auc_alvez']} en Álvez 2023 y {d['auc_hpa']} en el atlas del HPA. Los dos intervalos incluyen 0,50."),
        ("¿H7 es una prueba diagnóstica?", "No. Es un modelo de fase 1, fuera de recurso y fuera de proteína. El certificado lo dice."),
        ("¿Qué le falta a H7?", "ESM-2, TMbed, Papier 2024 y otras variables listadas en variables_faltantes. No se sustituyen."),
        ("¿H10 qué dice?", f"Estado {d['h10']}. Cuenta {d['especificas']} proteínas específicas. No hay umbral de sostiene o refuta."),
        ("¿Por qué H10 no se declara refutada?", "El preregistro no fijó ese umbral. Fijarlo después de ver el atlas sería elegirlo por el resultado."),
        ("¿El multiverso ya está completo?", "No. Lo calculado es la rejilla K4, 48 celdas. K5, el modelo y la normalización del tejido no se repitieron en cada celda."),
        ("¿Se puede publicar ya?", "No. El checklist de la sección 16 no está en verde. No hay push ni Firebase."),
        ("¿Se imputó el suero para el contraste?", "No. El cero de MaxQuant y la ausencia de DIA-NN se quedan como ausencia. La mediana del modelo solo se calcula dentro del pliegue de entrenamiento."),
        ("¿Se ejecutó DIA-NN otra vez?", "No. Se usa el report ya hecho, de 213 313 241 bytes."),
        ("¿Cuántas permutaciones tiene K5?", "10 000, en tres matrices. Fuente: k5_auc.json."),
        ("¿Cuántas permutaciones tiene el modelo?", "1 000 por modelo. La p no puede ser menor que 1/1001."),
        ("¿Qué semilla se usó?", "La base 20260914. Los desplazamientos de v3.3 están en la enmienda, no se eligieron después de ver el resultado."),
        ("¿La albúmina humana quedó marcada como contaminante?", "No. El cruce es por acceso y especie. La prueba de plasma clásico pasa."),
        ("¿Cuántas entradas tiene la biblioteca de contaminantes?", "381. Fuente: suero_controles.json, libreria_n_entradas."),
        ("¿Quién falla el control de sexo?", "CRC15 y CRC16. Fuente: qc_sexo_por_participante.tsv."),
        ("¿Quién falla el hígado?", "CRC13, CRC16 y CRC17. Fuente: qc_etiquetas_higado.tsv."),
        ("¿La mucosa de los 16 participantes es evaluable?", "Sí, y los 16 pasan. Fuente: qc_etiquetas_mucosa.tsv."),
        ("¿A quién saca la sensibilidad de etiquetas?", "A CRC16, el único que falla dos de tres controles. Fuente: k1_etiquetas.json."),
        ("¿v3 y v3.2 se reescribieron?", "No. Siguen congelados por sus manifiestos. v3.3 escribe solo en results/v3_3."),
        ("¿De dónde sale cada cifra de la página?", "De scripts/v3/exportar_sitio.py, que lee el certificado y el changelog."),
        ("¿De dónde sale cada cifra de la sección nueva de la tesis?", "De cifras.tex, escrito por scripts/v3/exportar_cifras_tesis.py."),
        ("¿Qué fase de Pepe es esta?", "La fase 1. Exploración en muestras ya recogidas. No hay una cohorte de cribado."),
        ("¿El intervalo de los pliegues es un intervalo de confianza?", "No. La dispersión entre pliegues se informa aparte. El intervalo de LeDell es el que lleva ese nombre en el JSON."),
        ("¿CHMP4A volvió a entrar?", "Está en la traza, acceso Q9BY43;Q9BY43-2, con Global.PG.Q 0,0374 y fracción de péptidos proteotípicos 0. No está entre las proteínas diferentes del contraste principal."),
        ("¿Qué no se debe decir en la sustentación?", "Biomarcador validado, prueba diagnóstica, y pólipo para la columna P."),
        ("¿Qué queda abierto a propósito?", "Los tres controles sin log2FC, el multiverso fuera de la rejilla K4, y los modelos de secuencia que no se ejecutaron."),
        ("¿El resumen de v3.2 sigue en la tesis?", "Sí, en el capítulo de resultados, marcado como corrida del 6 de noviembre. La sección vigente es la de v3.3."),
        ("¿Cuántas diapositivas hay?", "18, en tesis-ucaldas/sustentacion/diapositivas.tex, escritas por este mismo guion."),
    ]
    lineas = [
        "# Banco de preguntas, corrida v3.3",
        "",
        "Cada cifra sale de `results/v3_3/`. Este archivo lo escribe `scripts/v3/escribir_documentacion_v3_3.py`. El banco del 6 de noviembre describe v3.2 y no es la corrida vigente.",
        "",
    ]
    for i, (pregunta, respuesta) in enumerate(pares, start=1):
        lineas.append(f"{i}. **{pregunta}**")
        lineas.append(respuesta)
        lineas.append("")
    return "\n".join(lineas)


def main() -> None:
    datos = _datos()
    resumen = _resumen(datos)
    abstract = _abstract(datos)
    n_resumen = _palabras(resumen)
    n_abstract = _palabras(abstract)
    if n_resumen > 250 or n_abstract > 250:
        raise SystemExit(f"palabras resumen={n_resumen} abstract={n_abstract}")
    (RAIZ / "tesis-ucaldas/latex/chapters/03-resumen.tex").write_text(resumen, encoding="utf-8")
    (RAIZ / "tesis-ucaldas/latex/chapters/04-abstract.tex").write_text(abstract, encoding="utf-8")
    destino = RAIZ / "tesis-ucaldas/sustentacion"
    destino.mkdir(parents=True, exist_ok=True)
    (destino / "diapositivas.tex").write_text(_diapositivas(datos), encoding="utf-8")
    preguntas = _preguntas(datos)
    (RAIZ / "docs/sustentacion_v3_3.md").write_text(preguntas, encoding="utf-8")
    print(f"[v3.3] resumen {n_resumen} palabras, abstract {n_abstract}, preguntas {preguntas.count('**')}")


if __name__ == "__main__":
    main()
