"""D1. Determinantes del CEA sérico en TCGA COAD/READ.

Esta etapa no toca nada de las cohortes de proteómica. Lee la tabla limpia que
escribió ``scripts/download/tcga_coadread.py`` y responde dos preguntas
preregistradas: qué explica el CEA en sangre, y cómo son los tumores que
expresan mucho CEACAM5 y no lo elevan.

Si la descarga no se hizo o quedó marcada ``no_disponible``, la etapa lo dice y
no calcula nada. No hay valor por defecto ni dato de relleno.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import Settings
from oncogo_crc.utils.io import write_tsv
from oncogo_crc.v2 import tcga_cea as T
from oncogo_crc.v2.prereg import cargar


def _tarjeta(base: Path) -> dict | None:
    path = base / "estado.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def run(settings: Settings) -> None:
    prereg = cargar(settings.root)
    fuente = settings.resolve("external") / "tcga_coadread"
    salida = settings.resolve("results") / "v2" / "axes" / "tcga_cea"
    salida.mkdir(parents=True, exist_ok=True)

    tarjeta = _tarjeta(fuente)
    if tarjeta is None or tarjeta.get("estado") != "usada":
        motivo = (
            "no existe external/tcga_coadread/estado.json; hay que ejecutar "
            "scripts/download/tcga_coadread.py"
            if tarjeta is None
            else f"la fuente está marcada {tarjeta.get('estado')}: "
                 f"{tarjeta.get('motivo', 'sin motivo')}"
        )
        auditoria = {"version": "v2", "estado": "no_disponible", "motivo": motivo}
        (salida / "tcga_cea_audit_v2.json").write_text(
            json.dumps(auditoria, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        print(f"[v2-tcga] no_disponible: {motivo}")
        return

    tabla = pd.read_csv(fuente / tarjeta["tabla_limpia"], sep="\t", low_memory=False)
    d = T.preparar(tabla)
    print(f"[v2-tcga] {len(d)} muestras de tumor primario de TCGA COAD/READ.")

    censura = T.censura_del_ensayo(d["cea"])
    print(
        f"[v2-tcga] censura: {censura['n_en_el_valor_minimo']} valores en el mínimo "
        f"de {censura['n']}; Tobit "
        f"{'necesario' if censura['hace_falta_tobit'] else 'innecesario'}."
    )

    modelo = T.modelo_del_cea(d)
    if modelo["estado"] == "calculado":
        print(
            f"[v2-tcga] modelo: n={modelo['n_casos_completos']}, "
            f"R²={modelo['r2']:.4f}, R² solo CEACAM5={modelo['r2_solo_ceacam5']:.4f}, "
            f"aporte de CEACAM5={modelo['aporte_de_ceacam5_al_r2']:.4f}"
        )
        print(f"[v2-tcga] {modelo['lectura']}")
        write_tsv(pd.DataFrame(modelo["coeficientes"]), salida / "modelo_del_cea.tsv")
    else:
        print(f"[v2-tcga] modelo: {modelo['motivo']}")

    sensibilidad = T.sensibilidad_con_el_cero(
        d, suelo=float(censura.get("segundo_valor_mas_bajo") or 0.1)
    )

    no_secretor = T.fenotipo_no_secretor(d)
    if no_secretor["estado"] == "calculado":
        g = no_secretor["grupos"]
        print(
            f"[v2-tcga] fenotipo no secretor: {g['no_secretor']['n']} tumores con "
            f"CEACAM5 alto y CEA normal, frente a {g['secretor']['n']} que sí lo "
            "elevan."
        )
        print(f"[v2-tcga] {no_secretor['lectura']}")
        write_tsv(
            pd.DataFrame(
                [
                    {"contraste": k, **v}
                    for k, v in no_secretor["contrastes_no_secretor_vs_secretor"].items()
                ]
            ),
            salida / "fenotipo_no_secretor_contrastes.tsv",
        )
    else:
        print(f"[v2-tcga] fenotipo no secretor: {no_secretor['motivo']}")

    # La tabla por muestra con lo derivado, para que cada número del informe tenga
    # de dónde salir. Son identificadores públicos de TCGA, no datos de una persona
    # concreta del entorno del autor.
    columnas = [
        c
        for c in (
            "paciente", "muestra", "proyecto", "cea", "log10_cea",
            "log2_CEACAM5", "log2_CEACAM6", "log2_CEACAM1",
            "pureza_absolute", "fraccion_adn_tumoral", "msi", "msi_mantis",
            "lado", "estadio_grupo", "metastasis", "sexo", "edad",
            "subtipo_pancanatlas",
        )
        if c in d.columns
    ]
    write_tsv(d[columnas], salida / "tcga_cea_por_muestra.tsv")

    auditoria = {
        "version": "v2",
        "estado": "calculado",
        "fuente": {
            "nombre": tarjeta["nombre"],
            "fecha_de_acceso": tarjeta["fecha_de_acceso"],
            "cita": tarjeta["cita"],
            "licencia": tarjeta["licencia"],
            "sha256_de_la_tabla_limpia": tarjeta["sha256_de_la_tabla_limpia"],
            "cobertura": tarjeta["cobertura"],
            "cms": tarjeta["cms"],
            "discrepancias_con_el_enunciado": tarjeta["discrepancias_con_el_enunciado"],
        },
        "n_muestras": int(len(d)),
        "umbral_mantis": T.UMBRAL_MANTIS,
        "censura_del_ensayo": censura,
        "modelo_del_cea": modelo,
        "sensibilidad_imputando_el_cero": sensibilidad,
        "fenotipo_no_secretor": no_secretor,
        "relacion_con_c8": (
            "C8 encontró que la dirección del cambio en tejido no se traslada al "
            "plasma. Este eje mide lo mismo en el caso mejor estudiado que existe: "
            "el CEA y su propio gen. Si la expresión tisular explica poco del "
            "marcador en sangre, lo de C8 no es la rareza de un experimento pequeño."
        ),
        "no_entra_en_el_certificado": (
            "Ninguno de estos números es un subíndice. TCGA no tiene proteoma "
            "pareado en esta descarga, así que la pureza medida no se puede usar "
            "todavía para un R4 alternativo: eso depende de D2."
        ),
        "limitaciones": list(tarjeta.get("limitaciones", []))
        + [
            "El lado del tumor se deriva de códigos ICD-O-3 y una parte queda en "
            "«colon sin especificar». Esos casos entran como categoría propia en "
            "lugar de desaparecer.",
            "El fenotipo no secretor se define con un percentil de expresión y el "
            "umbral clínico de 5 ng/ml. Los contrastes que lo describen son "
            "exploratorios y no están corregidos por multiplicidad.",
            "Correlación no es causa: que el CEA suba con el estadio no dice que la "
            "masa tumoral sea la causa, solo que las dos cosas van juntas.",
        ],
        "advertencia_clinica": (
            "Ningún resultado de esta tesis sirve para decidir sobre un paciente."
        ),
        "semilla_base": int(prereg.get("semillas.base")),
        "nota_de_semilla": (
            "Nada de esta etapa es estocástico: son mínimos cuadrados y pruebas "
            "exactas. La semilla se registra por coherencia, no porque se use."
        ),
    }
    (salida / "tcga_cea_audit_v2.json").write_text(
        json.dumps(auditoria, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"[v2-tcga] escrito en {salida}")
