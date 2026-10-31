"""Reglas del cierre de tesis. No fusiona cohortes y no confirma péptidos mutados.

Una matriz de grupos proteicos dice si el gen se detectó y en qué dirección.
No dice si un péptido mutado estuvo en la sangre. Esa distinción es el límite
de lo que esta tesis puede afirmar sobre neoantígenos.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from oncogo_crc.utils.reproducibility_overlap import uniprot_base

CATALOGO = Path(__file__).resolve().parents[3] / "metadata" / "tesis" / "candidatas_literatura.tsv"

# Péptidos que la literatura de 2025-2026 pone en sangre o en el tumor.
# nivel_en_esta_tesis dice qué se puede hacer con las matrices ya cuantificadas.
NEOANTIGENOS: list[dict[str, str]] = [
    {
        "id": "rai2025_19_peptidos",
        "fuente": "Rai et al. J Proteome Res. 2025;24(8):4166-4180. doi:10.1021/acs.jproteome.5c00285",
        "verificacion": "verificado_2026-10-31",
        "que_reporta": "19 péptidos descritos como neoantígenos específicos del suero de cáncer colorrectal en PXD055218.",
        "que_se_puede_hacer_aqui": (
            "Auditar la proteína BAIAP3 y el diseño de tres grupos. "
            "No se puede confirmar ni rechazar los 19 péptidos: la matriz de grupos proteicos no conserva la secuencia mutada."
        ),
        "nivel_en_esta_tesis": "no_evaluable_en_matriz_proteica",
        "origen": "literatura",
    },
    {
        "id": "rai2025_prm_pxd065161",
        "fuente": "El mismo artículo. Validación PRM en PXD065161, n=10 por grupo.",
        "verificacion": "verificado_2026-10-31",
        "que_reporta": "El método PRM del suplemento cubre cuatro péptidos y uno de ellos es el neoantígeno. Los otros 18 quedan sin esa validación dirigida.",
        "que_se_puede_hacer_aqui": "Citar el límite 1 de 19. Rebuscar los crudos queda fuera de esta ejecución.",
        "nivel_en_esta_tesis": "cita_del_limite_no_reanalisis",
        "origen": "literatura",
    },
    {
        "id": "xiang2026_inmunopeptidoma_tumoral",
        "fuente": "Xiang et al. Cell Genomics. 2026;6:101062. doi:10.1016/j.xgen.2025.101062",
        "verificacion": "verificado_2026-10-31",
        "que_reporta": "En 10 pares tumor/normal, 96 neoepítopos mutados presentados por MHC-I. El 80,21 % sale de regiones no codificantes. Es inmunopeptidoma de tejido, no suero.",
        "que_se_puede_hacer_aqui": (
            "Usarlo como frontera: un neoantígeno de MHC en el tumor no es una proteína diferencial de suero. "
            "No se cruzan sus secuencias con CPTAC ni con PXD055218."
        ),
        "nivel_en_esta_tesis": "frontera_tejido_no_sangre",
        "origen": "literatura",
    },
    {
        "id": "hembrough2026_ev_plasma",
        "fuente": "Hembrough et al. Cancer Res. 2026;86(7 Suppl):Abstract 7710. doi:10.1158/1538-7445.am2026-7710",
        "verificacion": "verificado_2026-10-31",
        "que_reporta": "Resumen de congreso. De un paciente con cáncer colorrectal avanzado recuperan 5 péptidos de MHC-I en vesículas de 0,5 mL de plasma. Son péptidos de esa persona, no un panel diagnóstico.",
        "que_se_puede_hacer_aqui": "No se replica: no hay secuencias públicas ni cohorte. Sirve para decir que la detección en plasma es posible y específica de paciente.",
        "nivel_en_esta_tesis": "no_replicable_con_estos_datos",
        "origen": "literatura",
    },
]


def clasificar_corrida(sample_id: str) -> str | None:
    """Grupo clínico codificado en el nombre de corrida de PXD055218."""
    stem = str(sample_id).replace("\\", "/").split("/")[-1].split(" ")[0]
    if stem.startswith("CP_"):
        return "crc"
    if stem.startswith("CN_"):
        return "colonoscopy_negative"
    if stem.startswith("Normal_"):
        return "healthy"
    return None


def etiquetas_son_clinicas(metadata: pd.DataFrame | None) -> bool:
    """El corte alfabético de PXD036594 no es un diagnóstico."""
    if metadata is None or metadata.empty or "condition" not in metadata.columns:
        return False
    batch = metadata["batch"].astype(str) if "batch" in metadata.columns else pd.Series("", index=metadata.index)
    if batch.str.contains("decile_split|sin_etiqueta_clinica", regex=True).any():
        return False
    condiciones = set(metadata["condition"].astype(str))
    return bool(condiciones & {"crc", "healthy", "colonoscopy_negative"})


def grupos_de_muestras(columnas: list[str], metadata: pd.DataFrame | None = None) -> dict[str, list[str]]:
    """Prefiere el nombre CP_/CN_/Normal_. Si no está, usa metadata clínica."""
    grupos: dict[str, list[str]] = {"crc": [], "colonoscopy_negative": [], "healthy": []}
    for col in columnas:
        grupo = clasificar_corrida(col)
        if grupo is not None:
            grupos[grupo].append(col)
    if any(grupos.values()):
        return grupos
    if metadata is None or not etiquetas_son_clinicas(metadata):
        return grupos
    labels = metadata.set_index("sample_id")["condition"].astype(str)
    for col in columnas:
        grupo = labels.get(col)
        if grupo in grupos:
            grupos[grupo].append(col)
    return grupos


def a_escala_log2(bloque: np.ndarray) -> tuple[np.ndarray, str]:
    """Intensidad lineal pasa a log2. Un bloque que ya tiene negativos se deja."""
    valores = np.asarray(bloque, dtype=float)
    finitos = valores[np.isfinite(valores)]
    if finitos.size == 0:
        return valores, "vacio"
    if np.nanmin(finitos) < 0 or np.nanmedian(finitos) <= 64:
        return valores, "ya_en_log_o_razon"
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(valores > 0, np.log2(valores), np.nan), "log2_de_intensidad_lineal"


def auc_dirigido(positivo: np.ndarray, negativo: np.ndarray) -> float | None:
    """AUC con la primera clase como caso. None si falta un grupo o no hay valores."""
    a = np.asarray(positivo, dtype=float)
    b = np.asarray(negativo, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return None
    y = np.array([1] * len(a) + [0] * len(b))
    puntuacion = np.concatenate([a, b])
    if len(np.unique(puntuacion)) < 2:
        return 0.5
    return float(roc_auc_score(y, puntuacion))


def auc_discriminacion(positivo: np.ndarray, negativo: np.ndarray) -> float | None:
    """AUC sin fijar el signo: max(AUC, 1-AUC). Es el uso habitual en un marcador."""
    dirigido = auc_dirigido(positivo, negativo)
    if dirigido is None:
        return None
    return float(max(dirigido, 1.0 - dirigido))


def intervalo_auc(
    positivo: np.ndarray,
    negativo: np.ndarray,
    *,
    n_boot: int = 2000,
    semilla: int = 20261031,
) -> tuple[float | None, float | None, float | None]:
    """Bootstrap estratificado del AUC de discriminación. Devuelve punto, bajo y alto."""
    punto = auc_discriminacion(positivo, negativo)
    a = np.asarray(positivo, dtype=float)
    b = np.asarray(negativo, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if punto is None or n_boot < 1:
        return punto, None, None
    rng = np.random.default_rng(semilla)
    muestras = []
    for _ in range(n_boot):
        aa = a[rng.integers(0, len(a), len(a))]
        bb = b[rng.integers(0, len(b), len(b))]
        valor = auc_discriminacion(aa, bb)
        if valor is not None:
            muestras.append(valor)
    if len(muestras) < 20:
        return punto, None, None
    bajo, alto = np.percentile(muestras, [2.5, 97.5])
    return punto, float(bajo), float(alto)


# Tanaka et al., Cell Reports 2024, doi:10.1016/j.celrep.2024.113810.
# La cohorte 1 de PXD034575 son 16 tripletes: metástasis hepática, tumor
# primario y mucosa colónica benigna. El sufijo P no es un pólipo.
_LEYENDA_TANAKA = {
    "M": "metastasis_hepatica",
    "N": "mucosa_normal",
    "P": "tumor_primario",
}
_RE_TANAKA = re.compile(r"^LFQ intensity CRC(\d+)([MNP])$", re.IGNORECASE)


def etiqueta_tanaka(columna: str) -> tuple[str, str] | None:
    """Paciente y tejido real de una columna LFQ de PXD034575.

    Devuelve (paciente, tejido). None si la columna no es de ese depósito.
    """
    coincidencia = _RE_TANAKA.match(str(columna).strip())
    if coincidencia is None:
        return None
    return coincidencia.group(1), _LEYENDA_TANAKA[coincidencia.group(2).upper()]


def lectura_con_ruido(veredicto: str, lista_ruido: object) -> str:
    """La regla de AUC no se cambia. Una bandera previa impide llamarla específica.

    La bandera se escribió en el catálogo antes del contraste. Si el AUC
    cumple la regla y la proteína ya estaba en fase aguda, hemólisis o
    enfermedad inflamatoria, la lectura lo dice y no la presenta como marcador.
    """
    ruido = "" if lista_ruido is None or (isinstance(lista_ruido, float) and pd.isna(lista_ruido)) else str(lista_ruido).strip()
    if ruido.lower() in {"", "nan"}:
        ruido = ""
    if veredicto == "especifica" and ruido:
        return f"auc_pasa_pero_{ruido}"
    return veredicto


def puntuacion_enriquecimiento(puntajes_descendentes: np.ndarray, en_conjunto: np.ndarray) -> float:
    """Enrichment score ponderado de Subramanian, con el ranking ya ordenado."""
    puntajes = np.asarray(puntajes_descendentes, dtype=float)
    mascara = np.asarray(en_conjunto, dtype=bool)
    if puntajes.shape != mascara.shape:
        raise ValueError("El ranking y la máscara del término tienen que medir lo mismo.")
    n = int(puntajes.size)
    n_en = int(mascara.sum())
    if n_en < 10 or n_en >= n:
        return float("nan")
    pesos = np.where(mascara, np.abs(puntajes), 0.0)
    suma = float(pesos.sum())
    if suma <= 0:
        return float("nan")
    paso_dentro = pesos / suma
    paso_fuera = np.where(mascara, 0.0, 1.0 / (n - n_en))
    caminata = np.cumsum(paso_dentro - paso_fuera)
    return float(caminata[int(np.argmax(np.abs(caminata)))])


def veredicto_especificidad(
    auc_crc_cn: float | None,
    ci_bajo_crc_cn: float | None,
    auc_cn_sano: float | None,
) -> str:
    """Regla fijada antes de mirar el suero de tres grupos.

    Específica si separa el cáncer de la otra enfermedad intestinal, el límite
    inferior supera el azar y esa separación es mayor que la de la enfermedad
    frente a los sanos. Adaptada de Ögren 2026 al diseño de Rai 2025.
    """
    if auc_crc_cn is None or ci_bajo_crc_cn is None or auc_cn_sano is None:
        return "no_evaluable"
    if auc_crc_cn > 0.70 and ci_bajo_crc_cn > 0.5 and auc_crc_cn > auc_cn_sano:
        return "especifica"
    if auc_cn_sano >= auc_crc_cn and auc_cn_sano > 0.70:
        return "compartida_con_enfermedad"
    return "inespecifica"


def cargar_catalogo(ruta: Path | None = None) -> pd.DataFrame:
    path = ruta or CATALOGO
    tabla = pd.read_csv(path, sep="\t")
    if not (tabla["origen"] == "literatura").all():
        raise ValueError("El catálogo de candidatas solo admite origen=literatura.")
    return tabla


def _tokens_id(protein_id: str) -> set[str]:
    piezas = re.split(r"[;,\s|]+", str(protein_id).upper())
    tokens = set()
    for pieza in piezas:
        if not pieza:
            continue
        tokens.add(pieza.split("-")[0])
        base = uniprot_base(pieza)
        if base:
            tokens.add(base)
    return tokens


def localizar_candidata(identificadores: pd.Series, gen: str, uniprot: str | None) -> object | None:
    """Índice de la fila cuyo id es el gen o contiene su accesión UniProt."""
    gen_u = gen.upper()
    uni = uniprot.upper() if uniprot else ""
    for idx, protein_id in identificadores.items():
        tokens = _tokens_id(protein_id)
        if gen_u in tokens or (uni and uni in tokens):
            return idx
    return None


def cruzar_catalogo(
    catalogo: pd.DataFrame,
    contraste: pd.DataFrame,
    *,
    cohorte: str,
    mapa_uniprot: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Une cada candidata con un contraste. No inventa ausencias como especificidad."""
    tabla = contraste.copy()
    mapa = {k.upper(): v.upper() for k, v in (mapa_uniprot or {}).items()}
    filas = []
    for _, cand in catalogo.iterrows():
        gen = str(cand["gene"]).upper()
        idx = localizar_candidata(tabla["protein_id"], gen, mapa.get(gen))
        fila = {
            "gene": cand["gene"],
            "origen": cand["origen"],
            "fuente": cand["fuente"],
            "verificacion": cand["verificacion"],
            "lista_ruido": cand.get("lista_ruido", ""),
            "direccion_literatura": cand["direccion_literatura"],
            "cohorte": cohorte,
        }
        if idx is None:
            fila.update(
                {
                    "detectada": False,
                    "log2fc": np.nan,
                    "padj": np.nan,
                    "significativa": False,
                    "veredicto": "no_evaluable",
                }
            )
        else:
            r = tabla.loc[idx]
            log2fc = float(r["log2fc"]) if pd.notna(r["log2fc"]) else np.nan
            esperada = str(cand["direccion_literatura"])
            signo_ok = None
            if np.isfinite(log2fc) and log2fc != 0 and esperada in {"sube", "baja"}:
                signo_ok = (log2fc > 0 and esperada == "sube") or (log2fc < 0 and esperada == "baja")
            significativa = bool(r["significant"]) if "significant" in r.index else False
            if esperada == "opuesta":
                veredicto = "direccion_no_unica_en_la_fuente"
            elif not np.isfinite(log2fc):
                veredicto = "no_evaluable"
            elif significativa and signo_ok:
                veredicto = "replicada"
            elif significativa and signo_ok is False:
                veredicto = "significativa_en_sentido_opuesto"
            else:
                veredicto = "no_replicada"
            fila.update(
                {
                    "detectada": True,
                    "log2fc": log2fc,
                    "padj": float(r["padj"]) if pd.notna(r.get("padj", np.nan)) else np.nan,
                    "significativa": significativa,
                    "veredicto": veredicto,
                }
            )
        filas.append(fila)
    return pd.DataFrame(filas)


def tabla_neoantigenos() -> pd.DataFrame:
    return pd.DataFrame(NEOANTIGENOS)
