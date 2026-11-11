#!/usr/bin/env python3
"""Descarga la anotación de las proteínas de la celda y deja el manifiesto.

No rellena un campo si la petición falla: queda «no disponible» con el motivo.
Los modelos se alojan solo si la suma de BinaryCIF queda bajo el corte de
docs/sitio_v4_diseno.md (80 MB), fijado antes de esta descarga.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
V37 = REPO / "results" / "v3_7"
SALIDA = REPO / "sitio_fuente" / "fuentes" / "anotacion"
MANIFIESTO = REPO / "sitio_fuente" / "fuentes" / "anotacion_manifest.json"
CORTE_BYTES = 80 * 1024 * 1024
AGENTE = "OncoGO-CRC/4 (mailto:jancod@jancod.co)"
OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"


def base_de(token: str) -> str:
    parte, _, sufijo = token.rpartition("-")
    if sufijo.isdigit() and parte:
        return parte
    return token


def sha256(datos: bytes) -> str:
    return hashlib.sha256(datos).hexdigest()


def pedir(url: str, timeout: int = 90) -> tuple[int, dict, bytes]:
    peticion = urllib.request.Request(url, headers={"User-Agent": AGENTE, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(peticion, timeout=timeout) as respuesta:
            cuerpo = respuesta.read()
            cabeceras = {k: respuesta.headers[k] for k in respuesta.headers}
            return respuesta.status, cabeceras, cuerpo
    except urllib.error.HTTPError as error:
        cuerpo = error.read()
        return error.code, dict(error.headers or {}), cuerpo


def guardar(ruta: Path, cuerpo: bytes) -> str:
    ruta.parent.mkdir(parents=True, exist_ok=True)
    ruta.write_bytes(cuerpo)
    return sha256(cuerpo)


def accesos() -> tuple[list[str], dict[str, str], dict[str, list[str]]]:
    celda = pd.read_csv(V37 / "suero_k4_crc_cn.tsv", sep="\t")
    representados = {}
    grupos = {}
    unicos = []
    vistos = set()
    for protein_id in celda["protein_id"].astype(str):
        tokens = [token for token in protein_id.split(";") if token]
        grupos[protein_id] = tokens
        principal = tokens[0]
        representados[protein_id] = principal
        for token in tokens:
            if token not in vistos:
                vistos.add(token)
                unicos.append(token)
    return unicos, representados, grupos


def descargar_uniprot(unicos: list[str], entradas: list[dict]) -> dict[str, dict]:
    por_acceso: dict[str, dict] = {}
    release = ""
    for inicio in range(0, len(unicos), 20):
        lote = unicos[inicio:inicio + 20]
        consulta = " OR ".join(f"accession:{token}" for token in lote)
        url = "https://rest.uniprot.org/uniprotkb/stream?format=json&query=" + urllib.parse.quote(consulta)
        codigo, cabeceras, cuerpo = pedir(url, timeout=180)
        nombre = f"uniprot_lote_{inicio:04d}.json"
        digest = guardar(SALIDA / nombre, cuerpo)
        release = cabeceras.get("X-UniProt-Release", release)
        entradas.append({
            "id": nombre,
            "url": url,
            "fecha_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "version": release,
            "sha256": digest,
            "fichero": f"anotacion/{nombre}",
            "estado": "ok" if codigo == 200 else "no_disponible",
            "motivo": "" if codigo == 200 else f"HTTP {codigo}",
        })
        if codigo != 200:
            for token in lote:
                por_acceso[token] = {"estado": "no_disponible", "motivo": f"HTTP {codigo}"}
            continue
        try:
            resultados = json.loads(cuerpo).get("results") or []
        except json.JSONDecodeError:
            for token in lote:
                por_acceso[token] = {"estado": "no_disponible", "motivo": "JSON ilegible"}
            continue
        hallados = {item.get("primaryAccession"): item for item in resultados}
        for token in lote:
            item = hallados.get(token)
            if item is None:
                por_acceso[token] = {"estado": "no_disponible", "motivo": "UniProt no devolvió el acceso"}
                continue
            ruta = SALIDA / "uniprot" / f"{token}.json"
            crudo = json.dumps(item, ensure_ascii=False, sort_keys=True).encode()
            digest = guardar(ruta, crudo)
            entradas.append({
                "id": f"uniprot:{token}",
                "url": url,
                "fecha_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "version": release,
                "sha256": digest,
                "fichero": f"anotacion/uniprot/{token}.json",
                "estado": "ok",
                "motivo": "",
            })
            por_acceso[token] = {"estado": "ok", "entrada": item, "version": release}
        print("uniprot", inicio + len(lote), flush=True)
        time.sleep(0.2)
    return por_acceso


def descargar_alfafold(unicos: list[str], entradas: list[dict]) -> dict[str, dict]:
    bases = []
    vistos = set()
    for token in unicos:
        base = base_de(token)
        if base not in vistos:
            vistos.add(base)
            bases.append(base)
    por_base: dict[str, dict] = {}
    for base in bases:
        url = f"https://alphafold.ebi.ac.uk/api/prediction/{base}"
        ruta = SALIDA / "alfafold" / f"{base}.json"
        if ruta.is_file() and ruta.stat().st_size > 0:
            cuerpo = ruta.read_bytes()
            codigo = 200
        else:
            codigo, _cabeceras, cuerpo = pedir(url)
        digest = guardar(ruta, cuerpo)
        estado = "ok" if codigo == 200 else "no_disponible"
        motivo = "" if codigo == 200 else f"HTTP {codigo}"
        version = ""
        item = None
        if codigo == 200:
            try:
                lista = json.loads(cuerpo)
                item = lista[0] if isinstance(lista, list) and lista else None
                version = str((item or {}).get("latestVersion") or "")
            except json.JSONDecodeError:
                estado = "no_disponible"
                motivo = "JSON ilegible"
                item = None
        entradas.append({
            "id": f"alphafold:{base}",
            "url": url,
            "fecha_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "version": version,
            "sha256": digest,
            "fichero": f"anotacion/alfafold/{base}.json",
            "estado": estado,
            "motivo": motivo,
        })
        por_base[base] = {"estado": estado, "motivo": motivo, "item": item, "version": version}
        time.sleep(0.05)
    print("alphafold", len(bases), flush=True)
    return por_base


def alojar_modelos(por_base: dict[str, dict], entradas: list[dict]) -> dict:
    destino = SALIDA / "modelos"
    destino.mkdir(parents=True, exist_ok=True)
    total = 0
    ficheros = []
    for base, info in por_base.items():
        item = info.get("item") or {}
        url = item.get("bcifUrl") or ""
        if info.get("estado") != "ok" or not url:
            continue
        ruta = destino / f"AF-{base}-F1.bcif"
        if ruta.is_file() and ruta.stat().st_size > 0:
            cuerpo = ruta.read_bytes()
            codigo = 200
        else:
            codigo, _cabeceras, cuerpo = pedir(url, timeout=180)
        if codigo != 200 or not cuerpo:
            info["modelo_motivo"] = f"HTTP {codigo}"
            continue
        nombre = f"AF-{base}-F1.bcif"
        digest = guardar(ruta, cuerpo)
        total += len(cuerpo)
        ficheros.append((nombre, digest, url, len(cuerpo), info.get("version", "")))
        print("modelo", base, len(cuerpo), "suma", total, flush=True)
    alojados = total <= CORTE_BYTES and total > 0
    if not alojados:
        for nombre, *_resto in ficheros:
            ruta = destino / nombre
            if ruta.exists():
                ruta.unlink()
    else:
        for nombre, digest, url, tamano, version in ficheros:
            entradas.append({
                "id": f"modelo:{nombre}",
                "url": url,
                "fecha_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "version": version,
                "sha256": digest,
                "fichero": f"anotacion/modelos/{nombre}",
                "estado": "ok",
                "motivo": "",
                "bytes": tamano,
            })
    return {
        "bytes": total,
        "mb": round(total / (1024 * 1024), 3),
        "corte_bytes": CORTE_BYTES,
        "alojados": alojados,
        "n": len(ficheros),
        "motivo": "" if alojados else "La suma supera el corte de 80 MB fijado en docs/sitio_v4_diseno.md.",
    }


def descargar_obo(entradas: list[dict]) -> tuple[dict[str, list[str]], dict[str, str], str]:
    ruta_obo = SALIDA / "go-basic.obo"
    if ruta_obo.is_file() and ruta_obo.stat().st_size > 1_000_000:
        cuerpo = ruta_obo.read_bytes()
        codigo, cabeceras = 200, {}
    else:
        codigo, cabeceras, cuerpo = pedir(OBO_URL, timeout=180)
    digest = guardar(ruta_obo, cuerpo)
    texto = cuerpo.decode("utf-8", errors="replace")
    version = ""
    for linea in texto.splitlines():
        if linea.startswith("data-version:"):
            version = linea.split(":", 1)[1].strip()
            break
    entradas.append({
        "id": "go-basic.obo",
        "url": OBO_URL,
        "fecha_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "version": version,
        "sha256": digest,
        "fichero": "anotacion/go-basic.obo",
        "estado": "ok" if codigo == 200 else "no_disponible",
        "motivo": "" if codigo == 200 else f"HTTP {codigo}",
        "nota": "No es el obo de la corrida: results/v3_7 no guardó la ontología.",
    })
    padres: dict[str, list[str]] = {}
    nombres: dict[str, str] = {}
    aspectos: dict[str, str] = {}
    actual: dict | None = None

    def cerrar() -> None:
        nonlocal actual
        if not actual or not str(actual.get("id", "")).startswith("GO:"):
            actual = None
            return
        if actual.get("obsolete"):
            actual = None
            return
        ident = actual["id"]
        nombres[ident] = actual.get("name", ident)
        padres[ident] = list(actual.get("is_a", []))
        aspectos[ident] = {
            "biological_process": "BP",
            "molecular_function": "MF",
            "cellular_component": "CC",
        }.get(actual.get("namespace", ""), "")
        actual = None

    for linea in texto.splitlines():
        if linea == "[Term]":
            cerrar()
            actual = {}
            continue
        if linea.startswith("[") and linea != "[Term]":
            cerrar()
            continue
        if linea == "":
            cerrar()
            continue
        if actual is None or ": " not in linea:
            continue
        clave, valor = linea.split(": ", 1)
        if clave == "id":
            actual["id"] = valor.strip()
        elif clave == "name":
            actual["name"] = valor.strip()
        elif clave == "namespace":
            actual["namespace"] = valor.strip()
        elif clave == "is_a":
            actual.setdefault("is_a", []).append(valor.split()[0])
        elif clave == "is_obsolete" and valor.strip() == "true":
            actual["obsolete"] = True
    cerrar()
    print("obo", version, "terminos", len(nombres), "cabeceras", cabeceras.get("Content-Length"), flush=True)
    return padres, nombres, aspectos, version


def ancestros(ident: str, padres: dict[str, list[str]], cache: dict[str, set[str]]) -> set[str]:
    if ident in cache:
        return cache[ident]
    sale: set[str] = set()
    for padre in padres.get(ident, []):
        sale.add(padre)
        sale |= ancestros(padre, padres, cache)
    cache[ident] = sale
    return sale


def evidencia_directa(xref: dict) -> bool:
    codigos = []
    for prop in xref.get("properties") or []:
        if prop.get("key") == "GoEvidenceType":
            codigos.append(str(prop.get("value", "")).split(":", 1)[0])
    return any(codigo and codigo != "IEA" for codigo in codigos)


def aspecto_de(termino: str) -> str:
    prefijo = termino.split(":", 1)[0]
    return {"P": "BP", "F": "MF", "C": "CC"}.get(prefijo, "")


def propiedades(xref: dict) -> dict[str, str]:
    return {str(item.get("key")): str(item.get("value")) for item in xref.get("properties") or []}


def ficha_de(token: str, info: dict, padres, nombres, cache, alfafold: dict, modelos: dict, aspectos=None) -> dict:
    aspectos = aspectos or {}
    if info.get("estado") != "ok":
        return {"acceso": token, "estado": "no_disponible", "motivo": info.get("motivo", "")}
    entrada = info["entrada"]
    secuencia = (entrada.get("sequence") or {})
    nombre = ""
    descripcion = entrada.get("proteinDescription") or {}
    recomendado = (descripcion.get("recommendedName") or {}).get("fullName") or {}
    nombre = str(recomendado.get("value") or "")
    if not nombre:
        enviados = descripcion.get("submissionNames") or []
        if enviados:
            nombre = str((enviados[0].get("fullName") or {}).get("value") or "")
    genes = []
    for gen in entrada.get("genes") or []:
        gen_nombre = (gen.get("geneName") or {}).get("value")
        if gen_nombre:
            genes.append(str(gen_nombre))
    funcion, localizacion = [], []
    for comentario in entrada.get("comments") or []:
        textos = []
        for texto in comentario.get("texts") or []:
            if texto.get("value"):
                textos.append(str(texto["value"]))
        for sub in comentario.get("subcellularLocations") or []:
            sitio = ((sub.get("location") or {}).get("value"))
            if sitio:
                localizacion.append(str(sitio))
        if comentario.get("commentType") == "FUNCTION":
            funcion.extend(textos)
    go_directo = []
    vistos = set()
    for xref in entrada.get("uniProtKBCrossReferences") or []:
        if xref.get("database") != "GO":
            continue
        ident = str(xref.get("id") or "")
        if not ident or ident in vistos or not evidencia_directa(xref):
            continue
        vistos.add(ident)
        props = propiedades(xref)
        go_directo.append({
            "id": ident,
            "aspecto": aspecto_de(props.get("GoTerm", "")),
            "nombre": props.get("GoTerm", ident).split(":", 1)[-1],
            "evidencia": props.get("GoEvidenceType", ""),
        })
    propagado = set(item["id"] for item in go_directo)
    for item in go_directo:
        propagado |= ancestros(item["id"], padres, cache)
    go_propagado = []
    for ident in sorted(propagado):
        directo = next((item["aspecto"] for item in go_directo if item["id"] == ident), "")
        go_propagado.append({
            "id": ident,
            "nombre": nombres.get(ident, ident),
            "aspecto": aspectos.get(ident) or directo,
        })
    def cruces(base: str) -> list[dict]:
        salida = []
        for xref in entrada.get("uniProtKBCrossReferences") or []:
            if xref.get("database") != base:
                continue
            props = propiedades(xref)
            etiqueta = props.get("EntryName") or props.get("PathwayName") or props.get("ProteinName") or ""
            salida.append({"id": str(xref.get("id")), "nombre": etiqueta, "propiedades": props})
        return salida
    pdb = []
    for xref in entrada.get("uniProtKBCrossReferences") or []:
        if xref.get("database") != "PDB":
            continue
        props = propiedades(xref)
        pdb.append({
            "id": str(xref.get("id")),
            "metodo": props.get("Method", ""),
            "resolucion": props.get("Resolution", ""),
            "cadenas": props.get("Chains", ""),
        })
    rangos = []
    for feature in entrada.get("features") or []:
        if feature.get("type") not in {"Signal", "Domain"}:
            continue
        lugar = feature.get("location") or {}
        inicio = ((lugar.get("start") or {}).get("value"))
        fin = ((lugar.get("end") or {}).get("value"))
        if not isinstance(inicio, int) or not isinstance(fin, int):
            continue
        rangos.append({
            "tipo": "signal" if feature.get("type") == "Signal" else "domain",
            "inicio": inicio,
            "fin": fin,
            "nombre": str(feature.get("description") or feature.get("type")),
        })
    base = base_de(token)
    modelo = alfafold.get(base) or {"estado": "no_disponible", "motivo": "sin consulta"}
    item = modelo.get("item") or {}
    local = None
    if modelos.get("alojados") and item.get("bcifUrl"):
        ruta = SALIDA / "modelos" / f"AF-{base}-F1.bcif"
        if ruta.is_file():
            local = f"fuentes/modelos/AF-{base}-F1.bcif"
    return {
        "acceso": token,
        "estado": "ok",
        "motivo": "",
        "version_uniprot": info.get("version", ""),
        "revisada": "reviewed" in str(entrada.get("entryType", "")).lower(),
        "entry_type": entrada.get("entryType", ""),
        "annotation_score": entrada.get("annotationScore"),
        "nombre": nombre,
        "gen_uniprot": genes[0] if genes else "",
        "longitud": secuencia.get("length"),
        "masa": secuencia.get("molWeight"),
        "secuencia": secuencia.get("value") or "",
        "funcion": " ".join(funcion),
        "localizacion": localizacion,
        "palabras": [str(item.get("name")) for item in entrada.get("keywords") or [] if item.get("name")],
        "go_directo": go_directo,
        "go_propagado": go_propagado,
        "interpro": cruces("InterPro"),
        "pfam": cruces("Pfam"),
        "reactome": cruces("Reactome"),
        "pdb": pdb,
        "rangos": rangos,
        "alfafold": {
            "estado": modelo.get("estado"),
            "motivo": modelo.get("motivo", ""),
            "version": modelo.get("version", ""),
            "plddt_medio": item.get("globalMetricValue"),
            "muy_alto": item.get("fractionPlddtVeryHigh"),
            "fiable": item.get("fractionPlddtConfident"),
            "bajo": item.get("fractionPlddtLow"),
            "muy_bajo": item.get("fractionPlddtVeryLow"),
            "modelo_local": local,
            "modelo_es_canonico": token != base,
            "entrada": item.get("modelEntityId") or item.get("entryId") or "",
            "bcif_url": item.get("bcifUrl") or "",
        },
    }


def descargar_string(entradas: list[dict]) -> dict:
    prio = pd.read_csv(V37 / "priorizacion.tsv", sep="\t")
    genes = sorted(prio.loc[prio["ultimo_paso"] == 5, "gen"].astype(str))
    identificadores = "%0d".join(genes)
    url = (
        "https://string-db.org/api/json/network?species=9606&required_score=700"
        f"&caller_identity=OncoGO-CRC&identifiers={identificadores}"
    )
    codigo, cabeceras, cuerpo = pedir(url, timeout=120)
    digest = guardar(SALIDA / "string.json", cuerpo)
    version = cabeceras.get("X-String-Version") or cabeceras.get("X-StringDB-Version") or ""
    estado = "ok" if codigo == 200 else "no_disponible"
    motivo = "" if codigo == 200 else f"HTTP {codigo}"
    aristas = []
    if estado == "ok":
        try:
            aristas = json.loads(cuerpo)
            if not isinstance(aristas, list):
                estado = "no_disponible"
                motivo = "la respuesta no es una lista"
                aristas = []
        except json.JSONDecodeError:
            estado = "no_disponible"
            motivo = "JSON ilegible"
    entradas.append({
        "id": "string",
        "url": url,
        "fecha_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "version": version,
        "sha256": digest,
        "fichero": "anotacion/string.json",
        "estado": estado,
        "motivo": motivo,
    })
    print("string", estado, "aristas", len(aristas), flush=True)
    return {"estado": estado, "motivo": motivo, "version": version, "n_aristas": len(aristas), "genes": genes}


def main() -> None:
    SALIDA.mkdir(parents=True, exist_ok=True)
    unicos, representados, grupos = accesos()
    entradas: list[dict] = []
    print("accesos", len(unicos), flush=True)
    uniprot = descargar_uniprot(unicos, entradas)
    alfafold = descargar_alfafold(unicos, entradas)
    modelos = alojar_modelos(alfafold, entradas)
    print("modelos", modelos, flush=True)
    padres, nombres, aspectos, obo_version = descargar_obo(entradas)
    string = descargar_string(entradas)
    cache: dict[str, set[str]] = {}
    fichas = {}
    for token, info in uniprot.items():
        fichas[token] = ficha_de(token, info, padres, nombres, cache, alfafold, modelos, aspectos)
    documento = {
        "nota_go": (
            "Los términos y las evidencias salen de UniProt. "
            "go_gene_long.tsv está en el manifiesto de la corrida y no trae códigos de evidencia, "
            "y la corrida no guardó el obo. La propagación usa go-basic.obo descargado aparte. "
            "No es la ontología de la tesis."
        ),
        "obo_version": obo_version,
        "modelos": modelos,
        "string": {k: v for k, v in string.items() if k != "genes"},
        "genes_paso_5": string["genes"],
        "por_acceso": fichas,
        "por_proteina": representados,
        "grupos": grupos,
    }
    crudo = json.dumps(documento, ensure_ascii=False).encode()
    digest = guardar(SALIDA / "fichas.json", crudo)
    entradas.append({
        "id": "fichas",
        "url": "",
        "fecha_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "version": "",
        "sha256": digest,
        "fichero": "anotacion/fichas.json",
        "estado": "ok",
        "motivo": "Derivado de las respuestas guardadas. No es una respuesta de red.",
    })
    manifiesto = {
        "fecha_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "corte_modelos_bytes": CORTE_BYTES,
        "modelos": modelos,
        "obo_version": obo_version,
        "string": documento["string"],
        "entradas": entradas,
    }
    MANIFIESTO.write_text(json.dumps(manifiesto, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("manifiesto", MANIFIESTO, "fichas", len(fichas), flush=True)


def rehacer() -> None:
    """Completa accesos que el lote no devolvió y rehace fichas.json."""
    unicos, representados, grupos = accesos()
    previo = json.loads(MANIFIESTO.read_text(encoding="utf-8")) if MANIFIESTO.is_file() else {}
    entradas = [item for item in previo.get("entradas", []) if item.get("id") != "fichas"]
    for token in unicos:
        ruta = SALIDA / "uniprot" / f"{token}.json"
        if ruta.is_file() and ruta.stat().st_size > 2:
            continue
        url = f"https://rest.uniprot.org/uniprotkb/{urllib.parse.quote(token)}.json"
        codigo, cabeceras, cuerpo = pedir(url)
        if codigo != 200:
            print("isoforma", token, codigo, flush=True)
            continue
        item = json.loads(cuerpo)
        crudo = json.dumps(item, ensure_ascii=False, sort_keys=True).encode()
        digest = guardar(ruta, crudo)
        entradas.append({
            "id": f"uniprot:{token}",
            "url": url,
            "fecha_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "version": cabeceras.get("X-UniProt-Release", ""),
            "sha256": digest,
            "fichero": f"anotacion/uniprot/{token}.json",
            "estado": "ok",
            "motivo": "",
        })
        print("isoforma", token, flush=True)
        time.sleep(0.05)
    uniprot = {}
    for ruta in (SALIDA / "uniprot").glob("*.json"):
        uniprot[ruta.stem] = {
            "estado": "ok",
            "entrada": json.loads(ruta.read_text(encoding="utf-8")),
            "version": "2026_04",
        }
    for token in unicos:
        uniprot.setdefault(token, {"estado": "no_disponible", "motivo": "UniProt no devolvió el acceso"})
    alfafold = {}
    for ruta in (SALIDA / "alfafold").glob("*.json"):
        try:
            lista = json.loads(ruta.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            alfafold[ruta.stem] = {"estado": "no_disponible", "motivo": "JSON ilegible", "item": None, "version": ""}
            continue
        item = lista[0] if isinstance(lista, list) and lista else None
        alfafold[ruta.stem] = {
            "estado": "ok" if item else "no_disponible",
            "motivo": "" if item else "sin modelo",
            "item": item,
            "version": str((item or {}).get("latestVersion") or ""),
        }
    modelos = previo.get("modelos") or {"alojados": False, "bytes": 0, "mb": 0, "motivo": "sin medida", "n": 0}
    padres, nombres, aspectos, obo_version = descargar_obo(entradas)
    cache: dict[str, set[str]] = {}
    fichas = {
        token: ficha_de(token, info, padres, nombres, cache, alfafold, modelos, aspectos)
        for token, info in uniprot.items()
    }
    prio = pd.read_csv(V37 / "priorizacion.tsv", sep="\t")
    genes = sorted(prio.loc[prio["ultimo_paso"] == 5, "gen"].astype(str))
    documento = {
        "nota_go": (
            "Los términos y las evidencias salen de UniProt. "
            "go_gene_long.tsv está en el manifiesto de la corrida y no trae códigos de evidencia, "
            "y la corrida no guardó el obo. La propagación usa go-basic.obo descargado aparte. "
            "No es la ontología de la tesis."
        ),
        "obo_version": obo_version,
        "modelos": modelos,
        "string": previo.get("string") or {},
        "genes_paso_5": genes,
        "por_acceso": fichas,
        "por_proteina": representados,
        "grupos": grupos,
    }
    digest = guardar(SALIDA / "fichas.json", json.dumps(documento, ensure_ascii=False).encode())
    entradas = [item for item in entradas if item.get("id") != "go-basic.obo" or item.get("sha256") == hashlib.sha256((SALIDA / "go-basic.obo").read_bytes()).hexdigest()]
    entradas.append({
        "id": "fichas",
        "url": "",
        "fecha_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "version": "",
        "sha256": digest,
        "fichero": "anotacion/fichas.json",
        "estado": "ok",
        "motivo": "Derivado de las respuestas guardadas. No es una respuesta de red.",
    })
    previo["entradas"] = entradas
    previo["obo_version"] = obo_version
    previo["fecha_utc"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    MANIFIESTO.write_text(json.dumps(previo, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("rehecho", len(fichas), flush=True)


if __name__ == "__main__":
    if "--rehacer" in sys.argv:
        rehacer()
    else:
        main()
