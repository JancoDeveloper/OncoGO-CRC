"""ESM-2 fuera de la corrida. CPU, ventanas de 1022, reanudable. La corrida solo lee el TSV."""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from oncogo_crc.settings import load_settings
from oncogo_crc.v3.secuencia_features import (
    ESM_OMITIDAS,
    ESM_TSV,
    FASTA_UNIPROT,
    MODELO_ESM,
    media_ponderada,
    secuencias_por_gen,
    sha256_directorio,
    sha256_fichero,
    ventanas,
)

COLUMNAS = ["gen", "longitud", "n_ventanas", *[f"esm2_{i:03d}" for i in range(480)]]
ANCHO = 483


def dispositivo_de_calculo():
    import torch
    return torch.device("cpu")


def comprobar_cpu(dispositivo) -> None:
    if getattr(dispositivo, "type", None) != "cpu":
        raise RuntimeError("ESM-2 de esta enmienda solo corre en CPU")


def preparar(hilos: int):
    os.environ["OMP_NUM_THREADS"] = str(hilos)
    os.environ["MKL_NUM_THREADS"] = str(hilos)
    os.environ["OPENBLAS_NUM_THREADS"] = str(hilos)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(hilos)
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    import torch
    torch.set_num_threads(hilos)
    try:
        torch.set_num_interop_threads(1)
    except RuntimeError as exc:
        print(f"no se pudieron fijar los hilos de interoperación: {exc}", flush=True)
    torch.use_deterministic_algorithms(True)
    torch.manual_seed(0)
    return torch


def fabrica_atencion():
    """Cada ventana abre su propio contexto. Uno solo no se puede reentrar."""
    from contextlib import nullcontext
    try:
        from torch.nn.attention import SDPBackend, sdpa_kernel

        def abrir():
            return sdpa_kernel(SDPBackend.MATH)

        with abrir():
            pass
        return abrir
    except (ImportError, AttributeError, RuntimeError, TypeError) as exc:
        print(f"atención matemática no disponible: {type(exc).__name__}: {exc}", flush=True)

        def abrir():
            return nullcontext()

        return abrir


def cargar_modelo(ruta: Path):
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    dispositivo = dispositivo_de_calculo()
    comprobar_cpu(dispositivo)
    tokenizador = AutoTokenizer.from_pretrained(ruta, local_files_only=True)
    modelo = AutoModelForMaskedLM.from_pretrained(ruta, local_files_only=True)
    modelo.eval()
    modelo.to(dispositivo)
    return modelo, tokenizador


def vector_de_secuencia(modelo, tokenizador, secuencia: str, fabrica) -> tuple:
    import torch
    dispositivo = dispositivo_de_calculo()
    comprobar_cpu(dispositivo)
    trozos = ventanas(secuencia)
    if not trozos:
        raise ValueError("secuencia vacía")
    medias = []
    for trozo in trozos:
        tokens = tokenizador(trozo, return_tensors="pt")
        ids = tokens["input_ids"]
        if int(ids.shape[-1]) != len(trozo) + 2:
            raise ValueError(f"tokens {int(ids.shape[-1])} distintos de {len(trozo) + 2}")
        tokens = {clave: valor.to(dispositivo) for clave, valor in tokens.items()}
        with fabrica(), torch.inference_mode():
            salida = modelo(**tokens, output_hidden_states=True)
        oculto = salida.hidden_states[-1][0, 1:-1]
        if int(oculto.shape[0]) != len(trozo):
            raise ValueError(f"capa {int(oculto.shape[0])} distinta de {len(trozo)}")
        medias.append(oculto.mean(dim=0).detach().cpu().numpy())
    vector = media_ponderada(medias, [len(trozo) for trozo in trozos])
    if vector.shape != (480,):
        raise ValueError(f"dimensión {vector.shape}")
    return vector, len(trozos)


def procesar_gen(gen: str, secuencia: str | None, modelo, tokenizador, fabrica) -> tuple[str | None, str | None]:
    """Devuelve la línea del vector o la de la omisión. Nunca las dos, y nunca ninguna."""
    if not secuencia:
        return None, f"{gen}\t\tsin secuencia en el FASTA revisado"
    try:
        vector, n_ventanas = vector_de_secuencia(modelo, tokenizador, secuencia, fabrica)
    except Exception as exc:
        return None, f"{gen}\t{len(secuencia)}\t{motivo_de(exc)}"
    return texto_vector(gen, len(secuencia), n_ventanas, vector), None


def texto_vector(gen: str, longitud: int, n_ventanas: int, vector) -> str:
    numeros = [f"{float(valor):.9f}" for valor in vector]
    return "\t".join([gen, str(int(longitud)), str(int(n_ventanas)), *numeros])


def motivo_de(exc: BaseException) -> str:
    texto = f"{type(exc).__name__}: {exc}".replace("\t", " ").replace("\n", " ")
    return texto


def anexar(path: Path, linea: str, encabezado: str) -> None:
    nuevo = not path.exists() or path.stat().st_size == 0
    with path.open("a", encoding="utf-8") as salida:
        if nuevo:
            salida.write(encabezado + "\n")
        salida.write(linea + "\n")
        salida.flush()
        os.fsync(salida.fileno())


def leer_parcial(path: Path, ancho: int) -> dict[str, str]:
    if not path.is_file():
        return {}
    lineas = path.read_text(encoding="utf-8").splitlines()
    if not lineas:
        return {}
    buenos: dict[str, str] = {}
    incompletas = 0
    for linea in lineas[1:]:
        if not linea.strip():
            continue
        campos = linea.split("\t")
        if len(campos) != ancho:
            incompletas += 1
            continue
        gen = campos[0]
        if gen in buenos:
            raise RuntimeError(f"gen repetido en {path.name}: {gen}")
        buenos[gen] = linea
    if incompletas:
        print(f"{path.name}: se descartan {incompletas} líneas incompletas y se rehacen", flush=True)
        with path.open("w", encoding="utf-8") as salida:
            salida.write(lineas[0] + "\n")
            for linea in buenos.values():
                salida.write(linea + "\n")
    return buenos


def genes_pendientes(genes: list[str], hechos: set[str]) -> list[str]:
    vistos: set[str] = set()
    pendientes = []
    for gen in genes:
        if gen in hechos or gen in vistos:
            continue
        vistos.add(gen)
        pendientes.append(gen)
    return pendientes


def escribir_final(path: Path, encabezado: str, filas: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as salida:
        salida.write(encabezado + "\n")
        for gen in sorted(filas):
            salida.write(filas[gen] + "\n")


def residente() -> int:
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def memsize() -> int:
    return int(subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True).strip())


def vigilar_memoria(limite: int) -> None:
    rss = residente()
    if rss > limite:
        print(
            f"memoria residente {rss} pasa de la mitad de hw.memsize ({limite}). "
            "Se para: es un error de código, no de tamaño.",
            flush=True,
        )
        raise SystemExit(3)


def versiones() -> tuple[str, str]:
    import torch
    import transformers
    return torch.__version__, transformers.__version__


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", default="/Volumes/JANCODATA/OncoGO-CRC")
    parser.add_argument("--hilos", type=int, default=4)
    parser.add_argument("--tabla", default="")
    parser.add_argument("--desde-cero", action="store_true")
    args = parser.parse_args(argv)
    if args.hilos < 1:
        print("los hilos tienen que ser al menos 1", flush=True)
        return 2
    inicio = time.perf_counter()
    preparar(args.hilos)
    settings = load_settings(data_workspace=args.workspace)
    fasta = settings.resolve(FASTA_UNIPROT)
    modelo_dir = settings.resolve(MODELO_ESM)
    destino = settings.resolve(ESM_TSV)
    omitidas = settings.resolve(ESM_OMITIDAS)
    parcial = destino.with_name(destino.name + ".parcial")
    parcial_om = omitidas.with_name(omitidas.name + ".parcial")
    if args.desde_cero:
        for ruta in (parcial, parcial_om, destino, omitidas, Path(str(destino) + ".sello.json"), Path(str(omitidas) + ".sello.json")):
            if ruta.is_file():
                ruta.unlink()
    if not fasta.is_file() or not modelo_dir.is_dir():
        print(f"falta el FASTA o el modelo: {fasta} {modelo_dir}", flush=True)
        return 2
    tabla_ruta = Path(args.tabla) if args.tabla else settings.root / "results/v3_7/h7_tabla_proteinas.tsv"
    genes = sorted(pd.read_csv(tabla_ruta, sep="\t", usecols=["gen"])["gen"].astype(str).str.upper().unique())
    secuencias = secuencias_por_gen(fasta, set(genes))
    hechos = leer_parcial(parcial, ANCHO)
    ya_omitidos = leer_parcial(parcial_om, 3)
    pendientes = genes_pendientes(genes, set(hechos) | set(ya_omitidos))
    limite = memsize() // 2
    print(
        f"genes {len(genes)} con secuencia {len(secuencias)} pendientes {len(pendientes)} "
        f"hilos {args.hilos} limite_rss {limite}",
        flush=True,
    )
    vigilar_memoria(limite)
    contexto = fabrica_atencion()
    modelo, tokenizador = cargar_modelo(modelo_dir)
    vigilar_memoria(limite)
    encabezado = "\t".join(COLUMNAS)
    for i, gen in enumerate(pendientes, start=1):
        linea, omision = procesar_gen(gen, secuencias.get(gen), modelo, tokenizador, contexto)
        if omision is not None:
            anexar(parcial_om, omision, "gen\tlongitud\tmotivo")
            print(f"{i}/{len(pendientes)} {gen} omitida {omision.split(chr(9))[-1]}", flush=True)
        else:
            anexar(parcial, linea, encabezado)
            print(f"{i}/{len(pendientes)} {gen} ventanas {linea.split(chr(9))[2]} rss {residente()}", flush=True)
        vigilar_memoria(limite)
    filas = leer_parcial(parcial, ANCHO)
    omisiones = leer_parcial(parcial_om, 3)
    cubiertos = set(filas) | set(omisiones)
    if set(filas) & set(omisiones) or set(genes) - cubiertos or cubiertos - set(genes):
        print(f"el parcial no cuadra con la tabla: faltan {sorted(set(genes) - cubiertos)[:5]}", flush=True)
        return 4
    escribir_final(destino, encabezado, filas)
    escribir_final(omitidas, "gen\tlongitud\tmotivo", omisiones)
    hash_tsv = sha256_fichero(destino)
    hash_om = sha256_fichero(omitidas)
    torch_v, trans_v = versiones()
    n_ventanas = 0
    for linea in filas.values():
        n_ventanas += int(linea.split("\t")[2])
    sello = {
        "fecha_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "dispositivo": "cpu",
        "hilos": args.hilos,
        "torch": torch_v,
        "transformers": trans_v,
        "sha256_fasta": sha256_fichero(fasta),
        "sha256_modelo": sha256_directorio(modelo_dir),
        "sha256_tsv": hash_tsv,
        "n_calculadas": len(filas),
        "n_omitidas": len(omisiones),
        "n_ventanas": n_ventanas,
        "tiempo_s": round(time.perf_counter() - inicio, 3),
        "memoria_residente_maxima": residente(),
        "hw_memsize": memsize(),
        "longitud_maxima": 1022,
        "genes": genes,
    }
    Path(str(destino) + ".sello.json").write_text(json.dumps(sello, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    sello_om = {
        "fecha_utc": sello["fecha_utc"],
        "sha256_tsv": hash_om,
        "n_omitidas": len(omisiones),
        "hilos": args.hilos,
    }
    Path(str(omitidas) + ".sello.json").write_text(json.dumps(sello_om, indent=2) + "\n", encoding="utf-8")
    parcial.unlink(missing_ok=True)
    parcial_om.unlink(missing_ok=True)
    print(
        f"listo calculadas {len(filas)} omitidas {len(omisiones)} ventanas {n_ventanas} "
        f"sha256 {hash_tsv} rss {residente()} tiempo {sello['tiempo_s']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
