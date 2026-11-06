"""Borradores y página a partir de results/tesis_cierre. No recalcula contrastes."""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pandas as pd

DISCO = Path("/Volumes/JANCODATA/OncoGO-CRC/results/tesis_cierre")
REPO = Path(__file__).resolve().parents[1]
WEB = REPO / "docs" / "web"
BORRADOR = REPO / "docs" / "borradores"


def coma(valor, decimales: int = 2) -> str:
    if valor is None or (isinstance(valor, float) and pd.isna(valor)):
        return "—"
    numero = float(valor)
    if numero != 0 and abs(numero) < 10 ** (-decimales):
        return f"{numero:.1e}".replace(".", ",")
    return f"{numero:.{decimales}f}".replace(".", ",")


def nombres_go(identificadores: list[str]) -> dict[str, str]:
    nombres: dict[str, str] = {}
    unicos = list(dict.fromkeys(identificadores))
    for inicio in range(0, len(unicos), 40):
        lote = unicos[inicio : inicio + 40]
        url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/" + ",".join(lote)
        peticion = urllib.request.Request(url, headers={"Accept": "application/json"})
        try:
            with urllib.request.urlopen(peticion, timeout=30) as respuesta:
                datos = json.load(respuesta)
        except Exception:
            continue
        for termino in datos.get("results", []):
            nombres[termino.get("id", "")] = termino.get("name", "")
    return nombres


def terminos_legibles(texto: object, nombres: dict[str, str], cupo: int = 4) -> str:
    if not isinstance(texto, str) or not texto.strip():
        return "Sin término GO en el rango de 15 a 300 genes."
    piezas = []
    for trozo in texto.split(";"):
        if trozo.count("|") < 2:
            continue
        go_id, espacio, _n = trozo.split("|", 2)
        nombre = nombres.get(go_id, go_id)
        piezas.append(f"{nombre} ({espacio})")
        if len(piezas) >= cupo:
            break
    return "; ".join(piezas)


def barra(filas: list[tuple[str, float]], ancho: int = 560, alto_fila: int = 28) -> str:
    if not filas:
        return ""
    maximo = max(abs(v) for _, v in filas) or 1
    alto = 16 + alto_fila * len(filas)
    partes = [f'<svg viewBox="0 0 {ancho} {alto}" role="img" class="grafica">']
    for i, (nombre, valor) in enumerate(filas):
        y = 12 + i * alto_fila
        largo = abs(valor) / maximo * 180
        x = 250 - largo if valor < 0 else 250
        color = "#9a3b2f" if valor < 0 else "#0f6e6e"
        partes.append(f'<text x="8" y="{y + 14}" class="eje">{nombre}</text>')
        partes.append(f'<rect x="{x:.1f}" y="{y}" width="{largo:.1f}" height="16" fill="{color}"/>')
        partes.append(f'<text x="{250 + 188}" y="{y + 13}" class="eje">{coma(valor)}</text>')
    partes.append(f'<line x1="250" y1="6" x2="250" y2="{alto - 6}" stroke="#1c2430" stroke-width="1"/>')
    partes.append("</svg>")
    return "".join(partes)


def _guia_html() -> str:
    preguntas = [
        (
            "¿Cuál es el resultado, en una frase?",
            "El cambio de proteínas del tumor se repite entre tejidos. En este suero no aparece como una señal propia del cáncer: ninguna proteína separa el cáncer de la colonoscopia negativa con el umbral fijado de antemano.",
        ),
        (
            "¿Cero proteínas es un fracaso?",
            "Es el resultado. El umbral (FDR 0,05 y cambio de al menos 1,5 veces) se escribió antes de mirar el suero. Un cero con ese control responde la pregunta. Inventar un marcador después de relajar el umbral no la responde.",
        ),
        (
            "¿Por qué no basta con comparar contra personas sanas?",
            "Contra sanos suben SAA1 y LBP, proteínas secretadas de respuesta del organismo. La colonoscopia negativa es el grupo que decide si esa subida es del cáncer o de un intestino enfermo. Frente a ese grupo, no pasan el FDR.",
        ),
        (
            "¿Qué es la columna P?",
            "En PXD034575, P es tumor primario, M es metástasis hepática y N es mucosa normal. Son 16 pacientes, cada uno con las tres muestras. La versión 2 llamó pólipo a P. Esa etiqueta no se usa en esta corrida y los certificados viejos no se reescribieron.",
        ),
        (
            "¿0,83 y 0,104 son la misma correlación corregida?",
            "No. 0,83 es primario frente a mucosa, con ceros tratados como ausencia, en 874 genes, contra CPTAC. 0,104 es otro cálculo de la versión 2, sobre 1 993 proteínas y con la etiqueta antigua de columnas. Ninguna de las dos mide tejido contra sangre.",
        ),
        (
            "¿SAA1 sirve de marcador?",
            "Separa el cáncer de los sanos (log2FC 3,26, FDR 0,0002) y está secretada a la sangre. No separa el cáncer de la colonoscopia negativa. En la sustentación se presenta como respuesta del huésped.",
        ),
        (
            "¿Y BAIAP3?",
            "Es la proteína que Rai señaló. El AUC frente a la colonoscopia negativa es 0,82 y no tiene marca de fase aguda. El FDR es 0,20. Queda como candidata de otro estudio, con más personas y, si se puede, el CEA de la misma persona. No es un hallazgo de esta tesis.",
        ),
        (
            "¿CHMP4A, con un log2FC de 5,42?",
            "El cambio es grande y solo hay medición en 16 de 40 sueros. Un valor ausente no es un cero biológico. No se presenta como marcador.",
        ),
        (
            "¿Se confirmaron neoantígenos?",
            "No se pueden confirmar con esta matriz. Cada fila es un grupo de proteínas, no la secuencia de un péptido mutado. Rai describe 19 péptidos en suero, Xiang neoepítopos en el tumor y Hembrough cinco péptidos de una persona. Esos trabajos marcan el límite; no se reanalizan aquí.",
        ),
        (
            "¿Por qué no se juntaron los estudios en una sola tabla?",
            "CPTAC es tejido TMT ya en log2. Tanaka es tejido LFQ. Rai es suero DIA. Juntarlos mezclaría técnicas y personas distintas. Cada contraste se hizo dentro de su cohorte.",
        ),
    ]
    terminos = [
        ("log2FC", "Cuántas veces cambia la proteína, en logaritmo base 2. 1 significa el doble. −1 significa la mitad. El umbral 0,585 es log2(1,5): un cambio de al menos una vez y media."),
        ("FDR", "Entre las proteínas que se declaran distintas, la fracción esperada de falsos. Se usa Benjamini–Hochberg y el corte es 0,05. No es el p de una sola proteína."),
        ("t moderada, treat", "Prueba si el cambio supera 0,585, no si es distinto de cero. Así un cambio diminuto y muy preciso no entra en la lista. Smyth 2004 y McCarthy y Smyth 2009."),
        ("AUC", "Probabilidad de que una proteína ordene a una persona con cáncer por encima de una con la otra condición. 0,50 es el azar. Aquí se pidió más de 0,70, con el intervalo por encima de 0,50. El AUC no dice si la proteína sube o baja, y no sustituye al FDR."),
        ("Pareado", "En Tanaka, el tumor y la mucosa son del mismo paciente. El contraste resta dentro de la persona. En el suero de Rai los grupos son personas distintas: el contraste no es pareado."),
        ("DIA y grupo proteico", "DIA mide muchas proteínas en cada suero sin elegirlas de antemano. El informe las entrega en grupos. Un grupo no guarda la secuencia mutada de un péptido."),
        ("GO y NES", "Gene Ontology es el catálogo de procesos. NES dice si los genes de un proceso se acumulan en la cabeza del ranking. En este suero solo entran unos 12 procesos del tumor con al menos 10 genes, y no superan el azar."),
        ("Fase aguda", "Respuesta del organismo a inflamación, no una proteína exclusiva del tumor. C9, LRG1, SERPINA1, CRP y HP ya estaban en esa lista antes del contraste."),
    ]
    bloques_p = []
    for pregunta, respuesta in preguntas:
        bloques_p.append(
            f"<details><summary>{pregunta}</summary><p>{respuesta}</p></details>"
        )
    bloques_t = []
    for termino, decir in terminos:
        bloques_t.append(f"<article class=\"dato\"><h3>{termino}</h3><p>{decir}</p></article>")
    return f"""
  <h2 id="estudiar">Para sustentar: qué estudiar</h2>
  <p>El orden cabe en una tarde. Primero la frase del resultado, en voz alta y sin la hoja. Después el dibujo de las tres cohortes. Después dos proteínas: LBP, que pasa el FDR contra sanos y no contra la enfermedad, y BAIAP3, que pasa el AUC y no el FDR. Al final, la leyenda P, M y N, y la diferencia entre 0,83 y 0,104.</p>
  <ol>
    <li>Decir el resultado en una frase.</li>
    <li>Nombrar las tres cohortes y qué pregunta responde cada una.</li>
    <li>Explicar FDR y AUC con LBP y con BAIAP3.</li>
    <li>Explicar por qué P es tumor primario y por qué el cero de MaxQuant es una ausencia.</li>
    <li>Explicar por qué 0,83 y 0,104 son estadísticas distintas, y por qué ninguna es tejido contra sangre.</li>
    <li>Cerrar con el límite: 20 personas por grupo, sin el inmunoensayo de CEA, y sin secuencia de péptido mutado.</li>
  </ol>
  <h3>Palabras que hay que poder definir</h3>
  <div class="rejilla">{''.join(bloques_t)}</div>
  <h3>Si lo preguntan en la sustentación</h3>
  <p>Abre cada pregunta y responde antes de leer. La respuesta de abajo es la que sostienen estos números.</p>
  {''.join(bloques_p)}
  <h3>Frases que no se dicen</h3>
  <ul>
    <li>Que hay un marcador de sangre listo para la consulta.</li>
    <li>Que SAA1, LBP o CHMP4A distinguen el cáncer de otra enfermedad del intestino.</li>
    <li>Que la columna P es un pólipo o un adenoma.</li>
    <li>Que 0,83 corrige o sustituye a 0,104.</li>
    <li>Que esta matriz confirmó un neoantígeno.</li>
    <li>Que las vías del tumor se ven en este suero.</li>
  </ul>
"""


def _guia_md(primario, metastasis) -> str:
    return f"""# Para sustentar OncoGO-CRC

La página con gráficas y fichas está en [divulgacion.html](divulgacion.html). Esta nota es lo que hay que poder decir.

## El resultado, en una frase

Entre las 688 proteínas certificadas, el primario coincide en signo con CPTAC en 116 de 118. En el suero, el corte de 1,5 veces no separa el cáncer de la colonoscopia negativa: 145 proteínas tienen ese cambio descartado y 731 quedan indeterminadas.

Sirve como tesis, porque la pregunta queda respondida con un control de otra enfermedad intestinal. No sirve como prueba de un marcador para la consulta.

## Tres cohortes

| Cohorte | Qué es | Qué responde |
|---|---|---|
| PDC000116, Vasaikar, *Cell* 2019 | Tejido tumoral y mucosa, TMT, ya en log2 | El ancla. El tumor tiene un proteoma distinto. |
| PXD034575, Tanaka, *Cell Reports* 2024 | 16 pacientes. P primario, M metástasis hepática, N mucosa | Si ese proteoma se repite en otro tejido. En las 688 certificadas: 116 de 118 con el primario, 93 de 106 con la metástasis. El cierre del 31 de octubre, en otro universo, es {primario['mismo_signo']} de {primario['n']} y correlación {coma(primario['pearson'])}. |
| PXD055218, Rai, *J Proteome Res* 2025 | Suero DIA, 20 cáncer, 20 colonoscopia negativa, 20 sanos, 876 proteínas | Frente a la colonoscopia negativa el corte de 1,5 veces deja 0 proteínas. De 876, 145 son equivalentes y 731 indeterminadas. |

## Orden de una tarde

1. Decir el resultado en voz alta, sin la hoja.
2. Dibujar las tres cohortes y la pregunta de cada una.
3. Explicar LBP (pasa el FDR contra sanos, no contra la enfermedad) y BAIAP3 (AUC 0,82, FDR 0,20).
4. Explicar P, M y N, y por qué el cero de MaxQuant es una ausencia.
5. Explicar que 0,83 y 0,104 son estadísticas distintas. Ninguna mide tejido contra sangre.
6. Cerrar con el límite: 20 por grupo, sin el inmunoensayo de CEA, sin secuencia de péptido mutado.

## Palabras

- **log2FC.** Cambio en logaritmo base 2. 1 es el doble. El umbral 0,585 es log2(1,5).
- **FDR.** Fracción esperada de falsos entre lo que se declara distinto. Corte 0,05, Benjamini–Hochberg.
- **treat.** La t moderada pregunta si el cambio supera 0,585, no si es distinto de cero.
- **AUC.** Si la proteína ordena al cáncer por encima de la otra condición. 0,50 es azar. No dice si sube o baja y no sustituye al FDR.
- **Pareado.** En Tanaka, tumor y mucosa son del mismo paciente. En el suero, los grupos son personas distintas.
- **Grupo proteico.** La fila del DIA no guarda la secuencia mutada de un péptido.
- **Fase aguda.** Respuesta del organismo. C9, LRG1, SERPINA1, CRP y HP ya estaban marcadas así.

## Si lo preguntan

**¿Cero proteínas es un fracaso?** Es el resultado. El umbral se fijó antes de mirar el suero.

**¿Por qué no basta el sano?** Contra sanos suben SAA1 y LBP. Eso mide respuesta del organismo. La colonoscopia negativa es el control que decide la especificidad.

**¿SAA1 es un marcador?** Separa el cáncer de los sanos. No separa el cáncer de la otra enfermedad intestinal.

**¿BAIAP3?** Candidata de Rai. Pasa el AUC sin marca de fase aguda y no pasa el FDR. No es un hallazgo de esta tesis.

**¿CHMP4A?** log2FC alto con medición en 16 de 40 sueros. No se presenta como marcador.

**¿Neoantígenos?** Esta matriz no puede confirmarlos. Rai, Xiang y Hembrough describen péptidos. Aquí cada fila es un grupo de proteínas.

**¿Por qué no se fusionaron las cohortes?** Tejido TMT, tejido LFQ y suero DIA no son la misma medición ni las mismas personas.

## No se dice

- Que hay un marcador de sangre listo para la consulta.
- Que la columna P es un pólipo.
- Que 0,83 sustituye a 0,104.
- Que las vías del tumor se ven en este suero.
"""


def main() -> None:
    auditoria = json.loads((DISCO / "cierre_audit.json").read_text(encoding="utf-8"))
    fichas = pd.read_csv(DISCO / "fichas_proteinas.tsv", sep="\t")
    ids = []
    for texto in fichas["go"].dropna().astype(str):
        for trozo in texto.split(";"):
            if "|" in trozo:
                ids.append(trozo.split("|", 1)[0])
    nombres = nombres_go(ids)
    tejido = auditoria["pxd034575"]
    vias = auditoria["vias_suero"]["contrastes"]
    primario = tejido["concordancia_con_cptac"]["primario_vs_mucosa"]
    metastasis = tejido["concordancia_con_cptac"]["metastasis_vs_mucosa"]
    progresion = tejido["progresion_primario_hacia_metastasis"]
    cptac_sig = tejido["contrastes"]

    def fila(gen: str) -> pd.Series:
        return fichas.loc[fichas["gene"] == gen].iloc[0]

    saa1, lbp, chmp, apoa, baia = map(fila, ["SAA1", "LBP", "CHMP4A", "APOA4", "BAIAP3"])
    especificas = fichas.loc[fichas["veredicto"] == "especifica"].copy()

    WEB.mkdir(parents=True, exist_ok=True)
    BORRADOR.mkdir(parents=True, exist_ok=True)
    html = _html(auditoria, fichas, nombres, tejido, vias, primario, metastasis, progresion, cptac_sig, especificas)
    (WEB / "divulgacion.html").write_text(html, encoding="utf-8")
    sitio = REPO / "sitio"
    sitio.mkdir(parents=True, exist_ok=True)
    (sitio / "index.html").write_text(html, encoding="utf-8")
    (WEB / "divulgacion.md").write_text(_markdown(auditoria, primario, metastasis, progresion, vias, cptac_sig), encoding="utf-8")
    (WEB / "sustentacion.md").write_text(_guia_md(primario, metastasis), encoding="utf-8")
    (WEB / "diccionario.md").write_text(_diccionario_md(), encoding="utf-8")
    (BORRADOR / "tesis_borrador.tex").write_text(
        _latex(auditoria, primario, metastasis, progresion, vias, cptac_sig, saa1, lbp, chmp, apoa, baia),
        encoding="utf-8",
    )
    _docx(auditoria, primario, metastasis, progresion, vias, cptac_sig, saa1, lbp, chmp, apoa, baia, fichas, nombres)
    print(f"html {WEB / 'divulgacion.html'}")
    print(f"guia {WEB / 'sustentacion.md'}")
    print(f"tex {BORRADOR / 'tesis_borrador.tex'}")
    print(f"docx {BORRADOR / 'tesis_borrador.docx'}")


def _explicacion_html() -> str:
    return """
  <h2 id="teoria">La teoría detrás de cada paso</h2>
  <p>Cada número de arriba sale de un procedimiento con una regla escrita antes de mirar el resultado. Esta sección dice qué pregunta hace cada procedimiento y qué artículo lo sostiene.</p>

  <h3>Tres mediciones que no se suman</h3>
  <p>Una cohorte es un grupo de personas medidas con el mismo método. CPTAC es tejido marcado con TMT: varias muestras se comparan dentro del mismo experimento y el archivo ya llega en log2. Tanaka es tejido medido por LFQ: la intensidad de cada proteína se estima por separado en cada muestra. Rai es suero medido por DIA: se intentan ver muchas proteínas a la vez, sin elegirlas de antemano. Sumar las tres tablas mezclaría técnicas y personas. Por eso cada contraste ocurre dentro de su cohorte. <a href="https://doi.org/10.1371/journal.pcbi.1011828">Robles et al., 2024</a> ya integraron proteomas públicos de cáncer colorrectal; esta tesis no reclama ser la primera en juntar datos públicos.</p>

  <h3>El cero y el log2</h3>
  <p>En la matriz de Tanaka el cero de MaxQuant significa que esa proteína no se identificó en esa muestra, no que su cantidad sea cero. <a href="https://doi.org/10.1038/nbt.1511">Cox y Mann, 2008</a> describen ese programa. Antes del contraste, el cero se vuelve un hueco y solo la intensidad positiva pasa a log2. El log2 convierte una multiplicación en una resta: log2FC = 1 es el doble, log2FC = −1 es la mitad, y 0,585 es log2(1,5), un cambio de al menos una vez y media. Si el cero se dejara como número, la mediana de la matriz caería en cero y el programa creería que los datos ya están en logaritmo. Esa lectura produce cambios absurdos y no se usa.</p>

  <h3>Pareado o de grupos</h3>
  <p>En Tanaka, el primario, la metástasis y la mucosa son del mismo paciente. El contraste pareado resta dentro de la persona y así aparta lo que esa persona ya traía. En el suero de Rai los tres grupos son personas distintas, así que el contraste compara grupos. Pedir pareo donde no hay pares inventaría una relación que los datos no tienen.</p>

  <h3>La t moderada</h3>
  <p><a href="https://doi.org/10.2202/1544-6115.1027">Smyth, 2004</a> reparte la varianza entre las proteínas: una proteína con pocas mediciones no se declara distinta solo porque su desviación salió pequeña por azar. <a href="https://doi.org/10.1093/bioinformatics/btp053">McCarthy y Smyth, 2009</a> llaman treat a la versión que pregunta si el cambio supera un umbral, aquí 0,585, y no si es distinto de cero. Un cambio diminuto y muy preciso se queda fuera. No se imputan los huecos: una proteína ausente no se rellena con un promedio.</p>

  <h3>El FDR</h3>
  <p>Se prueban cientos de proteínas a la vez. Si cada una se juzgara con p &lt; 0,05, varias entrarían por azar. <a href="https://doi.org/10.1111/j.2517-6161.1995.tb02031.x">Benjamini y Hochberg, 1995</a> controlan la fracción esperada de falsos entre las que se declaran distintas. El corte es 0,05. El FDR no es el p de una sola proteína. En esta corrida, cáncer frente a colonoscopia negativa deja cero proteínas dentro de ese corte y del cambio mínimo. De las 876 proteínas del suero, 145 tienen un cambio menor que 1,5 veces demostrado y 731 quedan indeterminadas. La ausencia solo se afirma en las 145.</p>

  <h3>El AUC y la especificidad</h3>
  <p>El AUC dice si, al ordenar a las personas por una proteína, las del cáncer quedan por encima de las de la otra condición. 0,50 es el orden del azar. La regla de esta tesis pide tres cosas a la vez: AUC mayor que 0,70 al separar cáncer de colonoscopia negativa, el límite inferior del intervalo bootstrap por encima de 0,50, y ese AUC mayor que el de enfermedad frente a sanos. El intervalo usa 2 000 remuestreos y la semilla 20261031. El AUC no dice si la proteína sube o baja. Una proteína ya marcada como fase aguda, hemólisis o enfermedad inflamatoria conserva el número y cambia la lectura: el AUC puede pasar y la proteína seguir siendo una señal de inflamación. Esa marca se escribió en el catálogo antes del contraste.</p>

  <h3>Gene Ontology y las vías</h3>
  <p><a href="https://doi.org/10.1038/75556">Ashburner et al., 2000</a> arman Gene Ontology: un vocabulario de procesos, componentes y funciones. <a href="https://doi.org/10.1073/pnas.0506580102">Subramanian et al., 2005</a> proponen mirar si los genes de un proceso se acumulan al principio de un ranking, no solo si un gen pasa un corte. Aquí el ranking del suero es el signo del log2FC por menos el log10 del p. Se compara con el NES de los procesos que en el tejido de CPTAC tenían FDR de hasta 0,25. Entran solo los procesos con al menos 10 genes medidos en el suero. Son unos doce. Se permuta el ranking 200 veces. Esas correlaciones no superan el azar: este suero no contiene el grueso de los procesos que el tumor altera.</p>

  <h3>Atlas, ficha y estructura</h3>
  <p>La ficha anota la proteína para poder verificarla. <a href="https://doi.org/10.1126/science.1260419">Uhlén et al., 2015</a> publican el Human Protein Atlas: concentraciones en sangre de personas de referencia, no de estos pacientes. UniProt identifica la proteína. AlphaFold, <a href="https://doi.org/10.1038/s41586-021-03819-2">Jumper et al., 2021</a>, predice un pliegue. Un pliegue no demuestra que la cantidad haya cambiado en el suero. <a href="https://doi.org/10.1016/j.cell.2025.03.013">Malmström et al., 2025</a> publican un atlas de origen tisular del plasma; esta corrida usa la secreción del atlas de proteínas y los términos GO, y no descarga ese suplemento como una tabla aplicada paciente a paciente.</p>

  <h3>Por qué no hay neoantígeno en estas filas</h3>
  <p>Un neoantígeno es un péptido que el tumor presenta porque una mutación cambió su secuencia. La matriz de esta tesis guarda grupos de proteínas. Un grupo dice que esos péptidos se parecen a una proteína conocida. No guarda la letra cambiada. <a href="https://doi.org/10.1021/acs.jproteome.5c00285">Rai et al., 2025</a> describen 19 péptidos en suero y un PRM de cuatro péptidos, uno de ellos neoantigénico. <a href="https://doi.org/10.1016/j.xgen.2025.101062">Xiang et al., 2026</a> describen neoepítopos del tumor, no del suero de esta matriz. <a href="https://doi.org/10.1158/1538-7445.am2026-7710">Hembrough et al., 2026</a> recuperan cinco péptidos de una persona con cáncer colorrectal avanzado. Esos tres trabajos marcan el límite. No se reanalizan como si fueran filas de esta tabla.</p>

  <h2 id="datos">Diccionario de datos</h2>
  <p>Estos nombres aparecen en la página y en los archivos de <code>results/tesis_cierre</code> del disco. El diccionario dice qué es cada uno.</p>
  <table>
    <thead><tr><th>Nombre</th><th>Dónde está</th><th>Qué significa</th></tr></thead>
    <tbody>
      <tr><td>gene</td><td>Ficha y especificidad</td><td>Símbolo del gen humano con el que se emparejó la proteína. El emparejamiento usa el mapa del repositorio más el mapa revisado de UniProt, por símbolo o por acceso.</td></tr>
      <tr><td>uniprot</td><td>Ficha</td><td>Acceso de UniProt. BAIAP3 en esta corrida es O94812. El enlace de la ficha abre esa entrada.</td></tr>
      <tr><td>protein_id</td><td>Tablas de contraste</td><td>Identificador de la fila en la matriz. En el suero y en Tanaka es un acceso. En la tabla certificada de CPTAC es el símbolo del gen.</td></tr>
      <tr><td>log2fc</td><td>Contraste</td><td>Diferencia de medias en log2. Positivo: más abundante en el primer grupo del contraste. En tejido pareado, el primer grupo es el de la izquierda del nombre del archivo.</td></tr>
      <tr><td>padj</td><td>Contraste</td><td>p ajustado por Benjamini–Hochberg. El corte de la lista corta es 0,05, junto con |log2FC| mayor que 0,585.</td></tr>
      <tr><td>p_value</td><td>Contraste</td><td>p de esa proteína antes de corregir por las demás. No se usa solo para declararla distinta.</td></tr>
      <tr><td>n_obs</td><td>Contraste y ficha</td><td>Cuántas muestras tenían valor en ese contraste. CHMP4A tiene 16 de 40 en cáncer frente a sano. Un hueco no es un cero.</td></tr>
      <tr><td>significant</td><td>Contraste</td><td>Verdadero si pasa el FDR y el cambio mínimo. En cáncer frente a enfermedad, ninguna fila lo es.</td></tr>
      <tr><td>escala</td><td>Contraste de Tanaka</td><td><code>log2_de_intensidad_lineal</code> después de convertir los ceros en huecos. Si apareciera <code>ya_en_log_o_razon</code> en esa matriz, la corrida estaría mal leída.</td></tr>
      <tr><td>auc_crc_vs_enfermedad</td><td>Especificidad</td><td>AUC sin signo, el mayor entre el AUC y 1 menos el AUC, al separar cáncer de colonoscopia negativa.</td></tr>
      <tr><td>auc_crc_vs_enfermedad_ci_bajo</td><td>Especificidad</td><td>Límite inferior del intervalo bootstrap. La regla pide que supere 0,50.</td></tr>
      <tr><td>auc_enfermedad_vs_sano</td><td>Especificidad</td><td>AUC de la otra enfermedad frente a sanos. Si es mayor que el del cáncer, la proteína no se llama específica.</td></tr>
      <tr><td>veredicto</td><td>Especificidad</td><td><code>especifica</code>, <code>compartida_con_enfermedad</code>, <code>inespecifica</code> o <code>no_evaluable</code>. Sale solo del AUC. No relaja el FDR.</td></tr>
      <tr><td>lista_ruido</td><td>Catálogo</td><td><code>fase_aguda</code>, <code>hemolisis</code> o <code>eii</code>. Se escribió antes del contraste. No borra el AUC.</td></tr>
      <tr><td>lectura</td><td>Ficha</td><td>El veredicto más la marca. Si el AUC pasa y hay marca, la lectura es <code>auc_pasa_pero_fase_aguda</code> o la marca que corresponda.</td></tr>
      <tr><td>origen</td><td>Ficha</td><td><code>literatura</code> si la proteína venía del catálogo. <code>resultado_propio</code> si entró por pasar el FDR en esta corrida.</td></tr>
      <tr><td>go</td><td>Ficha</td><td>Hasta cinco términos biológicos y tres de componente celular, solo si el término tiene entre 15 y 300 genes. El nombre se leyó en QuickGO.</td></tr>
      <tr><td>hpa_localizacion</td><td>Ficha</td><td>Si el atlas anota la proteína como secretada a sangre. Es una anotación de referencia.</td></tr>
      <tr><td>hpa_ms_pg_l, hpa_inmuno_pg_l</td><td>Ficha</td><td>Concentración en sangre de personas de referencia, por espectrometría o por inmunoensayo, en pg/L. No es la concentración de estos pacientes.</td></tr>
      <tr><td>CP_, CN_, Normal_</td><td>Matriz de Rai</td><td>Cáncer, colonoscopia negativa y sano. Veinte columnas de cada prefijo. El orden alfabético no se usa como diagnóstico.</td></tr>
      <tr><td>LFQ intensity CRCnP, M, N</td><td>Matriz de Tanaka</td><td>P primario, M metástasis hepática, N mucosa. El número es el paciente. Hay 16 de cada letra.</td></tr>
      <tr><td>0,104</td><td>results/reproducibility</td><td>Correlación de una corrida anterior, 1 993 proteínas, con la etiqueta vieja de columnas. No es tejido contra sangre y no se sustituye por 0,83.</td></tr>
      <tr><td>0,83</td><td>cierre_audit.json</td><td>Correlación del log2FC del primario frente a mucosa con el log2FC certificado de CPTAC, en 874 genes, misma dirección en 687.</td></tr>
    </tbody>
  </table>

  <h2 id="procesos">Diccionario de procesos</h2>
  <p>El orden es el de la etapa <code>tesis_cierre</code>. No reescribe los certificados de la versión 1 ni de la versión 2. Los archivos quedan en el disco, bajo <code>results/tesis_cierre</code>.</p>
  <table>
    <thead><tr><th>Proceso</th><th>Qué hace</th><th>Qué deja escrito</th></tr></thead>
    <tbody>
      <tr><td>Grupos del suero</td><td>Lee la matriz DIA y separa las columnas CP, CN y Normal. Si esas columnas existen, no usa otra etiqueta clínica.</td><td>Tres contrastes: cáncer–sano, cáncer–enfermedad, enfermedad–sano.</td></tr>
      <tr><td>Contraste de grupos</td><td>t moderada de dos grupos, umbral 0,585, sin imputar. Luego FDR 0,05.</td><td><code>crc_vs_sano.tsv</code>, <code>crc_vs_enfermedad.tsv</code>, <code>enfermedad_vs_sano.tsv</code>.</td></tr>
      <tr><td>Especificidad</td><td>Para cada candidata del catálogo que esté medida, calcula el AUC y el intervalo. Aplica la regla de los tres requisitos.</td><td><code>especificidad_candidatas.tsv</code>. Ocho cumplen el AUC. Tres no tienen marca de ruido y aun así no pasan el FDR.</td></tr>
      <tr><td>Lectura con ruido</td><td>Si el veredicto es específico y la lista de ruido no está vacía, la lectura lo dice al lado. El AUC no se baja de 0,70 para compensar.</td><td>Columna <code>lectura</code>.</td></tr>
      <tr><td>Leyenda de Tanaka</td><td>Lee el sufijo P, M o N del nombre de la columna. No usa la etiqueta «pólipo» de la versión 2.</td><td>16 pacientes con primario, metástasis y mucosa.</td></tr>
      <tr><td>Cero y log2</td><td>Convierte intensidades no positivas en hueco y aplica log2 a las positivas. El contraste es pareado, con al menos 8 pares.</td><td><code>tanaka_primario_vs_mucosa.tsv</code>, <code>tanaka_metastasis_vs_mucosa.tsv</code>, <code>tanaka_metastasis_vs_primario.tsv</code>.</td></tr>
      <tr><td>Concordancia con CPTAC</td><td>Empareja por gen el log2FC de Tanaka con la tabla certificada de CPTAC y calcula el mismo signo y la correlación de Pearson.</td><td>Bloque <code>pxd034575</code> de <code>cierre_audit.json</code>.</td></tr>
      <tr><td>Vías</td><td>Ranking del suero contra los términos GO biológicos del tejido. Pide al menos 10 genes en el suero y 200 permutaciones.</td><td>Bloque <code>vias_suero</code>. Doce a trece términos. El nulo no se supera.</td></tr>
      <tr><td>Fichas</td><td>Reúne las proteínas que pasan un FDR o que tienen una lectura evaluable. Les añade GO, atlas, UniProt y AlphaFold.</td><td><code>fichas_proteinas.tsv</code> y las tarjetas de esta página.</td></tr>
      <tr><td>Péptidos</td><td>Escribe el límite: la matriz no confirma secuencia mutada. No busca espectros nuevos.</td><td><code>peptidos_neoantigenos.tsv</code>.</td></tr>
      <tr><td>PXD036594</td><td>No entra. <a href="https://doi.org/10.1021/acs.jproteome.2c00575">Lesur et al., 2023</a> miden plasma dirigido, y un orden alfabético de columnas no es un diagnóstico. El depósito está en <a href="https://www.ebi.ac.uk/pride/archive/projects/PXD036594">PXD036594</a>.</td><td>Esas columnas quedan como <code>sin_etiqueta_clinica</code>.</td></tr>
    </tbody>
  </table>

  <h2 id="estudios">Estudios enlazados y cómo se comprobaron</h2>
  <p>El 31 de octubre de 2026 cada DOI de esta lista se consultó en Crossref. La respuesta trae el título y el año. Coinciden con la cita. Donde la editorial respondió, también se abrió la página del artículo. El DOI <code>10.1016/j.xcrm.2025.101689</code> devolvió 404 y no se enlaza. Tampoco se enlaza una cita mal formada de ENCODE ni un DOI de memoria que no resolvió.</p>
  <table>
    <thead><tr><th>Para qué se usa</th><th>Estudio</th><th>Enlace comprobado</th></tr></thead>
    <tbody>
      <tr><td>Ancla de tejido</td><td>Vasaikar et al., Cell 2019. Proteogenomic analysis of human colon cancer.</td><td><a href="https://doi.org/10.1016/j.cell.2019.03.030">DOI</a> · <a href="https://pdc.cancer.gov/pdc/study/PDC000116">PDC000116</a></td></tr>
      <tr><td>Primario, metástasis y mucosa</td><td>Tanaka et al., Cell Reports 2024.</td><td><a href="https://doi.org/10.1016/j.celrep.2024.113810">DOI</a> · <a href="https://www.ebi.ac.uk/pride/archive/projects/PXD034575">PXD034575</a></td></tr>
      <tr><td>Suero con tres grupos</td><td>Rai et al., J Proteome Res 2025.</td><td><a href="https://doi.org/10.1021/acs.jproteome.5c00285">DOI</a> · <a href="https://www.ebi.ac.uk/pride/archive/projects/PXD055218">PXD055218</a></td></tr>
      <tr><td>t moderada</td><td>Smyth, Stat Appl Genet Mol Biol 2004.</td><td><a href="https://doi.org/10.2202/1544-6115.1027">DOI</a></td></tr>
      <tr><td>Umbral de cambio, treat</td><td>McCarthy y Smyth, Bioinformatics 2009.</td><td><a href="https://doi.org/10.1093/bioinformatics/btp053">DOI</a></td></tr>
      <tr><td>FDR</td><td>Benjamini y Hochberg, JRSS B 1995.</td><td><a href="https://doi.org/10.1111/j.2517-6161.1995.tb02031.x">DOI</a></td></tr>
      <tr><td>Cero de intensidad</td><td>Cox y Mann, Nat Biotechnol 2008. MaxQuant.</td><td><a href="https://doi.org/10.1038/nbt.1511">DOI</a></td></tr>
      <tr><td>Vocabulario GO</td><td>Ashburner et al., Nat Genet 2000.</td><td><a href="https://doi.org/10.1038/75556">DOI</a> · <a href="https://www.ebi.ac.uk/QuickGO/">QuickGO</a></td></tr>
      <tr><td>Enriquecimiento de vías</td><td>Subramanian et al., PNAS 2005.</td><td><a href="https://doi.org/10.1073/pnas.0506580102">DOI</a></td></tr>
      <tr><td>Concentración de referencia</td><td>Uhlén et al., Science 2015.</td><td><a href="https://doi.org/10.1126/science.1260419">DOI</a> · <a href="https://www.proteinatlas.org/">Human Protein Atlas</a></td></tr>
      <tr><td>Estructura predicha</td><td>Jumper et al., Nature 2021.</td><td><a href="https://doi.org/10.1038/s41586-021-03819-2">DOI</a> · <a href="https://alphafold.ebi.ac.uk/entry/O94812">AlphaFold de BAIAP3</a></td></tr>
      <tr><td>Identificador de proteína</td><td>UniProt.</td><td><a href="https://www.uniprot.org/uniprotkb/O94812">O94812</a></td></tr>
      <tr><td>Neoepítopos de tumor</td><td>Xiang et al., Cell Genomics 2026.</td><td><a href="https://doi.org/10.1016/j.xgen.2025.101062">DOI</a></td></tr>
      <tr><td>Péptidos de una persona</td><td>Hembrough et al., Cancer Res 2026, abstract 7710.</td><td><a href="https://doi.org/10.1158/1538-7445.am2026-7710">DOI</a></td></tr>
      <tr><td>Panel de plasma ya publicado</td><td>Ögren et al., Mol Cancer 2026. CEACAM5, FLT1, IL19 y ferritina.</td><td><a href="https://doi.org/10.1186/s12943-026-02654-1">DOI</a></td></tr>
      <tr><td>Origen del plasma</td><td>Malmström et al., Cell 2025.</td><td><a href="https://doi.org/10.1016/j.cell.2025.03.013">DOI</a></td></tr>
      <tr><td>Señal del huésped en vesículas</td><td>Seo et al., J Extracell Vesicles 2026.</td><td><a href="https://doi.org/10.1002/jev2.70278">DOI</a></td></tr>
      <tr><td>Candidatas que esta tesis no descubre</td><td>Yu et al., BMC Cancer 2026.</td><td><a href="https://doi.org/10.1186/s12885-026-15547-9">DOI</a></td></tr>
      <tr><td>Proteínas compartidas con enfermedad inflamatoria</td><td>Dhami et al., Proteomics Clin Appl 2026.</td><td><a href="https://doi.org/10.1002/prca.70041">DOI</a></td></tr>
      <tr><td>Fibronectina baja en plasma</td><td>Verathamjamras, Chantaraamporn et al., Clin Proteomics 2023.</td><td><a href="https://doi.org/10.1186/s12014-023-09407-y">DOI</a></td></tr>
      <tr><td>Plasma dirigido que no entra</td><td>Lesur et al., J Proteome Res 2023.</td><td><a href="https://doi.org/10.1021/acs.jproteome.2c00575">DOI</a></td></tr>
      <tr><td>Integración pública ya publicada</td><td>Robles, Prakash, Vizcaíno y Casal, PLOS Comput Biol 2024.</td><td><a href="https://doi.org/10.1371/journal.pcbi.1011828">DOI</a></td></tr>
      <tr><td>Otro proteoma de colon y recto, no PDC000116</td><td>Zhang et al., Nature 2014.</td><td><a href="https://doi.org/10.1038/nature13438">DOI</a></td></tr>
      <tr><td>Revisión narrativa, no proteómica</td><td>García-Foncillas, Rev Colomb Hematol Oncol 2026.</td><td><a href="https://doi.org/10.51643/22562915.853">DOI</a></td></tr>
      <tr><td>Plasma durante tratamiento, no cribado</td><td>Li et al., Nat Commun 2024.</td><td><a href="https://doi.org/10.1038/s41467-024-44911-1">DOI</a></td></tr>
    </tbody>
  </table>
"""


def _diccionario_md() -> str:
    return """# Diccionario y teoría

La explicación completa, con las tablas y los enlaces, está en [divulgacion.html](divulgacion.html), en las secciones Teoría, Datos, Procesos y Estudios.

Cada DOI de esa página se confirmó en Crossref el 31 de octubre de 2026: el título registrado coincide con la cita. El DOI 10.1016/j.xcrm.2025.101689 devolvió 404 y no se enlaza.

## Cómo está armada la página

1. El resultado, en una frase.
2. Las tres cohortes, cada una con el artículo y el depósito: [Vasaikar 2019](https://doi.org/10.1016/j.cell.2019.03.030) y [PDC000116](https://pdc.cancer.gov/pdc/study/PDC000116); [Tanaka 2024](https://doi.org/10.1016/j.celrep.2024.113810) y [PXD034575](https://www.ebi.ac.uk/pride/archive/projects/PXD034575); [Rai 2025](https://doi.org/10.1021/acs.jproteome.5c00285) y [PXD055218](https://www.ebi.ac.uk/pride/archive/projects/PXD055218).
3. Los números de tejido y de suero.
4. Las fichas, con [QuickGO](https://www.ebi.ac.uk/QuickGO/), el [Human Protein Atlas](https://www.proteinatlas.org/) y [AlphaFold](https://alphafold.ebi.ac.uk/entry/O94812).
5. La teoría de cada procedimiento: por qué no se suman las cohortes, qué es el cero de MaxQuant, el log2, el pareo, la t moderada, el FDR, el AUC, Gene Ontology y el límite de los neoantígenos.
6. El diccionario de datos y el de procesos.
7. Las preguntas para practicar la sustentación.

## Procedimientos, en corto

| Paso | Idea |
|---|---|
| Suero | Tres grupos de Rai. t moderada, umbral 0,585, FDR 0,05, sin imputar. |
| Especificidad | AUC mayor que 0,70, intervalo por encima de 0,50 y mayor que el AUC de enfermedad frente a sanos. La marca de fase aguda no baja el umbral. |
| Tanaka | P primario, M metástasis, N mucosa. El cero es un hueco. Luego log2 y contraste pareado. |
| Concordancia | Mismo signo y correlación contra la tabla certificada de CPTAC. 0,83 y 0,104 son estadísticas distintas. |
| Vías | Unos 12 términos GO del tejido tienen al menos 10 genes en el suero. No superan 200 permutaciones. |
| Fichas | GO, atlas de referencia, UniProt y AlphaFold. La estructura no valida la abundancia. |

Los métodos enlazados son [Smyth 2004](https://doi.org/10.2202/1544-6115.1027), [McCarthy y Smyth 2009](https://doi.org/10.1093/bioinformatics/btp053), [Benjamini y Hochberg 1995](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x), [Cox y Mann 2008](https://doi.org/10.1038/nbt.1511), [Ashburner et al. 2000](https://doi.org/10.1038/75556) y [Subramanian et al. 2005](https://doi.org/10.1073/pnas.0506580102).
"""


def _html(auditoria, fichas, nombres, tejido, vias, primario, metastasis, progresion, cptac_sig, especificas) -> str:
    barras_fc = barra(
        [
            ("SAA1, cáncer frente a sano", 3.2565411464710934),
            ("LBP, cáncer frente a sano", 1.254345833781091),
            ("CHMP4A, cáncer frente a sano", 5.419806415572392),
            ("APOA4, enfermedad frente a sano", -1.322996865727255),
        ]
    )
    ruido = {"fase_aguda": "fase aguda", "hemolisis": "hemólisis", "eii": "enfermedad inflamatoria"}
    barras_auc = barra(
        [
            (f"{g} · {ruido.get(r, r)}" if isinstance(r, str) and r.strip() else g, float(a))
            for g, r, a in especificas[["gene", "lista_ruido", "auc_crc_vs_enfermedad"]].itertuples(index=False)
        ]
    )
    tarjetas = []
    for _, f in fichas.sort_values("gene").iterrows():
        if f["lectura"] in {"no_evaluable", "inespecifica"} and f["origen"] != "resultado_propio":
            continue
        if f["veredicto"] == "inespecifica" and f["origen"] == "literatura":
            continue
        lectura = {
            "especifica": "AUC específico, sin marca de ruido",
            "auc_pasa_pero_fase_aguda": "el AUC pasa y ya estaba en fase aguda",
            "auc_pasa_pero_hemolisis": "el AUC pasa y ya estaba marcada por hemólisis",
            "auc_pasa_pero_eii": "el AUC pasa y ya estaba marcada en enfermedad inflamatoria",
            "pasa_fdr_en_algun_contraste": "pasa el FDR en un contraste",
            "inespecifica": "inespecífica",
        }.get(str(f["lectura"]), str(f["lectura"]))
        origen = "viene de la literatura" if f["origen"] == "literatura" else "sale de esta corrida"
        lineas = [
            f"Cáncer frente a sano: log2FC {coma(f.get('log2fc_crc_vs_sano'))}, FDR {coma(f.get('padj_crc_vs_sano'), 4)}, n {coma(f.get('n_obs_crc_vs_sano'), 0)}."
        ]
        if f.get("significativa_crc_vs_enfermedad") is True or (isinstance(f.get("auc_crc_vs_enfermedad"), float) and pd.notna(f.get("auc_crc_vs_enfermedad"))):
            extra = f"Cáncer frente a enfermedad: log2FC {coma(f.get('log2fc_crc_vs_enfermedad'))}, FDR {coma(f.get('padj_crc_vs_enfermedad'), 4)}"
            if pd.notna(f.get("auc_crc_vs_enfermedad")):
                extra += f", AUC {coma(f.get('auc_crc_vs_enfermedad'))}"
            lineas.append(extra + ".")
        if f.get("significativa_enfermedad_vs_sano") is True:
            lineas.append(
                f"Enfermedad frente a sano: log2FC {coma(f.get('log2fc_enfermedad_vs_sano'))}, FDR {coma(f.get('padj_enfermedad_vs_sano'), 4)}."
            )
        loc = f["hpa_localizacion"] if isinstance(f.get("hpa_localizacion"), str) else "Sin localización de secreción en el atlas"
        go = terminos_legibles(f.get("go"), nombres)
        tarjetas.append(
            f"""<article class="ficha">
            <h3>{f['gene']}</h3>
            <p class="meta">{lectura} · {origen}</p>
            <p>UniProt {f['uniprot'] or '—'}. {loc}.</p>
            {''.join(f'<p>{linea}</p>' for linea in lineas)}
            <p class="go">{go}</p>
            <p class="enlaces"><a href="{f['enlace_uniprot']}">UniProt</a> · <a href="{f['enlace_alphafold']}">AlphaFold</a></p>
            </article>"""
        )
    return f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>OncoGO-CRC · qué se ve del tumor en la sangre</title>
<style>
  :root {{ color-scheme: light; --ink: #1c2430; --paper: #f3efe6; --card: #fffdf8; --line: #d7d0c3; --teal: #0c5c5c; --rust: #8d3429; --nav: #142028; }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; font-family: "Iowan Old Style", Palatino, "Palatino Linotype", serif; background: var(--paper); color: var(--ink); line-height: 1.5; }}
  .avance {{ position: fixed; top: 0; left: 0; height: 3px; width: 100%; background: #e2b15a; transform: scaleX(0); transform-origin: 0 0; z-index: 30; }}
  .barra {{ position: sticky; top: 0; z-index: 20; display: flex; gap: 14px; align-items: center; background: var(--nav); color: #f4f1ea; padding: 10px 16px; overflow-x: auto; }}
  .barra a {{ color: #f4f1ea; text-decoration: none; font-family: "Avenir Next", "Segoe UI", sans-serif; font-size: 0.84rem; white-space: nowrap; }}
  .barra a:hover, .barra a:focus {{ color: #e2b15a; }}
  .marca {{ font-weight: 700; letter-spacing: 0.01em; }}
  #filtro {{ margin-left: auto; flex: 0 0 9.5rem; border: 0; border-radius: 999px; padding: 8px 12px; font-family: "Avenir Next", "Segoe UI", sans-serif; font-size: 0.86rem; }}
  header, main, footer {{ max-width: 1080px; margin: 0 auto; padding: 28px 20px; }}
  header {{ padding-bottom: 8px; }}
  h1 {{ font-size: clamp(1.9rem, 4vw, 2.7rem); line-height: 1.12; margin: 0 0 10px; letter-spacing: -0.02em; }}
  h2 {{ font-size: 1.45rem; margin: 2.2rem 0 0.7rem; padding-top: 1.2rem; border-top: 1px solid var(--line); scroll-margin-top: 4.2rem; }}
  h3 {{ scroll-margin-top: 4.2rem; }}
  .pregunta {{ font-size: 1.18rem; max-width: 42rem; }}
  .ojo {{ font-family: "Avenir Next", "Segoe UI", sans-serif; font-size: 0.82rem; letter-spacing: 0.04em; text-transform: uppercase; color: var(--teal); margin: 0 0 6px; }}
  .rejilla {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 12px; }}
  .dato, .aviso, .ficha, details {{ background: var(--card); border: 1px solid var(--line); border-radius: 12px; padding: 14px 16px; }}
  .dato strong {{ display: block; font-family: "Avenir Next", "Segoe UI", sans-serif; font-size: 1.7rem; letter-spacing: -0.03em; }}
  .aviso {{ border-left: 5px solid var(--rust); }}
  .grafica {{ width: 100%; height: auto; background: var(--card); border: 1px solid var(--line); border-radius: 12px; }}
  .eje {{ font-family: "Avenir Next", "Segoe UI", sans-serif; font-size: 12px; fill: var(--ink); }}
  .fichas {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 12px; }}
  .ficha h3 {{ margin: 0 0 6px; font-family: "Avenir Next", "Segoe UI", sans-serif; }}
  .meta, .go, .enlaces, footer, nav, .ojo {{ font-family: "Avenir Next", "Segoe UI", sans-serif; }}
  .meta, .go, .enlaces {{ font-size: 0.92rem; }}
  a {{ color: var(--teal); }}
  details {{ margin: 8px 0; }}
  summary {{ cursor: pointer; font-weight: 600; }}
  details p {{ margin: 0.6rem 0 0.2rem; }}
  ol li, ul li {{ margin: 0.45rem 0; }}
  .marco {{ overflow-x: auto; margin: 0.6rem 0 1.2rem; border: 1px solid var(--line); border-radius: 12px; background: var(--card); }}
  table {{ width: 100%; border-collapse: collapse; min-width: 640px; font-size: 0.92rem; }}
  th, td {{ border-bottom: 1px solid var(--line); padding: 8px 10px; text-align: left; vertical-align: top; }}
  th {{ background: #efe8dc; }}
  footer {{ font-size: 0.9rem; color: #3d4654; }}
  [hidden] {{ display: none !important; }}
</style>
</head>
<body>
<div class="avance" aria-hidden="true"></div>
<nav class="barra" aria-label="Secciones">
  <a class="marca" href="#resultado">OncoGO-CRC</a>
  <a href="#resultado">Resultado</a>
  <a href="#suero">Suero</a>
  <a href="#tejido">Tejido</a>
  <a href="#fichas">Fichas</a>
  <a href="#teoria">Teoría</a>
  <a href="#datos">Datos</a>
  <a href="#procesos">Procesos</a>
  <a href="#estudios">Estudios</a>
  <a href="#estudiar">Para sustentar</a>
  <input id="filtro" type="search" placeholder="Buscar" aria-label="Buscar proteína, técnica o estudio"/>
</nav>
<header>
  <p class="ojo">Universidad de Caldas · lectura pública · 31 de octubre de 2026</p>
  <h1>Qué parte del tumor colorrectal se deja ver en la sangre</h1>
  <p class="pregunta">La pregunta no es si hay proteínas distintas en el tumor. Esa parte ya está medida. La pregunta es cuáles de esas diferencias siguen ahí cuando la muestra es suero, y cuáles aparecen también en otra enfermedad del intestino.</p>
</header>
<main>
  <h2 id="resultado">Resultado de la tesis</h2>
  <div class="aviso">
    <p>El cambio de proteínas del tumor se repite entre tejidos. En este suero no aparece como una señal propia del cáncer. Sirve para responder la pregunta. No sirve para proponer un examen de sangre.</p>
  </div>
  <p>La página se lee de arriba abajo. Primero el resultado. Después las tres cohortes, cada una con su artículo y su depósito de datos. Luego los números, las fichas, la teoría de cada procedimiento, el diccionario de datos, el diccionario de procesos, la lista de estudios y las preguntas para practicar la sustentación. Cada DOI de esta página se confirmó en Crossref el 31 de octubre de 2026: el título registrado coincide con la cita. El identificador que no abre un artículo no está enlazado.</p>
  <section class="rejilla">
    <div class="dato"><strong>20 + 20 + 20</strong>Sueros de <a href="https://doi.org/10.1021/acs.jproteome.5c00285">Rai 2025</a>: cáncer, colonoscopia negativa y sanos. <a href="https://www.ebi.ac.uk/pride/archive/projects/PXD055218">PXD055218</a>.</div>
    <div class="dato"><strong>0</strong>Proteínas que separan el cáncer de la otra enfermedad intestinal con FDR 0,05 y |log2FC| &gt; 0,585. De 876, 145 tienen un cambio menor que 1,5 veces demostrado y 731 quedan indeterminadas.</div>
    <div class="dato"><strong>3</strong>Proteínas que separan el cáncer de los sanos con ese mismo umbral: SAA1, LBP y CHMP4A.</div>
    <div class="dato"><strong>{coma(primario['pearson'])}</strong>Correlación del primario de Tanaka con el tumor de CPTAC, en {primario['n']} genes. No es la cifra 0,104 de la versión 2.</div>
  </section>

  <h2>Tres lugares, tres preguntas</h2>
  <div class="rejilla">
    <article class="dato"><h3>Tejido de CPTAC</h3><p><a href="https://doi.org/10.1016/j.cell.2019.03.030">Vasaikar et al., <em>Cell</em> 2019</a>. Tumor y mucosa del mismo estudio, TMT. Es el ancla. Datos: <a href="https://pdc.cancer.gov/pdc/study/PDC000116">PDC000116</a>.</p></article>
    <article class="dato"><h3>Tejido de Tanaka</h3><p>16 pacientes. Columna P, tumor primario. Columna M, metástasis hepática. Columna N, mucosa normal. No hay pólipo. Artículo: <a href="https://doi.org/10.1016/j.celrep.2024.113810">Tanaka et al., <em>Cell Reports</em> 2024</a>. Datos: <a href="https://www.ebi.ac.uk/pride/archive/projects/PXD034575">PXD034575</a>.</p></article>
    <article class="dato"><h3>Suero de Rai</h3><p>DIA, tres grupos de 20. El grupo de colonoscopia negativa es el que decide si una proteína es del cáncer o de un intestino enfermo. Artículo: <a href="https://doi.org/10.1021/acs.jproteome.5c00285">Rai et al., <em>J Proteome Res</em> 2025</a>. Datos: <a href="https://www.ebi.ac.uk/pride/archive/projects/PXD055218">PXD055218</a>.</p></article>
  </div>

  <h2 id="tejido">La columna P no era un adenoma</h2>
  <div class="aviso">
    <p>La versión 2 llamó pólipo a la columna P y tumor a la columna M. El artículo de esas columnas describe tripletes de tumor primario, metástasis hepática y mucosa. Esta corrida usa esa leyenda. El cero de MaxQuant se trata como proteína no detectada y la intensidad se pasa a log2 antes del contraste pareado.</p>
    <p>Primario frente a mucosa: {cptac_sig['primario_vs_mucosa']['n_significativas']} proteínas de {cptac_sig['primario_vs_mucosa']['n_probadas']} pasan el umbral. Metástasis frente a mucosa: {cptac_sig['metastasis_vs_mucosa']['n_significativas']} de {cptac_sig['metastasis_vs_mucosa']['n_probadas']}. Metástasis frente a primario: {cptac_sig['metastasis_vs_primario']['n_significativas']} de {cptac_sig['metastasis_vs_primario']['n_probadas']}.</p>
    <p>Con CPTAC, el primario coincide de signo en {primario['mismo_signo']} de {primario['n']} genes ({coma(100 * primario['fraccion'], 1)} %) y la correlación es {coma(primario['pearson'])}. La metástasis coincide en {metastasis['mismo_signo']} de {metastasis['n']} ({coma(100 * metastasis['fraccion'], 1)} %), correlación {coma(metastasis['pearson'])}. Esa cuenta es el cierre del 31 de octubre, sobre todos los genes con efecto en las dos matrices.</p>
    <p>La corrida del 6 de noviembre, en <code>results/v3/k1_etiquetas.json</code>, usa las 688 proteínas certificadas. El primario frente a la mucosa coincide en 116 de 118 (intervalo 0,940 a 0,998; concordancia de Lin 0,734). La metástasis frente a la mucosa coincide en 93 de 106 (Lin 0,440). El 104 de 119 de la versión 2 tomó la metástasis como si fuera el tumor y se conserva con ese nombre.</p>
    <p>Entre primario y metástasis, el signo coincide en {progresion['mismo_signo_que_la_metastasis']} de {progresion['n']} genes. Solo el {coma(100 * progresion['fraccion_menor_entre_las_de_mismo_signo'], 1)} % de los que coinciden tiene magnitud menor en el primario. No queda un gradiente en el que el primario sea una versión pequeña de la metástasis.</p>
    <p>La correlación 0,104 y el 56 % de igual signo siguen en <code>results/reproducibility</code>. Salen de otro cálculo, sobre 1 993 proteínas, con la etiqueta antigua de columnas. No se sustituyen por estas cifras. Tampoco miden tejido contra sangre.</p>
  </div>

  <h2 id="suero">Lo que el suero separa</h2>
  <p>El eje es el log2FC. A la derecha, más abundante en el primer grupo del contraste. CHMP4A tiene 16 de 40 observaciones: el cambio es grande y el dato está a medias.</p>
  {barras_fc}
  <p>SAA1 (P0DJI8) log2FC {coma(3.2565411464710934)}, FDR {coma(0.0001997414304612, 4)}. LBP (P18428) log2FC {coma(1.254345833781091)}, FDR {coma(0.0069917556413199, 3)}. Las dos están anotadas como secretadas a la sangre. En la ficha de LBP aparece el término GO acute-phase response. La ficha de SAA1, con el mismo corte de términos, muestra quimiotaxis de neutrófilos y adhesión. Frente a la colonoscopia negativa no pasan el FDR. APOA4 baja en la enfermedad intestinal frente a los sanos (log2FC {coma(-1.322996865727255)}, FDR {coma(0.0222374286544479, 3)}) y no separa el cáncer de esa enfermedad.</p>

  <h2>AUC no es lo mismo que especificidad</h2>
  <p>Ocho candidatas de la literatura superan AUC 0,70 al separar cáncer de colonoscopia negativa, con el intervalo por encima de 0,50 y por encima del AUC de enfermedad frente a sano. Cinco ya estaban marcadas como fase aguda. Quedan tres sin esa marca: BAIAP3, FN1 y TPM4. Ninguna de las tres pasa el FDR. El AUC no dice si la proteína sube o baja. FN1 en este suero está un poco más baja en cáncer que en sanos (log2FC {coma(-0.2734827146432472)}).</p>
  {barras_auc}
  <p>BAIAP3, la proteína que Rai señaló en suero, tiene AUC {coma(0.8166666666666667)} y log2FC {coma(2.2032636469658655)} frente a la colonoscopia negativa, con FDR {coma(0.1991954563781721, 2)}. Es una señal de cribado, con 20 personas por grupo, no una validación.</p>

  <h2>Las vías del tumor casi no están en este suero</h2>
  <p>De los términos GO biológicos con FDR menor o igual a 0,25 en el tejido de CPTAC, solo entran al suero los que tienen al menos 10 genes medidos. Son {vias['crc_vs_sano']['n_terminos']} en cáncer frente a sano, {vias['crc_vs_enfermedad']['n_terminos']} en cáncer frente a enfermedad y {vias['enfermedad_vs_sano']['n_terminos']} en enfermedad frente a sano. La correlación de sus NES con el tejido es {coma(vias['crc_vs_sano']['pearson_nes'])}, {coma(vias['crc_vs_enfermedad']['pearson_nes'])} y {coma(vias['enfermedad_vs_sano']['pearson_nes'])}. El nulo, 200 permutaciones del ranking, da p {coma(vias['crc_vs_sano']['p_nulo'], 2)}, {coma(vias['crc_vs_enfermedad']['p_nulo'], 2)} y {coma(vias['enfermedad_vs_sano']['p_nulo'], 2)}. Es el cierre del 31 de octubre. La corrida del 6 de noviembre no repitió esa prueba con 10 000 permutaciones. No hay concordancia funcional que superar el azar. El suero de esta matriz no contiene el grueso de los procesos que el tumor altera.</p>

  <h2 id="fichas">Fichas</h2>
  <p>Cada término GO sale de la anotación humana ya usada por el pipeline, restringida a términos de 15 a 300 genes, y el nombre se leyó en <a href="https://www.ebi.ac.uk/QuickGO/">QuickGO</a> el 31 de octubre de 2026. Por ejemplo, <a href="https://www.ebi.ac.uk/QuickGO/term/GO:0006953">GO:0006953</a> es <em>acute-phase response</em>. La concentración es la del <a href="https://www.proteinatlas.org/">Human Protein Atlas</a> en personas de referencia, no en estos pacientes. <a href="https://alphafold.ebi.ac.uk/">AlphaFold</a> muestra una estructura predicha: no confirma que la abundancia haya cambiado. El artículo del atlas es <a href="https://doi.org/10.1126/science.1260419">Uhlén et al., <em>Science</em> 2015</a>. El de la estructura es <a href="https://doi.org/10.1038/s41586-021-03819-2">Jumper et al., <em>Nature</em> 2021</a>.</p>
  <div class="fichas">{''.join(tarjetas)}</div>

  {_explicacion_html()}
  {_guia_html()}

  <h2>Cómo leer un número de esta página</h2>
  <ol>
    <li>Si el FDR pasa y el |log2FC| pasa 0,585, la proteína entra en la lista corta. Hoy esa lista, en suero, tiene cuatro filas y un contraste distinto para APOA4.</li>
    <li>Si solo pasa el AUC, es una candidata de cribado. La bandera de fase aguda o de hemólisis se escribió antes del contraste y se muestra al lado.</li>
    <li>Un péptido mutado no está en estas filas. La matriz guarda grupos de proteínas. Rai describe 19 péptidos; Xiang, neoepítopos del tumor; Hembrough, cinco péptidos de una persona. Ninguno se confirma aquí.</li>
  </ol>
</main>
<footer>
  Corrida <code>tesis_cierre</code> del 31 de octubre de 2026. Semilla de los intervalos y del nulo: 20261031. Umbral treat 0,585, alfa 0,05. Los certificados de la versión 1 y de la versión 2 no se reescribieron. Los DOI de esta página se confirmaron ese día en Crossref.
</footer>
<script>
document.querySelectorAll("table").forEach((tabla) => {{
  const marco = document.createElement("div");
  marco.className = "marco";
  tabla.parentElement.insertBefore(marco, tabla);
  marco.appendChild(tabla);
}});
const filtro = document.querySelector("#filtro");
const avance = document.querySelector(".avance");
function aplicar() {{
  const t = (filtro.value || "").trim().toLowerCase();
  document.querySelectorAll("table tbody tr, .ficha, details").forEach((el) => {{
    el.hidden = t.length > 1 && !el.textContent.toLowerCase().includes(t);
  }});
}}
filtro.addEventListener("input", aplicar);
document.addEventListener("scroll", () => {{
  const h = document.documentElement;
  const p = h.scrollTop / (h.scrollHeight - h.clientHeight || 1);
  avance.style.transform = "scaleX(" + p + ")";
}}, {{ passive: true }});
document.querySelectorAll('a[href^="http"]').forEach((a) => {{
  a.target = "_blank";
  a.rel = "noopener noreferrer";
}});
</script>
</body>
</html>
"""


def _markdown(auditoria, primario, metastasis, progresion, vias, cptac_sig) -> str:
    return f"""# Qué se ve del tumor en la sangre

La página con las gráficas, las fichas, la teoría y los diccionarios está en [divulgacion.html](divulgacion.html). Lo que hay que estudiar está en [sustentacion.md](sustentacion.md). El mapa de la teoría y de los enlaces está en [diccionario.md](diccionario.md). Esta nota deja los números por escrito.

El suero de PXD055218 tiene tres grupos de 20. Con FDR 0,05 y |log2FC| mayor que 0,585, ninguna proteína separa el cáncer de la colonoscopia negativa. Frente a sanos pasan SAA1, LBP y CHMP4A. CHMP4A tiene 16 de 40 observaciones. APOA4 baja en la enfermedad intestinal frente a los sanos.

Ocho candidatas de la literatura superan la regla de AUC. Cinco están marcadas como fase aguda. BAIAP3, FN1 y TPM4 no llevan esa marca y tampoco pasan el FDR.

PXD034575, Tanaka et al. 2024, es tumor primario (P), metástasis hepática (M) y mucosa (N). No es un pólipo. Tras tratar el cero de MaxQuant como ausencia y pasar a log2, el primario coincide en signo con CPTAC en {primario['mismo_signo']} de {primario['n']} genes, correlación {coma(primario['pearson'])}. La metástasis correlaciona {coma(metastasis['pearson'])}. Esa cuenta es el cierre del 31 de octubre. La cifra 0,104 de la versión 2 es otro cálculo y sigue en los resultados de reproducibilidad.

La corrida del 6 de noviembre está en `results/v3/`. Entre las 688 proteínas certificadas, el primario frente a la mucosa coincide en signo en 116 de 118 (intervalo 0,940 a 0,998; Lin 0,734). La metástasis frente a la mucosa coincide en 93 de 106 (Lin 0,440). El suero tiene 876 proteínas. El corte de 1,5 veces sigue en cero entre cáncer y colonoscopia negativa. De esas 876, 145 tienen un cambio menor que 1,5 veces demostrado y 731 quedan indeterminadas. La ausencia solo se afirma en las 145.

De los términos GO del tejido, el suero solo comparte {vias['crc_vs_sano']['n_terminos']} con al menos 10 genes en cáncer frente a sano. Su correlación de NES no supera el nulo (p {coma(vias['crc_vs_sano']['p_nulo'], 2)}).

Primario frente a mucosa: {cptac_sig['primario_vs_mucosa']['n_significativas']} proteínas significativas de {cptac_sig['primario_vs_mucosa']['n_probadas']}. El signo entre primario y metástasis coincide en el {coma(100 * progresion['fraccion_mismo_signo'], 1)} % de {progresion['n']} genes.
"""


def _latex(auditoria, primario, metastasis, progresion, vias, cptac_sig, saa1, lbp, chmp, apoa, baia) -> str:
    return r"""\documentclass[11pt]{article}
\usepackage[spanish]{babel}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{geometry}
\usepackage{hyperref}
\geometry{margin=2.4cm}
\title{Borrador de tesis\\OncoGO-CRC: qué parte del proteoma tumoral se ve en el suero}
\author{Programa de Bioinformática y Biología Computacional\\Universidad de Caldas}
\date{31 de octubre de 2026}
\begin{document}
\maketitle
\begin{abstract}
Este borrador recoge solo cifras de la corrida \texttt{tesis\_cierre} del 31 de octubre de 2026. En el suero DIA de Rai y colaboradores (PXD055218; 20 cáncer, 20 colonoscopia negativa, 20 sanos) ninguna proteína separa el cáncer de la otra enfermedad intestinal con FDR 0,05 y cambio mínimo de 1,5 veces. Frente a sanos pasan SAA1, LBP y CHMP4A. CHMP4A se observa en 16 de 40 sueros. La columna P de PXD034575 es tumor primario y la M es metástasis hepática, no un pólipo. Con esa leyenda, el primario correlaciona """ + coma(primario["pearson"]) + r""" con CPTAC en """ + str(primario["n"]) + r""" genes. Las vías GO del tejido casi no están medidas en este suero, y las pocas que lo están no concuerdan por encima del azar.
\end{abstract}

\section{Pregunta}
El tumor colorrectal tiene un proteoma distinto del de la mucosa. La tesis pregunta cuánto de ese proteoma sigue siendo visible en sangre, y cuánto de lo visible en sangre es compartido con otra enfermedad intestinal. No propone un marcador para la consulta.

\section{Material y método}
Las cohortes no se fusionan. El contraste es una t moderada con umbral \texttt{treat} de 0,585 en escala log2, FDR de Benjamini--Hochberg 0,05, sin imputar ausencias \cite{smyth2004,mccarthy2009}. El suero usa los nombres de corrida CP, CN y Normal de PXD055218 \cite{rai2025}. PXD036594 no entra: un orden alfabético de columnas no es un diagnóstico \cite{lesur2023}.

PXD034575 se lee con la leyenda de Tanaka y colaboradores: P tumor primario, M metástasis hepática, N mucosa \cite{tanaka2024}. El cero de MaxQuant es ausencia. La intensidad positiva se transforma a log2 y el contraste es pareado por paciente. La versión 2 había llamado pólipo a P y tumor a M; esos certificados no se reescriben.

La especificidad por AUC exige separar cáncer de colonoscopia negativa con AUC mayor que 0,70, límite inferior del intervalo mayor que 0,50, y ese AUC mayor que el de enfermedad frente a sanos. Una proteína ya marcada como fase aguda, hemólisis o enfermedad inflamatoria no se presenta como específica aunque el AUC cumpla la regla. Las marcas salen del catálogo de literatura, no de este contraste.

Los términos GO se leen de la anotación humana del pipeline. La concordancia funcional toma los términos biológicos con FDR menor o igual a 0,25 en el GSEA de CPTAC y al menos 10 genes en el ranking de suero, y permuta ese ranking 200 veces.

\section{Resultados}
\subsection{Suero}
De unas 700 proteínas probadas, tres pasan el umbral en cáncer frente a sanos y una en enfermedad frente a sanos. Cáncer frente a colonoscopia negativa: cero. SAA1 sube (log2FC """ + coma(float(saa1.log2fc_crc_vs_sano)) + r""", FDR """ + coma(float(saa1.padj_crc_vs_sano), 4) + r"""). LBP sube (log2FC """ + coma(float(lbp.log2fc_crc_vs_sano)) + r""", FDR """ + coma(float(lbp.padj_crc_vs_sano), 3) + r"""). CHMP4A sube (log2FC """ + coma(float(chmp.log2fc_crc_vs_sano)) + r""", FDR """ + f"{float(chmp.padj_crc_vs_sano):.2e}".replace(".", ",") + r""") con """ + coma(float(chmp.n_obs_crc_vs_sano), 0) + r""" observaciones de 40. APOA4 baja en la enfermedad frente a los sanos (log2FC """ + coma(float(apoa.log2fc_enfermedad_vs_sano)) + r""", FDR """ + coma(float(apoa.padj_enfermedad_vs_sano), 3) + r""").

Ocho candidatas cumplen la regla de AUC. Cinco están en la lista de fase aguda (C9, LRG1, SERPINA1, HP, CRP). BAIAP3 queda sin esa lista: AUC """ + coma(float(baia.auc_crc_vs_enfermedad)) + r""", log2FC """ + coma(float(baia.log2fc_crc_vs_enfermedad)) + r""" frente a la enfermedad, FDR """ + coma(float(baia.padj_crc_vs_enfermedad), 2) + r""". FN1 y TPM4 también cumplen el AUC y no pasan el FDR. FN1 está ligeramente más baja en el cáncer. TPM4 se apoya en pocas observaciones del contraste contra la enfermedad.

\subsection{Tejido, con la leyenda corregida}
Primario frente a mucosa: """ + str(cptac_sig["primario_vs_mucosa"]["n_significativas"]) + r""" proteínas significativas de """ + str(cptac_sig["primario_vs_mucosa"]["n_probadas"]) + r""". Metástasis frente a mucosa: """ + str(cptac_sig["metastasis_vs_mucosa"]["n_significativas"]) + r""" de """ + str(cptac_sig["metastasis_vs_mucosa"]["n_probadas"]) + r""". Metástasis frente a primario: """ + str(cptac_sig["metastasis_vs_primario"]["n_significativas"]) + r""" de """ + str(cptac_sig["metastasis_vs_primario"]["n_probadas"]) + r""". CEACAM5 sube en el primario (log2FC 3,87, FDR 0,002, 13 pares).

El primario y CPTAC comparten signo en """ + str(primario["mismo_signo"]) + r""" de """ + str(primario["n"]) + r""" genes (""" + coma(100 * primario["fraccion"], 1) + r""" \%), correlación """ + coma(primario["pearson"]) + r""". La metástasis correlaciona """ + coma(metastasis["pearson"]) + r""" (""" + str(metastasis["mismo_signo"]) + r""" de """ + str(metastasis["n"]) + r"""). El signo del efecto primario y el de la metástasis coinciden en el """ + coma(100 * progresion["fraccion_mismo_signo"], 1) + r""" \% de """ + str(progresion["n"]) + r""" genes; entre esos, solo el """ + coma(100 * progresion["fraccion_menor_entre_las_de_mismo_signo"], 1) + r""" \% tiene magnitud menor en el primario.

La correlación 0,104 de la versión 2, sobre 1\,993 proteínas, permanece como está. Usó otra preparación y la etiqueta antigua de columnas. No mide la sangre.

\subsection{Vías}
Términos GO del tejido con al menos 10 genes en suero: """ + str(vias["crc_vs_sano"]["n_terminos"]) + r""", """ + str(vias["crc_vs_enfermedad"]["n_terminos"]) + r""" y """ + str(vias["enfermedad_vs_sano"]["n_terminos"]) + r""". Correlación de NES """ + coma(vias["crc_vs_sano"]["pearson_nes"]) + r""", """ + coma(vias["crc_vs_enfermedad"]["pearson_nes"]) + r""" y """ + coma(vias["enfermedad_vs_sano"]["pearson_nes"]) + r""". p del nulo """ + coma(vias["crc_vs_sano"]["p_nulo"], 2) + r""", """ + coma(vias["crc_vs_enfermedad"]["p_nulo"], 2) + r""" y """ + coma(vias["enfermedad_vs_sano"]["p_nulo"], 2) + r""". La concordancia funcional tejido--suero no se sostiene.

\section{Discusión}
Lo que separa el cáncer de una persona sana en este suero es, en dos de los tres casos, una proteína secretada de respuesta del huésped. Eso encaja con el atlas de origen plasmático de Malmström y con el hallazgo de Seo de que parte de la señal circulante viene del huésped, no de la célula tumoral \cite{malmstrom2025,seo2026}. No encaja con usar SAA1 o LBP como prueba de tumor.

BAIAP3 es la candidata de Rai que sobrevive a la regla de AUC sin estar en la lista de fase aguda, y no sobrevive al FDR. Con 20 por grupo el estudio es de cribado. Yu y colaboradores ya publicaron CA1 y SMPDL3A; aquí el suero no las confirma como específicas \cite{yu2026}.

El tejido, leído como primario y no como pólipo, sí se parece a CPTAC. La metástasis se parece menos. El relato de un adenoma intermedio no tiene apoyo en estas columnas.

Los 19 péptidos de Rai, los neoepítopos de Xiang y los cinco péptidos de Hembrough quedan fuera: una fila de grupo proteico no guarda la secuencia mutada \cite{xiang2026,hembrough2026}.

\section{Diez artículos que esta investigación puede sostener}
\begin{enumerate}
\item Especificidad del suero de cáncer colorrectal frente a colonoscopia negativa. El resultado central es el cero proteínas con FDR. \textit{Journal of Proteome Research} o \textit{Clinical Proteomics}.
\item Corrección de la leyenda de PXD034575 y concordancia del primario con CPTAC. La cifra 0,104 y la de esta corrida se publican como estadísticas distintas. \textit{PLOS Computational Biology}.
\item Cobertura de vías: el suero no contiene los términos GO del tumor, y los doce que contiene no concuerdan. \textit{Molecular \& Cellular Proteomics}.
\item SAA1 y LBP como respuesta del huésped cuando el control deja de ser una persona sana. Nota corta junto con el atlas de Malmström.
\item BAIAP3 entre el AUC y el FDR. Una proteína de la literatura, no un hallazgo de esta tesis.
\item FN1 y TPM4: AUC sin fase aguda y sin FDR, con el signo de FN1 hacia abajo, como en la proteómica de plasma de 2023 \cite{verathamjamras2023}.
\item CHMP4A y la ausencia: un log2FC alto con la mitad de los sueros sin medición no es un marcador.
\item CEACAM5 sube en el primario (log2FC 3,87) y el suero no dirigido no lo convierte en proteína específica. Dialoga con Ögren, cuyo panel de cáncer colorrectal es CEACAM5, FLT1, IL19 y ferritina \cite{ogren2026}.
\item Límite de los neoantígenos en matrices de grupos proteicos. Comentario metodológico, no un descubrimiento de péptidos.
\item Fichas públicas de las proteínas de esta corrida, con GO, concentración del Human Protein Atlas y enlace a AlphaFold. Nota de datos, no un biomarcador.
\end{enumerate}

\section{Qué hay que poder explicar en la sustentación}
El resultado, en una frase: el cambio de proteínas del tumor se repite entre tejidos y en este suero no aparece como señal propia del cáncer. La tesis responde la pregunta. No propone un examen de sangre.

El orden de estudio es decir esa frase, dibujar las tres cohortes, explicar LBP y BAIAP3, explicar la leyenda P, M y N, y separar la correlación """ + coma(primario["pearson"]) + r""" de la cifra 0,104. LBP pasa el FDR contra sanos y no contra la colonoscopia negativa. BAIAP3 tiene AUC """ + coma(float(baia.auc_crc_vs_enfermedad)) + r""" y FDR """ + coma(float(baia.padj_crc_vs_enfermedad), 2) + r""".

Hay que poder definir log2FC, el umbral 0,585 como log2 de 1,5, el FDR de Benjamini--Hochberg, la t moderada de \cite{smyth2004,mccarthy2009}, el AUC, el contraste pareado, el grupo proteico y la fase aguda. No se afirma un marcador clínico, ni que la columna P sea un pólipo, ni que 0,83 sustituya a 0,104, ni un neoantígeno, ni que las vías del tumor estén en este suero.

\section{Bibliografía citada en este borrador}
\begin{thebibliography}{99}
\bibitem{rai2025} Rai AB, et al. Data-independent acquisition analysis to identify serum biomarkers and neoantigens for differentiating colorectal cancer from inflammatory bowel diseases. \textit{J Proteome Res.} 2025;24:4166--4180. doi:10.1021/acs.jproteome.5c00285.
\bibitem{vasaikar2019} Vasaikar S, et al. Proteogenomic analysis of human colon cancer reveals new therapeutic opportunities. \textit{Cell.} 2019;177:1035--1049. doi:10.1016/j.cell.2019.03.030.
\bibitem{tanaka2024} Tanaka A, et al. Proteogenomic characterization of primary colorectal cancer and metastatic progression. \textit{Cell Reports.} 2024;43:113810. doi:10.1016/j.celrep.2024.113810.
\bibitem{zhang2014} Zhang B, et al. Proteogenomic characterization of human colon and rectal cancer. \textit{Nature.} 2014;513:382--387. doi:10.1038/nature13438.
\bibitem{robles2024} Robles J, Prakash A, Vizcaíno JA, Casal JI. Integrated meta-analysis of colorectal cancer public proteomic datasets. \textit{PLOS Comput Biol.} 2024;20:e1011828. doi:10.1371/journal.pcbi.1011828.
\bibitem{smyth2004} Smyth GK. Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. \textit{Stat Appl Genet Mol Biol.} 2004;3:Article 3.
\bibitem{mccarthy2009} McCarthy DJ, Smyth GK. Testing significance relative to a fold-change threshold is a TREAT. \textit{Bioinformatics.} 2009;25:765--771.
\bibitem{yu2026} Yu Y, et al. Proteome profiling reveals early diagnostic biomarker candidates for colorectal cancer. \textit{BMC Cancer.} 2026. doi:10.1186/s12885-026-15547-9.
\bibitem{seo2026} Seo Y, et al. Proteomic profiling of human extracellular vesicles reveals diagnostic biomarkers for colon adenocarcinoma. \textit{J Extracell Vesicles.} 2026. doi:10.1002/jev2.70278.
\bibitem{malmstrom2025} Malmström E, et al. Human proteome distribution atlas for tissue-specific plasma proteome dynamics. \textit{Cell.} 2025. doi:10.1016/j.cell.2025.03.013.
\bibitem{ogren2026} Ögren JÅ, et al. Composite proteomic and metabolomic plasma biomarkers for detection of colorectal, lung and ovarian cancers. \textit{Mol Cancer.} 2026. doi:10.1186/s12943-026-02654-1.
\bibitem{dhami2026} Dhami J, et al. A network-based association of IBD and colorectal cancer using proteomics data. \textit{Proteomics Clin Appl.} 2026. doi:10.1002/prca.70041.
\bibitem{xiang2026} Xiang H, et al. Predominant mutated non-canonical tumor-specific antigens in CRC. \textit{Cell Genomics.} 2026;6:101062. doi:10.1016/j.xgen.2025.101062.
\bibitem{hembrough2026} Hembrough T, et al. Extracellular vesicle immunopeptidomics recovers immunogenic neoantigen peptides from cancer patient plasma. \textit{Cancer Res.} 2026;86(7 Suppl):Abstract 7710. doi:10.1158/1538-7445.am2026-7710.
\bibitem{lesur2023} Lesur A, et al. Quantification of 782 plasma peptides by multiplexed targeted proteomics. \textit{J Proteome Res.} 2023;22:1630--1638. doi:10.1021/acs.jproteome.2c00575.
\bibitem{verathamjamras2023} Verathamjamras C, Chantaraamporn J, et al. Label-free quantitative proteomics of LRG, C9, FN, A1AT and AGP1 in colorectal cancer plasma. \textit{Clin Proteomics.} 2023. doi:10.1186/s12014-023-09407-y.
\bibitem{uhlen2015} Uhlén M, et al. Tissue-based map of the human proteome. \textit{Science.} 2015;347:1260419.
\bibitem{guinney2015} Guinney J, et al. The consensus molecular subtypes of colorectal cancer. \textit{Nat Med.} 2015;21:1350--1356.
\bibitem{garciafoncillas2026} García-Foncillas J. Biología molecular del cáncer colorrectal. \textit{Rev Colomb Hematol Oncol.} 2026. doi:10.51643/22562915.853.
\bibitem{subramanian2005} Subramanian A, et al. Gene set enrichment analysis. \textit{Proc Natl Acad Sci USA.} 2005;102:15545--15550.
\bibitem{jumper2021} Jumper J, et al. Highly accurate protein structure prediction with AlphaFold. \textit{Nature.} 2021;596:583--589.
\end{thebibliography}
\end{document}
"""


def _docx(auditoria, primario, metastasis, progresion, vias, cptac_sig, saa1, lbp, chmp, apoa, baia, fichas, nombres) -> None:
    try:
        from docx import Document
    except ImportError:
        import subprocess
        import sys

        subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "python-docx"])
        from docx import Document

    doc = Document()
    doc.add_heading("Borrador de tesis. OncoGO-CRC", 0)
    doc.add_paragraph("Qué parte del proteoma del tumor colorrectal se ve en el suero. 31 de octubre de 2026. Universidad de Caldas. Los números salen de la corrida tesis_cierre de esa fecha. Este archivo es un borrador para pasar a limpio, no el documento de grado firmado.")
    doc.add_heading("Resumen", level=1)
    doc.add_paragraph(
        "En el suero de Rai y colaboradores, tres grupos de 20, ninguna proteína separa el cáncer de la colonoscopia negativa "
        f"con FDR 0,05 y |log2FC| mayor que 0,585. Frente a sanos pasan SAA1, LBP y CHMP4A. "
        f"La columna P de PXD034575 es el tumor primario de Tanaka 2024, no un pólipo. "
        f"Ese primario correlaciona {coma(primario['pearson'])} con CPTAC en {primario['n']} genes. "
        "Las vías GO del tejido no se reproducen en este suero."
    )
    doc.add_heading("Pregunta", level=1)
    doc.add_paragraph("El tumor ya tiene un proteoma distinto de la mucosa. La tesis pregunta cuánto de esa diferencia llega a la sangre y cuánto de lo que llega a la sangre es compartido con otra enfermedad del intestino.")
    doc.add_heading("Método", level=1)
    for parrafo in [
        "Las cohortes no se fusionan. El contraste es una t moderada contra un cambio de 0,585 en log2, con FDR 0,05, sin imputar ausencias (Smyth 2004; McCarthy y Smyth 2009).",
        "El suero es PXD055218. CP es cáncer, CN es colonoscopia negativa, Normal es sano (Rai et al., J Proteome Res 2025, doi:10.1021/acs.jproteome.5c00285).",
        "PXD034575 se lee como Tanaka et al., Cell Reports 2024, doi:10.1016/j.celrep.2024.113810: P primario, M metástasis hepática, N mucosa. El cero de MaxQuant es ausencia. La versión 2 llamó pólipo a P; ese certificado no se reescribe.",
        "La regla de AUC pide más de 0,70 al separar cáncer de enfermedad, un intervalo cuyo límite inferior supere 0,50, y un AUC mayor que el de enfermedad frente a sanos. La marca de fase aguda, hemólisis o enfermedad inflamatoria se antepone a esa regla.",
        "La concordancia de vías usa términos GO biológicos con FDR de hasta 0,25 en CPTAC y al menos 10 genes en el suero, con 200 permutaciones del ranking.",
    ]:
        doc.add_paragraph(parrafo)
    doc.add_heading("Resultados", level=1)
    doc.add_heading("Suero", level=2)
    doc.add_paragraph(
        f"SAA1 sube frente a sanos (log2FC {coma(float(saa1.log2fc_crc_vs_sano))}, FDR {coma(float(saa1.padj_crc_vs_sano), 4)}, 40 observaciones). "
        f"LBP sube (log2FC {coma(float(lbp.log2fc_crc_vs_sano))}, FDR {coma(float(lbp.padj_crc_vs_sano), 3)}). "
        f"CHMP4A sube (log2FC {coma(float(chmp.log2fc_crc_vs_sano))}, FDR {float(chmp.padj_crc_vs_sano):.2e}) con {coma(float(chmp.n_obs_crc_vs_sano), 0)} de 40 observaciones. "
        f"APOA4 baja en la enfermedad frente a sanos (log2FC {coma(float(apoa.log2fc_enfermedad_vs_sano))}, FDR {coma(float(apoa.padj_enfermedad_vs_sano), 3)}). "
        "Cáncer frente a colonoscopia negativa: ninguna proteína pasa el umbral."
    )
    doc.add_paragraph(
        f"BAIAP3 tiene AUC {coma(float(baia.auc_crc_vs_enfermedad))} y log2FC {coma(float(baia.log2fc_crc_vs_enfermedad))} frente a la enfermedad, FDR {coma(float(baia.padj_crc_vs_enfermedad), 2)}. "
        "C9, LRG1, SERPINA1, HP y CRP cumplen el AUC y ya estaban marcadas como fase aguda. FN1 y TPM4 cumplen el AUC, no pasan el FDR, y FN1 apunta hacia abajo."
    )
    doc.add_heading("Tejido", level=2)
    doc.add_paragraph(
        f"Primario frente a mucosa: {cptac_sig['primario_vs_mucosa']['n_significativas']} de {cptac_sig['primario_vs_mucosa']['n_probadas']} proteínas. "
        f"Metástasis frente a mucosa: {cptac_sig['metastasis_vs_mucosa']['n_significativas']} de {cptac_sig['metastasis_vs_mucosa']['n_probadas']}. "
        f"Metástasis frente a primario: {cptac_sig['metastasis_vs_primario']['n_significativas']} de {cptac_sig['metastasis_vs_primario']['n_probadas']}. "
        "CEACAM5 sube en el primario (log2FC 3,87, FDR 0,002)."
    )
    doc.add_paragraph(
        f"Concordancia con CPTAC, primario: {primario['mismo_signo']} de {primario['n']} genes con el mismo signo, correlación {coma(primario['pearson'])}. "
        f"Metástasis: correlación {coma(metastasis['pearson'])}. "
        f"El signo entre primario y metástasis coincide en el {coma(100 * progresion['fraccion_mismo_signo'], 1)} %; "
        f"la magnitud es menor en el primario solo en el {coma(100 * progresion['fraccion_menor_entre_las_de_mismo_signo'], 1)} % de esos genes. "
        "La correlación 0,104 de la versión 2 es otra estadística y no describe la sangre."
    )
    doc.add_heading("Vías y fichas", level=2)
    doc.add_paragraph(
        f"Términos GO compartidos con al menos 10 genes: {vias['crc_vs_sano']['n_terminos']}, {vias['crc_vs_enfermedad']['n_terminos']} y {vias['enfermedad_vs_sano']['n_terminos']}. "
        f"Correlación de NES {coma(vias['crc_vs_sano']['pearson_nes'])}, {coma(vias['crc_vs_enfermedad']['pearson_nes'])} y {coma(vias['enfermedad_vs_sano']['pearson_nes'])}. "
        f"p del nulo {coma(vias['crc_vs_sano']['p_nulo'], 2)}, {coma(vias['crc_vs_enfermedad']['p_nulo'], 2)} y {coma(vias['enfermedad_vs_sano']['p_nulo'], 2)}. "
        "No hay concordancia funcional por encima del azar."
    )
    doc.add_paragraph("Las fichas de proteína, con el nombre del término GO leído en QuickGO, están en la página docs/web/divulgacion.html. La estructura de AlphaFold no valida una abundancia.")
    doc.add_heading("Diez artículos posibles", level=1)
    articulos = [
        "Especificidad del suero frente a colonoscopia negativa. El resultado es cero proteínas con FDR. Journal of Proteome Research o Clinical Proteomics.",
        "Leyenda de PXD034575: primario y metástasis, y la concordancia con CPTAC junto a la cifra histórica 0,104, explicadas como estadísticas distintas. PLOS Computational Biology.",
        "El suero no cubre las vías GO del tumor, y las que cubre no concuerdan con el nulo. Molecular & Cellular Proteomics.",
        "SAA1 y LBP como respuesta del huésped cuando el control es una persona sana. Nota con el atlas de Malmström 2025.",
        "BAIAP3 entre el AUC y el FDR. Auditoría de Rai 2025, no un descubrimiento.",
        "FN1 y TPM4 sin fase aguda y sin FDR. FN1 hacia abajo, en la dirección ya publicada en plasma en 2023.",
        "CHMP4A y los valores ausentes: un cambio grande con 16 de 40 observaciones no se presenta como marcador.",
        "CEACAM5 sube en el tejido primario y el suero no dirigido no lo hace específico. Diálogo con el panel de Ögren 2026 (CEACAM5, FLT1, IL19 y ferritina).",
        "Qué no puede decir una matriz de grupos proteicos sobre los 19 péptidos de Rai, los neoepítopos de Xiang y el abstract de Hembrough.",
        "Fichas públicas con GO, Human Protein Atlas y AlphaFold como capa de verificación. Nota de datos.",
    ]
    for texto in articulos:
        doc.add_paragraph(texto, style="List Number")
    doc.add_heading("Para la sustentación", level=1)
    doc.add_paragraph(
        "Resultado en una frase: el cambio de proteínas del tumor se repite entre tejidos y en este suero no aparece como señal propia del cáncer. "
        "Sirve para responder la pregunta. No sirve para proponer un examen de sangre."
    )
    doc.add_paragraph("Orden de una tarde de estudio:")
    for texto in [
        "Decir el resultado en voz alta, sin la hoja.",
        "Dibujar las tres cohortes: CPTAC como ancla, Tanaka como repetición en tejido, Rai como la pregunta de la sangre.",
        "Explicar LBP (FDR contra sanos, no contra la colonoscopia negativa) y BAIAP3 (AUC 0,82, FDR 0,20).",
        "Explicar que P es tumor primario, M metástasis hepática y N mucosa, y que el cero de MaxQuant es una ausencia.",
        f"Separar la correlación {coma(primario['pearson'])} de la cifra 0,104. Son estadísticas distintas. Ninguna mide tejido contra sangre.",
        "Cerrar con el límite: 20 personas por grupo, sin el inmunoensayo de CEA y sin secuencia de péptido mutado.",
    ]:
        doc.add_paragraph(texto, style="List Number")
    doc.add_heading("Palabras que hay que definir", level=2)
    for texto in [
        "log2FC: cambio en logaritmo base 2. 1 es el doble. El umbral 0,585 es log2(1,5).",
        "FDR: fracción esperada de falsos entre lo declarado distinto. Corte 0,05, Benjamini–Hochberg.",
        "treat: la t moderada pregunta si el cambio supera 0,585, no si es distinto de cero.",
        "AUC: si la proteína ordena al cáncer por encima de la otra condición. 0,50 es azar. No sustituye al FDR.",
        "Pareado: en Tanaka, tumor y mucosa son del mismo paciente. En el suero, los grupos son personas distintas.",
        "Grupo proteico: la fila del DIA no guarda la secuencia mutada de un péptido.",
        "Fase aguda: respuesta del organismo. C9, LRG1, SERPINA1, CRP y HP ya estaban marcadas así.",
    ]:
        doc.add_paragraph(texto, style="List Bullet")
    doc.add_heading("Respuestas cortas", level=2)
    for texto in [
        "Cero proteínas no es un fracaso del análisis. Es el resultado con el umbral escrito de antemano.",
        "Comparar solo contra sanos enciende SAA1 y LBP. Eso mide respuesta del organismo.",
        "SAA1 no distingue el cáncer de la otra enfermedad intestinal.",
        "BAIAP3 es una candidata de Rai, no un hallazgo de esta tesis.",
        "CHMP4A tiene un cambio grande y medición en 16 de 40 sueros.",
        "Los neoantígenos no se confirman con una matriz de grupos proteicos.",
        "Las cohortes no se fusionan porque TMT, LFQ y DIA no son la misma medición.",
    ]:
        doc.add_paragraph(texto, style="List Bullet")
    doc.add_paragraph("No se dice que haya un marcador de sangre listo, que la columna P sea un pólipo, que 0,83 sustituya a 0,104, ni que las vías del tumor se vean en este suero.")
    doc.add_heading("Bibliografía de este borrador", level=1)
    refs = [
        "Rai AB, et al. J Proteome Res. 2025;24:4166-4180. doi:10.1021/acs.jproteome.5c00285.",
        "Vasaikar S, et al. Cell. 2019;177:1035-1049. doi:10.1016/j.cell.2019.03.030. Ficha de PDC000116. No es el artículo de Cancer Cell 2023.",
        "Tanaka A, et al. Cell Reports. 2024;43:113810. doi:10.1016/j.celrep.2024.113810.",
        "Zhang B, et al. Nature. 2014;513:382-387. doi:10.1038/nature13438.",
        "Robles J, Prakash A, Vizcaíno JA, Casal JI. PLOS Comput Biol. 2024;20:e1011828. doi:10.1371/journal.pcbi.1011828.",
        "Smyth GK. Stat Appl Genet Mol Biol. 2004;3:Article 3.",
        "McCarthy DJ, Smyth GK. Bioinformatics. 2009;25:765-771.",
        "Yu Y, et al. BMC Cancer. 2026. doi:10.1186/s12885-026-15547-9.",
        "Seo Y, et al. J Extracell Vesicles. 2026. doi:10.1002/jev2.70278.",
        "Malmström E, et al. Cell. 2025. doi:10.1016/j.cell.2025.03.013.",
        "Ögren JÅ, et al. Mol Cancer. 2026. doi:10.1186/s12943-026-02654-1. Panel de CCR: CEACAM5, FLT1, IL19 y ferritina, AUC 0,89.",
        "Dhami J, et al. Proteomics Clin Appl. 2026. doi:10.1002/prca.70041.",
        "Xiang H, et al. Cell Genomics. 2026;6:101062. doi:10.1016/j.xgen.2025.101062.",
        "Hembrough T, et al. Cancer Res. 2026;86(7 Suppl):Abstract 7710. doi:10.1158/1538-7445.am2026-7710.",
        "Lesur A, et al. J Proteome Res. 2023;22:1630-1638. doi:10.1021/acs.jproteome.2c00575.",
        "Verathamjamras C, Chantaraamporn J, et al. Clin Proteomics. 2023. doi:10.1186/s12014-023-09407-y.",
        "Uhlén M, et al. Science. 2015;347:1260419.",
        "Guinney J, et al. Nat Med. 2015;21:1350-1356.",
        "García-Foncillas J. Rev Colomb Hematol Oncol. 2026. doi:10.51643/22562915.853. Revisión narrativa, no proteómica.",
        "Subramanian A, et al. Proc Natl Acad Sci USA. 2005;102:15545-15550.",
        "Jumper J, et al. Nature. 2021;596:583-589.",
        "Li Y, et al. Nat Commun. 2024. doi:10.1038/s41467-024-44911-1. Plasma durante cetuximab, no cribado.",
        "Franzetti F, et al. Proteomes. 2026. doi:10.3390/proteomes14020016.",
        "Zhang X, et al. J Natl Cancer Inst. 2026. doi:10.1093/jnci/djag217. Plasma prediagnóstico.",
    ]
    for texto in refs:
        doc.add_paragraph(texto, style="List Bullet")
    doc.add_paragraph("No entran en esta bibliografía, porque el identificador no abrió un artículo el 31 de octubre de 2026: doi:10.1016/j.xcrm.2025.101689, la cita «ENCODE, Nature 12564», y el DOI 10.64898/2026.08.30.26361622 de la memoria del 30 de octubre. Los preprints de Lee, Tan y Ramsbottom pueden ir a un anexo marcados como no revisados por pares.")
    doc.save(BORRADOR / "tesis_borrador.docx")


if __name__ == "__main__":
    main()
