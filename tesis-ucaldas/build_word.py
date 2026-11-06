"""Construye el .docx a partir de la plantilla oficial y de los capítulos LaTeX."""

from __future__ import annotations

import re
import shutil
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Emu, Pt, RGBColor

RAIZ = Path(__file__).resolve().parent
LATEX = RAIZ / "latex"
PLANTILLA = RAIZ / "assets" / "guia-presentacion-tesis.docx"
SALIDA = RAIZ / "word" / "tesis-maestria-ucaldas.docx"
BIB = LATEX / "references.bib"

ARCHIVOS = [
    "chapters/03-resumen.tex",
    "chapters/04-abstract.tex",
    "chapters/05-preliminares.tex",
    "chapters/06-introduccion.tex",
    "chapters/07-marco-teorico.tex",
    "chapters/08-metodos.tex",
    "chapters/09-resultados.tex",
    "chapters/90-conclusiones.tex",
    "chapters/91-anexos.tex",
]


def leer_bib(ruta: Path) -> dict[str, dict[str, str]]:
    entradas: dict[str, dict[str, str]] = {}
    clave = None
    campos: dict[str, str] = {}
    for linea in ruta.read_text(encoding="utf-8").splitlines():
        inicio = re.match(r"@\w+\{(\w+),", linea)
        if inicio:
            clave = inicio.group(1)
            campos = {}
            continue
        if linea.strip() == "}" and clave:
            entradas[clave] = campos
            clave = None
            continue
        campo = re.match(r"\s*(\w+)\s*=\s*\{(.*)\},?\s*$", linea)
        if campo and clave:
            valor = campo.group(2).replace(r"\_", "_")
            while valor.startswith("{") and valor.endswith("}"):
                valor = valor[1:-1]
            campos[campo.group(1)] = valor
    return entradas


def apellidos(author: str) -> list[str]:
    personas = [p.strip() for p in author.split(" and ") if p.strip()]
    salida = []
    for persona in personas:
        salida.append(persona.split(",")[0].strip() if "," in persona else persona)
    return salida


def cita(entrada: dict[str, str], textual: bool) -> str:
    apes = apellidos(entrada.get("author", ""))
    if len(apes) == 1:
        nombre = apes[0]
    elif len(apes) == 2:
        nombre = f"{apes[0]} y {apes[1]}"
    else:
        nombre = f"{apes[0]} et al."
    anio = entrada.get("year", "s. f.")
    if textual:
        return f"{nombre} ({anio})"
    return f"({nombre}, {anio})"


def persona_apa(persona: str) -> str:
    persona = persona.strip()
    if "," not in persona:
        return persona
    apellido, nombres = persona.split(",", 1)
    iniciales = " ".join(f"{parte[0].upper()}." for parte in nombres.split() if parte)
    return f"{apellido.strip()}, {iniciales}".strip().rstrip(",")


def referencia_apa(campos: dict[str, str]) -> str:
    personas = [p.strip() for p in campos.get("author", "").split(" and ") if p.strip()]
    if len(personas) > 20:
        formateadas = [persona_apa(p) for p in personas[:19]] + ["..."] + [persona_apa(personas[-1])]
        autores = ", ".join(formateadas)
    elif len(personas) == 1:
        autores = persona_apa(personas[0])
    else:
        autores = ", ".join(persona_apa(p) for p in personas[:-1]) + " y " + persona_apa(personas[-1])
    anio = campos.get("year", "s. f.")
    titulo = campos.get("title", "")
    revista = campos.get("journal") or campos.get("booktitle") or ""
    volumen = campos.get("volume", "")
    numero = campos.get("number", "")
    paginas = campos.get("pages", "")
    doi = campos.get("doi", "")
    pieza = f"{autores} ({anio}). {titulo}."
    if revista:
        pieza += f" {revista}"
        if volumen:
            pieza += f", {volumen}"
        if numero:
            pieza += f"({numero})"
        if paginas:
            pieza += f", {paginas}"
        pieza += "."
    if doi:
        pieza += f" https://doi.org/{doi}"
    return pieza


def rfonts(run, nombre: str) -> None:
    run.font.name = nombre
    rpr = run._element.get_or_add_rPr()
    fuentes = rpr.find(qn("w:rFonts"))
    if fuentes is None:
        fuentes = OxmlElement("w:rFonts")
        rpr.append(fuentes)
    for atributo in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        fuentes.set(qn(atributo), nombre)


def campo(parrafo, instruccion: str, reserva: str = "") -> None:
    run = parrafo.add_run()
    rfonts(run, "Arial")
    run.font.size = Pt(11)

    def fld(tipo: str):
        nodo = OxmlElement("w:fldChar")
        nodo.set(qn("w:fldCharType"), tipo)
        return nodo

    run._r.append(fld("begin"))
    texto = OxmlElement("w:instrText")
    texto.set(qn("xml:space"), "preserve")
    texto.text = instruccion
    run._r.append(texto)
    run._r.append(fld("separate"))
    visible = parrafo.add_run(reserva)
    rfonts(visible, "Arial")
    visible.font.size = Pt(11)
    cierre = parrafo.add_run()
    rfonts(cierre, "Arial")
    cierre._r.append(fld("end"))


def poner_numero(header, alineacion) -> None:
    header.is_linked_to_previous = False
    parrafo = header.paragraphs[0]
    parrafo.clear()
    parrafo.alignment = alineacion
    campo(parrafo, " PAGE ", "1")


def vaciar_header(header) -> None:
    header.is_linked_to_previous = False
    header.paragraphs[0].clear()


def numeracion(section, fmt: str, inicio: int | None) -> None:
    sect = section._sectPr
    nodo = sect.find(qn("w:pgNumType"))
    if nodo is None:
        nodo = OxmlElement("w:pgNumType")
        sect.append(nodo)
    nodo.set(qn("w:fmt"), fmt)
    if inicio is None:
        if qn("w:start") in nodo.attrib:
            del nodo.attrib[qn("w:start")]
    else:
        nodo.set(qn("w:start"), str(inicio))


def margenes(section) -> None:
    section.page_width = Cm(21.59)
    section.page_height = Cm(27.94)
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)
    section.left_margin = Cm(3.60)
    section.right_margin = Cm(2.54)
    section.header_distance = Cm(1.25)
    section.footer_distance = Cm(1.25)
    section.different_first_page_header_footer = False


def bordes(tabla) -> None:
    tbl = tabla._tbl
    tblPr = tbl.tblPr if tbl.tblPr is not None else OxmlElement("w:tblPr")
    existentes = tblPr.find(qn("w:tblBorders"))
    if existentes is not None:
        tblPr.remove(existentes)
    marco = OxmlElement("w:tblBorders")
    for lado in ("top", "left", "bottom", "right", "insideH", "insideV"):
        nodo = OxmlElement(f"w:{lado}")
        nodo.set(qn("w:val"), "single")
        nodo.set(qn("w:sz"), "4")
        nodo.set(qn("w:space"), "0")
        nodo.set(qn("w:color"), "000000")
        marco.append(nodo)
    tblPr.append(marco)


def activar_pares(doc: Document) -> None:
    ajustes = doc.settings.element
    for etiqueta in ("w:evenAndOddHeaders", "w:mirrorMargins", "w:updateFields"):
        if ajustes.find(qn(etiqueta)) is None:
            nodo = OxmlElement(etiqueta)
            if etiqueta == "w:updateFields":
                nodo.set(qn("w:val"), "true")
            ajustes.append(nodo)


def preparar_estilos(doc: Document) -> None:
    normal = doc.styles["Normal"]
    normal.font.name = "Arial"
    normal.font.size = Pt(11)
    normal.font.color.rgb = RGBColor(0, 0, 0)
    pf = normal.paragraph_format
    pf.line_spacing = 1.5
    pf.space_before = Pt(0)
    pf.space_after = Pt(16)
    pf.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    especificacion = {
        "Heading 1": (20, True, 100, 24, 0),
        "Heading 2": (16, False, 20, 24, 1),
        "Heading 3": (14, True, 12, 10, 2),
    }
    for nombre, (tamano, negrita, antes, despues, nivel) in especificacion.items():
        estilo = doc.styles[nombre]
        estilo.font.name = "Arial"
        estilo.font.size = Pt(tamano)
        estilo.font.bold = negrita
        estilo.font.color.rgb = RGBColor(0, 0, 0)
        formato = estilo.paragraph_format
        formato.space_before = Pt(antes)
        formato.space_after = Pt(despues)
        formato.line_spacing_rule = WD_LINE_SPACING.SINGLE
        formato.line_spacing = 1.0
        formato.page_break_before = False
        formato.keep_with_next = True
        ppr = estilo.element.get_or_add_pPr()
        outline = ppr.find(qn("w:outlineLvl"))
        if outline is None:
            outline = OxmlElement("w:outlineLvl")
            ppr.append(outline)
        outline.set(qn("w:val"), str(nivel))


class Constructor:
    def __init__(self, doc: Document, entradas: dict[str, dict[str, str]]) -> None:
        self.doc = doc
        self.entradas = entradas
        self.fase = "preliminares"
        self.capitulo = 0
        self.seccion = 0
        self.tabla = 0
        self.figura = 0
        self.ecuacion = 0
        self.en_lista = False
        self.en_ingles = False

    def salto_impar(self, reiniciar: int | None, fmt: str) -> None:
        seccion = self.doc.add_section(WD_SECTION.ODD_PAGE)
        margenes(seccion)
        poner_numero(seccion.header, WD_ALIGN_PARAGRAPH.RIGHT)
        poner_numero(seccion.even_page_header, WD_ALIGN_PARAGRAPH.LEFT)
        vaciar_header(seccion.footer)
        vaciar_header(seccion.even_page_footer)
        numeracion(seccion, fmt, reiniciar)

    def portada(self) -> None:
        self._centrado("Qué parte del tumor colorrectal se deja ver en la sangre", 16, True, antes=80)
        self._centrado("Juan Diego Polanco Franco", 14, False, antes=28)
        self._centrado("Universidad de Caldas", 14, False, antes=28)
        self._centrado("Facultad de Inteligencia Artificial e Ingenierías", 12, False, antes=12)
        self._centrado("Manizales, Colombia", 12, False, antes=28)
        self._centrado("2026", 12, False, antes=6)
        self.doc.add_page_break()
        self._centrado("Qué parte del tumor colorrectal se deja ver en la sangre", 16, True, antes=36)
        self._centrado("Juan Diego Polanco Franco", 14, False, antes=18)
        parrafo = self.doc.add_paragraph()
        parrafo.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = parrafo.add_run(
            "Tesis o trabajo de investigación presentada(o) como requisito parcial para optar al título de:"
        )
        rfonts(run, "Arial")
        run.font.size = Pt(11)
        self._centrado("Magíster en Bioinformática y Biología Computacional", 12, True, antes=8)
        nota = self.doc.add_paragraph()
        nota.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = nota.add_run(
            "Programa de investigación, 60 créditos, SNIES 103917. Convenio SUMA con la Universidad de Manizales, "
            "la Universidad Autónoma de Manizales y la Universidad Católica de Manizales."
        )
        rfonts(run, "Arial")
        run.font.size = Pt(10)
        self._etiqueta("Director (a):", "[título académico y nombre completo del director]")
        self._etiqueta("Codirector (a):", "[título académico y nombre completo del codirector]")
        self._etiqueta("Grupo de investigación:", "[nombre del grupo de investigación]")
        self._etiqueta("Línea de investigación:", "[nombre de la línea de investigación]")
        self._centrado("Universidad de Caldas", 12, False, antes=18)
        self._centrado("Facultad de Inteligencia Artificial e Ingenierías", 12, False)
        self._centrado("Manizales, Colombia", 12, False)
        self._centrado("2026", 12, False)

    def _centrado(self, texto: str, tamano: int, negrita: bool, antes: int = 0) -> None:
        parrafo = self.doc.add_paragraph()
        parrafo.alignment = WD_ALIGN_PARAGRAPH.CENTER
        parrafo.paragraph_format.space_before = Pt(antes)
        parrafo.paragraph_format.space_after = Pt(0)
        run = parrafo.add_run(texto)
        rfonts(run, "Arial")
        run.font.size = Pt(tamano)
        run.bold = negrita

    def _etiqueta(self, etiqueta: str, hueco: str) -> None:
        parrafo = self.doc.add_paragraph()
        parrafo.alignment = WD_ALIGN_PARAGRAPH.LEFT
        parrafo.paragraph_format.space_after = Pt(2)
        run = parrafo.add_run(etiqueta + " ")
        rfonts(run, "Arial")
        run.font.size = Pt(11)
        run.bold = True
        vacio = parrafo.add_run(hueco)
        rfonts(vacio, "Arial")
        vacio.font.size = Pt(11)
        vacio.italic = True
        vacio.font.color.rgb = RGBColor(0x80, 0x80, 0x80)

    def titulo(self, texto: str, nivel: int, impar: bool = False) -> None:
        if impar and self.fase == "cuerpo":
            self.salto_impar(None, "decimal")
        elif impar and self.fase == "preliminares":
            self.salto_impar(None, "upperRoman")
        self.doc.add_heading(texto, level=nivel)

    def procesar(self, ruta: Path) -> None:
        lineas = [self._sin_comentario(linea) for linea in ruta.read_text(encoding="utf-8").splitlines()]
        buffer: list[str] = []
        i = 0
        while i < len(lineas):
            linea = lineas[i].strip()
            if linea.startswith(r"\begin{"):
                self._volcar(buffer)
                buffer = []
                entorno = re.match(r"\\begin\{(\w+)\}", linea).group(1)
                bloque = [linea]
                i += 1
                while i < len(lineas) and not lineas[i].strip().startswith(f"\\end{{{entorno}}}"):
                    bloque.append(lineas[i])
                    i += 1
                if i < len(lineas):
                    bloque.append(lineas[i])
                self._entorno(entorno, bloque)
            elif linea.startswith(r"\chapter"):
                self._volcar(buffer)
                buffer = []
                self._capitulo(linea)
            elif linea.startswith(r"\section"):
                self._volcar(buffer)
                buffer = []
                self._seccion(linea, 2)
            elif linea.startswith(r"\subsection"):
                self._volcar(buffer)
                buffer = []
                self._seccion(linea, 3)
            elif linea in {r"\clearpage", r"\cleardoublepage"}:
                self._volcar(buffer)
                buffer = []
            elif linea == r"\mainmatter":
                self._volcar(buffer)
                buffer = []
                self.fase = "cuerpo"
                self.salto_impar(1, "decimal")
            elif linea == r"\tableofcontents":
                self._volcar(buffer)
                buffer = []
                self.titulo("Tabla de contenido", 1, impar=True)
                campo(self.doc.add_paragraph(), r' TOC \o "1-3" \h \z \u ', "Actualice los campos para ver el contenido.")
            elif linea == r"\listoftables":
                self._volcar(buffer)
                buffer = []
                self.titulo("Lista de tablas", 1, impar=True)
                campo(self.doc.add_paragraph(), r' TOC \h \z \c "Tabla" ', "Actualice los campos para ver la lista de tablas.")
            elif linea == r"\listoffigures":
                self._volcar(buffer)
                buffer = []
                self.titulo("Lista de figuras", 1, impar=True)
                campo(self.doc.add_paragraph(), r' TOC \h \z \c "Figura" ', "Actualice los campos para ver la lista de figuras.")
            elif linea.startswith(r"\addcontentsline") or linea.startswith(r"\label") or not linea:
                if not linea:
                    self._volcar(buffer)
                    buffer = []
            else:
                buffer.append(linea)
            i += 1
        self._volcar(buffer)

    def _sin_comentario(self, linea: str) -> str:
        salida = []
        escapado = False
        for caracter in linea:
            if caracter == "%" and not escapado:
                break
            salida.append(caracter)
            escapado = caracter == "\\" and not escapado
            if caracter != "\\":
                escapado = False
        return "".join(salida).rstrip()

    def _capitulo(self, linea: str) -> None:
        estrella = linea.startswith(r"\chapter*")
        titulo = re.search(r"\{(.*)\}", linea).group(1)
        if estrella:
            self.titulo(titulo, 1, impar=True)
        else:
            self.capitulo += 1
            self.seccion = 0
            self.titulo(f"{self.capitulo}. {titulo}", 1, impar=True)

    def _seccion(self, linea: str, nivel: int) -> None:
        titulo = re.search(r"\{(.*)\}", linea).group(1)
        if "\\section*{" in linea or "\\subsection*{" in linea:
            self.titulo(titulo, nivel, impar=False)
            return
        if nivel == 2:
            self.seccion += 1
            self.titulo(f"{self.capitulo}.{self.seccion} {titulo}", nivel)
        else:
            self.titulo(titulo, nivel)

    def _volcar(self, buffer: list[str]) -> None:
        texto = " ".join(parte.strip() for parte in buffer if parte.strip())
        if texto:
            parrafo = self.doc.add_paragraph()
            parrafo.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            self._runs(parrafo, texto)

    def _runs(self, parrafo, texto: str, tamano: int = 11) -> None:
        patron = re.compile(
            r"\\parencite\{(\w+)\}|\\textcite\{(\w+)\}|\\textbf\{([^{}]*)\}|\\placeholder\{([^{}]*)\}|"
            r"\\eqref\{[^}]+\}|\\ref\{[^}]+\}|\\fuente\{([^{}]*)\}"
        )
        cursor = 0
        for hallazgo in patron.finditer(texto):
            if hallazgo.start() > cursor:
                self._plano(parrafo, texto[cursor:hallazgo.start()], tamano)
            if hallazgo.group(1):
                self._plano(parrafo, cita(self.entradas[hallazgo.group(1)], False), tamano)
            elif hallazgo.group(2):
                self._plano(parrafo, cita(self.entradas[hallazgo.group(2)], True), tamano)
            elif hallazgo.group(3) is not None:
                self._plano(parrafo, hallazgo.group(3), tamano, negrita=True)
            elif hallazgo.group(4) is not None:
                self._plano(parrafo, f"[{hallazgo.group(4)}]", tamano, cursiva=True, gris=True)
            elif hallazgo.group(0).startswith(r"\eqref"):
                self._plano(parrafo, f"({self.ecuacion or 1})", tamano)
            elif hallazgo.group(0).startswith(r"\ref"):
                self._plano(parrafo, str(self.figura or self.tabla or ""), tamano)
            elif hallazgo.group(5) is not None:
                self._plano(parrafo, "Fuente: " + hallazgo.group(5), 10, cursiva=True)
            cursor = hallazgo.end()
        if cursor < len(texto):
            self._plano(parrafo, texto[cursor:], tamano)

    def _plano(self, parrafo, texto: str, tamano: int, negrita: bool = False, cursiva: bool = False, gris: bool = False) -> None:
        texto = (
            texto.replace(r"\%", "%")
            .replace(r"\,", "")
            .replace("~", " ")
            .replace("---", "—")
            .replace("--", "–")
            .replace(r"\log_{2}", "log2")
            .replace("{,}", ",")
            .replace("$", "")
            .replace(r"\&", "&")
        )
        texto = re.sub(r"\\[a-zA-Z]+\*?", "", texto)
        texto = texto.replace("{", "").replace("}", "")
        if not texto:
            return
        run = parrafo.add_run(texto)
        rfonts(run, "Arial")
        run.font.size = Pt(tamano)
        run.bold = negrita
        run.italic = cursiva
        if gris:
            run.font.color.rgb = RGBColor(0x80, 0x80, 0x80)

    def _entorno(self, nombre: str, bloque: list[str]) -> None:
        cuerpo = "\n".join(bloque)
        if nombre == "itemize":
            for item in re.findall(r"\\item\s+(.*)", cuerpo):
                parrafo = self.doc.add_paragraph()
                parrafo.paragraph_format.left_indent = Cm(1.0)
                parrafo.paragraph_format.space_after = Pt(4)
                self._runs(parrafo, "• " + item.strip())
        elif nombre == "equation":
            self.ecuacion += 1
            formula = re.search(r"\\begin\{equation\}(.*)\\label", cuerpo, re.S)
            texto = formula.group(1) if formula else cuerpo
            texto = texto.replace(r"\log_{2}", "log2").replace("{,}", ",").replace("$", "")
            texto = re.sub(r"\\[a-zA-Z]+", "", texto).replace("{", "").replace("}", "").strip()
            tabla = self.doc.add_table(rows=1, cols=2)
            tabla.autofit = True
            izquierda = tabla.cell(0, 0).paragraphs[0]
            izquierda.alignment = WD_ALIGN_PARAGRAPH.CENTER
            run = izquierda.add_run(texto)
            rfonts(run, "Arial")
            run.font.size = Pt(11)
            run.italic = True
            derecha = tabla.cell(0, 1).paragraphs[0]
            derecha.alignment = WD_ALIGN_PARAGRAPH.RIGHT
            numero = derecha.add_run(f"({self.ecuacion})")
            rfonts(numero, "Arial")
            numero.font.size = Pt(11)
            self.doc.add_paragraph()
        elif nombre in {"table", "longtable"}:
            self._tabla(cuerpo, nombre == "longtable")
        elif nombre == "figure":
            self._figura(cuerpo)
        elif nombre == "tabular":
            self._filas(cuerpo, None)
        elif nombre == "otherlanguage":
            interno = re.sub(r"\\begin\{otherlanguage\}\{english\}", "", cuerpo)
            interno = interno.replace(r"\end{otherlanguage}", "")
            temporal = LATEX / "_abstract_tmp.tex"
            temporal.write_text(interno, encoding="utf-8")
            self.procesar(temporal)
            temporal.unlink()
        elif nombre in {"center", "minipage"}:
            texto = re.sub(r"\\begin\{[^}]+\}(\[[^\]]*\])?(\{[^}]*\})?", "", cuerpo)
            texto = texto.replace(r"\end{center}", "").replace(r"\end{minipage}", "")
            for parte in texto.split(r"\par"):
                limpio = parte.strip()
                if limpio:
                    parrafo = self.doc.add_paragraph()
                    parrafo.alignment = WD_ALIGN_PARAGRAPH.CENTER
                    self._runs(parrafo, limpio)

    def _tabla(self, cuerpo: str, larga: bool) -> None:
        caption = re.search(r"\\caption\{([^{}]*)\}", cuerpo)
        if caption and not caption.group(0).startswith(r"\caption*"):
            self.tabla += 1
            self._leyenda("Tabla", self.tabla, caption.group(1))
        self._filas(cuerpo, "Tabla")
        fuente = re.search(r"\\fuente\{([^{}]*)\}", cuerpo)
        if fuente:
            parrafo = self.doc.add_paragraph()
            self._plano(parrafo, "Fuente: " + fuente.group(1), 10, cursiva=True)

    def _filas(self, cuerpo: str, tipo: str | None) -> None:
        inicio = cuerpo.find(r"\toprule")
        if inicio < 0:
            return
        resto = cuerpo[inicio:]
        resto = re.sub(r"\\(toprule|midrule|bottomrule|endfirsthead|endhead|endfoot|caption\*?\{[^{}]*\}|label\{[^}]*\})", "\n", resto)
        filas = []
        for trozo in re.split(r"\\\\", resto):
            celdas = [self._celda(parte) for parte in trozo.split("&")]
            if any(celdas) and (not filas or celdas != filas[-1]):
                filas.append(celdas)
        if not filas:
            return
        ancho = max(len(fila) for fila in filas)
        tabla = self.doc.add_table(rows=1, cols=ancho)
        bordes(tabla)
        for indice, valor in enumerate(filas[0]):
            tabla.rows[0].cells[indice].text = ""
            parrafo = tabla.rows[0].cells[indice].paragraphs[0]
            run = parrafo.add_run(valor)
            rfonts(run, "Arial")
            run.bold = True
            run.font.size = Pt(10)
        encabezado = tabla.rows[0]._tr.get_or_add_trPr()
        marca = OxmlElement("w:tblHeader")
        encabezado.append(marca)
        for fila in filas[1:]:
            celdas = tabla.add_row().cells
            for indice in range(ancho):
                celdas[indice].text = ""
                parrafo = celdas[indice].paragraphs[0]
                run = parrafo.add_run(fila[indice] if indice < len(fila) else "")
                rfonts(run, "Arial")
                run.font.size = Pt(10)
        self.doc.add_paragraph()

    def _celda(self, texto: str) -> str:
        texto = re.sub(r"\\fuente\{[^{}]*\}", "", texto)
        texto = texto.replace(r"\%", "%").replace(r"\,", "").replace("{,}", ",").replace("$", "")
        texto = re.sub(r"\\[a-zA-Z]+\*?(?:\[[^\]]*\])?(?:\{[^{}]*\})?", "", texto)
        return " ".join(texto.replace("{", "").replace("}", "").split())

    def _figura(self, cuerpo: str) -> None:
        ruta = re.search(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", cuerpo)
        caption = re.search(r"\\caption\{([^{}]*)\}", cuerpo)
        self.figura += 1
        if ruta:
            imagen = LATEX / ruta.group(1)
            parrafo = self.doc.add_paragraph()
            parrafo.alignment = WD_ALIGN_PARAGRAPH.CENTER
            parrafo.add_run().add_picture(str(imagen), width=Cm(14.5))
        if caption:
            self._leyenda("Figura", self.figura, caption.group(1))
        fuente = re.search(r"\\fuente\{([^{}]*)\}", cuerpo)
        if fuente:
            parrafo = self.doc.add_paragraph()
            parrafo.alignment = WD_ALIGN_PARAGRAPH.CENTER
            self._plano(parrafo, "Fuente: " + fuente.group(1), 10, cursiva=True)

    def _leyenda(self, tipo: str, numero: int, texto: str) -> None:
        parrafo = self.doc.add_paragraph()
        parrafo.alignment = WD_ALIGN_PARAGRAPH.LEFT if tipo == "Tabla" else WD_ALIGN_PARAGRAPH.CENTER
        parrafo.paragraph_format.space_after = Pt(6)
        etiqueta = parrafo.add_run(f"{tipo} ")
        rfonts(etiqueta, "Arial")
        etiqueta.bold = True
        etiqueta.font.size = Pt(11)
        campo(parrafo, f" SEQ {tipo} \\* ARABIC ", str(numero))
        resto = parrafo.add_run(f": {texto}")
        rfonts(resto, "Arial")
        resto.font.size = Pt(11)

    def referencias(self) -> None:
        self.salto_impar(None, "decimal")
        self.doc.add_heading("Referencias bibliográficas", level=1)
        nota = self.doc.add_paragraph()
        self._plano(
            nota,
            "Las referencias siguen la 7.ª edición del Manual de publicaciones de la APA, en orden alfabético. "
            "La guía de la Biblioteca no impone un estilo. Cada DOI se confirmó en Crossref el 31 de octubre de 2026. "
            "Para volver a generarlas con Zotero, Mendeley o JabRef, siga word/INSTRUCCIONES-WORD.md.",
            11,
        )
        ordenadas = sorted(self.entradas.items(), key=lambda par: (apellidos(par[1].get("author", ""))[:1] or [""], par[1].get("year", "")))
        for _, campos in ordenadas:
            parrafo = self.doc.add_paragraph()
            parrafo.paragraph_format.left_indent = Cm(1.27)
            parrafo.paragraph_format.first_line_indent = Cm(-1.27)
            parrafo.paragraph_format.space_after = Pt(8)
            parrafo.alignment = WD_ALIGN_PARAGRAPH.LEFT
            self._plano(parrafo, referencia_apa(campos), 11)


def main() -> None:
    shutil.copy(PLANTILLA, SALIDA)
    doc = Document(SALIDA)
    cuerpo = doc.element.body
    for hijo in list(cuerpo):
        if hijo.tag != qn("w:sectPr"):
            cuerpo.remove(hijo)
    preparar_estilos(doc)
    activar_pares(doc)
    seccion = doc.sections[0]
    margenes(seccion)
    vaciar_header(seccion.header)
    vaciar_header(seccion.even_page_header)
    vaciar_header(seccion.first_page_header)
    vaciar_header(seccion.footer)
    entradas = leer_bib(BIB)
    constructor = Constructor(doc, entradas)
    constructor.portada()
    for relativo in ARCHIVOS:
        constructor.procesar(LATEX / relativo)
    constructor.referencias()
    propiedades = doc.core_properties
    propiedades.author = "Polanco Franco, Juan Diego"
    propiedades.title = "Qué parte del tumor colorrectal se deja ver en la sangre"
    propiedades.subject = "Tesis para optar al título de Magíster en Bioinformática y Biología Computacional"
    propiedades.keywords = "cáncer colorrectal; proteómica; suero; especificidad; Gene Ontology; neoantígeno"
    propiedades.category = "Trabajo de grado – Maestría"
    doc.save(SALIDA)
    print(SALIDA)


if __name__ == "__main__":
    main()
