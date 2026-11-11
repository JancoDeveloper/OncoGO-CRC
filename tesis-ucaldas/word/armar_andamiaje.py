"""Arma el .docx de andamiaje sobre los estilos del archivo ya basado en la guía."""
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.text import WD_LINE_SPACING
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

ORIGEN = "tesis-maestria-ucaldas.docx"


def limpiar(doc: Document) -> None:
    body = doc.element.body
    for child in list(body):
        if child.tag != qn("w:sectPr"):
            body.remove(child)


def todo(doc: Document, texto: str) -> None:
    p = doc.add_paragraph()
    run = p.add_run("TODO Juan Diego: " + texto)
    run.italic = True
    run.font.color.rgb = RGBColor(0x80, 0x80, 0x80)
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE


def cuerpo(doc: Document, texto: str) -> None:
    p = doc.add_paragraph(texto)
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE


def viñeta(doc: Document, texto: str) -> None:
    p = doc.add_paragraph(texto, style="List Bullet")
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE


def titulo(doc: Document, texto: str, nivel: int, salto: bool) -> None:
    if salto:
        doc.add_section(WD_SECTION.ODD_PAGE)
    doc.add_heading(texto, level=nivel)


def margenes(section) -> None:
    section.page_width = Cm(21.59)
    section.page_height = Cm(27.94)
    section.left_margin = Cm(3.60)
    section.right_margin = Cm(2.54)
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)


def main() -> None:
    doc = Document(ORIGEN)
    limpiar(doc)
    for section in doc.sections:
        margenes(section)

    cuerpo(doc, "Qué parte del tumor colorrectal se deja ver en la sangre")
    cuerpo(doc, "Juan Diego Polanco Franco")
    cuerpo(doc, "Universidad de Caldas")
    cuerpo(doc, "Facultad de Inteligencia Artificial e Ingenierías")
    cuerpo(doc, "Manizales, Colombia")
    cuerpo(doc, "2026")
    titulo(doc, "Portada", 1, True)
    cuerpo(doc, "Qué parte del tumor colorrectal se deja ver en la sangre")
    cuerpo(doc, "Juan Diego Polanco Franco")
    cuerpo(doc, "Tesis o trabajo de investigación presentada(o) como requisito parcial para optar al título de:")
    cuerpo(doc, "Magíster en Bioinformática y Biología Computacional")
    cuerpo(doc, "Programa de investigación, 60 créditos, SNIES 103917. Convenio SUMA con la Universidad de Manizales, la Universidad Autónoma de Manizales y la Universidad Católica de Manizales.")
    cuerpo(doc, "Director: Doctor en Ingeniería, Jorge Alberto Jaramillo Garzón")
    cuerpo(doc, "Grupo de investigación: Grupo de Investigación en Tecnologías de la Información y Redes (GITIR)")
    cuerpo(doc, "Línea de investigación: Modelos Computacionales y Bioinformática")
    cuerpo(doc, "Universidad de Caldas. Facultad de Inteligencia Artificial e Ingenierías. Manizales, Colombia. 2026.")
    cuerpo(doc, "Andamiaje, 11 de noviembre de 2026. No es la prosa final.")

    bloques = [
        ("Declaración de uso de herramientas de inteligencia artificial",
         "Cursor 3.28.7 y SciSpace, consultado entre 2026-03 y 2026-11. No las cites como fuente de un hecho científico. Puedes reescribir el alcance."),
        ("Resumen", "Redacta como máximo 250 palabras. Usa 0 diferentes, 166 equivalentes y 118 indeterminadas, de 284; 147 de 147; AUC 0,880, 0,915 y 0,647 (p 0,12). No copies las viñetas."),
        ("Abstract", "Traduce tu resumen al inglés. No añadas un resultado."),
        ("Tabla de contenido", "Clic derecho, Actualizar campos, toda la tabla."),
        ("Lista de tablas", "Clic derecho, Actualizar campos, toda la tabla."),
        ("Lista de figuras", "Clic derecho, Actualizar campos, toda la tabla."),
        ("Lista de abreviaturas", "Revisa CN, CP, HC, TOST y TREAT. HC no es el contraste que decide la pregunta. CN no son personas sanas."),
        ("1. Introducción", "De 2 a 4 páginas: problema, justificación, pregunta y objetivos. No prometas un examen."),
        ("2. Marco teórico y estado del arte", "Define solo lo que luego lee un resultado. Equivalente no es identidad biológica."),
        ("3. Materiales y métodos", "Reescribe con tus palabras los parámetros del preregistro. No cambies alfa 0,05 ni el margen 0,585. El párrafo de Cursor y SciSpace está en la guía."),
        ("4. Resultados", "Escribe la interpretación de las figuras f01 a f09. Los captions técnicos están en el LaTeX. El cero no lleva un segmento inventado."),
        ("5. Discusión", "Incluye las limitaciones: trabajo in silico, sin examen clínico. Las desviaciones van en esta sección, no en un capítulo 7."),
        ("6. Conclusiones y recomendaciones", "Una conclusión por objetivo, con la cifra y con lo que no demuestra. No anticipes un número futuro."),
        ("Anexos", "Elige las tablas suplementarias. La palabra cumplido del TSV de correspondencia no es un biomarcador."),
        ("Referencias bibliográficas", "APA 7. Actualiza los campos si usas el gestor. La biblioteca verificada está en latex/references.bib."),
    ]
    for nombre, nota in bloques:
        titulo(doc, nombre, 1, True)
        todo(doc, nota)
        cuerpo(doc, "El detalle, las cifras y las preguntas guía están en GUIA-DE-REDACCION.md y en el PDF compilado desde latex/main.tex.")

    for section in doc.sections:
        margenes(section)
    doc.core_properties.author = "Polanco Franco, Juan Diego"
    doc.core_properties.title = "Qué parte del tumor colorrectal se deja ver en la sangre"
    doc.core_properties.subject = "Tesis para optar al título de Magíster en Bioinformática y Biología Computacional"
    doc.core_properties.keywords = "cáncer colorrectal; proteómica; suero; especificidad; Gene Ontology"
    doc.core_properties.category = "Trabajo de grado – Maestría"
    doc.save(ORIGEN)
    print("secciones", len(doc.sections))


if __name__ == "__main__":
    main()
