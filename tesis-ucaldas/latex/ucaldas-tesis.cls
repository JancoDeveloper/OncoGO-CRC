% Clase para imitar la guía de presentación de tesis de la Biblioteca
% de la Universidad de Caldas (2022). No existe plantilla LaTeX oficial.
% Compilar con LuaLaTeX + Biber.
%
% Si un capítulo se publicara como artículo, la carta de autorización del
% repositorio puede pedir embargo. Eso no es un requisito de esta clase.
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ucaldas-tesis}[2026/11/06 Tesis de maestría, Universidad de Caldas]
\LoadClass[11pt,letterpaper,twoside,openright]{book}

\RequirePackage{fontspec}
% TeX Gyre Heros es la familia libre de métrica afín a Arial, disponible en Overleaf.
% Arial propietario no está en TeX Live. Desviación documentada en VALIDACION-FORMATO.md.
\setmainfont{TeX Gyre Heros}
\setsansfont{TeX Gyre Heros}

\RequirePackage[english,spanish]{babel}
\RequirePackage{csquotes}
\RequirePackage{geometry}
\geometry{
  paperwidth=21.59cm,
  paperheight=27.94cm,
  inner=3.60cm,
  outer=2.54cm,
  top=2.54cm,
  bottom=2.54cm,
  headheight=14pt,
  headsep=0.6cm
}
\RequirePackage{setspace}
\onehalfspacing
\setlength{\parindent}{0pt}
% Equivale al Enter adicional de la guía, sin insertar párrafos vacíos.
\setlength{\parskip}{0.9\baselineskip}

\RequirePackage{titlesec}
\titleformat{\chapter}[hang]
  {\normalfont\fontsize{20}{24}\selectfont\bfseries\singlespacing}
  {\thechapter.}{0.6em}{}
\titleformat{name=\chapter,numberless}[hang]
  {\normalfont\fontsize{20}{24}\selectfont\bfseries\singlespacing}
  {}{0pt}{}
\titlespacing*{\chapter}{0pt}{100pt}{24pt}

% El ejemplo de la guía: subtítulo 16 pt sin negrita; el de 14 pt va en negrita.
\titleformat{\section}
  {\normalfont\fontsize{16}{19}\selectfont\singlespacing}
  {\thesection}{0.6em}{}
\titlespacing*{\section}{0pt}{20pt}{24pt}
\titleformat{\subsection}
  {\normalfont\fontsize{14}{17}\selectfont\bfseries\singlespacing}
  {\thesubsection}{0.6em}{}
\titlespacing*{\subsection}{0pt}{12pt}{10pt}

\RequirePackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[LE]{\fontsize{11}{13}\selectfont\thepage}
\fancyhead[RO]{\fontsize{11}{13}\selectfont\thepage}
\fancyhead[C]{\fontsize{8}{10}\selectfont\textit{Andamiaje, 11 de noviembre de 2026}}
\renewcommand{\headrulewidth}{0pt}
\fancypagestyle{plain}{%
  \fancyhf{}%
  \fancyhead[LE]{\fontsize{11}{13}\selectfont\thepage}%
  \fancyhead[RO]{\fontsize{11}{13}\selectfont\thepage}%
  \fancyhead[C]{\fontsize{8}{10}\selectfont\textit{Andamiaje, 11 de noviembre de 2026}}%
  \renewcommand{\headrulewidth}{0pt}%
}

\RequirePackage{graphicx}
\RequirePackage{longtable}
\RequirePackage{booktabs}
\RequirePackage{array}
\RequirePackage{caption}
\captionsetup{
  font={normalsize,sf},
  labelfont={bf},
  labelsep=colon,
  justification=justified,
  singlelinecheck=false
}
\captionsetup[table]{position=top}
\captionsetup[figure]{position=bottom}
\renewcommand{\listfigurename}{Lista de figuras}
\renewcommand{\listtablename}{Lista de tablas}
\addto\captionsspanish{\renewcommand{\contentsname}{Tabla de contenido}}

\RequirePackage{amsmath}
\RequirePackage{xcolor}
\newcommand{\placeholder}[1]{\textcolor{gray}{\textit{[#1]}}}
\newcommand{\todojd}[1]{\par\noindent\textcolor{gray}{\textit{\% TODO Juan Diego: #1}}\par}
\newcommand{\objetivo}[1]{\par\noindent\textbf{Objetivo de la sección.} #1\par}

\RequirePackage[
  backend=biber,
  style=apa,
  sortcites=true,
  doi=true,
  isbn=false,
  url=true
]{biblatex}
% Para Vancouver, cambie style=apa por style=vancouver y vuelva a compilar con Biber.
% La guía de 2022 no impone un estilo: APA 7 es el valor por defecto de este proyecto.
\DeclareLanguageMapping{spanish}{spanish-apa}

\RequirePackage[
  colorlinks=true,
  linkcolor=black,
  citecolor=black,
  urlcolor=blue,
  pdftitle={Qué parte del tumor colorrectal se deja ver en la sangre},
  pdfauthor={Polanco Franco, Juan Diego},
  pdfsubject={Tesis para optar al título de Magíster en Bioinformática y Biología Computacional},
  pdfkeywords={cáncer colorrectal; proteómica; suero; especificidad; Gene Ontology; neoantígeno}
]{hyperref}

\newcommand{\fuente}[1]{\par{\fontsize{10}{12}\selectfont Fuente: #1\par}}

\AtBeginDocument{\pagenumbering{Roman}}
