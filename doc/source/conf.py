# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Rayleigh'
copyright = '2019, Nick Featherstone'
author = 'Nick Featherstone'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
]

master_doc = 'index'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
'papersize': 'letterpaper',
'babel': '\usepackage[english]{babel}',

# The font size ('10pt', '11pt' or '12pt').
'pointsize': '11pt',

# Additional stuff for the LaTeX preamble.
'preamble': r'''
\usepackage{graphicx}
\usepackage{amsmath}
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt}
\usepackage{textpos}

% use a larger page size; otherwise, it is difficult to have complete
% code listings and output on a single page
\usepackage{fullpage}


% use the listings package for code snippets. define keywords for prm files
% and for gnuplot
\usepackage{listings}
\lstset{frame=tb,
  language=Fortran,
  aboveskip=3mm,
  belowskip=3mm,
  showstringspaces=false,
  columns=flexible,
  basicstyle={\small\ttfamily},
  numbers=none,
  numberstyle=\tiny\color{gray},
  keywordstyle=\color{blue},
  commentstyle=\color{dkgreen},
  stringstyle=\color{mauve},
  breaklines=true,
  breakatwhitespace=true,
  tabsize=3
}

\definecolor{dkgreen}{rgb}{0,0.6,0}
\definecolor{gray}{rgb}{0.5,0.5,0.5}
\definecolor{mauve}{rgb}{0.58,0,0.82}

\newcommand{\rayleigh}{\textsc{Rayleigh}}
  ''',

'maketitle': r'''
\definecolor{dark_grey}{gray}{0.3}

%LINE 1%
{
\renewcommand{\familydefault}{\sfdefault}

\pagenumbering{gobble}
\begin{center}
\resizebox{\textwidth}{!}{\textcolor{dark_grey}{\fontfamily{\sfdefault}\selectfont
COMPUTATIONAL INFRASTRUCTURE FOR GEODYNAMICS (CIG)
}}

%\hrule

%LINE 2%
\color{dark_grey}
\rule{\textwidth}{2pt}

%LINE 3%
%\color{dark_grey}
% FILL: additional organizations
% e.g.: {\Large Organization 1\\Organization 2}
%{\Large }
\end{center}

%COLOR AND CODENAME BLOCK%
\begin{center}
\resizebox{\textwidth}{!}{\colorbox{orange}{\fontfamily{\rmdefault}\selectfont \textcolor{white} {\Large
\hspace{0.2in}\rayleigh{}\hspace{0.1in}} }}
\end{center}

%MAIN PICTURE%
\begin{textblock*}{0in}(0.0in,0.0in)
% FILL: image height
% e.g. height=6.5in
\begin{center}
\vspace{3em}
\includegraphics[height=4.0in]
% FILL: image file name of your software
% e.g. cover_image.png
{../../src/user_guide/rayleigh_manual_image_300dpi.jpeg}
\hspace{1em}
\end{center}
\end{textblock*}

%USER MANUAL%
\color{dark_grey}
\vspace{0.5in}
\hfill{\Huge \fontfamily{\sfdefault}\selectfont User Manual \\
\raggedleft \huge \fontfamily{\sfdefault}\selectfont Version
% keep the following line as is so that we can replace this using a script:
0.9.1 %VERSION-INFO%
\\\large(generated \today)\\
{\Large Nicholas Featherstone\\}
}
%AUTHOR(S) & WEBSITE%
\null
\vspace{17em}
\color{dark_grey}
{\fontfamily{\sfdefault}\selectfont
\large
\vspace{0.58in}
\noindent with contributions by: \\
    Kyle Augustson, Wolfgang Bangerth, Rene Gassm\"oller, Sebastian Glane, Brad Hindman, Lorraine Hwang, Hiro Matsui, Ryan Orvedahl, Krista Soderlund, Cian Wilson, Maria Weber, Rakesh Yadav
    \\
\vspace{1.0em}

{\noindent
{\href{https://geodynamics.org}{geodynamics.org}}
}
}


%LINE%
{\noindent
\color{dark_grey}
\rule{\textwidth}{2pt}
}

%COPYRIGHT STATEMENT
\textcopyright Copyright 2018, Regents of the University of California

}

\pagebreak
'''
}
