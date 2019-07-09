\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ICASP13Paper}[2018/08/07 ICASP13Paper]
% Written by:
% Joseph Nagel, Stefano Marelli
% ETH Zurich - Chair of Risk, Safety and Uncertainty Quantification
% September 2014

% Modified by:
% Junho Song
% SNU Seoul - Chair of Organizing Committee
% August 2018

%% Article options
\DeclareOption{10pt}{
	\PassOptionsToClass{\CurrentOption}{article}
}

%% Fallback
\DeclareOption*{
	\ClassWarning{myclass}{Unknown option '\CurrentOption'}
}


%% Custom package options
\DeclareOption{sansserif}{
	\PassOptionsToPackage{\CurrentOption}{paxcommands}
}
\DeclareOption{neverindent}{
	\PassOptionsToPackage{\CurrentOption}{paxcommands}
}


%% Execute default options
\ExecuteOptions{10pt}

%% Process given options
\ProcessOptions\relax

%% Load base class (paper in a4 format)
\LoadClass[a4paper,12pt,twocolumn]{article}

% Page layout
\RequirePackage[left=1.7cm,right=1.7cm,top=2.95cm,bottom=2.95cm]{geometry}
\RequirePackage{graphicx}

% Amsmath package
\RequirePackage{amsmath}
\RequirePackage{amstext}
\RequirePackage{amsfonts}
\RequirePackage{amssymb}

% Text encoding
\RequirePackage[utf8x]{inputenc}
\RequirePackage[T1]{fontenc}

% Fonts
\RequirePackage{mathptmx}% Times Roman font

% Bibliography
\RequirePackage{natbib}

% Caption formatting
\RequirePackage[font={small,it},justification=RaggedRight,singlelinecheck=false]{caption}

% Paralist
\RequirePackage{paralist}

% hyperlinks
\RequirePackage[colorlinks=true, urlcolor=black, linkcolor=black, 
citecolor=black]{hyperref}


%% Headers/footers: fancyhdr config
\RequirePackage{fancyhdr}
%\pagestyle{empty} 

\pagestyle{fancy} \fancyhead{}
\fancyfoot[C]{\thepage} \fancyhead[R]{ \footnotesize \textit{
  13th International Conference on Applications of Statistics and Probability 
  in Civil Engineering, ICASP13 \\
  Seoul, South Korea, May 26-30, 2019} }
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}

% Titlesec formatting
\RequirePackage{titlesec}
\titlelabel{\thetitle.\quad}
\titleformat*{\section}{\normalfont\scshape}
\titleformat*{\subsection}{\normalfont\itshape}
\titleformat*{\subsubsection}{\normalfont\itshape}
\titlespacing*{\section}{0pt}{12pt plus 4pt minus 2pt}{0pt plus 2pt minus 2pt}
\titlespacing*{\subsection}{0pt}{12pt plus 4pt minus 2pt}{0pt plus 2pt minus 2pt}
\titlespacing*{\subsubsection}{0pt}{12pt plus 4pt minus 2pt}{0pt plus 2pt minus 
2pt}
\renewcommand{\bibsection}{}


%% HEADER MACROS SECTION: REQUIRED TO COMPLETE THE TEMPLATE
% bibtex file to be used to create the references
\newcommand{\referencefile}[1]{\newcommand{\bibfile}{#1}}

% short abstract
\newcommand{\shortabstract}[1]{\newcommand{\shortabs}{#1}}

% Expanding author list with the \addauthor command
\def \authlist {}
\newcommand{\addauthor}[2]{
\expandafter\def\expandafter\authlist\expandafter{\authlist{}
    \fontsize{14pt}{16pt}\selectfont #1 \\
    \normalfont \textit{#2} \\ \vspace{1ex} }
}

% Patch to allow for text spanning on two pages columns at the beginning of the 
%document
\RequirePackage{etoolbox}
\patchcmd[\long]%   prefix
  {\@topnewpage}%   cmd
  {\@nodocument}{}% search & replace
  {}{}%             success & failure
  


%% Title page (header, title and authors)
\AtBeginDocument{
\twocolumn[{
\begin{@twocolumnfalse}
  \begin{flushleft}
    \fontsize{18pt}{22pt}\selectfont
    \@title
  \end{flushleft}
  \vspace{1ex}
  \begin{flushleft}
  	% use the \authlist expanded by the \addauthor directives
	\authlist
  \end{flushleft}
  \vspace{0ex}
  \begin{flushleft}
    ABSTRACT: \shortabs
  \end{flushleft}
  \vspace{3ex}
\end{@twocolumnfalse}
}]
\noindent
}


%% Add the reference section at the end of the document
%\AtEndDocument{
%
%}
\endinput
