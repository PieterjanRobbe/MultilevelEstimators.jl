## tex.jl : tex utilities for making reports

const LINESTYLES = ["dashed","dotted","dashdotted"]
const PGFPLOTSMARKERS = ["*","square*","triangle*","diamond*","pentagon*","x"]

tikz_header(x_label,y_label,optional,name) = "%!TEX root = ../$(name)\n
\\begin{tikzpicture}[trim axis left,trim axis right]
\\begin{axis}[
width=\\figurewidth,
height=\\figureheight,
scale only axis,
xlabel={\\small $(x_label)},
every x tick label/.append style={font=\\scriptsize},
every x label/.append style={font=\\scriptsize},
ylabel={\\small $(y_label)},
every y tick label/.append style={font=\\scriptsize},
every y label/.append style={font=\\scriptsize},
legend style={draw=none,font=\\scriptsize,at={(0.03,0.03)},anchor=south west,fill=none},
$(isempty(optional) ? "" : optional)"*"]\n
"

tikz_footer() = "
\\end{axis}
\\end{tikzpicture}
"

tikz_add_plot(color_name,line_style,marker_type,name,legend_name,table_options) = "
\\addplot [color=$(color_name),$(line_style),mark=$(marker_type),mark options={solid,fill=$(color_name)},line width=0.75pt,mark size=1.2,line cap=round$(isempty(legend_name) ? ", forget plot" : "")]
table[$(table_options)]{data/$(name).txt};
$(isempty(legend_name) ? "" : "\\addlegendentry{$(legend_name)};")
"

figure_header() = "
\\begin{figure}[h]
\\centering
\\setlength{\\figureheight}{0.3\\textwidth}
\\setlength{\\figurewidth}{0.3\\textwidth}
"

figure_footer() = "
\\end{figure}
"

file_contents(title) = "
\\documentclass[11pt, oneside]{article}

\\usepackage[margin=2cm]{geometry}
\\usepackage{graphicx}
\\usepackage{pgfplots}
\\newlength\\figureheight
\\newlength\\figurewidth
\\usepackage{tikz}

\\title{Report: $(title)}
\\date{}

\\begin{document}
\\maketitle\n"*
figure_header()*
"\\input{figures/E.tex}\n \\caption{\\label{fig:E}Decay of the expected value \$E[\\Delta Q]\$.}"*
figure_footer()*"\n"*
figure_header()*
"\\input{figures/V.tex}\n \\caption{\\label{fig:V}Decay of the variance \$V[\\Delta Q]\$.}"*
figure_footer()*"\n"*
#figure_header()*
#"\\input{figures/samples.tex}\n \\caption{\\label{fig:samples}Total number of samples \$N_\\ell\$ taken on each level.}"*
#figure_footer()*"\n"*
figure_header()*
"\\input{figures/time.tex}\n \\caption{\\label{fig:time}Total simulation run time.}"*
figure_footer()*"\n"*
figure_header()*
"\\input{figures/cost.tex}\n \\caption{\\label{fig:cost}Total simulation standard cost.}"*
figure_footer()*
"\\end{document}"
