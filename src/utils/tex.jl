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
legend style={draw=none,font=\\scriptsize,at={(0.03,0.03)},anchor=south west,fill=none, legend cell align={left}},
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

preamble(d) = "$(d == 2 ? square() : d == 3 ? cube() : "\n\n")"

square() = "\n\\usepackage{booktabs}\n\n\\newcommand{\\drawsquare}[3]{
\\edef\\temp{
\\noexpand\\addplot[line width=3pt,white,fill=#3,forget plot,shift={(#1,#2)}]
}
\\temp
table[row sep=crcr] {%
x	y\\\\
0	0\\\\
1	0\\\\
1	1\\\\
0	1\\\\
}--cycle;
}\n"

cube() = "\n\\usepackage{booktabs}\n\n\\newcommand{\\drawcube}[4]{
{	
\\edef\\temp{
\\noexpand\\addplot3[area legend,solid,fill=#4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
}
\\temp
table[row sep=crcr] {%
x	y	z\\\\
1	0	0\\\\
1	0	1\\\\
1	1	1\\\\
1	1	0\\\\
}--cycle;

\\edef\\temp{
\\noexpand\\addplot3[area legend,solid,fill=#4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
}
\\temp
table[row sep=crcr] {%
x	y	z\\\\
0	1	0\\\\
1	1	0\\\\
1	1	1\\\\
0	1	1\\\\
}--cycle;

\\edef\\temp{
\\noexpand\\addplot3[area legend,solid,fill=#4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
}
\temp
table[row sep=crcr] {%
x	y	z\\\\
0	0	1\\\\
1	0	1\\\\
1	1	1\\\\
0	1	1\\\\
}--cycle;
}\n"

file_contents(title,d) = "\\documentclass[11pt, oneside]{article}

\\usepackage[margin=2cm]{geometry}
\\usepackage{graphicx}
\\usepackage{pgfplots}
    \\pgfplotsset{compat=newest}
    \\newlength\\figureheight
    \\newlength\\figurewidth
\\usepackage{tikz}"*
preamble(d)*"
\\title{Report: $(title)}
\\date{}

\\begin{document}
\\maketitle\n"*
figure_header()*
"\\input{figures/E.tex}\n \\caption{\\label{fig:E}Decay of the expected value.}"*
figure_footer()*"\n"*
figure_header()*
"\\input{figures/V.tex}\n \\caption{\\label{fig:V}Decay of the variance.}"*
figure_footer()*"\n"*
"$(d == 1 ? figure_header()*"\\input{figures/samples.tex}\n \\caption{\\label{fig:samples}Total number of samples \$N_\\ell\$ taken on each level.}"*figure_footer()*"\n" : d < 4 ? "\\input{figures/index_set.tex}\n" : "" )"*"\n"*
figure_header()*
"\\input{figures/runtime.tex}\n \\caption{\\label{fig:time}Total simulation run time.}"*
figure_footer()*"\n"*
figure_header()*
"\\input{figures/cost.tex}\n \\caption{\\label{fig:cost}Total simulation standard cost.}"*
figure_footer()*
"\\end{document}"

tikz_index_set_2d(name,i,max_level) = "%!TEX root = ../$(name)\n
\\pgfplotstableread[header=false]{data/index_set_$(i).txt}\\indexset
\\pgfplotstablegetrowsof{\\indexset}
\\pgfmathsetmacro{\\rows}{\\pgfplotsretval-1}\n
\\begin{tikzpicture}[trim axis left,trim axis right]
\\begin{axis}[
width=\\figurewidth,
height=\\figureheight,
scale only axis,
xmin=-0.1,
xmax=$(max_level+1).1,
xticklabel={},
xmajorticks=false,
ymin=-0.1,
ymax=$(max_level+1).1,
yticklabel={},
ymajorticks=false,
axis line style={ultra thin, draw opacity=0}
]\n
\\foreach \\j in {0,...,\\rows} {
    \\pgfplotstablegetelem{\\j}{0}\\of\\indexset
        \\pgfmathsetmacro{\\a}{\\pgfplotsretval}
        \\pgfplotstablegetelem{\\j}{1}\\of\\indexset
        \\pgfmathsetmacro{\\b}{\\pgfplotsretval}
	\\drawsquare{\\a}{\\b}{white!90!black}
}
\\end{axis}\n
\\end{tikzpicture}
"

tikz_index_set_3d(name,i,max_level) = "%!TEX root = ../$(name)\n
\\pgfplotstableread[header=false]{data/index_set_$(i).txt}\\indexset
\\pgfplotstablegetrowsof{\\indexset}
\\pgfmathsetmacro{\\rows}{\\pgfplotsretval-1}\n
\\begin{tikzpicture}[trim axis left,trim axis right]
\\begin{axis}[
width=\\figurewidth,
height=\\figureheight,
scale only axis,
view={120}{17},
xmin=-0.1,
xmax=$(max_level+1).1,
xticklabel={},
xmajorticks=false,
ymin=-0.1,
ymax=$(max_level+1).1,
yticklabel={},
ymajorticks=false,
zmin=-0.1,
zmax=$(max_level+1).1,
zticklabel={},
zmajorticks=false,
axis line style={ultra thin, draw opacity=0}
]\n
\\foreach \\j in {0,1,...,\\rows} {
	\\pgfplotstablegetelem{\\j}{2}\\of\\indexset
   \\pgfmathsetmacro{\\a}{\\pgfplotsretval}
	\\pgfplotstablegetelem{\\j}{0}\\of\\indexset
   \\pgfmathsetmacro{\\b}{\\pgfplotsretval}
	\\pgfplotstablegetelem{\\j}{1}\\of\\indexset
   \\pgfmathsetmacro{\\c}{\\pgfplotsretval}
   \\drawcube{\\a}{\\b}{\\c}{white!90!black}
}
\\end{axis}\n
\\end{tikzpicture}
"

function tikz_index_set_table(name,tols)
    str = "%!TEX root = ../$(name)\n"#\n\\begin{table}[t]\n\\centering\n"
    for tabnum = 1:div(length(tols)-1,9)+1
        str = string(str,figure_header())
        str = string(str,"\\begin{tabular}{ccc}")
        # TODO newpage ?
        for col in 1:min(3,div(length(tols)-(tabnum-1)*9-1,3)+1)
            # header
            str = string(str," \\toprule\n\$\\epsilon = \\pgfmathprintnumber[/pgf/number format/sci,precision=3,fixed zerofill]\{", @sprintf("%7.3e",tols[(tabnum-1)*9+(col-1)*3+1]), "\}\$ \&")
            str = length(tols) > (tabnum-1)*9+(col-1)*3+1 ? string(str," \$\\epsilon = \\pgfmathprintnumber[/pgf/number format/sci,precision=3,fixed zerofill]\{", @sprintf("%7.3e",tols[(tabnum-1)*9+(col-1)*3+2]), "\}\$") : str
            str = string(str," \&")
            str = length(tols) > (tabnum-1)*9+(col-1)*3+2 ? string(str," \$\\epsilon = \\pgfmathprintnumber[/pgf/number format/sci,precision=3,fixed zerofill]\{", @sprintf("%7.3e",tols[(tabnum-1)*9+(col-1)*3+3]), "\}\$") : str
            str = string(str,"\\\\ \\midrule \n")
            str = string(str," \\input\{figures/index_set_",(tabnum-1)*9+(col-1)*3+1, ".tex\} \&")
            str = length(tols) > (tabnum-1)*9+(col-1)*3+1 ? string(str," \\input\{figures/index_set_",(tabnum-1)*9+(col-1)*3+2, ".tex\}") : str
            str = string(str," \&")
            str = length(tols) > (tabnum-1)*9+(col-1)*3+2 ? string(str," \\input\{figures/index_set_",(tabnum-1)*9+(col-1)*3+3, ".tex\}") : str
            str = string(str,"\\\\ \n")
        end
        str = string(str,"\n\\end{tabular}\\\\\n\n")
        str = string(str,"\\caption{\\label{fig:index_set}Shape of the index set for different tolerances on the RMSE ($(tabnum)/$(div(length(tols)-1,9)+1)).}")
        str = string(str,figure_footer())
    end
    str
end
