## tex.jl : utilities for making ".tex" files
#
# Utility functions for creating ".tex" files (for report.jl).
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## Plot ##
abstract type Plot end

## LinePlot ##
struct LinePlot <: Plot
    color::String
    line_style::String
    mark::String
    legend::String
    data_file::String
    table_options::String
end

function generate(plot::LinePlot)
    str = "\\addplot ["
    str = string(str, "color=", plot.color, ", ")
    str = string(str, plot.line_style, ", ")
    str = string(str, "mark=", plot.mark, ", ")
    str = string(str, "mark options={solid, fill=", plot.color, "}, ")
    str = string(str, "line width=0.75pt, ")
    str = string(str, "mark size=1.2, ")
    str = string(str, "line cap=round, ")
    str = string(str, isempty(plot.legend) ? "forget plot, " : "", "]")
    str = string(str, "table[", plot.table_options, "]")
    str = string(str, "{data/", plot.data_file, ".txt};\n")
    if !isempty(plot.legend)
        str = string(str, "\\addlegendentry{", plot.legend, "};")
    end
    return str
end

## TikzPicture ##
struct TikzPicture
    x_label
    y_label
    x_min
    x_max
    x_tick
    x_mode
    y_mode
    legend_at
    legend_anchor
    plots::Vector{Plot}
end

function generate(pic::TikzPicture)
    str = String[]
    push!(str, shebang(pic))
    push!(str, begintikzpicture(pic))
    push!(str, beginaxis(pic))
    push!(str, width(pic))
    push!(str, height(pic))
    push!(str, scale_only_axis(pic))
    push!(str, x_min(pic))
    push!(str, x_max(pic))
    push!(str, x_tick(pic))
    push!(str, x_mode(pic))
    push!(str, x_label(pic))
    push!(str, x_axis_style(pic))
    push!(str, y_label(pic))
    push!(str, y_mode(pic))
    push!(str, y_axis_style(pic))
    push!(str, legend_style(pic))
    push!(str, endaxisoptions(pic))
    push!(str, generate.(pic.plots)...)
    push!(str, endaxis(pic))
    push!(str, endtikzpicture(pic))
    str = str[map(!isempty, str)]
    join(str, "\n")
end

shebang(pic::TikzPicture) = string("%!TEX root = ../report.tex")
begintikzpicture(pic::TikzPicture) = string("\\begin{tikzpicture}[trim axis left,trim axis right]")
endtikzpicture(pic::TikzPicture) = "\\end{tikzpicture}"
beginaxis(pic::TikzPicture) = "\\begin{axis}["
endaxisoptions(pic::TikzPicture) = "]"
endaxis(pic::TikzPicture) = "\\end{axis}"
width(pic::TikzPicture) = string("width=\\figurewidth,")
height(pic::TikzPicture) = "height=\\figureheight,"
scale_only_axis(pic::TikzPicture) = "scale only axis,"
x_label(pic::TikzPicture) = string("xlabel={\\small ", pic.x_label, "},")
x_min(pic::TikzPicture) = isempty(pic.x_min) ? "" : string("xmin=", pic.x_min, ",")
x_max(pic::TikzPicture) = isempty(pic.x_max) ? "" : string("xmax=", pic.x_max, ",")
x_tick(pic::TikzPicture) = isempty(pic.x_tick) ? "" : string("xtick={", pic.x_tick, "},")
y_label(pic::TikzPicture) = string("ylabel={\\small ", pic.y_label, "},")
y_mode(pic::TikzPicture) = isempty(pic.y_mode) ? "" : string("ymode=", pic.y_mode, ",")
x_mode(pic::TikzPicture) = isempty(pic.x_mode) ? "" : string("xmode=", pic.x_mode, ",")
axis_style(pic::TikzPicture, dir::String) = string("every ", dir, " tick label/.append style={font=\\scriptsize},\nevery ", dir, " label/.append style={font=\\scriptsize},")
x_axis_style(pic::TikzPicture) = axis_style(pic, "x")
y_axis_style(pic::TikzPicture) = axis_style(pic, "y")
additional_options(pic::TikzPicture) = pic.additional_options
legend_at(pic::TikzPicture) = isempty(pic.legend_at) ? "(0.03,0.03)" : pic.legend_at
legend_anchor(pic::TikzPicture) = isempty(pic.legend_anchor) ? "south west" : pic.legend_anchor
legend_style(pic::TikzPicture) = string("legend style={draw=none, font=\\scriptsize, at={", legend_at(pic), "}, anchor=", legend_anchor(pic), ", fill=none, legend cell align=left}")

## Figure ##
struct Figure
    name::String
    label::String        
    caption::String
end

begin_figure(fig::Figure) = "\\begin{figure}[h]\n\\centering\n\\setlength{\\figureheight}{0.3\\textwidth}\n\\setlength{\\figurewidth}{0.3\\textwidth}"
end_figure(fig::Figure) = "\\end{figure}"
caption(fig::Figure) = string("\\caption{\\label{fig:", fig.label, "}", fig.caption, "}")
input(fig::Figure) = string("\\input{figures/", fig.name, ".tex}")

function generate(fig::Figure)
    str = String[]
    push!(str, begin_figure(fig))
    push!(str, input(fig))
    push!(str, caption(fig))
    push!(str, end_figure(fig))
    join(str, "\n")
end 

## Document ##
struct Document
    packages::Vector{String}
    preamble::Vector{String}
    title::String
    contents::Vector{Figure}
end

function Document(h::History, folder::String)
    
    packages = String[]
    push!(packages, "\\usepackage[margin=2cm]{geometry}")
    push!(packages, "\\usepackage{amsfonts}")
    push!(packages, "\\usepackage{graphicx}")
    push!(packages, "\\usepackage{pgfplots}")
    push!(packages, "\\pgfplotsset{compat=newest}")
    push!(packages, "\\newlength\\figureheight")
    push!(packages, "\\newlength\\figurewidth")
    push!(packages, "\\usepackage{tikz}")

    preamble = String[]

    contents = Figure[]
    write_data_rates(h, folder)
    write_figures_rates(h, folder)
    push!(contents, Figure("E", "E", "Decay of the expected value."))
    push!(contents, Figure("V", "V", "Decay of the variance."))
    push!(contents, Figure("W", "W", "Increase of the cost."))

    d = h[:ndims]
    if d == 1
        write_data_samples(h, folder)
        write_figure_samples(h, folder)
        push!(contents, Figure("nb_of_samples", "N", "Total number of samples."))
    end

    write_data_runtime(h, folder)
    write_figure_runtime(h, folder)
    push!(contents, Figure("runtime", "runtime", "Runtime versus accuracy \$\\epsilon\$."))
    write_data_cost(h, folder)
    write_figure_cost(h, folder)
    push!(contents, Figure("cost", "cost", "Standard cost versus accuracy \$\\epsilon\$."))

    Document(packages, preamble, string("Report: ", h[:name]), contents)
end

function generate(doc::Document)
    str = String[]
    push!(str, document_header(doc))
    push!(str, join(doc.packages, "\n"))
    push!(str, join(doc.preamble, "\n"))
    push!(str, title(doc))
    push!(str, date(doc))
    push!(str, begin_document(doc))
    push!(str, maketitle(doc))
    push!(str, generate.(doc.contents)...)
    push!(str, end_document(doc))
    join(str, "\n")
end

document_header(doc::Document) = "\\documentclass[11pt, oneside]{article}"
begin_document(doc::Document) = "\\begin{document}"
end_document(doc::Document) = "\\end{document}"
title(doc::Document) = string("\\title{", doc.title[1:end-4], "}")
date(doc::Document) = "\\date{}"
maketitle(doc::Document) = "\\maketitle"

## Colors, Linestyles, Markers ##
# See https://github.com/JuliaAttic/Color.jl/issues/75
function jet(n)
    RGB{Float64}[RGB(
        clamp(min(4x-1.5, -4x+4.5), 0, 1),
        clamp(min(4x-0.5, -4x+3.5), 0, 1),
        clamp(min(4x+0.5, -4x+2.5), 0, 1))
    for x in range(0, stop=1 , length=n)]
end

function linecolor(i, n=8)
    n = max(i, n)
    if n == 2
        return i == 1 ? "red" : "blue"
    elseif n == 3
        return i == 1 ? "red" : i == 2 ? "blue" : "green"
    else
        cm = jet(n)
        color = cm[i]
        string("{rgb,1:red,", color.r, "; green," , color.g, "; blue,", color.b, "}")
    end
end

function linestyle(i)
    LINESTYLES = ["solid", "dashed", "dotted", "dashdotted"]
    n = length(LINESTYLES)
    LINESTYLES[mod(i-1, n)+1]
end

function marker(i)
    MARKERS = ["*", "square*", "triangle*", "diamond*", "pentagon*", "x"]
    n = length(MARKERS)
    MARKERS[mod(i-1, n)+1]
end




#=
const LINESTYLES = ["dashed","dotted","dashdotted"]
const PGFPLOTSMARKERS = ["*","square*","triangle*","diamond*","pentagon*","x"]

tikz_header(x_label,y_label,optional,name,scale_factor,trim) = "%!TEX root = ../$(name)\n
\\begin{tikzpicture}$(trim ? "[trim axis left,trim axis right]" : "")
\\begin{axis}[
width=$(scale_factor)\\figurewidth,
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

tikz_add_bar_plot(name,color,legend) = "
\\addplot[area legend,draw opacity=0,ybar,bar width=0.6,draw=none,fill=$(color)] plot table[] {data/$(name).txt};
$(isempty(legend) ? "" : "\\addlegendentry{$(legend)};")
"

tikz_add_cuboid(x,y,z,z_height,color) = string("\\cuboid{",@sprintf("%2.1f",x-0.5),"}{",@sprintf("%2.1f",y-0.5),"}{",@sprintf("%7.5f",z),"}{",color,"}{0.8}{",@sprintf("%7.5f",z_height),"}\n")

figure_header() = "
\\begin{figure}[h]
\\centering
\\setlength{\\figureheight}{0.3\\textwidth}
\\setlength{\\figurewidth}{0.3\\textwidth}
"

figure_footer() = "
\\end{figure}
"

preamble(d,is_multigrid) = "$(d == 2 ? square()*(is_multigrid && d < 3 ? cuboid() : "") : d == 3 ? cube() : "\n\n")"

square() = "\n\n\\newcommand{\\drawsquare}[3]{
\\edef\\temp{
\\noexpand\\addplot[line width=3pt,white,fill= #3,forget plot,shift={(#1,#2)}]
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

cube() = "\n\n\\newcommand{\\drawcube}[4]{
\\edef\\temp{
\\noexpand\\addplot3[area legend,solid,fill= #4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
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
\\noexpand\\addplot3[area legend,solid,fill= #4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
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
\\noexpand\\addplot3[area legend,solid,fill= #4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
}
\\temp
table[row sep=crcr] {%
x	y	z\\\\
0	0	1\\\\
1	0	1\\\\
1	1	1\\\\
0	1	1\\\\
}--cycle;
}\n"

cuboid() = "\n\n\\newcommand{\\cuboid}[6]{
\\addplot3[area legend,solid,fill= #4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
table[row sep=crcr] {%
x  y  z\\\\
#5 0  0 \\\\
#5 0  #6\\\\
#5 #5 #6\\\\
#5 #5 0 \\\\
}--cycle;
\\addplot3[area legend,solid,fill= #4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
table[row sep=crcr] {%
x  y  z\\\\
0  #5 0\\\\
#5 #5 0\\\\
#5 #5 #6\\\\
0  #5 #6\\\\
}--cycle;
\\addplot3[area legend,solid,fill= #4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
table[row sep=crcr] {%
x  y  z\\\\
0  0  #6\\\\
#5 0  #6\\\\
#5 #5 #6\\\\
0  #5 #6\\\\
}--cycle;
}\n"

file_contents(title,d,is_adaptive,is_multigrid) = "\\documentclass[11pt, oneside]{article}

\\usepackage[margin=2cm]{geometry}
\\usepackage{graphicx}
\\usepackage{pgfplots}
\\pgfplotsset{compat=newest}
\\newlength\\figureheight
\\newlength\\figurewidth
\\usepackage{tikz}\n"*
"$((is_multigrid && d ==1) || d > 1 ? "\\usepackage{booktabs}\n\\usepackage{pgfplotstable}\n" : "") "*
preamble(d,is_multigrid)*"
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
"$(is_adaptive ? d < 4 ? "\\input{figures/adaptive_index_set.tex}\n" : "" : "" )"*
figure_header()*
"\\input{figures/runtime.tex}\n \\caption{\\label{fig:time}Total simulation run time.}"*
figure_footer()*"\n"*
figure_header()*
"\\input{figures/cost.tex}\n \\caption{\\label{fig:cost}Total simulation standard cost.}"*
figure_footer()*
"$(is_multigrid && d < 3 ? "\n\\input{figures/sample_reuse.tex}\n\n" : "")"*
"\\end{document}"

tikz_index_set_2d(name,i,max_level,is_adaptive,has_active,has_old,has_maximum,scaling,mode) = "%!TEX root = ../$(name)\n"*
tikz_load_index_set_tables(is_adaptive,i,has_active,has_old,has_maximum)*
"\\begin{tikzpicture}[trim axis left,trim axis right]
\\begin{axis}[
width=\\figurewidth$(scaling != 1 ? "/$(scaling)" : ""),
height=\\figureheight$(scaling != 1 ? "/$(scaling)" : ""),
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
]\n"*
tikz_draw_index_sets_2d_from_table(is_adaptive,has_active,has_old,has_maximum,mode)*
"\\end{axis}\n
\\end{tikzpicture}
"

function tikz_draw_index_sets_2d_from_table(is_adaptive,has_active,has_old,has_maximum,mode)
    if is_adaptive
        str = ""
        str = has_old ? string(str,tikz_draw_index_set_2d_from_table("oldset","white!90!black")) : str
        str = has_active ? string(str,tikz_draw_index_set_2d_from_table("activeset","orange!50!white")) : str
        str = has_maximum ? string(str,tikz_draw_index_set_2d_from_table("maximumindex","blue!50!white")) : str
        return str
    else
        if mode == "active"
            return tikz_draw_index_set_2d_from_table("indexset","orange!50!white")
        elseif mode == "maximum"
            return tikz_draw_index_set_2d_from_table("indexset","blue!50!white")
        else
            return tikz_draw_index_set_2d_from_table("indexset","white!90!black")
        end
    end
end

tikz_draw_index_set_2d_from_table(name,color) = 
"\\foreach \\j in {0,...,\\$(name)rows} {
\\pgfplotstablegetelem{\\j}{0}\\of\\$(name)
\\pgfmathsetmacro{\\a}{\\pgfplotsretval}
\\pgfplotstablegetelem{\\j}{1}\\of\\$(name)
\\pgfmathsetmacro{\\b}{\\pgfplotsretval}
\\drawsquare{\\a}{\\b}{$(color)}\n}\n"

function tikz_load_index_set_tables(is_adaptive,i,has_active,has_old,has_maximum)
    if is_adaptive
        str = ""
        str = has_active ? string(str,tikz_load_table("activeset","level_$(i)_active.txt")) : str
        str = has_old ? string(str,tikz_load_table("oldset","level_$(i)_old.txt")) : str
        str = has_maximum ? string(str,tikz_load_table("maximumindex","level_$(i)_maximum.txt")) : str
        return str
    else
        return tikz_load_table("indexset","index_set_$(i).txt")
    end
end

tikz_load_table(name,filename) = "\\pgfplotstableread[header=false]{data/$(filename)}\\$(name)%
\\pgfplotstablegetrowsof{\\$(name)}%
\\pgfmathsetmacro{\\$(name)rows}{\\pgfplotsretval-1}%\n"

tikz_index_set_3d(name,i,max_level,is_adaptive,has_active,has_old,has_maximum,scaling,mode) = "%!TEX root = ../$(name)\n"*
tikz_load_index_set_tables(is_adaptive,i,has_active,has_old,has_maximum)*
"\\begin{tikzpicture}[trim axis left,trim axis right]
\\begin{axis}[
width=\\figurewidth$(scaling != 1 ? "/$(scaling)" : ""),
height=\\figureheight$(scaling != 1 ? "/$(scaling)" : ""),
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
]\n"*
tikz_draw_index_sets_3d_from_table(is_adaptive,has_active,has_old,has_maximum,mode)*
"\\end{axis}\n
\\end{tikzpicture}
"

function tikz_draw_index_sets_3d_from_table(is_adaptive,has_active,has_old,has_maximum,mode)
    if is_adaptive
        str = ""
        str = has_old ? string(str,tikz_draw_index_set_3d_from_table("oldset","white!90!black")) : str
        str = has_active ? string(str,tikz_draw_index_set_3d_from_table("activeset","orange!50!white")) : str
        str = has_maximum ? string(str,tikz_draw_index_set_3d_from_table("maximumindex","blue!50!white")) : str
        return str
    else
        if mode == "active"
            return tikz_draw_index_set_3d_from_table("indexset","orange!50!white")
        elseif mode == "maximum"
            return tikz_draw_index_set_3d_from_table("indexset","blue!50!white")
        else
            return tikz_draw_index_set_3d_from_table("indexset","white!90!black")
        end
    end
end

tikz_draw_index_set_3d_from_table(name,color) = 
"\\foreach \\j in {0,...,\\$(name)rows} {
\\pgfplotstablegetelem{\\j}{2}\\of\\$(name)
\\pgfmathsetmacro{\\a}{\\pgfplotsretval}
\\pgfplotstablegetelem{\\j}{0}\\of\\$(name)
\\pgfmathsetmacro{\\b}{\\pgfplotsretval}
\\pgfplotstablegetelem{\\j}{1}\\of\\$(name)
\\pgfmathsetmacro{\\c}{\\pgfplotsretval}
\\drawcube{\\a}{\\b}{\\c}{$(color)}\n}\n"

tikz_samples_table(name,tols) = tikz_index_set_table_internal(name,"sample_reuse_",tols,"\$\\epsilon = \\pgfmathprintnumber[/pgf/number format/sci,precision=3,sci zerofill]",false,"\\label{fig:samples_reused}Total number of samples and percentage of reused samples for different tolerances on the RMSE") 

tikz_index_set_table(name,tols) = tikz_index_set_table_internal(name,"index_set_",tols,"\$\\epsilon = \\pgfmathprintnumber[/pgf/number format/sci,precision=3,sci zerofill]",false,"\\label{fig:index_set}Shape of the index set for different tolerances on the RMSE") 

tikz_adaptive_index_set_table(name,levels) = tikz_index_set_table_internal(name,"adaptive_index_set_",levels,"\$L = \\pgfmathprintnumber[/pgf/number format/fixed]",true,"\\label{fig:adaptive_index_set}Shape of the index set for different level parameters in the adaptive algorithm") 

function tikz_index_set_table_internal(name,fname,iters,headerstring,is_adaptive,legend)
    str = "%!TEX root = ../$(name)\n"
    for tabnum = 1:div(length(iters)-1,9)+1
        str = string(str,figure_header())
        str = string(str,"\\begin{tabular}{ccc}")
        # add adaptive legend
        if is_adaptive
            str = string(str,"\\multicolumn{1}{l}{\\raisebox{-0.025\\textwidth}{\\input{figures/index_set_legend_old.tex}} = old set} &\n")
            str = string(str,"\\multicolumn{1}{l}{\\raisebox{-0.025\\textwidth}{\\input{figures/index_set_legend_active.tex}} = active set} &\n")
            str = string(str,"\\multicolumn{1}{l}{\\raisebox{-0.025\\textwidth}{\\input{figures/index_set_legend_maximum.tex}} = maximum profit} \\\\\n")
        end
        for col in 1:min(3,div(length(iters)-(tabnum-1)*9-1,3)+1)
            # header
            str = string(str," \\toprule\n",headerstring,"{", @sprintf("%7.3e",iters[(tabnum-1)*9+(col-1)*3+1]), "}\$ &")
            str = length(iters) > (tabnum-1)*9+(col-1)*3+1 ? string(str," ",headerstring,"{", @sprintf("%7.3e",iters[(tabnum-1)*9+(col-1)*3+2]), "}\$") : str
            str = string(str," &")
            str = length(iters) > (tabnum-1)*9+(col-1)*3+2 ? string(str," ",headerstring,"{", @sprintf("%7.3e",iters[(tabnum-1)*9+(col-1)*3+3]), "}\$") : str
            str = string(str,"\\\\ \\midrule \n")
            str = string(str," \\input{figures/$(fname)",(tabnum-1)*9+(col-1)*3+1, ".tex} &")
            str = length(iters) > (tabnum-1)*9+(col-1)*3+1 ? string(str," \\input{figures/$(fname)",(tabnum-1)*9+(col-1)*3+2, ".tex}") : str
            str = string(str," &")
            str = length(iters) > (tabnum-1)*9+(col-1)*3+2 ? string(str," \\input{figures/$(fname)",(tabnum-1)*9+(col-1)*3+3, ".tex}") : str
            str = string(str,"\\\\ \n")
        end
        str = string(str,"\n\\end{tabular}\\\\\n\n")
        str = string(str,"\\caption{$(legend) ($(tabnum)/$(div(length(iters)-1,9)+1)).}")
        str = string(str,figure_footer())
    end
    str
end

tikz_read_fractions(i) = "
\\pgfplotstableread[header=false]{data/fractions_$(i).txt}\\fractions
\\pgfplotstablegetrowsof{\\fractions}
\\pgfmathsetmacro{\\rows}{\\pgfplotsretval-1}
"

tikz_add_fractions_1d(i) = tikz_read_fractions(i)*"
\\pgfplotsset{
after end axis/.code={
\\foreach \\j in {0,...,\\rows} {
\\pgfplotstablegetelem{\\j}{0}\\of\\fractions
\\pgfmathsetmacro{\\a}{\\pgfplotsretval}
\\pgfplotstablegetelem{\\j}{1}\\of\\fractions
\\pgfmathsetmacro{\\b}{\\pgfplotsretval}
\\pgfplotstablegetelem{\\j}{2}\\of\\fractions
\\pgfmathsetmacro{\\c}{\\pgfplotsretval}
\\node[right, align=left, text=black, anchor=south] at (axis cs:\\a,\\b) {\\tiny \\c\\%};
}
}
}
"

tikz_bar_plot_1d_header(L,M,fname) = tikz_header("level \$\\ell\$","number of samples \$N_\\ell\$",
                                                 string("ymode=log,\n",
                                                        "xmin=-.3,\n",
                                                        "xmax=",@sprintf("%2.1f",L+0.3),",\n",
                                                        "ymin=1,\n",
                                                        "ymax=10^",@sprintf("%i",ceil(log10(M))),",\n",
                                                        "xtick={0,1,...,$(L)},\n",
                                                        "axis x line*=left,\n",
                                                        "x axis line style={draw opacity=0},\n",
                                                        "xtick style={draw=none},\n",
                                                        "axis y line*=left,\n",
                                                        "y axis line style={draw opacity=0},\n",
                                                        "ytick style={draw=none},\n",
                                                        "ymajorgrids,\n",
                                                        "grid style={line width=1pt,white},\n",
                                                        "axis on top,\n",
                                                        "legend style={legend cell align=left,align=left,font=\\tiny,",
                                                        "draw=none,at={(1.03,1.03)},anchor=north east}\n"),
                                                 fname,
                                                 @sprintf("%2.1f",0.8), 
                                                 false
                                                 )

tikz_bar_plot_2d_header(L,M,fname) = tikz_header("\$\\ell_x\$","\$\\ell_y\$",
                                                 string("view={120}{17},\n",
                                                        "every x tick/.append style={draw=none},\n",
                                                        "xmin=-.3,\n",
                                                        "xmax=",@sprintf("%2.1f",L+1+0.3),",\n",
                                                        "xtick={0,1,...,",@sprintf("%i",L),"},",
                                                        "every y tick/.append style={draw=none},\n",
                                                        "ymin=-.3,\n",
                                                        "ymax=",@sprintf("%2.1f",L+1),",\n",
                                                        "ytick={0,1,...,",@sprintf("%i",L),"},\n",
                                                        "every z tick label/.append style=",
                                                        "{font=\\color{white!15!black}\\scriptsize},\n",
                                                        "every z label/.append style=",
                                                        "{font=\\color{white!15!black}\\scriptsize},\n",
                                                        "zmin=0,\n",
                                                        "zmax=",@sprintf("%i",ceil(M)),",\n",
                                                        "every z tick/.append style={draw=none},\n",
                                                        "ztick={0,1,...,",@sprintf("%i",ceil(M)),"},\n",
                                                        "zticklabels={",
                                                        string(join([@sprintf("\$10^%i\$",i) for i in 0:ceil(M)],",")),
                                                        "},\n",
                                                        "zmajorgrids=true,\n",
                                                        "set layers,\n",
                                                        "zlabel=number of samples \$N_\\ell\$,\n",
                                                        "legend style={legend cell align=left,align=left,font=\\tiny,",
                                                        "draw=none,at={(1.03,0.9)},anchor=north east}"),
                                                 fname, 
                                                 @sprintf("%2.1f",0.8), 
                                                 false
                                                 )

tikz_3d_bar_plot_legend() = "
\\addlegendimage{area legend,solid,fill=red,draw opacity=0}
\\addlegendentry{original}
\\addlegendimage{area legend,solid,fill=blue,draw opacity=0}
\\addlegendentry{reused}
"
=#
