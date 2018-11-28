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

function Document(h::History, folder::AbstractString)
    
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
