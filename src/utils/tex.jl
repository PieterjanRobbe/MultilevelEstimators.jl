## tex.jl : tex utilities for making reports

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

preamble(d,is_multigrid) = "$(d == 2 ? square()*(is_multigrid && d < 3? cuboid() : "") : d == 3 ? cube() : "\n\n")"

square() = "\n\n\\newcommand{\\drawsquare}[3]{
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

cube() = "\n\n\\newcommand{\\drawcube}[4]{
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

cuboid() = "\n\n\\newcommand{\\cuboid}[6]{
\\addplot3[area legend,solid,fill=#4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
table[row sep=crcr] {%
x  y  z\\\\
#5 0  0 \\\\
#5 0  #6\\\\
#5 #5 #6\\\\
#5 #5 0 \\\\
}--cycle;
\\addplot3[area legend,solid,fill=#4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
table[row sep=crcr] {%
x  y  z\\\\
0  #5 0\\\\
#5 #5 0\\\\
#5 #5 #6\\\\
0  #5 #6\\\\
}--cycle;
\\addplot3[area legend,solid,fill=#4,draw=black,rounded corners=0.2pt,forget plot,shift={(#1,#2,#3)}]
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

tikz_load_table(name,filename) = "\\pgfplotstableread[header=false]{data/$(filename)}\\$(name)\%
\\pgfplotstablegetrowsof{\\$(name)}\%
\\pgfmathsetmacro{\\$(name)rows}{\\pgfplotsretval-1}\%\n"

tikz_index_set_3d(name,i,max_level,is_adaptive) = "%!TEX root = ../$(name)\n
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
			str = string(str,"\\multicolumn{1}{l}{\\raisebox{-0.025\\textwidth}{\\input{figures/index_set_legend_old.tex}} = old set} \&\n")
			str = string(str,"\\multicolumn{1}{l}{\\raisebox{-0.025\\textwidth}{\\input{figures/index_set_legend_active.tex}} = active set} \&\n")
			str = string(str,"\\multicolumn{1}{l}{\\raisebox{-0.025\\textwidth}{\\input{figures/index_set_legend_maximum.tex}} = maximum profit} \\\\\n")
		end
		for col in 1:min(3,div(length(iters)-(tabnum-1)*9-1,3)+1)
            # header
            str = string(str," \\toprule\n",headerstring,"\{", @sprintf("%7.3e",iters[(tabnum-1)*9+(col-1)*3+1]), "\}\$ \&")
            str = length(iters) > (tabnum-1)*9+(col-1)*3+1 ? string(str," ",headerstring,"\{", @sprintf("%7.3e",iters[(tabnum-1)*9+(col-1)*3+2]), "\}\$") : str
            str = string(str," \&")
            str = length(iters) > (tabnum-1)*9+(col-1)*3+2 ? string(str," ",headerstring,"\{", @sprintf("%7.3e",iters[(tabnum-1)*9+(col-1)*3+3]), "\}\$") : str
            str = string(str,"\\\\ \\midrule \n")
			str = string(str," \\input\{figures/$(fname)",(tabnum-1)*9+(col-1)*3+1, ".tex\} \&")
			str = length(iters) > (tabnum-1)*9+(col-1)*3+1 ? string(str," \\input\{figures/$(fname)",(tabnum-1)*9+(col-1)*3+2, ".tex\}") : str
            str = string(str," \&")
			str = length(iters) > (tabnum-1)*9+(col-1)*3+2 ? string(str," \\input\{figures/$(fname)",(tabnum-1)*9+(col-1)*3+3, ".tex\}") : str
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
