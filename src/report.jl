## report.jl : loads a history file and make a report with the main results

function report(h::History;folder="report $(h[:name])"::AbstractString)

    # problem file name
    fname = isempty(h[:name]) ? "report.tex" : "report_$(h[:name]).tex"

    # make the required directories
    !isdir(folder) && mkdir(folder)
    !isdir(joinpath(folder,"data")) && mkdir(joinpath(folder,"data"))
    !isdir(joinpath(folder,"figures")) && mkdir(joinpath(folder,"figures"))

    # make all data files
    write_E(h,folder)
    write_dE(h,folder)
    write_V(h,folder)
    write_dV(h,folder)
    write_W(h,folder)
    write_samples(h,folder)
    write_time(h,folder)
    write_cost(h,folder)

    # make the tikz files
    write_figure_E(h,folder,fname)
    write_figure_V(h,folder,fname)
    write_figure_samples(h,folder,fname)
    write_figure_cost(h,folder,fname)
    write_figure_time(h,folder,fname)

    # make the report file
    write_main_file(h,folder,fname)

end

function write_E(h::History, folder::AbstractString)
    open(joinpath(folder,"data","E.txt"), "w") do f
        L = h[:max_level]
        x = 0:L
        y = log2.(abs.(h[:E]))
        for i in 1:L+1
            write(f, @sprintf("%i %12.5f\n",x[i],y[i]))
        end
    end
end

function write_dE(h::History, folder::AbstractString)
    open(joinpath(folder,"data","dE.txt"), "w") do f
        L = h[:max_level]
        x = 1:L
        y = log2.(abs.(h[:dE]))[2:end]
        for i in 1:L
            write(f, @sprintf("%i %12.5f\n",x[i],y[i]))
        end
    end
end

function write_V(h::History, folder::AbstractString)
    open(joinpath(folder,"data","V.txt"), "w") do f
        L = h[:max_level]
        x = 0:L
        y = log2.(h[:V])
        for i in 1:L+1
            write(f, @sprintf("%i %12.5f\n",x[i],y[i]))
        end
    end
end

function write_dV(h::History, folder::AbstractString)
    open(joinpath(folder,"data","dV.txt"), "w") do f
        L = h[:max_level]
        x = 1:L
        y = log2.(h[:dV])[2:end]
        for i in 1:L
            write(f, @sprintf("%i %12.5f\n",x[i],y[i]))
        end
    end
end

function write_W(h::History, folder::AbstractString)
    open(joinpath(folder,"data","W.txt"), "w") do f
        L = h[:max_level]
        x = 0:L
        y = log2.(h[:W])
        for i in 1:L+1
            write(f, @sprintf("%i %12.5e\n",x[i],y[i]))
        end
    end
end

function write_samples(h::History, folder::AbstractString)
    open(joinpath(folder,"data","samples.txt"), "w") do f
        m = h[:max_level] + 1
        n = h.iter
        S = zeros(m,n)
        for j = 1:n
            s = h[j][:nsamples]
            S[1:length(s),j] = s
        end
        for i in 1:m
            str = @sprintf("%i",i-1)
            for j in 1:n
                str = string(str,@sprintf(" %i",S[i,j]))
            end
            write(f,string(str,"\n"))
        end
    end
end

function write_time(h::History, folder::AbstractString)
    open(joinpath(folder,"data","time.txt"), "w") do f
        x = [h[i][:tol] for i in 1:h.iter]
        y = cumsum([h[i][:runtime] for i in 1:h.iter])
        for i in 1:length(x)
            write(f, @sprintf("%12.5f %12.5e\n",x[i],y[i]))
        end
    end
end

function write_cost(h::History, folder::AbstractString)
    open(joinpath(folder,"data","cost.txt"), "w") do f
        x = [h[i][:tol] for i in 1:h.iter]
        y = [h[i][:cost] for i in 1:h.iter]
        for i in 1:length(x)
            write(f, @sprintf("%12.5f %12.5e\n",x[i],y[i]))
        end
    end
end

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

tikz_add_plot(color_name,line_style,name,legend_name,table_options) = "
\\addplot [color=$(color_name),$(line_style),mark=*,mark options={solid,fill=$(color_name)},line width=0.75pt,mark size=1.2]
table[$(table_options)]{data/$(name).txt};
$(isempty(legend_name) ? "" : "\\addlegendentry{$(legend_name)};")
"

function write_figure_E(h::History, folder::AbstractString, fname::AbstractString)
    open(joinpath(folder,"figures","E.tex"), "w") do f
        L = h[:max_level]
        str = tikz_header("level \$\\ell\$",raw"$\log_2(|E[\;\cdot\;]|)$","xmin=0,\n xmax=$(L),\n xtick={0,1,...,$(L)},\n",fname)
        str = string(str,tikz_add_plot("red","solid","E",raw"$Q_\ell$",""))
        str = string(str,tikz_add_plot("red","dashed","dE",raw"$\Delta Q_\ell$",""))
        str = string(str,tikz_footer())
        write(f, str)
    end
end

function write_figure_V(h::History, folder::AbstractString, fname::AbstractString)
    open(joinpath(folder,"figures","V.tex"), "w") do f
        L = h[:max_level]
        str = tikz_header("level \$\\ell\$",raw"$\log_2(V[\;\cdot\;])$","xmin=0,\n xmax=$(L),\n xtick={0,1,...,$(L)},\n",fname)
        str = string(str,tikz_add_plot("blue","solid","V",raw"$Q_\ell$",""))
        str = string(str,tikz_add_plot("blue","dashed","dV",raw"$\Delta Q_\ell$",""))
        str = string(str,tikz_footer())
        write(f, str)
    end
end

function write_figure_samples(h::History, folder::AbstractString, fname::AbstractString)
    open(joinpath(folder,"figures","samples.tex"), "w") do f
        L = h[:max_level]
        str = tikz_header("level \$\\ell\$","number of samples \$N_\\ell\$","ymode=log,\n xmin=0,\n xmax=$(L),\n xtick={0,1,...,$(L)},\n",fname)
        cm = ColorMap("jet")
        v = linspace(0,1,h.iter)
        for i in 1:h.iter
            color = cm(v[i])
            color_name = "{rgb,1:red,$(color[1]); green,$(color[2]); blue,$(color[3])}"
            str = string(str,tikz_add_plot(color_name,"solid","samples","","x index = {0}, y index = {$(i)}"))
        end
        str = string(str,tikz_footer())
        write(f, str)
    end
end

function write_figure_time(h::History, folder::AbstractString, fname::AbstractString)
    open(joinpath(folder,"figures","time.tex"), "w") do f
        str = tikz_header("accuracy","run time","xmode = log,\n ymode = log",fname)
        str = string(str,tikz_add_plot("blue","solid","time","",""))
        str = string(str,tikz_footer())
        write(f, str)
    end
end

function write_figure_cost(h::History, folder::AbstractString, fname::AbstractString)
    open(joinpath(folder,"figures","cost.tex"), "w") do f
        str = tikz_header("accuracy","standard cost","xmode = log,\n ymode = log",fname)
        str = string(str,tikz_add_plot("red","solid","cost","",""))
        str = string(str,tikz_footer())
        write(f, str)
    end
end

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
"\\input{figures/E.tex}\n \\caption{\\label{fig:E}Decay of the expected value \$E[\\Delta Q_\\ell]\$.}"*
figure_footer()*"\n"*
figure_header()*
"\\input{figures/V.tex}\n \\caption{\\label{fig:V}Decay of the variance \$V[\\Delta Q_\\ell]\$.}"*
figure_footer()*"\n"*
figure_header()*
"\\input{figures/samples.tex}\n \\caption{\\label{fig:samples}Total number of samples \$N_\\ell\$ taken on each level.}"*
figure_footer()*"\n"*
figure_header()*
"\\input{figures/time.tex}\n \\caption{\\label{fig:time}Total simulation run time.}"*
figure_footer()*"\n"*
figure_header()*
"\\input{figures/cost.tex}\n \\caption{\\label{fig:cost}Total simulation standard cost.}"*
figure_footer()*
"\\end{document}"

function write_main_file(h::History,folder::AbstractString,fname::AbstractString)
    open(joinpath(folder,fname), "w") do f
        write(f, file_contents(h[:name]))
    end
end
