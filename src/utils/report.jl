## report.jl : create a summary of the estimator
#
# Create a report that contains useful information about the estimator.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## report ##
function report(h::History, folder::String=h[:name][1:end-4])

    # make the required directories
    !isdir(folder) && mkdir(folder)
    !isdir(joinpath(folder,"data")) && mkdir(joinpath(folder,"data"))
    !isdir(joinpath(folder,"figures")) && mkdir(joinpath(folder,"figures"))

    # create report
    doc = Document(h, folder)
    open(joinpath(folder, "report.tex"), "w") do f
        write(f, generate(doc))
    end
    
end

## rates ##
function write_data_rates(h::History, folder::String)
    d = h[:ndims]
    for name in ["E", "dE", "V", "dV", "W"]
        for idx in Base.Iterators.drop(CartesianIndices(tuple(fill(0:1, d)...)), 1)
            cntr = first(name) == 'd' ? 1 : 0
            x = Int64[]
            y = Float64[]
            while (cntr*idx).I ∈ h[:index_set]
                push!(x, cntr)
                push!(y, log2(abs(h[Symbol(name)][(cntr*idx).I])))
                cntr += 1
            end
            open(joinpath(folder, "data", string(name, join(string.(idx.I)), ".txt")), "w") do f
                for n in 1:length(x)
                    write(f, @sprintf("%i %12.5e\n", x[n], y[n]))
                end
            end
        end
    end
end

## samples ##
function write_data_samples(h::History, folder::String)
    d = h[:ndims]
    if d == 1
        n = first(maximum(h[:index_set]))
        m = length(h)
        S = fill(0, n+1, m+1)
        S[:,1] = 0:n
        for i in 1:n+1, j in 2:m+1
            S[i,j] = get(h[j-1][:nb_of_samples], Level(i-1), 0)
        end
        writedlm(joinpath(folder, "data", "nb_of_samples.txt"), S)
    else
        error("not implemented yet")
    end
end             

## complexity ##
function write_data_runtime(h::History, folder::String)
    x = [h[i][:tol] for i in 1:length(h)]
    y = cumsum([Dates.value(h[i][:time_stamp]-h.t_start)/1000 for i in 1:length(h)])
    open(joinpath(folder, "data", "runtime.txt"), "w") do f
        for i in 1:length(x)
            write(f, @sprintf("%12.5e %12.5e\n", x[i], y[i]))
        end
    end
end

function write_data_cost(h::History, folder::String)
    x = [h[i][:tol] for i in 1:length(h)]
    y = [sum(h[:W][index]*h[i][:nb_of_samples][index] for index in h[i][:index_set]) for i in 1:length(h)]
    open(joinpath(folder, "data", "cost.txt"), "w") do f
        for i in 1:length(x)
            write(f, @sprintf("%12.5e %12.5e\n", x[i], y[i]))
        end
    end
end

# TODO maybe use @generated function with dispatch in h[:type]
function write_figures_rates(h::History, folder::String)
    d = h[:ndims]
    if d == 1
        L = first(maximum(h[:index_set]))
        # E and V
        for (i, name) in enumerate(["E", "V"])
            plots = LinePlot[]
            plot1 = LinePlot(linecolor(i, 3), linestyle(1), marker(1), "\$Q_\\ell\$", string(name, 1), "")
            plot2 = LinePlot(linecolor(i, 3), linestyle(2), marker(1), "\$Q_\\ell-Q_{\\ell-1}\$", string("d", name, 1), "") 
            push!(plots, plot1)
            push!(plots, plot2)
            bar = i == 1 ? "|" : ""
            y_label = string("\$\\log_2(", bar, "\\mathbb{", name, "}[\\;\\cdot\\;]", bar, ")\$")
            pic = TikzPicture("level \$\\ell\$", y_label, 0, L, string("0,1,...,", L), "", "", "", "", plots)
            open(joinpath(folder, "figures", string(name, ".tex")), "w") do f
                write(f, generate(pic))
            end
        end
        # W
        plots = LinePlot[]
        plot = LinePlot(linecolor(3, 3), linestyle(1), marker(1), "\$Q_\\ell-Q_{\\ell-1}\$", "W1", "")
        push!(plots, plot)
        pic = TikzPicture("level \$\\ell\$", "\$\\log_2(\\mathbb{W}[\\;\\cdot\\;])\$", 0, L, string("0,1,...,", L), "", "", "(0.03,0.97)", "north west", plots)
        open(joinpath(folder, "figures", "W.tex"), "w") do f
            write(f, generate(pic))
        end
    else
        error("not implemented yet")
    end
end

function write_figure_samples(h::History, folder::String)
    d = h[:ndims]
    if d == 1
        L = first(maximum(h[:index_set]))
        plots = LinePlot[]
        for i in length(h):-1:1
            plot = LinePlot(linecolor(i, max(2, length(h))), linestyle(1), marker(1), string("\$\\epsilon = \$", shorte(h[i][:tol])) , "nb_of_samples", string("x index = {0}, y index = {", i, "}"))
            push!(plots, plot)
        end
        pic = TikzPicture("level \$\\ell\$", "number of samples \$N_\\ell\$", 0, L, string("0,1,...,", L), "", "log", "(1.03,0.5)", "west", plots)
        open(joinpath(folder, "figures", "nb_of_samples.tex"), "w") do f
            write(f, generate(pic))
        end
    else
        error("not implemented yet")
    end
end

function write_figure_runtime(h::History, folder::String)
    plots = LinePlot[]
    plot = LinePlot(linecolor(1, 2), linestyle(1), marker(1), "" , "runtime", "")
    push!(plots, plot)
    pic = TikzPicture("accuracy \$\\epsilon\$", "runtime [s]", "", "", "", "log", "log", "", "", plots)
    open(joinpath(folder, "figures", "runtime.tex"), "w") do f
        write(f, generate(pic))
    end
end

function write_figure_cost(h::History, folder::String)
    plots = LinePlot[]
    plot = LinePlot(linecolor(2, 2), linestyle(1), marker(1), "" , "cost", "")
    push!(plots, plot)
    pic = TikzPicture("accuracy \$\\epsilon\$", "standard cost", "", "", "", "log", "log", "", "", plots)
    open(joinpath(folder, "figures", "cost.tex"), "w") do f
        write(f, generate(pic))
    end
end

#=

for (name,color,ylabel) in zip(["E","V"],["red","blue"],["|E[\\;\\cdot\\;]|","V[\\;\\cdot\\;]"])
    ex = :(
           function $(Symbol("write_figure_",name))(h::History, folder::AbstractString, fname::AbstractString)
               open(joinpath(folder,"figures",string($(name),".tex")), "w") do f
                   L = maximum(maximum.(h[:index_set]))
                   d = length(h[:index_set][1])
                   str = tikz_header("level \$\\ell\$",string("\$\\log_2(",$(ylabel),")\$"),
                                     "xmin=0,\n xmax=$(L),\n xtick={0,1,...,$(L)},\n",fname,"",true)
                   if d > 1
                       for i in 2:2^d
                           proto_idx = ind2sub(tuple([2 for i = 1:d]...),i).-1
                           marker_idx = mod(i-2,length(PGFPLOTSMARKERS))+1
                           linestyle_idx = mod(sum(proto_idx)-1,length(LINESTYLES))+1
                           name = string($(name),join(string.(proto_idx)))
                           str = string(str,tikz_add_plot($(color),"solid",PGFPLOTSMARKERS[marker_idx],name,"",""))
                           legend_name = string("\$\\Delta Q_{",replace(string(proto_idx),"1","\\ell"),"}\$")
                           str = string(str,tikz_add_plot($(color),LINESTYLES[linestyle_idx],
                                                          PGFPLOTSMARKERS[marker_idx],
                                                          string("d",name),legend_name,"")
                                       )
                       end
                   else
                       str = string(str,tikz_add_plot($(color),"solid",PGFPLOTSMARKERS[1],
                                                      string($(name),"1"),"\$Q_\\ell\$",""))
                       str = string(str,tikz_add_plot($(color),LINESTYLES[1],PGFPLOTSMARKERS[1],
                                                      string("d",$(name),"1"),"\$Q_\\ell-Q_{\\ell-1}\$",""))
                   end
                   str = string(str,tikz_footer())
                   write(f, str)
               end
           end
          )
    eval(ex)
end
=#











#=
function report(h::History, folder::AbstractString=h[:name])

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
    write_index_set(h,folder)
    write_adaptive_index_set(h,folder)
    write_runtime(h,folder)
    write_runtime_vs_rmse(h,folder)
    write_cost(h,folder)
    write_sample_reuse(h,folder)

    # make the tikz files
    write_figure_E(h,folder,fname)
    write_figure_V(h,folder,fname)
    write_figure_W(h,folder,fname)
    write_figure_samples(h,folder,fname)
    write_figure_index_set(h,folder,fname)
    write_figure_adaptive_index_set(h,folder,fname)
    write_figure_runtime(h,folder,fname)
    write_figure_runtime_vs_rmse(h,folder,fname)
    write_figure_cost(h,folder,fname)
    write_figure_sample_reuse(h,folder,fname)

    # make the report file
    write_main_file(h,folder,fname)
end

for (name,start) in zip(["E","dE","V","dV","W"],[1,2,1,2,1])
    ex = :(
           function $(Symbol("write_",name))(h::History, folder::AbstractString)
               d = h[:ndims]
               for i in 2:2^d
                   proto_idx = ind2sub(tuple([2 for i = 1:d]...),i).-1
                   cntr = 0
                   x = Float64[]
                   y = Float64[]
                   while cntr.*proto_idx ∈ h[:index_set]
                       idx = findfirst([cntr.*proto_idx].==h[:index_set])
                       push!(x,cntr)
                       push!(y,log2(abs(h[Symbol(Symbol($(name)))][idx])))
                       cntr += 1
                   end
                   open(joinpath(folder,"data",string($(name),join(string.(proto_idx)),".txt")), "w") do f
                       for j in $(start):length(x)
                           write(f, @sprintf("%i %12.5e\n",x[j],y[j]))
                       end
                   end
               end
           end
          )
    eval(ex)
end

function get_samples(h::History)
    d = h[:ndims]
    idx_set = h[:index_set]
    m = maximum.([getindex.(idx_set,i) for i in 1:d]) .+ 1
    n = h.iter
    S_ = zeros(Int,m...,n)
    h[:index_set]
    h[:nsamples]
    for j = 1:n
        s = h[j][:nsamples]
        for (i,idx) in enumerate(h[j][:index_set])
            S_[idx.+1...,j] = s[i]
        end
    end
    S = copy(S_)
    for idx in Iterators.product(colon.(1,m)...)
        for i in Iterators.product(colon.(1,idx)...)
            if i != idx
                for j in 1:size(S_,d+1)
                    S[i...,j] += S_[idx...,j]
                end
            end
        end
    end
    return S,S_
end

function write_samples(h::History, folder::AbstractString)
    d = h[:ndims]
    if d == 1
        S,S_ = get_samples(h)
        open(joinpath(folder,"data","samples.txt"), "w") do f
            for i in 1:size(S_,1)
                str = @sprintf("%i",i-1)
                for j in 1:size(S_,2)
                    str = string(str,@sprintf(" %i",S_[i,j]))
                end
                write(f,string(str,"\n"))
            end
        end
    end
end

function write_sample_reuse(h::History, folder::AbstractString)
    d = h[:ndims]
    if h[:method] isa MG
        S,S_ = get_samples(h)
        for j in 1:h.iter
            open(joinpath(folder,"data","samples_total_$(j).txt"), "w") do f
                for idx in Iterators.product(colon.(1,map(i->size(S,i),1:d))...)
                    if S[idx...,j] > 0
                        for i in 1:d
                            write(f,@sprintf("%i ",idx[i]-1))
                        end
                        write(f,@sprintf("%i\n",S[idx...,j]))
                    end
                end
            end
            open(joinpath(folder,"data","samples_reused_$(j).txt"), "w") do f
                for idx in Iterators.product(colon.(1,map(i->size(S,i),1:d))...)
                    if S[idx...,j] > 0
                        for i in 1:d
                            write(f,@sprintf("%i ",idx[i]-1))
                        end
                        write(f,@sprintf("%i\n",S[idx...,j]-S_[idx...,j]))
                    end
                end
            end
        end
        write_fractions(h, folder)
    end
end

function write_fractions(h::History, folder::AbstractString)
    S,S_ = get_samples(h)
    F = (S.-S_)./S*100
    for j in 1:h.iter
        open(joinpath(folder,"data","fractions_$(j).txt"), "w") do f
            for idx in Iterators.product(colon.(1,map(i->size(S,i),1:h[:ndims]))...)
                if S[idx...,j] > 0
                    for i in 1:h[:ndims]
                        write(f,@sprintf("%i ",idx[i]-1))
                    end
                    write(f,@sprintf("%i %i\n",S[idx...,j],ceil(F[idx...,j])))
                end
            end
        end
    end
end

function write_index_set(h::History, folder::AbstractString)
    d = h[:ndims]
    tols = [h[i][:tol] for i in 1:h.iter]
    for i in 1:length(tols)
        open(joinpath(folder,"data",string("index_set_",i,".txt")), "w") do f
            index_set = h[i][:index_set]
            for idx in index_set
                for j in 1:d
                    write(f, @sprintf("%i ",idx[j]))
                end
                write(f, "\n")
            end
        end
    end
end

function write_adaptive_index_set(h::History, folder::AbstractString)
    d = h[:ndims]
    dicts = h[:adaptive_index_set]
    for i in 1:length(dicts)
        for (val,name,mode) in zip(0:3,["old","active","active","maximum"],["w" "w" "a" "w"])
            open(joinpath(folder,"data",string("level_",i,"_",name,".txt")), mode) do f
                index_set = collect(keys(dicts[i]))[collect(values(dicts[i])).==val] 
                for idx in index_set
                    for j in 1:d
                        write(f, @sprintf("%i ",idx[j]))
                    end
                    write(f, "\n")
                end
            end # open
        end
    end
end

function write_runtime_vs_rmse(h::History, folder::AbstractString)
    open(joinpath(folder,"data","runtime_vs_rmse.txt"), "w") do f
        x = [h[i][:rmse] for i in 1:h.iter]
        y = cumsum([h[i][:runtime] for i in 1:h.iter])
        for i in 1:length(x)
            write(f, @sprintf("%12.5e %12.5e\n",x[i],y[i]))
        end
    end
end

function write_runtime(h::History, folder::AbstractString)
    open(joinpath(folder,"data","runtime.txt"), "w") do f
        x = [h[i][:tol] for i in 1:h.iter]
        y = cumsum([h[i][:runtime] for i in 1:h.iter])
        for i in 1:length(x)
            write(f, @sprintf("%12.5e %12.5e\n",x[i],y[i]))
        end
    end
end

function write_cost(h::History, folder::AbstractString)
    open(joinpath(folder,"data","cost.txt"), "w") do f
        x = [h[i][:tol] for i in 1:h.iter]
        y = [h[i][:cost] for i in 1:h.iter]
        for i in 1:length(x)
            write(f, @sprintf("%12.5e %12.5e\n",x[i],y[i]))
        end
    end
end

for (name,color,ylabel) in zip(["E","V"],["red","blue"],["|E[\\;\\cdot\\;]|","V[\\;\\cdot\\;]"])
    ex = :(
           function $(Symbol("write_figure_",name))(h::History, folder::AbstractString, fname::AbstractString)
               open(joinpath(folder,"figures",string($(name),".tex")), "w") do f
                   L = maximum(maximum.(h[:index_set]))
                   d = length(h[:index_set][1])
                   str = tikz_header("level \$\\ell\$",string("\$\\log_2(",$(ylabel),")\$"),
                                     "xmin=0,\n xmax=$(L),\n xtick={0,1,...,$(L)},\n",fname,"",true)
                   if d > 1
                       for i in 2:2^d
                           proto_idx = ind2sub(tuple([2 for i = 1:d]...),i).-1
                           marker_idx = mod(i-2,length(PGFPLOTSMARKERS))+1
                           linestyle_idx = mod(sum(proto_idx)-1,length(LINESTYLES))+1
                           name = string($(name),join(string.(proto_idx)))
                           str = string(str,tikz_add_plot($(color),"solid",PGFPLOTSMARKERS[marker_idx],name,"",""))
                           legend_name = string("\$\\Delta Q_{",replace(string(proto_idx),"1","\\ell"),"}\$")
                           str = string(str,tikz_add_plot($(color),LINESTYLES[linestyle_idx],
                                                          PGFPLOTSMARKERS[marker_idx],
                                                          string("d",name),legend_name,"")
                                       )
                       end
                   else
                       str = string(str,tikz_add_plot($(color),"solid",PGFPLOTSMARKERS[1],
                                                      string($(name),"1"),"\$Q_\\ell\$",""))
                       str = string(str,tikz_add_plot($(color),LINESTYLES[1],PGFPLOTSMARKERS[1],
                                                      string("d",$(name),"1"),"\$Q_\\ell-Q_{\\ell-1}\$",""))
                   end
                   str = string(str,tikz_footer())
                   write(f, str)
               end
           end
          )
    eval(ex)
end

function write_figure_W(h::History, folder::AbstractString, fname::AbstractString)
    open(joinpath(folder,"figures","W.tex"), "w") do f
        L = maximum(maximum.(h[:index_set]))
        d = h[:ndims]
        str = tikz_header("level \$\\ell\$",string("\$\\log_2(C_\\ell)\$"),
                          "xmin=0,\n xmax=$(L),\n xtick={0,1,...,$(L)},\n",fname,"",true)
        for i in 2:2^d
            proto_idx = ind2sub(tuple([2 for i = 1:d]...),i).-1
            marker_idx = mod(i-2,length(PGFPLOTSMARKERS))+1
            linestyle_idx = mod(sum(proto_idx)-1,length(LINESTYLES))+1
            name = string("W",join(string.(proto_idx)))
            legend_name = string("\$W_{",replace(string(proto_idx),"1","\\ell"),"}\$")
            str = string(str,tikz_add_plot("green","solid",PGFPLOTSMARKERS[marker_idx],name,legend_name,""))
        end
        str = string(str,tikz_footer())
        write(f, str)
    end
end

function write_figure_samples(h::History, folder::AbstractString, fname::AbstractString)
    d = h[:ndims]
    if d == 1
        open(joinpath(folder,"figures","samples.tex"), "w") do f
            L = maximum(maximum.(h[:index_set]))
            str = tikz_header("level \$\\ell\$","number of samples \$N_\\ell\$",
                              "ymode=log,\n xmin=0,\n xmax=$(L),\n xtick={0,1,...,$(L)},\n",fname,"",true)
            cm = ColorMap("jet")
            v = linspace(0,1,h.iter)
            for i in 1:h.iter
                color = cm(v[i])
                color_name = "{rgb,1:red,$(color[1]); green,$(color[2]); blue,$(color[3])}"
                str = string(str,tikz_add_plot(color_name,"solid","*","samples","","x index = {0}, y index = {$(i)}"))
            end
            str = string(str,tikz_footer())
            write(f, str)
        end
    end
end

function write_figure_sample_reuse(h::History, folder::AbstractString, fname::AbstractString)
    if h[:method] isa MG
        d = h[:ndims]
        if d == 1
            write_figure_sample_reuse_1d(h, folder, fname)
        elseif d == 2
            write_figure_sample_reuse_2d(h, folder, fname)
        end
    end
end

function write_figure_sample_reuse_1d(h::History, folder::AbstractString, fname::AbstractString)
    for i in 1:h.iter
        open(joinpath(folder,"figures","sample_reuse_$(i).tex"), "w") do f
            L = maximum(maximum.(h[:index_set]))
            S,S_ = get_samples(h)
            M = maximum(S)
            str = tikz_bar_plot_1d_header(L,M,fname)
            str = string(str,tikz_add_bar_plot("samples_total_$(i)","red","original"))
            str = string(str,tikz_add_bar_plot("samples_reused_$(i)","blue","reused"))
            str = string(str,tikz_add_fractions_1d(i))
            str = string(str,tikz_footer())
            write(f, str)
        end
    end
    open(joinpath(folder,"figures","sample_reuse.tex"), "w") do f
        tols = [h[i][:tol] for i in 1:h.iter]
        str = tikz_samples_table(fname,tols)
        write(f, str)
    end
end

function write_figure_sample_reuse_2d(h::History, folder::AbstractString, fname::AbstractString)
    for i in 1:h.iter
        open(joinpath(folder,"figures","sample_reuse_$(i).tex"), "w") do f
            L = maximum(maximum.(h[:index_set]))
            S,S_ = get_samples(h)
            M = maximum(log10.(S))
            str = tikz_bar_plot_2d_header(L,M,fname)
            for idx in Iterators.product(colon.(1,map(i->size(S,i),1:h[:ndims]))...)
                if S[idx...,i] - S_[idx...,i] != 0
                    reused = log10(S[idx...,i] - S_[idx...,i])
                    str = string(str,tikz_add_cuboid(idx[1]-1,idx[2]-1,0,reused,"blue"))
                else
                    reused = 0
                end
                if S_[idx...,i] != 0
                    original = log10(S[idx...,i]) - reused
                    str = string(str,tikz_add_cuboid(idx[1]-1,idx[2]-1,reused,original,"red"))
                end
            end
            str = string(str,tikz_3d_bar_plot_legend())
            str = string(str,tikz_footer())
            write(f, str)
        end
    end
    open(joinpath(folder,"figures","sample_reuse.tex"), "w") do f
        tols = [h[i][:tol] for i in 1:h.iter]
        str = tikz_samples_table(fname,tols)
        write(f, str)
    end
end


function write_figure_index_set(h::History, folder::AbstractString, fname::AbstractString)
    d = h[:ndims]
    if d == 2
        write_figure_index_set_2(h, folder, fname)
    elseif d == 3
        write_figure_index_set_3(h, folder, fname)
    end
end

for d in ["2" "3"]
    ex = :(function $(Symbol("write_figure_index_set_",d))(h::History, folder::AbstractString, fname::AbstractString)
               for i = 1:h.iter
                   open(joinpath(folder,"figures","index_set_$(i).tex"), "w") do f
                       max_level = maximum(maximum.(h[:index_set]))
                       str = $(Symbol("tikz_index_set_",d,"d"))(fname,i,max_level,false,false,false,false,1,"")
                       write(f, str)
                   end
               end
               open(joinpath(folder,"figures","index_set.tex"), "w") do f
                   tols = [h[i][:tol] for i in 1:h.iter]
                   str = tikz_index_set_table(fname,tols)
                   write(f, str)
               end
           end
          )
    eval(ex)
end

function write_figure_adaptive_index_set(h::History, folder::AbstractString, fname::AbstractString)
    d = h[:ndims]
    if h[:method] isa AD || h[:method] isa MG{d,<:AD} where {d}
        if d == 2
            write_figure_adaptive_index_set_2(h, folder, fname)
        elseif d == 3
            write_figure_adaptive_index_set_3(h, folder, fname)
        end
    end
end

for d in ["2" "3"]
    ex = :(function $(Symbol("write_figure_adaptive_index_set_",d))(h::History, folder::AbstractString, fname::AbstractString)
               dicts = h[:adaptive_index_set]
               for i = 1:length(dicts)
                   open(joinpath(folder,"figures","adaptive_index_set_$(i).tex"), "w") do f
                       max_level = maximum(maximum.(keys(dicts[end])))+1
                       v = collect(values(dicts[i]))
                       has_active = any(v.==1) || any(v.==2)
                       has_old = any(v.==0)
                       has_maximum = any(v.==3)
                       str = $(Symbol("tikz_index_set_",d,"d"))(fname,i,max_level,true,has_active,has_old,has_maximum,1,"")
                       write(f, str)
                   end
               end
               open(joinpath(folder,"figures","adaptive_index_set.tex"), "w") do f
                   str = tikz_adaptive_index_set_table(fname,0:length(dicts)-1)
                   write(f, str)
               end
               for mode in ["active","old","maximum"]
                   open(joinpath(folder,"data","index_set_legend_$(mode).txt"), "w") do f
                       dim = h[:ndims]
                       for j in 1:dim
                           write(f, @sprintf("%i ",0))
                       end
                   end
                   open(joinpath(folder,"figures","index_set_legend_$(mode).tex"), "w") do f
                       max_level = maximum(maximum.(keys(dicts[end])))+1
                       scaling = (max_level+1.2)/1.2
                       write(f, $(Symbol("tikz_index_set_",d,"d"))(fname,string("legend_",mode),0,false,false,false,false,scaling,mode))
                   end
               end
           end
          )
    eval(ex)
end

for (name,color,ylabel) in zip(["runtime","cost","runtime_vs_rmse"],["blue","red","blue"],["run time","standard cost","run time"])
    ex = :(
           function $(Symbol("write_figure_",name))(h::History, folder::AbstractString, fname::AbstractString)
               open(joinpath(folder,"figures",string($(name),".tex")), "w") do f
                   str = tikz_header("accuracy",$(ylabel),"xmode = log,\n ymode = log",fname,"",true)
                   str = string(str,tikz_add_plot($(color),"solid","*",$(name),"",""))
                   str = string(str,tikz_footer())
                   write(f, str)
               end
           end
          )
    eval(ex)
end

function write_main_file(h::History,folder::AbstractString,fname::AbstractString)
    open(joinpath(folder,fname), "w") do f
        write(f, file_contents(h[:name],length(h[:index_set][1]),h[:method] isa AD || h[:method] isa MG{d,<:AD} where {d},isa(h[:method],MG)))
    end
end
=#
