## report.jl : create a summary of the estimator
#
# Create a report that contains useful information about the estimator.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## report ##
function report(h::History, folder::AbstractString=first(split(h[:name], ".")))

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
function write_data_rates(h::History, folder::AbstractString)
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
function write_data_samples(h::History, folder::AbstractString)
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
function write_data_runtime(h::History, folder::AbstractString)
    x = [h[i][:tol] for i in 1:length(h)]
    y = cumsum([Dates.value(h[i][:time_stamp]-h.t_start)/1000 for i in 1:length(h)])
    open(joinpath(folder, "data", "runtime.txt"), "w") do f
        for i in 1:length(x)
            write(f, @sprintf("%12.5e %12.5e\n", x[i], y[i]))
        end
    end
end

function write_data_cost(h::History, folder::AbstractString)
    x = [h[i][:tol] for i in 1:length(h)]
    y = [sum(h[:W][index]*h[i][:nb_of_samples][index] for index in h[i][:index_set]) for i in 1:length(h)]
    open(joinpath(folder, "data", "cost.txt"), "w") do f
        for i in 1:length(x)
            write(f, @sprintf("%12.5e %12.5e\n", x[i], y[i]))
        end
    end
end

# TODO maybe use @generated function with dispatch in h[:type]
function write_figures_rates(h::History, folder::AbstractString)
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

function write_figure_samples(h::History, folder::AbstractString)
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

function write_figure_runtime(h::History, folder::AbstractString)
    plots = LinePlot[]
    plot = LinePlot(linecolor(1, 2), linestyle(1), marker(1), "" , "runtime", "")
    push!(plots, plot)
    pic = TikzPicture("accuracy \$\\epsilon\$", "runtime [s]", "", "", "", "log", "log", "", "", plots)
    open(joinpath(folder, "figures", "runtime.tex"), "w") do f
        write(f, generate(pic))
    end
end

function write_figure_cost(h::History, folder::AbstractString)
    plots = LinePlot[]
    plot = LinePlot(linecolor(2, 2), linestyle(1), marker(1), "" , "cost", "")
    push!(plots, plot)
    pic = TikzPicture("accuracy \$\\epsilon\$", "standard cost", "", "", "", "log", "log", "", "", plots)
    open(joinpath(folder, "figures", "cost.tex"), "w") do f
        write(f, generate(pic))
    end
end
