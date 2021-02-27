## print.jl : useful print statements
#
# A collection of useful print statements for Estimators.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2021

#
# convenience functions
#
long(num) = sprintf1("%12.5e", num)

short(num) = sprintf1("%5.3f", num)

shorte(num) = sprintf1("%5.3e", num)

get_crayon() = crayon"blue bold"

highlight_cols(k) = Highlighter((v, i, j) -> any(j .== k), get_crayon())

#
# header and footer
#
function print_head_foot(strs...)
    pre = ["***" for _ in strs]
    data = hcat(pre, vcat(strs...), pre)
    pretty_table(data, noheader=true, vlines=[:begin, :end], columns_width=[3, 71, 3],
                 highlighters=highlight_cols([1, 3]), alignment=:l)
end

function print_header(estimator::Estimator, ϵ::Real)
    str1 = string("MultilevelEstimators.jl @", now()) 
    str2 = string("This is an ", estimator)
    str3 = string("Simulating ", first(split(estimator[:name], ".")))
    str4 = string("Tolerance on RMSE ϵ = ", shorte(ϵ))
    print_head_foot(str1, str2, str3, str4)
end

function print_footer()
    str1 = string("MultilevelEstimators.jl @", now())
    str2 = "Successfull termination"
    print_head_foot(str1, str2)
end

#
# index or level
#
print_elname(::Estimator{<:SL}) = "level"

print_elname(::Estimator{<:AbstractML}) = "level"

print_elname(::Estimator{<:AbstractMI}) = "index"

#
# status
#
function print_status(estimator::Estimator)
    # header
    header = [print_elname(estimator), "|Eℓ|", "|ΔEℓ|", "Vℓ", "ΔVℓ", "Nℓ", "Wℓ"]
    # data
    indices = collect(keys(estimator))
    E = [abs.(mean0(estimator, index)) for index in indices]
    dE = [abs.(mean(estimator, index)) for index in indices]
    V = [var0(estimator, index) for index in indices]
    dV = [var(estimator, index) for index in indices]
    N = [nb_of_samples(estimator, index) for index in indices]
    W = [cost(estimator, index) for index in indices]
    data = hcat(indices, E, dE, V, dV, N, W)
    # formatters
    e_fmt = ft_printf("%5.3e", [2, 3, 4, 5, 7])
    # print table
    pretty_table(data, header, header_crayon=get_crayon(), header_alignment=:c, 
                 formatters = e_fmt, highlighters=highlight_cols(1),
                 equal_columns_width=true)
end

#
# optimal number of samples
#
function print_optimal_nb_of_samples(estimator::Estimator, samples)
    # header
    header = [print_elname(estimator), "Nℓ"]
    # data
    indices = [key for key in keys(estimator) if samples[key] > 0]
    N = [print_nb_of_samples(estimator, index, samples[index]) for index in indices]
    data = hcat(indices, N)
    # formatters
    s_fmt = (v, i, j) -> sprintf1("%9s", v)
    # print table
    println("Samples will be updated ", estimator isa U ? "with" : "to", ":")
    pretty_table(data, header, header_crayon=get_crayon(), header_alignment=:c,
                 formatters = s_fmt, highlighters=highlight_cols(1))
end

print_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index, n::Integer) = string(n)

print_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index::Index, n::Integer) = string(estimator[:nb_of_shifts](index), " x ", n)

print_nb_of_samples(estimator::Estimator, index::Index) = print_nb_of_samples(estimator, index, nb_of_samples(estimator, index))

#
# profits
#
function print_largest_profit(estimator, max_index, max_profit, indices, profits)
    # header
    header = [print_elname(estimator), "profit"]
    # data
    data = hcat(indices, profits)
    # print table
    println("Profit indicators:")
    pretty_table(data, header, header_crayon=get_crayon(), header_alignment=:c,
                 formatters = ft_printf("%9s"), highlighters=highlight_cols(1))
    println("Index with max profit is ", max_index, " (value = ", long(max_profit), ").")
end

#
# probabilities
#
function print_pmf(estimator::Estimator)
    P = pmf(estimator)
    # header
    header = [print_elname(estimator), "P"]
    # data
    indices = sort(collect(keys(P)))
    data = hcat(indices, P[indices])
    # print table
    println("Using approximate probability mass function:")
    pretty_table(data, header, header_crayon=get_crayon(), header_alignment=:c,
                 formatters = ft_printf("%9s"), highlighters=highlight_cols(1))
end

#
# convergence
#
function print_convergence(estimator::Estimator, converged::Bool)
    print_status(estimator)
    converged && println(string("Convergence reached. RMSE ≈", long(rmse(estimator)), "."))
    print_footer()
end

function print_mse_analysis(estimator::Estimator, ϵ::Real, θ::Real)
    println("Checking convergence...")
    print_rates(estimator)
    println(string("  ==> Variance of the estimator ≈", long(varest(estimator)), "."))
    println(string("  ==> Bias of the estimator ≈", long(bias(estimator)), "."))
    println(string("  ==> MSE splitting parameter ≈ ", short(θ), "."))
    if !converged(estimator, ϵ, θ) > ϵ
        println(string("No convergence yet. RMSE ≈", long(rmse(estimator)), " > ", shorte(ϵ), "."))
        println("Adding an extra level...")
    end
end

function print_qmc_convergence(estimator::Estimator{<:AbstractIndexSet, <:QMC}, ϵ::Real, θ::Real)
    println(string("Checking if variance of estimator is larger than θ*ϵ^2..."))
    converged = !(varest(estimator) > θ*ϵ^2)
    println(string("  ==> ", converged ? "no!" : "yes", " (", long(varest(estimator))[2:end], converged ? " < " : " > ", short(θ), " *", long(ϵ^2), ")"))
end

function print_unbiased_convergence(estimator::Estimator{<:U}, ϵ::Real)
    println(string("Checking if variance of estimator is larger than ϵ^2..."))
    converged = !(varest(estimator) > ϵ^2)
    println(string("  ==> ", converged ? "no!" : "yes", " (", long(varest(estimator))[2:end], converged ? " < " : " > ", long(ϵ^2)[2:end], ")"))
    print_rates(estimator)
end

#
# taking ...
#
print_sample!_header(estimator::Estimator, index::Index, n::Integer, warm_up::Bool) = print(string("Taking ", print_nb_of_samples(estimator, index, n), print_warm_up(Val(warm_up)), " sample", print_with_s(n), " at ", print_elname(estimator), " ", index, "..."))

print_warm_up(::Val{true}) = " warm-up"

print_warm_up(::Val{false}) = " additional"

print_with_s(n::Integer) = n > 1 ? "s" : ""

print_sample!_footer() = println(" done")

#
# restart
#
function print_restart(index, restart, samples_dir, wildcard)
    if any(values(restart))
        param_file = joinpath(samples_dir, join(index.I, "_"), wildcard, "params.dat")
        println("\nWritten parameter values to ", param_file)
        print("for sample numbers ")
        tf = findall([restart[i] for i in sort(collect(keys(restart)))])
        if length(tf) < 11
            println(join(tf, ","), ".")
        else
            println(join(tf[1:5], ","), ", ..., ", join(tf[end-4:end], ","), ".")
        end
    end
end

#
# rates
#
function print_rates(estimator::Estimator)
    str = string("  ==> Rates: α ≈ ", print_rate(estimator, α))
    str = string(str, ", β ≈ ", print_rate(estimator, β))
    str = string(str, ", γ ≈ ", print_rate(estimator, γ), ".")
    println(str)
end

print_rate(estimator::Estimator{<:AbstractML}, f::Function) = short(f(estimator)[1])

print_rate(estimator::Estimator{<:AbstractMI}, f::Function) = string("(", join(short.(f(estimator)), ", "), ")")

#
# max level/index
#
warn_max_level(estimator::Estimator) = @warn string("Maximum ", _warn_max_level_name(estimator), " L = ", estimator[:max_index_set_param], " reached, no convergence yet.")

_warn_max_level_name(estimator::Estimator{<:AbstractML}) = "level"

_warn_max_level_name(estimator::Estimator{<:AbstractMI}) = "index set parameter"

warn_max_index(estimator::Estimator{<:AD}, index) = @warn string("Trying to add index ", index, " outside search space to index set, ignoring...")

#
# level parameter
#
print_level(estimator::Estimator{<:SL}, L) = println(string("Currently running on finest level."))

print_level(estimator::Estimator{<:ML}, L) = println(string("Currently running on level ", L, "."))

print_level(estimator::Estimator, L) = println(string("Currently running with L = ", L, "."))

#
# index set
#
print_index_set(estimator::Estimator, index_set) = nothing

print_index_set(estimator::Estimator{<:AbstractMI}, index_set) = begin
    if ndims(estimator) == 2
        println("Shape of the index set:")
        print(union(keys(estimator), index_set))
    end
end

print_index_set(estimator::Estimator{<:U}, index_set) = begin
    if ndims(estimator) == 2
        println("Shape of the index set:")
		print(filter(i -> has_samples_at_index(estimator, i), keys(estimator)))
    end
end

print_index_set(estimator::Estimator{<:AD}, index_set) = begin
    if ndims(estimator) == 2
        println("Shape of the index set (\u25A3 = active set):")
        print(collect(old_set(estimator)), collect(active_set(estimator)))
    end
end
